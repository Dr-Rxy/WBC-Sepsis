# Shared functions for the MIMIC-IV / eICU / NWICU analysis pipelines.

suppressPackageStartupMessages({
  library(data.table)
  library(cluster)   # silhouette(), clusGap()
  library(mclust)    # Mclust(), adjustedRandIndex()
  library(survival)  # coxph(), cox.zph()
  # MASS is loaded lazily inside .lr_confint() for profile-likelihood CIs.
})

# Shared constants

PALETTE <- c(C1 = "#C0392B", C2 = "#2471A3")  # subgroup colours
PALETTE_LAC <- c("FALSE" = "#2471A3", "TRUE" = "#C0392B")  # low lactate (blue), high lactate (red)

`%||%` <- function(x, y) if (is.null(x)) y else x

# WBC time windows (hours from time zero).
#
# The four windows are DISJOINT: together they partition [-8, 96) with no shared
# boundary point and no overlap, so every measurement contributes to exactly one
# clustering feature. Windows are named by their START hour, not by their centre:
# wbc_baseline starts 8 h before time zero and the other three start at 24, 48
# and 72 h.
#
# This differs from the published grid, which used closed intervals and a 32-h
# baseline window (-8, 24) that shared the 12-24 h range with wbc_t24h (12, 36).
# See docs/MODIFICATIONS.md for the full specification difference.
WBC_WINDOWS <- list(
  wbc_baseline = c(-8, 24), wbc_t24h = c(24, 48),
  wbc_t48h     = c(48, 72), wbc_t72h = c(72, 96)
)
CLUSTER_FEATURES <- names(WBC_WINDOWS)

LABELS <- list(
  c1_lab      = "C1 high-WBC subgroup",
  c2_lab      = "C2 normal-WBC subgroup"
)

# Logging and robustness

open_log <- function(path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  con <- file(path, open = "wt")
  sink(con, split = TRUE)          # stdout -> console + file
  cat(sprintf("=== run log opened: %s ===\n", path))
  con
}

close_log <- function(con) {
  cat("=== run log closed ===\n")
  sink()
  try(close(con), silent = TRUE)
}

# Registry of analysis blocks that raised an error, read by safe_all_ok()/safe_failed().
.safe_env <- new.env(parent = emptyenv())
.safe_env$failed <- character()

safe <- function(expr, label = "block") {
  tryCatch(
    expr,
    error = function(e) {
      .safe_env$failed <- c(.safe_env$failed, label)
      cat(sprintf("\n[SAFE] '%s' FAILED: %s\n", label, conditionMessage(e)))
      NULL
    }
  )
}

safe_all_ok <- function() length(.safe_env$failed) == 0L
safe_failed <- function() .safe_env$failed

# Data cleaning

#' In-place type coercion: bool→0/1 int, num_pattern→numeric; other char cols →numeric only if <5% new NAs.
coerce_types <- function(df, num_pattern, bool_cols = character()) {
  num_cols <- grep(num_pattern, names(df), value = TRUE)
  for (col in names(df)) {
    x <- df[[col]]
    if (!is.character(x)) next
    if (col %in% bool_cols) {
      set(df, j = col, value = as.integer(tolower(x) %in% c("t", "true", "1")))
    } else if (col %in% num_cols) {
      set(df, j = col, value = suppressWarnings(as.numeric(x)))
    } else {
      test <- suppressWarnings(as.numeric(x))
      if (sum(is.na(test)) - sum(is.na(x)) < nrow(df) * 0.05) set(df, j = col, value = test)
    }
  }
  for (col in intersect(c(num_cols, bool_cols), names(df))) {
    if (!is.numeric(df[[col]])) set(df, j = col, value = as.numeric(df[[col]]))
  }
  invisible(df)
}

build_wbc_windows <- function(long_dt, id_col, windows = WBC_WINDOWS) {
  ids <- sort(unique(long_dt[[id_col]]))
  out <- data.table(id = ids); setnames(out, "id", id_col)
  for (wn in names(windows)) {
    lo <- windows[[wn]][1]; hi <- windows[[wn]][2]
    w <- long_dt[hrs >= lo & hrs < hi,
                 .(val = median(wbc, na.rm = TRUE)), by = c(id_col)]
    setnames(w, "val", wn)
    out <- merge(out, w, by = id_col, all.x = TRUE)
  }
  out[]
}

missingness_report <- function(df, vars, label = "", warn_at = 0.20) {
  vars <- intersect(vars, names(df))
  rep <- data.table(
    Variable = vars,
    n_missing = sapply(vars, function(v) sum(is.na(df[[v]]))),
    pct_missing = sapply(vars, function(v) mean(is.na(df[[v]])) * 100)
  )
  if (nzchar(label)) cat(sprintf("\n--- Missingness: %s ---\n", label))
  print(rep)
  hi <- rep[pct_missing > warn_at * 100]
  if (nrow(hi) > 0) {
    warning(sprintf("High missingness (>%.0f%%): %s",
                    warn_at * 100, paste(hi$Variable, collapse = ", ")))
  }
  rep[]
}

# Imputation — MICE primary + complete-case sensitivity

#' MICE (PMM) imputation of the WBC windows; `aux` predictors are not clustered on.
impute_wbc_mice <- function(df, features = CLUSTER_FEATURES, aux = character(),
                            m = 5, seed = 2026, maxit = 5) {
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("Package 'mice' is required for the primary analysis. install.packages('mice')")
  }
  aux <- intersect(aux, names(df))
  sub <- as.data.frame(df[, c(features, aux), with = FALSE])
  mids <- mice::mice(sub, m = m, method = "pmm", maxit = maxit,
                     seed = seed, printFlag = FALSE)
  per_imp <- lapply(seq_len(m), function(i)
    as.matrix(mice::complete(mids, action = i)[, features, drop = FALSE]))
  Xmean <- Reduce(`+`, per_imp) / m
  colnames(Xmean) <- features
  list(X = Xmean, mids = mids, per_imp = per_imp)
}

complete_case_matrix <- function(df, features = CLUSTER_FEATURES) {
  ok <- stats::complete.cases(as.data.frame(df[, features, with = FALSE]))
  list(X = as.matrix(df[ok, features, with = FALSE]), rows = which(ok))
}

# Standardization — batch z-scores, store center/scale for reuse

standardize <- function(df, vars, suffix = "_z") {
  vars <- intersect(vars, names(df))
  params <- list()
  for (v in vars) {
    s <- scale(as.numeric(df[[v]]))
    set(df, j = paste0(v, suffix), value = as.numeric(s))
    params[[v]] <- c(center = attr(s, "scaled:center"), scale = attr(s, "scaled:scale"))
  }
  invisible(params)
}

apply_scaling <- function(X, center, scale) {
  Xs <- sweep(sweep(X, 2, center, "-"), 2, scale, "/")
  Xs[!is.finite(Xs)] <- 0
  Xs
}

# Clustering & external-prediction interface

#' K selection (within-SS, silhouette, Gap); K is taken from the silhouette maximum.
select_k <- function(X_scaled, kmax = 6, seed = 42, sub_n = 5000, gap_B = 50) {
  set.seed(seed)
  idx <- sample(seq_len(nrow(X_scaled)), min(sub_n, nrow(X_scaled)))
  d_sub <- dist(X_scaled[idx, ])
  ks <- 2:kmax
  wss <- sil <- numeric(length(ks))
  for (i in seq_along(ks)) {
    km <- kmeans(X_scaled, centers = ks[i], nstart = 25, iter.max = 50)
    wss[i] <- km$tot.withinss
    sil[i] <- summary(silhouette(km$cluster[idx], d_sub))$avg.width
  }
  gap_col <- rep(NA_real_, length(ks))
  gres <- safe({
    cg <- clusGap(X_scaled[idx, ], FUN = function(x, k)
      list(cluster = kmeans(x, k, nstart = 25, iter.max = 50)$cluster),
      K.max = kmax, B = gap_B, verbose = FALSE)
    cg$Tab[ks, "gap"]
  }, "clusGap")
  if (!is.null(gres)) gap_col <- gres
  data.table(K = ks, WSS = wss, Silhouette = sil, Gap = gap_col)
}

fit_kmeans <- function(X_scaled, k, nstart = 50, iter.max = 100, seed = 2026) {
  set.seed(seed)
  km <- kmeans(X_scaled, centers = k, nstart = nstart, iter.max = iter.max)
  list(km = km, centers = km$centers)
}

cluster_stability <- function(X_scaled, k, reference, nrep = 25, seed = 2026) {
  set.seed(seed)
  aris <- vapply(seq_len(nrep), function(i) {
    cl <- kmeans(X_scaled, centers = k, nstart = 1, iter.max = 100)$cluster
    adjustedRandIndex(cl, reference)
  }, numeric(1))
  list(mean = mean(aris), min = min(aris), aris = aris)
}

#' Relabel clusters so that the highest-baseline-WBC subgroup is C1 (=1).
align_labels_by_wbc <- function(cluster, wbc_baseline) {
  ord <- data.table(cl = cluster, wbc = wbc_baseline)[
    , .(m = mean(wbc, na.rm = TRUE)), by = cl][order(-m)]
  map <- setNames(seq_len(nrow(ord)), as.character(ord$cl))
  list(labels = as.integer(map[as.character(cluster)]),
       map = map,
       wbc_gap = if (nrow(ord) >= 2) abs(ord$m[1] - ord$m[2]) else NA_real_)
}

predict_cluster <- function(newX_raw, model) {
  X <- as.matrix(as.data.frame(newX_raw)[, model$features, drop = FALSE])
  mode(X) <- "numeric"
  for (j in seq_len(ncol(X))) {                       # median-fill any residual gaps
    na_idx <- is.na(X[, j])
    if (any(na_idx)) X[na_idx, j] <- median(X[, j], na.rm = TRUE)
  }
  Xs <- apply_scaling(X, model$center, model$scale)
  d <- as.matrix(proxy_dist(Xs, model$centers))
  apply(d, 1, which.min)
}

proxy_dist <- function(A, B) {
  a2 <- rowSums(A^2); b2 <- rowSums(B^2)
  sqrt(pmax(outer(a2, b2, "+") - 2 * A %*% t(B), 0))
}

export_cluster_model <- function(features, center, scale, centers, label_map,
                                 ref_stats, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  model <- list(features = features, center = center, scale = scale,
                centers = centers, label_map = label_map,
                ref_stats = ref_stats)
  saveRDS(model, path)
  model
}

# Modeling helpers

.lr_confint <- function(m, term) {
  ci <- safe(suppressMessages(stats::confint(m))[term, ], "confint(LR)")
  if (is.null(ci) || any(!is.finite(ci))) ci <- stats::confint.default(m)[term, ]
  ci
}

fit_glm_set <- function(df, formulas, term = "km_class") {
  out <- data.table()
  for (nm in names(formulas)) {
    res <- safe({
      m <- glm(as.formula(formulas[[nm]]), data = df, family = binomial)
      ci <- .lr_confint(m, term)
      data.table(Model = nm, OR = unname(exp(coef(m)[term])),
                 lower = unname(exp(ci[1])), upper = unname(exp(ci[2])),
                 p = summary(m)$coefficients[term, 4])
    }, paste0("glm:", nm))
    if (!is.null(res)) out <- rbind(out, res)
  }
  out
}

fit_cox_set <- function(df, formulas, term = "km_class") {
  out <- data.table()
  for (nm in names(formulas)) {
    res <- safe({
      fit <- coxph(formulas[[nm]], data = df, model = TRUE)  # model=TRUE so cox.zph
      s <- summary(fit)                                      # reuses the stored frame
      zph <- safe(cox.zph(fit), paste0("cox.zph:", nm))
      ph_p <- if (!is.null(zph) && term %in% rownames(zph$table))
        zph$table[term, "p"] else NA_real_
      ph_global <- if (!is.null(zph)) zph$table["GLOBAL", "p"] else NA_real_
      data.table(Model = nm, HR = s$coefficients[term, "exp(coef)"],
                 lower = s$conf.int[term, "lower .95"],
                 upper = s$conf.int[term, "upper .95"],
                 p = s$coefficients[term, "Pr(>|z|)"],
                 C_index = s$concordance["C"],
                 PH_p_term = ph_p, PH_p_global = ph_global)
    }, paste0("cox:", nm))
    if (!is.null(res)) out <- rbind(out, res)
  }
  out
}

fit_cox_timevarying <- function(df, term, covars = character(),
                                cutpoints = c(7, 14), time = "surv_time",
                                event = "event") {
  d <- as.data.frame(df)
  keep <- c(time, event, term, covars)
  d <- d[stats::complete.cases(d[, keep, drop = FALSE]), keep, drop = FALSE]
  sp <- survival::survSplit(
    stats::as.formula(sprintf("Surv(%s, %s) ~ .", time, event)),
    data = d, cut = cutpoints, episode = "tgroup", start = "tstart")
  cv <- if (length(covars)) paste("+", paste(covars, collapse = " + ")) else ""
  f <- stats::as.formula(sprintf("Surv(tstart, %s, %s) ~ %s:strata(tgroup)%s",
                                 time, event, term, cv))
  fit <- coxph(f, data = sp)
  s <- summary(fit)
  rows <- grep(sprintf("^%s:strata\\(tgroup\\)", term),
               rownames(s$coefficients), value = TRUE)
  bounds <- c(0, cutpoints, Inf)
  labs <- vapply(seq_along(rows), function(i)
    if (is.finite(bounds[i + 1])) sprintf("%g-%g d", bounds[i], bounds[i + 1])
    else sprintf(">%g d", bounds[i]), character(1))
  out <- data.table(
    Interval = labs,
    HR = s$coefficients[rows, "exp(coef)"],
    lower = s$conf.int[rows, "lower .95"],
    upper = s$conf.int[rows, "upper .95"],
    p = s$coefficients[rows, "Pr(>|z|)"])
  attr(out, "global_ph_p") <- safe(cox.zph(fit)$table["GLOBAL", "p"], "zph_tv")
  out[]
}

mediation_multi <- function(dm, treat, mediator, outcome, adjustment_sets,
                            nboot = 500, seed = 2026) {
  out <- data.table()
  for (nm in names(adjustment_sets)) {
    cv <- adjustment_sets[[nm]]
    res <- safe({
      f_a <- as.formula(sprintf("%s ~ %s%s", mediator, treat, cv))
      f_b <- as.formula(sprintf("%s ~ %s + %s%s", outcome, treat, mediator, cv))
      a_eff <- coef(lm(f_a, data = dm))[treat]
      mb    <- glm(f_b, data = dm, family = binomial)
      b_eff <- coef(mb)[mediator]; dir_eff <- coef(mb)[treat]
      ind_eff <- a_eff * b_eff
      tot_eff <- dir_eff + ind_eff
      prop <- ind_eff / tot_eff * 100
      # Per-sample ratio bootstrap: indirect/total is recomputed within each resample.
      set.seed(seed)
      boot_ratio <- replicate(nboot, {
        bd <- dm[sample(.N, .N, replace = TRUE)]
        tryCatch({
          aa <- coef(lm(f_a, data = bd))[treat]
          mbb <- glm(f_b, data = bd, family = binomial)
          bb <- coef(mbb)[mediator]; cc <- coef(mbb)[treat]
          ind <- aa * bb
          tot <- ind + cc
          ind / tot * 100
        }, error = function(e) NA_real_)
      })
      boot_ratio <- boot_ratio[is.finite(boot_ratio)]
      ci_prop <- quantile(boot_ratio, c(0.025, 0.975), na.rm = TRUE)
      data.table(Adjustment = nm, Total = tot_eff, Indirect = ind_eff,
                 Prop_mediated_pct = prop,
                 CI_lower_pct = ci_prop[1],
                 CI_upper_pct = ci_prop[2])
    }, paste0("mediation:", nm))
    if (!is.null(res)) out <- rbind(out, res)
  }
  out
}

#' Causal mediation (natural indirect/direct effects) on the risk-difference scale.
causal_mediation <- function(df, treat, mediator, outcome, confounders,
                             sims = 1000, seed = 2026) {
  if (!requireNamespace("mediation", quietly = TRUE)) {
    stop("Package 'mediation' required for causal mediation.")
  }
  cov_str <- paste(confounders, collapse = " + ")
  m_f <- sprintf("%s ~ %s + %s", mediator, treat, cov_str)
  y_f <- sprintf("%s ~ %s + %s + %s", outcome, treat, mediator, cov_str)
  # eval(parse()) embeds the formula inline so mediate() can re-evaluate the call
  med.fit <- eval(parse(text = sprintf("lm(%s, data = df)", m_f)))
  out.fit <- eval(parse(text = sprintf("glm(%s, data = df, family = binomial)", y_f)))
  set.seed(seed)
  s <- summary(mediation::mediate(med.fit, out.fit, treat = treat, mediator = mediator,
                                  boot = TRUE, sims = sims))
  data.table(
    ACME = s$d.avg, ACME_lo = s$d.avg.ci[1], ACME_hi = s$d.avg.ci[2], ACME_p = s$d.avg.p,
    ADE = s$z.avg, ADE_lo = s$z.avg.ci[1], ADE_hi = s$z.avg.ci[2], ADE_p = s$z.avg.p,
    Total = s$tau.coef, Total_lo = s$tau.ci[1], Total_hi = s$tau.ci[2],
    Prop_mediated_pct = s$n.avg * 100, CI_lower_pct = s$n.avg.ci[1] * 100,
    CI_upper_pct = s$n.avg.ci[2] * 100
  )
}

stratified_effects <- function(df, outcome, treat, strata_var, covars,
                               interaction_term) {
  cv <- if (length(covars)) paste("+", paste(covars, collapse = " + ")) else ""
  out <- data.table()
  for (lv in sort(unique(na.omit(df[[strata_var]])))) {
    sub <- df[get(strata_var) == lv]
    res <- safe({
      m <- glm(as.formula(sprintf("%s ~ %s%s", outcome, treat, cv)),
               data = sub, family = binomial)
      ci <- .lr_confint(m, treat)
      data.table(Stratum = sprintf("%s = %s", strata_var, lv), n = nrow(sub),
                 OR = unname(exp(coef(m)[treat])),
                 lower = unname(exp(ci[1])), upper = unname(exp(ci[2])),
                 p = summary(m)$coefficients[treat, 4])
    }, paste0("stratum:", lv))
    if (!is.null(res)) out <- rbind(out, res)
  }
  int_p <- safe({
    mi <- glm(as.formula(sprintf("%s ~ %s * %s%s", outcome, treat, strata_var, cv)),
              data = df, family = binomial)
    summary(mi)$coefficients[interaction_term, 4]
  }, "interaction_p")
  list(strata = out, interaction_p = if (is.null(int_p)) NA_real_ else int_p)
}

#' Additive-scale interaction (RERI) for two binary exposures; reference = both unexposed.
compute_reri <- function(df, outcome, expo_a, expo_b, covars = character(),
                         nboot = 2000, seed = 2026) {
  cv <- if (length(covars)) paste("+", paste(covars, collapse = " + ")) else ""
  ab <- paste0(expo_a, ":", expo_b)
  f <- as.formula(sprintf("%s ~ %s + %s + %s%s", outcome, expo_a, expo_b, ab, cv))
  point <- function(dat) {
    m <- glm(f, data = dat, family = binomial); b <- coef(m)
    or10 <- exp(b[expo_a]); or01 <- exp(b[expo_b])
    or11 <- exp(b[expo_a] + b[expo_b] + b[ab])
    list(or10 = unname(or10), or01 = unname(or01), or11 = unname(or11),
         reri = unname(or11 - or10 - or01 + 1),
         mult_or = unname(exp(b[ab])), mult_p = summary(m)$coefficients[ab, 4])
  }
  p0 <- point(df)
  set.seed(seed)
  bt <- replicate(nboot, {
    bd <- df[sample(.N, .N, replace = TRUE)]
    tryCatch(point(bd)$reri, error = function(e) NA_real_)
  })
  bt_valid <- bt[is.finite(bt)]
  ci <- if (length(bt_valid) >= 10) {
    quantile(bt_valid, c(0.025, 0.975))
  } else rep(NA_real_, 2)
  data.table(
    Measure = c("OR (A only)", "OR (B only)", "OR (A and B)",
                "Multiplicative interaction OR", "Multiplicative p",
                "RERI (additive)", "RERI CI lower", "RERI CI upper"),
    Value = c(p0$or10, p0$or01, p0$or11, p0$mult_or, p0$mult_p,
              p0$reri, ci[[1]], ci[[2]]))
}

# Table 1 — Two-group baseline characteristics

make_table1 <- function(df, group, vars_spec,
                        c1_name = "C1 high-WBC", c2_name = "C2 normal-WBC") {
  g <- df[[group]]
  t1 <- data.table()
  fmt_p <- function(p) ifelse(is.na(p), "-", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
# vars_spec entries: list(var, label, type[, digits[, how]])
  for (v in vars_spec) {
    var <- v[[1]]; label <- v[[2]]; vtype <- v[[3]]
    if (!var %in% names(df)) next
    dig <- if (length(v) >= 4) v[[4]] else 1L
    how <- if (length(v) >= 5) v[[5]] else "mean"
    label <- sub(", (mean \\(SD\\)|median \\(IQR\\))$", "", label)
    if (vtype == "cat") {
      for (lvl in na.omit(unique(df[[var]]))) {
        r1 <- mean(df[g == 1, get(var)] == lvl, na.rm = TRUE) * 100
        r2 <- mean(df[g == 2, get(var)] == lvl, na.rm = TRUE) * 100
        tbl <- table(g, df[[var]] == lvl)
        pval <- if (all(dim(tbl) == 2)) safe(chisq.test(tbl)$p.value, "chisq") else NA_real_
        t1 <- rbind(t1, setNames(data.table(paste0("  ", as.character(lvl)),
          sprintf("%.1f%%", r1), sprintf("%.1f%%", r2), fmt_p(pval)),
          c("Characteristic", c1_name, c2_name, "p value")))
      }
    } else {
      c1 <- na.omit(df[g == 1, get(var)]); c2 <- na.omit(df[g == 2, get(var)])
      if (vtype == "cont") {
        pval <- safe(wilcox.test(c1, c2)$p.value, "wilcox")
        fmt_num <- function(x) if (how == "mean") {
          sprintf(sprintf("%%.%df (%%.1f)", dig), mean(x), sd(x))
        } else {
          sprintf(sprintf("%%.%df (%%.%df\U2013%%.%df)", dig, dig, dig),
                  median(x), quantile(x, .25), quantile(x, .75))
        }
        t1 <- rbind(t1, setNames(data.table(
          sprintf("%s, %s", label, if (how == "mean") "mean (SD)" else "median (IQR)"),
          fmt_num(c1), fmt_num(c2), fmt_p(pval)),
          c("Characteristic", c1_name, c2_name, "p value")))
      } else {  # bin
        tbl <- matrix(c(sum(c1), length(c1) - sum(c1),
                        sum(c2), length(c2) - sum(c2)), nrow = 2)
        pval <- safe(chisq.test(tbl)$p.value, "chisq")
        t1 <- rbind(t1, setNames(data.table(label,
          sprintf("%d (%.1f%%)", as.integer(sum(c1)), mean(c1) * 100),
          sprintf("%d (%.1f%%)", as.integer(sum(c2)), mean(c2) * 100), fmt_p(pval)),
          c("Characteristic", c1_name, c2_name, "p value")))
      }
    }
  }
  t1
}

check_min_n <- function(group_vec, min_n = 50, what = "subgroup") {
  tab <- table(group_vec)
  if (any(tab < min_n)) {
    warning(sprintf("Skipping %s: group(s) below n=%d (%s)", what, min_n,
                    paste(sprintf("%s=%d", names(tab), tab), collapse = ", ")))
    return(FALSE)
  }
  TRUE
}
