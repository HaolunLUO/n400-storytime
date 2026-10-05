#!/usr/bin/env Rscript
# Baseline as covariate (Alday 2019) vs baseline subtraction, on the revision
# core model. Same rows, same fixed terms, same random effects as core.rds.
#   subtract:  dv_<roi> ~ core                  (baseline coefficient fixed at 1)
#   covariate: raw_<roi> ~ core + base_<roi>_c  (coefficient estimated)
#   none:      raw_<roi> ~ core                 (coefficient fixed at 0)
# ROI from the first argument: cp (default) or fz. Fz reuses the CP core
# fixed and random structure with dv_fz / raw_fz / base_fz.

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
})

if (!nzchar(Sys.getenv("SLURM_JOB_ID"))) {
  stop("Refusing to fit outside Slurm. Submit run_basecov.sbatch.")
}

args <- commandArgs(trailingOnly = TRUE)
roi <- if (length(args) >= 1) args[[1]] else "cp"
if (!(roi %in% c("cp", "fz"))) stop(paste("unknown roi", roi))
pool <- "/orcd/pool/005/haolun52/n400_storytime_v1"
in_csv <- file.path(pool, "baseline_cov", "lmm_input_basecov.csv")
core_rds <- file.path(pool, "revision", "lmm", "core.rds")
out_dir <- file.path(pool, "baseline_cov", if (roi == "cp") "lmm" else paste0("lmm_", roi))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ctrl <- lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 1000000L),
  check.conv.singular = .makeCC(action = "ignore", tol = 1e-4)
)

# Same filtering, factors, contrasts and age centre as fit_revision_lmm.R.
d <- read.csv(in_csv, stringsAsFactors = FALSE)
d <- d[d$keep == 1 & is.finite(d$dv_cp) & is.finite(d$dv_fz), ]
d$subject <- factor(d$subject)
d$item <- factor(d$item_id)
d$section <- factor(d$section)
d$group <- factor(d$group)
age_subjects <- tapply(d$age, d$subject, function(x) x[[1]])
age_center <- mean(age_subjects[is.finite(age_subjects)])
g1 <- c(TD = 2 / 3, dyslexia_normal_CAP = -1 / 3, dyslexia_atypical_CAP = -1 / 3)
g2 <- c(TD = 0, dyslexia_normal_CAP = 0.5, dyslexia_atypical_CAP = -0.5)
d$g1 <- unname(g1[as.character(d$group)])
d$g2 <- unname(g2[as.character(d$group)])
d$age_c <- d$age - age_center
d <- droplevels(d[is.finite(d$age), ])
# Raw units, grand-mean centred: the coefficient is directly the fraction of
# baseline that gets subtracted (1 = classic correction).
dv <- paste0("dv_", roi)
raw <- paste0("raw_", roi)
base_c <- paste0("base_", roi, "_c")
d[[base_c]] <- d[[paste0("base_", roi)]] - mean(d[[paste0("base_", roi)]])
cat("rows", nrow(d), "subjects", nlevels(d$subject), "items", nlevels(d$item),
    "age_center", age_center, "\n")

core <- readRDS(core_rds)
f_core <- formula(core)
if (nrow(model.frame(core)) != nrow(d)) {
  stop(sprintf("row mismatch with core.rds: %d vs %d", nrow(model.frame(core)), nrow(d)))
}
f_core <- update(f_core, as.formula(paste(dv, "~ .")))
cat("roi", roi, "formula:", deparse1(f_core), "\n")

specs <- list(
  subtract = f_core,
  covariate = update(f_core, as.formula(paste(raw, "~ . +", base_c))),
  none = update(f_core, as.formula(paste(raw, "~ .")))
)

key <- c("sur_z", "sur_z:age_c", "sur_z:g1", "sur_z:g2", base_c)
rows <- list()
fits <- list()
for (nm in names(specs)) {
  t0 <- proc.time()[["elapsed"]]
  m <- lmer(specs[[nm]], data = d, REML = TRUE, control = ctrl)
  el <- proc.time()[["elapsed"]] - t0
  fits[[nm]] <- m
  saveRDS(m, file.path(out_dir, paste0(nm, ".rds")))
  sink(file.path(out_dir, paste0(nm, "_summary.txt")))
  print(summary(m))
  sink()
  sm <- coef(summary(m))
  vc <- as.data.frame(VarCorr(m))
  msgs <- m@optinfo$conv$lme4$messages
  for (t in intersect(key, rownames(sm))) {
    rows[[paste(nm, t)]] <- data.frame(
      model = nm, term = t,
      estimate = sm[t, "Estimate"], se = sm[t, "Std. Error"],
      df = sm[t, "df"], p = sm[t, "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
  }
  rows[[paste(nm, "fit")]] <- data.frame(
    model = nm, term = c("residual_sd", "subject_sur_slope_sd", "singular", "elapsed_s"),
    estimate = c(
      sigma(m),
      vc$sdcor[vc$grp %in% c("subject", "subject.1") & vc$var1 %in% "sur_z"][1],
      as.numeric(isSingular(m, tol = 1e-4)),
      el
    ),
    se = NA, df = NA, p = NA, stringsAsFactors = FALSE
  )
  if (length(msgs)) cat(nm, "convergence:", paste(msgs, collapse = " | "), "\n")
  cat(nm, "done", round(el, 1), "s sigma", sigma(m), "\n")
}
res <- do.call(rbind, rows)
rownames(res) <- NULL
write.csv(res, file.path(out_dir, "basecov_comparison.csv"), row.names = FALSE)

# subtract vs covariate are nested (beta_base = 1 vs free) only after an offset,
# so compare the raw-DV models by ML LRT: none (beta=0) vs covariate (free),
# and covariate vs offset(base_<roi>_c) (beta=1, same as subtraction).
ml_cov <- lmer(specs$covariate, data = d, REML = FALSE, control = ctrl)
ml_none <- lmer(specs$none, data = d, REML = FALSE, control = ctrl)
ml_sub <- lmer(update(specs$none, as.formula(paste(". ~ . + offset(", base_c, ")"))),
               data = d, REML = FALSE, control = ctrl)
lrt <- data.frame(
  test = c("free_vs_no_correction", "free_vs_full_subtraction"),
  chisq = c(2 * (logLik(ml_cov) - logLik(ml_none)), 2 * (logLik(ml_cov) - logLik(ml_sub))),
  df = 1
)
lrt$p <- pchisq(lrt$chisq, 1, lower.tail = FALSE)
write.csv(lrt, file.path(out_dir, "basecov_lrt.csv"), row.names = FALSE)

print(res, digits = 4)
print(lrt, digits = 4)
cat("basecov done\n")
