# 嵌套模型似然比检验：证明"加白蛋白没有增量"不是只看 C 指数
suppressPackageStartupMessages({library(survival)})
# Resolve the repository root whether this file is run with Rscript or source()d
.this_file <- tryCatch(
  normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
  error = function(e) NA_character_)
if (is.na(.this_file)) .this_file <- normalizePath(sys.frame(1)$ofile)
B <- normalizePath(file.path(dirname(.this_file), ".."))
d <- read.csv(file.path(B, "data/cohort_all.csv"), fileEncoding = "UTF-8-BOM")
d <- d[d$win72 == 1, ]
COV <- c("admission_age","male","race_white","myocardial_infarct","congestive_heart_failure",
         "peripheral_vascular_disease","cerebrovascular_disease","diabetes","renal_disease",
         "malignant_cancer","severe_liver_disease","hypertensive","sapsii_idx1","crrt_flag")
for (v in c("crrt_flag")) { d[[v]] <- as.numeric(d[[v]]); d[[v]][is.na(d[[v]])] <- 0 }

z <- function(x) as.numeric(scale(x))
d$PNI_z <- z(d$PNI); d$ALB_z <- z(d$alb); d$LYM_z <- z(-d$lym)  # 淋巴越低越危险，取负号统一方向
d$NLR_z <- z(d$NLR)

f <- function(y, extra) as.formula(paste("Surv(t90, death90) ~", paste(c(extra, COV), collapse = " + ")))

cmp <- function(label, base_terms, add_terms) {
  dd <- d[complete.cases(d[, c("t90","death90", COV, base_terms, add_terms)]), ]
  m0 <- coxph(f(dd, base_terms), data = dd)
  m1 <- coxph(f(dd, c(base_terms, add_terms)), data = dd)
  a <- anova(m0, m1)
  c0 <- summary(m0)$concordance[1]; c1 <- summary(m1)$concordance[1]
  cat(sprintf("%-34s n=%d  C: %.4f -> %.4f (ΔC=%+.4f)  LRT chi2=%.3f df=%d  P=%.4f\n",
              label, nrow(dd), c0, c1, c1 - c0,
              a$Chisq[2], a$Df[2], a$`Pr(>|Chi|)`[2]))
  invisible(c(c0, c1, a$`Pr(>|Chi|)`[2]))
}

cat("=== 嵌套模型似然比检验（校正同一组协变量）===\n")
cmp("Lymphocyte -> + Albumin",   "LYM_z", "ALB_z")
cmp("PNI -> + NLR",              "PNI_z", "NLR_z")
cmp("NLR -> + PNI",              "NLR_z", "PNI_z")
cmp("Albumin -> + Lymphocyte",   "ALB_z", "LYM_z")
