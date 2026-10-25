# 协变量共线性诊断（VIF）
suppressPackageStartupMessages({library(survival)})
B <- "C:/Users/stay1/Desktop/机器学习/脓毒症+谵妄+COPD/重建"
d <- read.csv(file.path(B, "data/cohort_all.csv"), fileEncoding = "UTF-8-BOM")
d <- d[d$win72 == 1, ]
COV <- c("admission_age","male","race_white","myocardial_infarct","congestive_heart_failure",
         "peripheral_vascular_disease","cerebrovascular_disease","diabetes","renal_disease",
         "malignant_cancer","severe_liver_disease","hypertensive","sapsii_idx1","crrt_flag")
d$crrt_flag <- as.numeric(d$crrt_flag); d$crrt_flag[is.na(d$crrt_flag)] <- 0
d$PNI_z <- as.numeric(scale(d$PNI))

dd <- d[complete.cases(d[, c("t90","death90", COV, "PNI_z")]), ]
X <- as.data.frame(model.matrix(~ ., data = dd[, c("PNI_z", COV)])[, -1])

res <- c()
for (v in names(X)) {
  others <- setdiff(names(X), v)
  fit <- lm(X[[v]] ~ ., data = X[, others, drop = FALSE])
  r2 <- summary(fit)$r.squared
  res[v] <- 1 / (1 - r2)
}
res <- sort(res, decreasing = TRUE)
print(round(head(res, 10), 2))
cat("最大 VIF =", round(max(res), 2), "  (>5 提示共线, >10 严重)\n")
cat("样本:", nrow(dd), " 事件:", sum(dd$death90), "\n")
