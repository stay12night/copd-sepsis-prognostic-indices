# 模板文章里有、我们没有的统计诊断：PH 假定检验 + VIF + EPV
suppressPackageStartupMessages({library(survival)})
B <- "C:/Users/stay1/Desktop/机器学习/脓毒症+谵妄+COPD/重建"
d <- read.csv(file.path(B, "data/cohort_all.csv"), fileEncoding = "UTF-8-BOM")
d <- d[d$win72 == 1, ]
COV <- c("admission_age","male","race_white","myocardial_infarct","congestive_heart_failure",
         "peripheral_vascular_disease","cerebrovascular_disease","diabetes","renal_disease",
         "malignant_cancer","severe_liver_disease","hypertensive","sapsii_idx1","crrt_flag")
d$crrt_flag <- as.numeric(d$crrt_flag); d$crrt_flag[is.na(d$crrt_flag)] <- 0
z <- function(x) as.numeric(scale(x))
d$PNI_z <- z(d$PNI); d$ALB_z <- z(d$alb); d$LYM_z <- z(-d$lym); d$NLR_z <- z(d$NLR)

cat("=== 事件数与 EPV ===\n")
cat("事件:", sum(d$death90, na.rm=TRUE), " 参数:", 1+length(COV), " EPV =",
    round(sum(d$death90, na.rm=TRUE)/(1+length(COV)), 1), "\n\n")

cat("=== 比例风险假定（cox.zph，全局检验）===\n")
for (v in c("PNI_z","ALB_z","LYM_z","NLR_z")) {
  dd <- d[complete.cases(d[, c("t90","death90", COV, v)]), ]
  m <- coxph(as.formula(paste("Surv(t90, death90) ~", paste(c(v, COV), collapse="+"))), data=dd)
  zz <- cox.zph(m)
  cat(sprintf("%-8s n=%d  global chisq=%.2f  df=%-2d P=%.4f\n",
              v, nrow(dd), zz$table["GLOBAL","chisq"], zz$table["GLOBAL","df"], zz$table["GLOBAL","p"]))
  # 打印单个变量里 P<0.05 的
  tt <- zz$table[-nrow(zz$table), , drop=FALSE]
  bad <- rownames(tt)[tt[,"p"] < 0.05]
  if (length(bad)) cat("   违反假定:", paste(bad, collapse=", "), "\n")
}

cat("\n=== 协变量共线性（VIF，以 PNI 模型为例）===\n")
dd <- d[complete.cases(d[, c("t90","death90", COV, "PNI_z")]), ]
m <- coxph(as.formula(paste("Surv(t90, death90) ~", paste(c("PNI_z", COV), collapse="+"))), data=dd)
X <- model.matrix(~ ., data = dd[, c("PNI_z", COV)])
vif <- sapply(colnames(X)[-1], function(v) {
  r2 <- summary(lm(as.formula(paste("`"+v+"` ~ .")), data = as.data.frame(X[, -1, drop=FALSE])))$r.squared
  1/(1-r2)
})
vif <- sort(vif, decreasing = TRUE)
print(round(head(vif, 8), 2))
cat("最大 VIF =", round(max(vif),2), "（>5 提示共线，>10 严重）\n")
