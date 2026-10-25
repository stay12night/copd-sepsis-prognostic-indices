## ============ 重建版：全套统计表 ============
suppressPackageStartupMessages({library(survival); library(rms)})
options(warn=-1)
B <- "C:/Users/stay1/Desktop/机器学习/脓毒症+谵妄+COPD/重建"
setwd(B); OUT <- file.path(B, "tables"); dir.create(OUT, showWarnings=FALSE)
d <- read.csv("data/cohort_all.csv", encoding="UTF-8", fileEncoding="UTF-8")
d$crrt_flag <- as.numeric(d$crrt_flag); d$crrt_flag[is.na(d$crrt_flag)] <- 0
d$t <- d$t90; d$e <- d$death90
D <- d[d$win72 == 1, ]                      # 主分析集

COV <- c("admission_age","male","race_white","myocardial_infarct","congestive_heart_failure",
         "peripheral_vascular_disease","cerebrovascular_disease","diabetes","renal_disease",
         "severe_liver_disease","hypertensive","sapsii_idx1","crrt_flag")
IDX <- c(PNI="PNI", RAR="RAR", Albumin="alb", RDW="rdw", Lymphocyte="lym",
         NLR="NLR", SIRI="SIRI", ACAG="ACAG", SII="SII", PLR="PLR")
rawf <- function(y,x){ok<-is.finite(x);y<-y[ok];x<-x[ok];r<-rank(x)
  (sum(r[y==1])-sum(y)*(sum(y)+1)/2)/(sum(y)*(length(y)-sum(y)))}
dirf <- function(y,x){a<-rawf(y,x); ifelse(a>0.5,a,1-a)}

## ---------- Table 1 ----------
lab <- c(admission_age="Age, y", BMI="BMI, kg/m2", heart_rate_idx1="Heart rate, bpm",
  sbp_idx1="Systolic BP, mmHg", mbp_idx1="Mean BP, mmHg", resp_rate_idx1="Respiratory rate, /min",
  temperature_idx1="Temperature, C", spo2_idx1="SpO2, %", wbc_idx1="White blood cells, x10^9/L",
  neut="Neutrophils, x10^9/L", lym="Lymphocytes, x10^9/L", mono="Monocytes, x10^9/L",
  plt="Platelets, x10^9/L", rdw="RDW, %", alb="Albumin, g/dL", bun_idx1="BUN, mg/dL",
  creatinine_idx1="Creatinine, mg/dL", glucose_idx1="Glucose, mg/dL", sodium_idx1="Sodium, mmol/L",
  potassium_idx1="Potassium, mmol/L", chloride_idx1="Chloride, mmol/L",
  bicarbonate_idx1="Bicarbonate, mmol/L", aniongap_idx1="Anion gap, mmol/L",
  inr_idx1="INR", ph_idx1="pH", pco2_idx1="PaCO2, mmHg", po2_idx1="PaO2, mmHg",
  ld_ldh_idx1="LDH, U/L", hemoglobin_idx1="Hemoglobin, g/dL",
  sofa_score="SOFA score", sapsii_idx1="SAPS II", apsiii_min="APS III",
  los_hospital="Hospital LOS, d", los_icu="ICU LOS, d",
  PNI="PNI", RAR="RDW-to-albumin ratio", ACAG="Albumin-corrected anion gap", NLR="NLR",
  SII="SII", SIRI="SIRI", PLR="PLR")
iqr <- function(z){m<-median(z,na.rm=TRUE); dg<-ifelse(m>=100,0,ifelse(m>=10,1,2))
  sprintf(paste0("%.",dg,"f (%.",dg,"f-%.",dg,"f)"), m, quantile(z,.25,na.rm=TRUE), quantile(z,.75,na.rm=TRUE))}
smdc <- function(x,g){m<-tapply(x,g,mean,na.rm=TRUE); s<-tapply(x,g,sd,na.rm=TRUE)
  mm<-mean(x,na.rm=TRUE); ss<-sd(x,na.rm=TRUE)
  sqrt(mean((m-mm)^2+(s-ss)^2,na.rm=TRUE))/ss}
t1 <- NULL
for (v in names(lab)){
  x <- D[[v]]; if (is.null(x)) next
  x[!is.finite(x)] <- NA; if (all(is.na(x))) next
  p <- wilcox.test(x ~ D$e)$p.value
  t1 <- rbind(t1, data.frame(Variable=lab[[v]], Survivors=iqr(x[D$e==0]), NonSurvivors=iqr(x[D$e==1]),
    P=ifelse(p<0.001,"<0.001",sprintf("%.3f",p)), SMD=sprintf("%.2f", smdc(x, D$e))))
}
addb <- function(name, v){
  v <- as.integer(v); p <- suppressWarnings(tryCatch(chisq.test(table(v, D$e))$p.value, error=function(x) NA))
  f <- function(z) sprintf("%d (%.1f)", sum(z==1,na.rm=TRUE), 100*mean(z==1,na.rm=TRUE))
  pp <- tapply(v==1, D$e, mean, na.rm=TRUE); po <- mean(v==1, na.rm=TRUE)
  data.frame(Variable=name, Survivors=f(v[D$e==0]), NonSurvivors=f(v[D$e==1]),
    P=ifelse(is.na(p),"",ifelse(p<0.001,"<0.001",sprintf("%.3f",p))),
    SMD=sprintf("%.2f", sqrt(mean((pp-po)^2,na.rm=TRUE)/(po*(1-po)))))
}
t1 <- rbind(t1, addb("Male sex, n (%)", D$male))
for (x in list(c("Myocardial infarction","myocardial_infarct"),c("Congestive heart failure","congestive_heart_failure"),
  c("Peripheral vascular disease","peripheral_vascular_disease"),c("Cerebrovascular disease","cerebrovascular_disease"),
  c("Diabetes","diabetes"),c("Renal disease","renal_disease"),c("Malignant cancer","malignant_cancer"),
  c("Severe liver disease","severe_liver_disease"),c("Hypertension","hypertensive"),
  c("Invasive mechanical ventilation","ventilator_flag"),c("CRRT","crrt_flag"),
  c("Acute kidney injury","aki_flag"),c("Delirium","delirium_flag"))){
  if (x[2] %in% names(D)) t1 <- rbind(t1, addb(paste0(x[1],", n (%)"), D[[x[2]]]))
}
t1 <- rbind(data.frame(Variable="No. of patients", Survivors=as.character(sum(D$e==0)),
  NonSurvivors=as.character(sum(D$e==1)), P="", SMD=""), t1)
write.csv(t1, file.path(OUT,"T1_baseline.csv"), row.names=FALSE, fileEncoding="UTF-8")

## ---------- Table 2: AUC ----------
set.seed(2026); Bt <- 2000
t2 <- NULL
for (k in names(IDX)){
  x <- D[[IDX[[k]]]]
  a <- dirf(D$e, x)
  bs <- replicate(Bt, {i <- sample(nrow(D), replace=TRUE); dirf(D$e[i], x[i])})
  t2 <- rbind(t2, data.frame(Index=k, n=sum(is.finite(x)), AUC=round(a,3),
    lcl=round(quantile(bs,.025),3), ucl=round(quantile(bs,.975),3),
    Direction=ifelse(rawf(D$e,x) > 0.5, "Higher is riskier", "Lower is riskier")))
}
write.csv(t2, file.path(OUT,"T2_auc.csv"), row.names=FALSE, fileEncoding="UTF-8")

## ---------- Table 3: Cox ----------
t3 <- NULL
for (k in names(IDX)){
  x <- D[[IDX[[k]]]]; sgn <- ifelse(rawf(D$e,x) > 0.5, 1, -1); x[!is.finite(x)] <- NA
  D$z <- as.numeric(scale(x)) * sgn
  mu <- coxph(Surv(t,e)~z, data=D); su <- summary(mu)
  ma <- coxph(as.formula(paste("Surv(t,e)~z+", paste(COV, collapse="+"))), data=D); sa <- summary(ma)
  t3 <- rbind(t3, data.frame(Index=k,
    HR_u=sprintf("%.2f (%.2f-%.2f)", su$conf.int[1,1],su$conf.int[1,3],su$conf.int[1,4]),
    P_u=ifelse(su$coefficients[1,5]<0.001,"<0.001",sprintf("%.3f",su$coefficients[1,5])),
    HR_a=sprintf("%.2f (%.2f-%.2f)", sa$conf.int[1,1],sa$conf.int[1,3],sa$conf.int[1,4]),
    P_a=ifelse(sa$coefficients[1,5]<0.001,"<0.001",sprintf("%.3f",sa$coefficients[1,5])),
    C_u=round(concordance(mu)$concordance,3), C_a=round(concordance(ma)$concordance,3)))
}
write.csv(t3, file.path(OUT,"T3_cox.csv"), row.names=FALSE, fileEncoding="UTF-8")

## ---------- Table 4: 竞争模型 ----------
t4 <- NULL
for (w in c(72, 24)){
  dd <- d[d[[paste0("win",w)]] == 1, ]
  dd$PNI_z <- -as.numeric(scale(dd$PNI)); dd$NLR_z <- as.numeric(scale(dd$NLR))
  dd$ALB_z <- -as.numeric(scale(dd$alb)); dd$LYM_z <- -as.numeric(scale(dd$lym))
  specs <- list("PNI alone"=c("PNI_z"), "NLR alone"=c("NLR_z"), "PNI + NLR"=c("PNI_z","NLR_z"),
                "Albumin alone"=c("ALB_z"), "Albumin + NLR"=c("ALB_z","NLR_z"),
                "Lymphocyte alone"=c("LYM_z"))
  for (nm in names(specs)){
    m <- coxph(as.formula(paste("Surv(t,e)~", paste(c(specs[[nm]], COV), collapse="+"))), data=dd)
    s <- summary(m)
    i <- which(rownames(s$conf.int) %in% specs[[nm]])[1]
    t4 <- rbind(t4, data.frame(Window=paste0(w," h"), N=nrow(dd), Events=sum(dd$e), Model=nm,
      Term=rownames(s$conf.int)[i],
      HR=sprintf("%.2f (%.2f-%.2f)", s$conf.int[i,1],s$conf.int[i,3],s$conf.int[i,4]),
      P=ifelse(s$coefficients[i,5]<0.001,"<0.001",sprintf("%.3f",s$coefficients[i,5])),
      C=round(concordance(m)$concordance,3)))
  }
}
write.csv(t4, file.path(OUT,"T4_compete.csv"), row.names=FALSE, fileEncoding="UTF-8")

## ---------- Table 5: 窗口敏感性 ----------
t5 <- NULL
for (w in c(72, 48, 24)){
  dd <- d[d[[paste0("win",w)]] == 1, ]
  for (k in c("PNI","RAR","Albumin","Lymphocyte","NLR","ACAG")){
    x <- dd[[IDX[[k]]]]; sgn <- ifelse(rawf(dd$e,x) > 0.5, 1, -1); x[!is.finite(x)] <- NA
    dd$z <- as.numeric(scale(x)) * sgn
    s <- summary(coxph(as.formula(paste("Surv(t,e)~z+", paste(COV, collapse="+"))), data=dd))
    t5 <- rbind(t5, data.frame(Window=paste0(w," h"), Index=k,
      HR=sprintf("%.2f (%.2f-%.2f)", s$conf.int[1,1],s$conf.int[1,3],s$conf.int[1,4]),
      P=ifelse(s$coefficients[1,5]<0.001,"<0.001",sprintf("%.3f",s$coefficients[1,5]))))
  }
}
write.csv(t5, file.path(OUT,"T5_windows.csv"), row.names=FALSE, fileEncoding="UTF-8")

## ---------- Table S1: 样条 ----------
dat <- datadist(D); options(datadist="dat")
ts1 <- NULL
for (k in c("PNI","RAR","NLR","alb")){
  f <- as.formula(paste0("Surv(t,e)~rcs(", k, ",4)+", paste(COV, collapse="+")))
  m <- cph(f, data=D, x=TRUE, y=TRUE); a <- anova(m)
  ts1 <- rbind(ts1, data.frame(Index=k, P_overall=signif(a[1,"P"],3),
    P_nonlinear=signif(a[grep("Nonlinear",rownames(a)),"P"],3)))
}
write.csv(ts1, file.path(OUT,"S1_rcs.csv"), row.names=FALSE, fileEncoding="UTF-8")

## ---------- 两两 AUC ----------
set.seed(11); Bp <- 1500
pw <- NULL
for (p in list(c("PNI","Albumin"), c("PNI","NLR"), c("PNI","ACAG"), c("PNI","RAR"),
               c("Lymphocyte","Albumin"), c("NLR","Albumin"), c("NLR","ACAG"), c("RAR","ACAG"))){
  x1 <- D[[IDX[[p[1]]]]]; x2 <- D[[IDX[[p[2]]]]]
  dif <- dirf(D$e,x1) - dirf(D$e,x2)
  bs <- replicate(Bp, {i<-sample(nrow(D),replace=TRUE); dirf(D$e[i],x1[i]) - dirf(D$e[i],x2[i])})
  pv <- 2*min(mean(bs<=0), mean(bs>=0))
  pw <- rbind(pw, data.frame(Comparison=paste(p[1],"vs",p[2]), Diff=round(dif,3),
    lcl=round(quantile(bs,.025),3), ucl=round(quantile(bs,.975),3),
    P=ifelse(pv<0.001,"<0.001",sprintf("%.3f",pv))))
}
write.csv(pw, file.path(OUT,"T6_auc_pairs.csv"), row.names=FALSE, fileEncoding="UTF-8")

cat("=== T2 AUC ===\n"); print(t2, row.names=FALSE)
cat("\n=== T3 Cox (per SD) ===\n"); print(t3, row.names=FALSE)
cat("\n=== T4 竞争 ===\n"); print(t4[,c("Window","Model","Term","HR","P")], row.names=FALSE)
cat("\n=== T5 窗口 ===\n"); print(t5, row.names=FALSE)
cat("\n=== S1 样条 ===\n"); print(ts1, row.names=FALSE)
cat("\n=== 两两 AUC ===\n"); print(pw, row.names=FALSE)
cat("\nN =", nrow(D), " events =", sum(D$e), "\n")
