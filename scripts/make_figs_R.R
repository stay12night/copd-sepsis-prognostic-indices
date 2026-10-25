# ============================================================
# R 版作图（ggplot2 + survminer），替换原先 Python 版本
# 输出：600 dpi PNG / TIFF / PDF / EPS，单栏 85mm、双栏 170mm
# ============================================================
suppressPackageStartupMessages({
  library(ggplot2); library(survival); library(survminer)
  library(dplyr); library(grid); library(patchwork)
})

B   <- "C:/Users/stay1/Desktop/机器学习/脓毒症+谵妄+COPD/重建"
FIG <- file.path(B, "figures")
TAB <- file.path(B, "tables")
dir.create(FIG, showWarnings = FALSE)

MM <- 1/25.4
W1 <- 85 * MM
W2 <- 170 * MM

# 配色：NPG 风格（与模板一致）
COL <- c(npg_red = "#E64B35", npg_teal = "#00A087", npg_blue = "#3C5488",
         npg_green = "#4DBBD5", npg_orange = "#F39B7F", grey = "#8491B4",
         dark = "#000000", light = "#BEBEBE")

theme_set(theme_classic(base_size = 8, base_family = "sans") +
          theme(axis.line = element_line(linewidth = 0.4),
                axis.ticks = element_line(linewidth = 0.4),
                plot.title = element_text(size = 9, face = "bold", hjust = 0),
                legend.key.size = unit(3, "mm")))

save_all <- function(p, name, w, h) {
  for (ext in c("png", "tif", "pdf", "eps")) {
    f <- file.path(FIG, paste0(name, ".", ext))
    if (ext == "png") ggsave(f, p, width = w, height = h, dpi = 600, bg = "white")
    if (ext == "tif") ggsave(f, p, width = w, height = h, dpi = 600, bg = "white",
                             device = "tiff", compression = "lzw")
    if (ext == "pdf") ggsave(f, p, width = w, height = h, device = cairo_pdf, bg = "white")
    if (ext == "eps") ggsave(f, p, width = w, height = h, device = cairo_ps,
                             bg = "white", family = "sans")
  }
  cat("  ->", name, "\n")
}

# ---------------- 数据 ----------------
d <- read.csv(file.path(B, "data", "cohort_all.csv"), fileEncoding = "UTF-8-BOM")
d <- d[d$win72 == 1, ]
cat("n =", nrow(d), " 死亡 =", sum(d$death90), "\n")

auc_val <- function(x, y) {          # y=1 死亡
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]
  w <- wilcox.test(x[y == 1], x[y == 0])$statistic
  w / (sum(y == 1) * sum(y == 0))
}
roc_df <- function(x, y, flip = FALSE, lab = "") {
  if (flip) x <- -x
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]
  th <- sort(unique(x), decreasing = TRUE)
  fpr <- sapply(th, function(t) mean(x[y == 0] >= t))
  tpr <- sapply(th, function(t) mean(x[y == 1] >= t))
  data.frame(fpr = c(0, fpr, 1), tpr = c(0, tpr, 1), index = lab)
}

# ============================================================
# 图 1  流程图
# ============================================================
fig1 <- function() {
  bx <- data.frame(
    y  = c(4, 3, 2, 1),
    xmin = 0.5, xmax = 5.0, ymin = c(3.65, 2.65, 1.65, 0.65), ymax = c(4.35, 3.35, 2.35, 1.35),
    lab = c("MIMIC-IV v3.1\nICU admissions\nn = 94,458",
            "Sepsis-3 at the first ICU admission\nICU stay > 24 h\nn = 14,952",
            "COPD\n(Charlson ICD mapping)\nn = 1,553",
            "First albumin and lymphocyte count\nwithin 72 h of ICU admission\nn = 531")
  )
  ex <- data.frame(
    xmin = 5.9, xmax = 10.3, ymin = c(2.70, 1.70, 0.70), ymax = c(3.55, 2.55, 1.55),
    lab = c("Excluded, n = 79,506\nnot a first ICU admission,\nno Sepsis-3, or ICU stay \u2264 24 h",
            "Excluded, n = 13,399\nno COPD",
            "Excluded, n = 1,022\nlaboratory values not\navailable within 72 h")
  )
  ar <- data.frame(x = c(2.75, 2.75, 2.75, 5.0, 5.0, 5.0),
                   xend = c(2.75, 2.75, 2.75, 5.9, 5.9, 5.9),
                   y = c(3.65, 2.65, 1.65, 3.1, 2.1, 1.1),
                   yend = c(3.35, 2.35, 1.35, 3.1, 2.1, 1.1))
  p <- ggplot() +
    geom_rect(data = bx, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "white", colour = "black", linewidth = 0.35) +
    geom_rect(data = ex, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "#FBEEEA", colour = COL["npg_red"], linewidth = 0.35) +
    geom_segment(data = ar, aes(x = x, xend = xend, y = y, yend = yend),
                 linewidth = 0.35, colour = "black",
                 arrow = arrow(length = unit(1.2, "mm"), type = "closed")) +
    geom_text(data = bx, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = lab),
              size = 2.1, lineheight = 0.95) +
    geom_text(data = ex, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = lab),
              size = 1.85, lineheight = 0.95) +
    coord_cartesian(xlim = c(0.2, 10.6), ylim = c(0.5, 4.5), expand = FALSE) +
    theme_void()
  save_all(p, "Fig1_flow", 94 * MM, 128 * MM)
}

# ============================================================
# 图 2  判别力
# ============================================================
fig2 <- function() {
  specs <- list(
    list(x = d$NLR, flip = FALSE, lab = "NLR", col = COL[["npg_red"]]),
    list(x = d$PNI, flip = TRUE,  lab = "PNI", col = COL[["npg_orange"]]),
    list(x = d$alb, flip = TRUE,  lab = "Albumin", col = COL[["grey"]]),
    list(x = d$sapsii_idx1, flip = FALSE, lab = "SAPS II", col = "black")
  )
  rr <- do.call(rbind, lapply(specs, function(s) {
    dd <- roc_df(s$x, d$death90, s$flip, s$lab)
    dd$auc <- auc_val(if (s$flip) -s$x else s$x, d$death90)
    dd$col <- s$col
    dd
  }))
  rr$label <- sprintf("%s (%.3f)", rr$index, rr$auc)
  rr$label <- factor(rr$label, levels = unique(rr$label))
  pA <- ggplot(rr, aes(fpr, tpr, colour = label)) +
    geom_abline(intercept = 0, slope = 1, linetype = "22", linewidth = 0.3, colour = COL[["light"]]) +
    geom_step(linewidth = 0.55) +
    scale_colour_manual(values = setNames(unique(rr$col), unique(rr$label)), name = "AUC") +
    scale_x_continuous("1 \u2212 specificity", limits = c(0, 1), expand = c(0.01, 0)) +
    scale_y_continuous("Sensitivity", limits = c(0, 1), expand = c(0.01, 0)) +
    coord_equal() +
    theme(legend.position = c(0.98, 0.03), legend.justification = c(1, 0),
          legend.background = element_blank(), legend.title = element_text(size = 7),
          legend.text = element_text(size = 6.5)) +
    ggtitle("A  ROC for 90-day mortality")

  t2 <- read.csv(file.path(TAB, "T2_auc.csv"), fileEncoding = "UTF-8-BOM")
  t2 <- t2[order(t2$AUC), ]
  t2$lab <- factor(t2$Index, levels = t2$Index)
  t2$grp <- ifelse(t2$Index == "SAPS II", "SAPS II", "Index")
  if (!"SAPS II" %in% t2$Index)
    t2 <- rbind(t2, data.frame(Index = "SAPS II", n = 531, AUC = 0.599,
                               lcl = 0.547, ucl = 0.650, Direction = "Higher is riskier",
                               lab = NA, grp = "SAPS II"))
  t2 <- t2[order(t2$AUC), ]; t2$lab <- factor(t2$Index, levels = t2$Index)
  pB <- ggplot(t2, aes(y = lab)) +
    geom_vline(xintercept = 0.5, linetype = "22", linewidth = 0.3) +
    geom_vline(xintercept = 0.599, linetype = "12", linewidth = 0.3, colour = COL[["grey"]]) +
    geom_errorbarh(aes(xmin = lcl, xmax = ucl, colour = grp), height = 0, linewidth = 0.45) +
    geom_point(aes(x = AUC, colour = grp), size = 1.2) +
    scale_colour_manual(values = c(Index = COL[["npg_blue"]], `SAPS II` = "black"), guide = "none") +
    scale_x_continuous("AUC (95% CI)", limits = c(0.46, 0.66), breaks = seq(0.5, 0.65, 0.05)) +
    ylab(NULL) +
    ggtitle("B  Discrimination of all indices")

  p <- pA + pB + plot_layout(widths = c(1, 1.1))
  save_all(p, "Fig2_ROC", W2, 66 * MM)
}

# ============================================================
# 图 3  森林图
# ============================================================
parse_hr <- function(s) {
  h  <- as.numeric(sub(" .*", "", s))
  ci <- gsub("^[^(]*\\(([^)]*)\\)$", "\\1", s)
  ci <- gsub("\u2013", "-", ci)
  lo <- suppressWarnings(as.numeric(sub("-.*", "", ci)))
  hi <- suppressWarnings(as.numeric(sub(".*-", "", ci)))
  data.frame(h = h, lo = lo, hi = hi)
}
fig3 <- function() {
  t3 <- read.csv(file.path(TAB, "T3_cox.csv"), fileEncoding = "UTF-8-BOM")
  u <- parse_hr(t3$HR_u); a <- parse_hr(t3$HR_a)
  df <- data.frame(index = t3$Index, hu = u$h, lu = u$lo, uu = u$hi,
                   ha = a$h, la = a$lo, ua = a$hi, p = as.numeric(t3$P_a), c = t3$C_a)
  df <- df[order(df$ha), ]
  df$index <- factor(df$index, levels = df$index)
  df$txt <- sprintf("%.2f (%.2f\u2013%.2f)", df$ha, df$la, df$ua)
  df$ptxt <- ifelse(df$p < 0.001, "<0.001", sprintf("%.3f", df$p))
  df$ctxt <- sprintf("%.3f", df$c)

  p <- ggplot(df, aes(y = index)) +
    geom_vline(xintercept = 1, linetype = "22", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = lu, xmax = uu), height = 0, linewidth = 0.45,
                   colour = COL[["light"]], position = position_nudge(y = 0.16)) +
    geom_point(aes(x = hu), size = 0.9, shape = 1, colour = COL[["grey"]],
               position = position_nudge(y = 0.16)) +
    geom_errorbarh(aes(xmin = la, xmax = ua), height = 0, linewidth = 0.55,
                   colour = COL[["npg_blue"]], position = position_nudge(y = -0.16)) +
    geom_point(aes(x = ha), size = 1.5, shape = 15, colour = COL[["npg_blue"]],
               position = position_nudge(y = -0.16)) +
    geom_text(aes(x = 1.74, label = txt), size = 2.1, hjust = 0) +
    geom_text(aes(x = 2.16, label = ctxt), size = 2.1, hjust = 0) +
    geom_text(aes(x = 2.38, label = ptxt, fontface = ifelse(df$p < 0.05, "bold", "plain")),
              size = 2.1, hjust = 0) +
    scale_x_log10(breaks = c(0.8, 1.0, 1.25, 1.5), limits = c(0.8, 2.75),
                  labels = c("0.8", "1.0", "1.25", "1.5")) +
    annotation_logticks(sides = "b", linewidth = 0.3, size = 1.5) +
    xlab("Hazard ratio per 1-SD (95% CI)") + ylab(NULL) +
    annotate("text", x = 1.74, y = length(levels(df$index)) + 0.75, label = "HR (95% CI)",
             size = 2.4, fontface = "bold", hjust = 0) +
    annotate("text", x = 2.16, y = length(levels(df$index)) + 0.75, label = "C-index",
             size = 2.4, fontface = "bold", hjust = 0) +
    annotate("text", x = 2.38, y = length(levels(df$index)) + 0.75, label = "P",
             size = 2.4, fontface = "bold", hjust = 0) +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(4, 52, 4, 4))
  save_all(p, "Fig3_forest", W2, 62 * MM)
}

# ============================================================
# 图 4  营养项 vs 炎症项
# ============================================================
fig4 <- function() {
  t4 <- read.csv(file.path(TAB, "T4_compete.csv"), fileEncoding = "UTF-8-BOM")
  ord <- c("Lymphocyte alone", "Albumin + NLR", "Albumin alone",
           "PNI alone", "PNI + NLR", "NLR alone")
  short <- c("Lymphocyte alone" = "Lymphocyte", "Albumin + NLR" = "Albumin + NLR*",
             "Albumin alone" = "Albumin", "PNI alone" = "PNI",
             "PNI + NLR" = "PNI + NLR\u2020", "NLR alone" = "NLR")
  pl <- lapply(c("72 h", "24 h"), function(w) {
    s <- t4[t4$Window == w, ]
    h <- parse_hr(s$HR)
    df <- data.frame(model = s$Model, h = h$h, lo = h$lo, hi = h$hi, p = as.numeric(s$P))
    df$model <- factor(df$model, levels = rev(ord))
    df$sig <- ifelse(df$p < 0.05, "P < 0.05", "P \u2265 0.05")
    df$txt <- sprintf("%.2f (%.2f\u2013%.2f)", df$h, df$lo, df$hi)
    ttl <- sprintf("%s   n = %d, deaths = %d", ifelse(w == "72 h", "A  Within 72 h", "B  Within 24 h"),
                   s$N[1], s$Events[1])
    ggplot(df, aes(y = model)) +
      geom_vline(xintercept = 1, linetype = "22", linewidth = 0.3) +
      geom_errorbarh(aes(xmin = lo, xmax = hi, colour = sig), height = 0, linewidth = 0.5) +
      geom_point(aes(x = h, colour = sig), size = 1.5) +
      geom_text(aes(x = 1.78, label = txt), size = 2.0, hjust = 0) +
      scale_colour_manual(values = c("P < 0.05" = COL[["npg_red"]], "P \u2265 0.05" = COL[["grey"]]),
                          name = NULL) +
      scale_x_log10(breaks = c(0.8, 1.0, 1.25, 1.5), limits = c(0.75, 2.35),
                    labels = c("0.8", "1.0", "1.25", "1.5")) +
      xlab("Hazard ratio (95% CI)") + ylab(NULL) +
      ggtitle(ttl) + coord_cartesian(clip = "off") +
      theme(legend.position = "none", plot.margin = margin(4, 34, 4, 4))
  })
  pl[[1]] <- pl[[1]] + scale_y_discrete(labels = function(x) short[x])
  pl[[2]] <- pl[[2]] + scale_y_discrete(labels = function(x) short[x]) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p <- pl[[1]] + pl[[2]] + plot_layout(widths = c(1, 1))
  save_all(p, "Fig4_compete", W2, 60 * MM)
}

# ============================================================
# 图 S1  KM（survminer，带风险表与置信带）
# ============================================================
figS1 <- function() {
  s <- d[!is.na(d$NLR), ]
  s$tert <- cut(s$NLR, breaks = quantile(s$NLR, c(0, 1/3, 2/3, 1)),
                labels = c("T1 (lowest)", "T2", "T3 (highest)"), include.lowest = TRUE)
  fit <- survfit(Surv(t90, death90) ~ tert, data = s)

  # 配色：深蓝 / 橙 / 砖红，三者明度与色相都拉开
  PAL <- c("#2B6CB0", "#E8A33D", "#C0392B")

  p <- ggsurvplot(
    fit, data = s,
    conf.int = FALSE,               # 去掉置信带，避免三层阴影叠成一片
    risk.table = TRUE,
    risk.table.height = 0.22,
    risk.table.title = "No. at risk",
    censor = TRUE, censor.shape = "|", censor.size = 1.8, censor.alpha = 0.5,
    palette = PAL,
    size = 0.75,                    # 线条加粗
    linetype = c("solid", "dashed", "dotted"),   # 线型+颜色双重区分，黑白打印也分得开
    legend.title = "",
    legend.labs = c("T1 (lowest)", "T2", "T3 (highest)"),
    legend = c(0.78, 0.92),         # 图例放右上，避开曲线
    xlab = "Days since ICU admission",
    ylab = "Survival probability",
    break.time.by = 30, xlim = c(0, 90), ylim = c(0.4, 1),
    break.y.by = 0.1,
    pval = FALSE,
    ggtheme = theme_classic(base_size = 8, base_family = "sans") +
      theme(
        axis.line    = element_line(linewidth = 0.45, colour = "black"),
        axis.ticks   = element_line(linewidth = 0.45, colour = "black"),
        axis.text    = element_text(colour = "black", size = 7),
        axis.title   = element_text(colour = "black", size = 8),
        legend.text  = element_text(size = 7),
        legend.key.width = unit(5, "mm"),
        legend.background = element_blank(),
        panel.grid = element_blank()
      ),
    tables.theme = theme_cleantable(base_size = 7, base_family = "sans") +
      theme(plot.title = element_text(size = 7, hjust = 0, face = "plain")),
    tables.y.text = TRUE,
    risk.table.y.text.col = TRUE,
    risk.table.y.text = TRUE,
    fontsize = 3,
    font.family = "sans"
  )

  W <- 96 * MM; H <- 118 * MM
  png(file.path(FIG, "FigS1_KM.png"), width = W, height = H, units = "in", res = 600)
  print(p, newpage = FALSE); dev.off()
  tiff(file.path(FIG, "FigS1_KM.tif"), width = W, height = H, units = "in",
       res = 600, compression = "lzw"); print(p, newpage = FALSE); dev.off()
  cairo_pdf(file.path(FIG, "FigS1_KM.pdf"), width = W, height = H)
  print(p, newpage = FALSE); dev.off()
  setEPS(); postscript(file.path(FIG, "FigS1_KM.eps"), width = W, height = H,
                       family = "sans", colormodel = "cmyk")
  print(p, newpage = FALSE); dev.off()
  cat("  -> FigS1_KM\n")
}

fig1(); fig2(); fig3(); fig4(); figS1()
cat("全部完成\n")
