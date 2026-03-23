# ==============================================================================
# NHANES 真实数据四场景对比分析 - 完整修正版
# ==============================================================================
library(osDesign); library(plyr); library(MESS); library(mvtnorm)
source("C:/Users/良/Desktop/文章/数据/真实数据-0217/calculate.R") # 确保包含补全后的函数

set.seed(2026)
t_val <- 0.32 

# 1. 数据准备
nhanes_raw <- read.csv("C:/Users/良/Desktop/文章/数据/真实数据-0217/nhanes_clinical_dataset.csv")
dat_final <- na.omit(data.frame(
  id = as.character(nhanes_raw$SEQN),
  y  = as.numeric(nhanes_raw$Diabetes),
  x1 = nhanes_raw$Age, x2 = nhanes_raw$Gender, x3 = nhanes_raw$BMI,
  x4 = nhanes_raw$Smoking, x5 = nhanes_raw$Hypertension,
  z1 = nhanes_raw$SII, z2 = nhanes_raw$TyG_Index
))
n_total <- nrow(dat_final)

# 预计算分层 G (用于 R-balanced 抽样)
fit_pre <- glm(y ~ x1+x2+x3+x4+x5, data = dat_final, family = binomial)
probs_pre <- predict(fit_pre, type = "response")
dat_final$G <- findInterval(probs_pre, c(-Inf, quantile(probs_pre[dat_final$y == 1], c(0.25, 0.5, 0.75)), Inf))

# ==============================================================================
# 2. 执行纵向四列计算
# ==============================================================================

# --- 第一列：Full (X-only) 基准 ---
cat(">>> 计算场景 1: Full (X-only)\n")
fit_c1 <- Dat_mleFit_general(dat_final, y ~ x1+x2+x3+x4+x5)
res_auc_c1 <- Dat_auc_mle_old(dat_final, fit_c1, model.matrix(fit_c1$model))
res_snb_c1 <- Dat_snb_mle_old(t_val, dat_final, fit_c1, model.matrix(fit_c1$model))

col1 <- c(res_auc_c1$auc, res_auc_c1$auc, res_auc_c1$sq_var_auc, NA, NA, NA, NA, res_snb_c1$snb, res_snb_c1$snb, sqrt(res_snb_c1$var_snb))

# --- 第二列：Benchmark (X+Z) 全量对比 ---
cat(">>> 计算场景 2: Benchmark (X+Z)\n")
fit_c2 <- Dat_mleFit_general(dat_final, y ~ x1+x2+x3+x4+x5+z1+z2)
res_auc_c2 <- Dat_auc_mle_old(dat_final, fit_c2, model.matrix(fit_c2$model))
res_snb_c2 <- Dat_snb_mle_old(t_val, dat_final, fit_c2, model.matrix(fit_c2$model))
res_idi_c2 <- Dat_idi_mle_final_full(dat_final, fit_c2, fit_c1) 
# 如果你没写全量版 NRI，这里暂时用 corrected 逻辑但传入全量格式
res_nri_c2 <- Dat_nri_mle_final_full(dat_final, fit_c2, fit_c1)

col2 <- c(res_auc_c1$auc, res_auc_c2$auc, res_auc_c2$sq_var_auc, res_idi_c2$idi, sqrt(res_idi_c2$var_idi), res_nri_c2$nri, sqrt(res_nri_c2$var_nri), res_snb_c1$snb, res_snb_c2$snb, sqrt(res_snb_c2$var_snb))

# --- 第三列：Case-Control 二阶段模拟 ---
cat(">>> 计算场景 3: Case-Control (1:1)\n")
# 1. 1:1 随机抽样对照
dat_cc_sample <- rbind(dat_final[dat_final$y == 1, ], 
                       dat_final[sample(which(dat_final$y == 0), sum(dat_final$y == 1)), ])
# 2. 【核心修正】Case-Control 视为不分层，强制将所有人 G 设为 1 以修复“下标出界”报错
dat_cc_temp <- dat_final; dat_cc_temp$G <- 1
dat_cc_sample$G <- 1

mle_list_cc <- Dat_format(dat_cc_temp, dat_cc_sample, "MLE")
fit_cc_new  <- Dat_mleFit_new(mle_list_cc)

res_auc_c3 <- Dat_auc_mle_new(mle_list_cc, fit_cc_new, n_total)
res_snb_c3 <- Dat_snb_mle_new(t_val, mle_list_cc, fit_cc_new)
res_idi_c3 <- Dat_idi_mle_final(mle_list_cc, fit_cc_new, fit_c1, dat_cc_temp) # 注意传入 dat_cc_temp
res_nri_c3 <- Dat_nri_mle_final_corrected(mle_list_cc, fit_cc_new, fit_c1, dat_cc_temp)

col3 <- c(res_auc_c1$auc, res_auc_c3$auc, sqrt(res_auc_c3$var_auc), res_idi_c3$idi, sqrt(res_idi_c3$var_idi), res_nri_c3$nri, sqrt(res_nri_c3$var_nri), res_snb_c1$snb, res_snb_c3$snb, sqrt(res_snb_c3$var_snb))

# --- 第四列：R-balanced 风险均衡抽样 ---
cat(">>> 计算场景 4: R-balanced\n")
# 计算抽样概率 p0_rb
case_counts <- xtabs(~y + G, data = dat_final)[2, ]
ctrl_counts <- xtabs(~y + G, data = dat_final)[1, ]
p0_rb <- pmin(1, round((case_counts * 1) / ctrl_counts, 3))

dat_rb_sample <- Dat_genR(dat_final, n_total, NULL, 4, p0_rb, 1)
mle_list_rb <- Dat_format(dat_final, dat_rb_sample, "MLE")
fit_rb_new  <- Dat_mleFit_new(mle_list_rb)

res_auc_c4 <- Dat_auc_mle_new(mle_list_rb, fit_rb_new, n_total)
res_snb_c4 <- Dat_snb_mle_new(t_val, mle_list_rb, fit_rb_new)
res_idi_c4 <- Dat_idi_mle_final(mle_list_rb, fit_rb_new, fit_c1, dat_final)
res_nri_c4 <- Dat_nri_mle_final_corrected(mle_list_rb, fit_rb_new, fit_c1, dat_final)

col4 <- c(res_auc_c1$auc, res_auc_c4$auc, sqrt(res_auc_c4$var_auc), res_idi_c4$idi, sqrt(res_idi_c4$var_idi), res_nri_c4$nri, sqrt(res_nri_c4$var_nri), res_snb_c1$snb, res_snb_c4$snb, sqrt(res_snb_c4$var_snb))

# ==============================================================================
# 3. 最终汇总与展示
# ==============================================================================
indicators <- c("AUC_old", "AUC_new", "AUC_new_se", "IDI", "IDI_se", "NRI", "NRI_se", "SNB_old", "SNB_new", "SNB_new_se")
final_table <- data.frame(Indicator = indicators, 
                          Full_X_Only = col1, 
                          Benchmark = col2, 
                          CaseControl = col3, 
                          R_balanced = col4)

cat("\n--- NHANES 糖尿病预测模型：四场景指标对比表 ---\n")
print(final_table, digits = 4, row.names = FALSE)

# ==============================================================================
# 1. 核心数据准备：解决“找不到对象”报错
# ==============================================================================
# 提取概率
prob_baseline <- predict(fit_c1$model, type = "response")
X_full_with_z <- as.matrix(cbind(1, dat_final[, c("x1","x2","x3","x4","x5","z1","z2")]))
prob_rb <- as.vector(1 / (1 + exp(-X_full_with_z %*% fit_rb_new$beta.mle)))

# 构造绘图所需的长格式数据框
calib_df <- rbind(
  data.frame(pred = prob_baseline, obs = dat_final$y, 模型 = "基础模型 (仅 X)"),
  data.frame(pred = prob_rb, obs = dat_final$y, 模型 = "优化模型 (X+Z)")
)

# ==============================================================================
# 2. 糖尿病预测模型校准曲线图 (保持原边界)
# ==============================================================================
library(ggplot2)



p1 <- ggplot(calib_df, aes(x = pred, y = obs, linetype = 模型, size = 模型)) +
  # 理想参考线
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey40") +
  # LOESS 平滑拟合
  geom_smooth(method = "loess", se = FALSE, color = "grey20") +
  # 线型与粗细映射
  scale_linetype_manual(values = c("基础模型 (仅 X)" = "longdash", "优化模型 (X+Z)" = "solid")) +
  scale_size_manual(values = c("基础模型 (仅 X)" = 0.8, "优化模型 (X+Z)" = 1.2)) +
  # 【设定 0-1 完整范围】
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "",
       x = "平均预测风险", y = "实际发生比例",
       linetype = "模型类型", size = "模型类型") +
  theme_classic() +
  # 优化图例：合并线型与粗细，使其在图例中显示清晰
  guides(linetype = guide_legend(override.aes = list(size = c(0.8, 1.2), linetype = c("longdash", "solid"))),
         size = "none") +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 14), # 标题居中
        text = element_text(family = "STHeiti"))

ggsave("校准曲线图_临床版.png", plot = p1, width = 8, height = 6, dpi = 300, bg = "white")

# ==============================================================================
# 3. 糖尿病预测决策曲线分析 (保持 0.5/0.4 边界)
# ==============================================================================
# 恢复为你要求的计算范围：仅到 0.5
thresholds <- seq(0.05, 0.5, by = 0.01)

dca_base <- sapply(thresholds, function(t) Dat_snb_mle_old(t, dat_final, fit_c1, model.matrix(fit_c1$model))$snb)
dca_rb <- sapply(thresholds, function(t) Dat_snb_mle_new(t, mle_list_rb, fit_rb_new)$snb)
pi_y_val <- mean(dat_final$y)
dca_all <- pi_y_val - (1 - pi_y_val) * (thresholds / (1 - thresholds))

dca_plot_df <- data.frame(阈值 = thresholds, 基础模型 = dca_base, 优化模型 = dca_rb, 全部干预 = dca_all, 不干预 = 0)



p2 <- ggplot(dca_plot_df, aes(x = 阈值)) +
  geom_line(aes(y = 全部干预, linetype = "全部干预"), color = "grey70", size = 0.6) +
  geom_line(aes(y = 不干预, linetype = "不干预"), color = "black", size = 0.6) +
  geom_line(aes(y = 基础模型, linetype = "基础模型 (仅 X)"), color = "grey50", size = 0.9) +
  geom_line(aes(y = 优化模型, linetype = "优化模型 (X+Z)"), color = "grey20", size = 1.3) +
  scale_linetype_manual(name = "决策策略", 
                        values = c("全部干预" = "dotted", "不干预" = "solid", 
                                   "基础模型 (仅 X)" = "longdash", "优化模型 (X+Z)" = "solid")) +
  # 严格保持你要求的边界：横轴到 0.5，纵轴到 0.4
  coord_cartesian(xlim = c(0.05, 0.5), ylim = c(-0.05, 0.4)) + 
  labs(title = "",
       x = "风险阈值概率", y = "标准化净获益") +
  theme_bw() +
  guides(linetype = guide_legend(override.aes = list(
    size = c(0.9, 1.3, 0.6, 0.6),
    color = c("grey50", "grey20", "grey70", "black"),
    linetype = c("longdash", "solid", "dotted", "solid")
  ))) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5), # 标题居中
        text = element_text(family = "STHeiti"))

ggsave("决策曲线分析_临床版.png", plot = p2, width = 8, height = 6, dpi = 300, bg = "white")

# 导出校准图数据 (平滑曲线所需的全量预测点)
write.csv(calib_df, "calibration_data.csv", row.names = FALSE)

# 导出决策曲线数据 (Threshold, Baseline, R-balanced, Treat_All)
write.csv(dca_plot_df, "dca_data.csv", row.names = FALSE)