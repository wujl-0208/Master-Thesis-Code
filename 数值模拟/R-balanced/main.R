# 加载所需包
require(osDesign)
require(plyr)
require(MESS)
require(mvtnorm)
require(ggplot2)
require(tidyr)

# 加载外部函数文件 (请确保路径正确)
source("C:/Users/良/Desktop/文章/代码/R-Banlanced/Data_sampling.R")
source("C:/Users/良/Desktop/文章/代码/R-Banlanced/calculate.R")

# ==============================================================================
# 3. 模拟主程序 (E-balanced Design + SNB 计算)
# ==============================================================================

# --- 3.1 模拟参数设定 ---
set.seed(0)
Prev <- "Rare"
Rho <- "Independent"
ratio <- 1
n <- 3000
K <- 100 # 模拟次数

# 设定指标阈值
q <- 0.2     # PCF 阈值
p <- 0.9     # PNF 阈值
t_snb <- 0.08 # 【新增】SNB 阈值 (根据患病率调整，Rare建议 0.1-0.2)

# 根据患病率设置真实的beta系数
if(Prev == "Rare") {
  beta_true <- log(c(0.03, 0.6, 1.6, 0.6, 1.5))
}
if(Prev == "Moder") {
  beta_true <- log(c(0.13, 0.6, 1.6, 0.6, 1.5))
}
numbeta_new <- length(beta_true)
numbeta_old <- length(beta_true) - 1
numG <- 4 # 定义分层数量
epsilon_beta <- abs(beta_true)/5

# --- 3.2 计算固定的抽样比率 p0 (E-balanced) ---
# 生成一次性的初始大数据集
x1_init <- rnorm(n)
x2_init <- runif(n)
x3_init <- rbinom(n, size=1, prob=0.2)
z_init <- rnorm(n)
X_init <- cbind(rep(1,n), x1_init, x2_init, x3_init, z_init)
y_init <- rbinom(n, size=1, prob=exp(X_init %*% beta_true)/(1+exp(X_init %*% beta_true)))

dat_init_temp <- data.frame(y=y_init, x1=x1_init, x2=x2_init, x3=x3_init, z=z_init)
glm_fit_init <- glm(y ~ 1 + x1 + x2 + x3, data = dat_init_temp, family = binomial(link = "logit"))
X_old_init <- model.matrix(y ~ 1 + x1 + x2 + x3, data = dat_init_temp)
prob_init <- exp(X_old_init %*% glm_fit_init$coefficients) / (1 + exp(X_old_init %*% glm_fit_init$coefficients))
Quan_init <- unname(quantile(prob_init[dat_init_temp$y == 1], c(0.25, 0.5, 0.75)))

G_init <- sapply(prob_init, function(x) {
  if (x <= Quan_init[1]) 1
  else if (x > Quan_init[1] && x <= Quan_init[2]) 2
  else if (x > Quan_init[2] && x <= Quan_init[3]) 3
  else 4
})
dat_init <- data.frame(y=y_init, G=G_init)

x0 <- xtabs(~y+G, data=dat_init)[2,]*ratio
p0 <- round(x0/xtabs(~y+G, data=dat_init)[1,], 3)
p0[is.nan(p0)] <- 0 

# --- 3.3 初始化结果存储对象 ---
beta_names_new <- paste0("beta_new_", 0:(numbeta_new - 1))
beta_names_old <- paste0("beta_old_", 0:(numbeta_old - 1))
beta_new_se_names <- paste0(beta_names_new, "_se")
beta_old_se_names <- paste0(beta_names_old, "_se")

# 【修改】增加 SNB 相关的列名
col_names <- c(
  "Iteration",
  "AUC_New", "AUC_New_se",
  "AUC_Old", "AUC_Old_se",
  "PCF_New", "PCF_New_se",
  "PCF_Old", "PCF_Old_se",
  "PNF_New", "PNF_New_se",
  "PNF_Old", "PNF_Old_se",
  "SNB_New", "SNB_New_se", # 【新增】
  "SNB_Old", "SNB_Old_se", # 【新增】
  "IDI", "IDI_se",
  "NRI", "NRI_se",
  beta_names_new, beta_new_se_names,
  beta_names_old, beta_old_se_names
)

all_results_df <- data.frame(matrix(NA, nrow = K, ncol = length(col_names)))
colnames(all_results_df) <- col_names

# --- 3.4 开始模拟循环 ---
cat("Starting E-balanced simulation with K =", K, "iterations...\n")
for(k in 1:K)
{
  tryCatch({
    # --- A. 数据生成 ---
    generated_data <- Dat_gen(n, beta_true, Rho)
    dat <- generated_data$dat
    dat_ol <- generated_data$dat_ol
    
    # --- B. E-balanced 分层 ---
    glm.fit_k <- glm(y ~ 1 + x1 + x2 + x3, data = dat, family = binomial(link = "logit"))
    X_old_k <- model.matrix(y ~ 1 + x1 + x2 + x3, data = dat)
    prob_k <- exp(X_old_k %*% glm.fit_k$coefficients) / (1 + exp(X_old_k %*% glm.fit_k$coefficients))
    Quan_k <- unname(quantile(prob_k[dat$y == 1], c(0.25, 0.5, 0.75), na.rm = TRUE)) 
    
    dat$G <- sapply(prob_k, function(x) {
      if (is.na(x) || x <= Quan_k[1]) 1 
      else if (x > Quan_k[1] && x <= Quan_k[2]) 2
      else if (x > Quan_k[2] && x <= Quan_k[3]) 3
      else 4
    })
    
    # --- C. 第二阶段抽样 ---
    dat_balanced_k <- Dat_genR(dat, n, beta_true, numG, p0, ratio)
    dat_mlelist <- Dat_format(dat, dat_balanced_k, "MLE")
    
    # --- D. 新模型拟合与计算 ---
    mleFit_list_new <- Dat_mleFit_new(dat_mlelist) 
    
    # 指标计算
    auc_list_mle_new <- Dat_auc_mle_new(dat_mlelist, mleFit_list_new, n)
    pcf_list_mle <- Dat_pcf_mle(q, dat_mlelist, mleFit_list_new)
    pnf_list_mle <- Dat_pnf_mle(p, dat_mlelist, mleFit_list_new)
    snb_list_mle <- Dat_snb_mle_new(t_snb, dat_mlelist, mleFit_list_new) # 【新增】
    
    # --- E. 旧模型拟合与计算 ---
    mleFit_list_old <- Dat_mleFit_old(dat_ol)
    X_old <- model.matrix(y ~ x1 + x2 + x3, data = dat_ol)
    
    # 指标计算
    auc_list_mle_old <- Dat_auc_mle_old(dat_ol, mleFit_list_old, X_old)
    pcf_list_mle_old <- Dat_pcf_mle_old(q, dat_ol, mleFit_list_old, X_old)
    pnf_list_mle_old <- Dat_pnf_mle_old(p, dat_ol, mleFit_list_old, X_old)
    snb_list_mle_old <- Dat_snb_mle_old(t_snb, dat_ol, mleFit_list_old, X_old) # 【新增】
    
    # --- F. 增量指标 IDI / NRI ---
    idi_list <- Dat_idi_mle_final(dat_mlelist, mleFit_list_new, mleFit_list_old, dat)
    nri_list <- Dat_nri_mle_final_corrected(dat_mlelist, mleFit_list_new, mleFit_list_old, dat)
    
    # --- G. 结果存储 ---
    all_results_df[k, "Iteration"] <- k
    
    # New Model
    all_results_df[k, "AUC_New"] <- auc_list_mle_new$auc
    all_results_df[k, "AUC_New_se"] <- sqrt(auc_list_mle_new$var_auc)
    all_results_df[k, "PCF_New"] <- pcf_list_mle$pcf
    all_results_df[k, "PCF_New_se"] <- sqrt(pcf_list_mle$var_pcf)
    all_results_df[k, "PNF_New"] <- pnf_list_mle$pnf
    all_results_df[k, "PNF_New_se"] <- sqrt(pnf_list_mle$var_pnf)
    all_results_df[k, "SNB_New"] <- snb_list_mle$snb        # 【新增】
    all_results_df[k, "SNB_New_se"] <- sqrt(snb_list_mle$var_snb) # 【新增】
    
    # Old Model
    all_results_df[k, "AUC_Old"] <- auc_list_mle_old$auc
    all_results_df[k, "AUC_Old_se"] <- auc_list_mle_old$sq_var_auc 
    all_results_df[k, "PCF_Old"] <- pcf_list_mle_old$pcf
    all_results_df[k, "PCF_Old_se"] <- pcf_list_mle_old$sq_var_pcf 
    all_results_df[k, "PNF_Old"] <- pnf_list_mle_old$pnf
    all_results_df[k, "PNF_Old_se"] <- pnf_list_mle_old$sq_var_pnf 
    all_results_df[k, "SNB_Old"] <- snb_list_mle_old$snb        # 【新增】
    all_results_df[k, "SNB_Old_se"] <- sqrt(snb_list_mle_old$var_snb) # 【新增】
    
    # IDI & NRI
    all_results_df[k, "IDI"] <- idi_list$idi
    all_results_df[k, "IDI_se"] <- sqrt(idi_list$var_idi)
    all_results_df[k, "NRI"] <- nri_list$nri
    all_results_df[k, "NRI_se"] <- sqrt(nri_list$var_nri)
    
    # Beta
    all_results_df[k, beta_names_new] <- mleFit_list_new$beta.mle
    all_results_df[k, beta_names_old] <- mleFit_list_old$beta.mle
    all_results_df[k, beta_new_se_names] <- sqrt(mleFit_list_new$beta.var.mle)
    all_results_df[k, beta_old_se_names] <- sqrt(mleFit_list_old$beta.var.full)     
    
    if(k %% 10 == 0) cat("Iteration:", k, "completed.\n")
    
  }, error = function(e) {
    cat("Error in iteration", k, ":", e$message, "\n")
  })
}

cat("Simulation finished.\n\n")

# --- 4. 汇总与分析结果 ---
cat("--- Final Summary Table ---\n")

get_summary_stats <- function(metric_name) {
  est_col <- all_results_df[[metric_name]]
  se_col <- all_results_df[[paste0(metric_name, "_se")]]
  
  mean_est <- mean(est_col, na.rm = TRUE)
  mean_se <- mean(se_col, na.rm = TRUE)
  empirical_sd <- sd(est_col, na.rm = TRUE)
  
  return(c(Mean_Estimate = mean_est, Mean_Theoretical_SD = mean_se, Empirical_SD = empirical_sd))
}

summary_auc_new <- get_summary_stats("AUC_New")
summary_auc_old <- get_summary_stats("AUC_Old")
summary_pcf_new <- get_summary_stats("PCF_New")
summary_pcf_old <- get_summary_stats("PCF_Old")
summary_pnf_new <- get_summary_stats("PNF_New")
summary_pnf_old <- get_summary_stats("PNF_Old")
summary_snb_new <- get_summary_stats("SNB_New") # 【新增】
summary_snb_old <- get_summary_stats("SNB_Old") # 【新增】
summary_idi <- get_summary_stats("IDI")
summary_nri <- get_summary_stats("NRI")

summary_df <- data.frame(
  Metric = c("AUC (New)", "AUC (Old)",
             "PCF (New)", "PCF (Old)",
             "PNF (New)", "PNF (Old)",
             "SNB (New)", "SNB (Old)", # 【新增】
             "IDI", "NRI"),
  rbind(summary_auc_new, summary_auc_old,
        summary_pcf_new, summary_pcf_old,
        summary_pnf_new, summary_pnf_old,
        summary_snb_new, summary_snb_old,
        summary_idi, summary_nri)
)

print(summary_df, row.names = FALSE, digits = 4)

# 保存结果
# require(writexl)
# write_xlsx(all_results_df, "simulation_results_with_SNB.xlsx")