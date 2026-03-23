require(osDesign)
require(plyr)
require(MESS)
require(mvtnorm)
require(ggplot2) # Added for plotting
require(tidyr)
source("C:/Users/良/Desktop/文章/代码/case-control/Data_sampling.R")
source("C:/Users/良/Desktop/文章/代码/case-control/calculate.R")
set.seed(0)
# 2. Case-control Design
Prev <- "Rare"
Rho <- 0.3#"Independent"
ratio <- 1
n <- 3000
if(Prev == "Rare")
{beta_true <- log(c(0.03,0.6,1.6,0.6,1.5))}
if(Prev == "Moder")
{beta_true <- log(c(0.13,0.6,1.6,0.6,1.5))}

epsilon_beta <- abs(beta_true)/5
if (any(epsilon_beta == 0)) {
  epsilon_beta[epsilon_beta == 0] <- 0.001
}
numbeta_new <- length(beta_true)
numbeta_old <- length(beta_true) - 1

x1_init <- rnorm(n)
x2_init <- runif(n)
x3_init <- rbinom(n, size=1, prob=0.2)
z_init <- rnorm(n)
X_init <- cbind(rep(1,n), x1_init, x2_init, x3_init, z_init)
y_init <- rbinom(n, size=1, prob=exp(X_init %*% beta_true)/(1+exp(X_init %*% beta_true)))
dat_init <- data.frame(cbind(y=y_init,x1=x1_init,x2=x2_init,x3=x3_init,z=z_init))

x0 <- sum(dat_init$y)*ratio
p0 <- round(x0/sum(1-dat_init$y), 3)
# --- 为PCF和PNF设置阈值 ---
q <- 0.2 # PCF 阈值
p <- 0.9 # PNF 阈值

K=1000
dat_cc <- list(NULL)

# --- 初始化统一的结果存储对象 ---
# 定义所有beta参数的名称
beta_names_new <- paste0("beta_new_", 0:(numbeta_new - 1))
beta_names_old <- paste0("beta_old_", 0:(numbeta_old - 1))
beta_new_se_names <- paste0(beta_names_new, "_se") # 新模型beta标准差的列名
beta_old_se_names <- paste0(beta_names_old, "_se") # 旧模型beta标准差的列名
# 定义统一数据框的所有列名
col_names <- c(
  "Iteration",
  "AUC_New", "AUC_New_se",
  "AUC_Old", "AUC_Old_se",
  "PCF_New", "PCF_New_se",
  "PCF_Old", "PCF_Old_se",
  "PNF_New", "PNF_New_se",
  "PNF_Old", "PNF_Old_se",
  "IDI", "IDI_se",
  "NRI", "NRI_se",
  beta_names_new, beta_new_se_names,
  beta_names_old, beta_old_se_names
)

# 创建空的 all_results_df 数据框
all_results_df <- data.frame(matrix(NA, nrow = K, ncol = length(col_names)))
colnames(all_results_df) <- col_names

# --- 开始模拟循环 ---
cat("Starting simulation with K =", K, "iterations...\n")
for(k in 1:K)
{
  tryCatch({
    # --- 1. 数据生成与准备 ---
    generated_data <- Dat_gen(n, beta_true, Rho)
    dat <- generated_data$dat
    dat_ol <- generated_data$dat_ol
    #dat$id <- 1:n # 强烈建议加入唯一ID，用于未来更可靠的对齐
    dat_cc_k <- Dat_gencc(dat, n, beta_true, p0, ratio)
    dat_mlelist <- Dat_format(dat, dat_cc_k, "MLE")
    #dat_ol$id <- 1:n
    
    # --- 2. 拟合与计算AUC ---
    mleFit_list_new <- Dat_mleFit_new(dat_mlelist)
    auc_list_mle_new <- Dat_auc_mle_new(dat_mlelist, mleFit_list_new, n)
    
    mleFit_list_old <- Dat_mleFit_old(dat_ol)
    X_old <- model.matrix(y ~ x1 + x2 + x3, data = dat_ol)
    auc_list_mle_old <- Dat_auc_mle_old(dat_ol, mleFit_list_old, X_old)
    pcf_list_mle_old <- Dat_pcf_mle_old(q, dat_ol, mleFit_list_old, X_old)
    pnf_list_mle_old <- Dat_pnf_mle_old(p, dat_ol, mleFit_list_old, X_old)

    
    # --- 3. 拟合与计算IDI ---
    
    idi_list <- Dat_idi_mle_final(dat_mlelist, mleFit_list_new, mleFit_list_old, dat)

    
    # --- 4. 拟合与计算NRI ---
    nri_list <- Dat_nri_mle_final_corrected(
      dat_mlelist, 
      mleFit_list_new, 
      mleFit_list_old,
      dat  # <-- 使用 dat, 而不是 dat_ol
    )
    
    # --- 4. 计算 PCF 和 PNF (基于新模型) ---
    pcf_list_mle <- Dat_pcf_mle(q, dat_mlelist, mleFit_list_new)
    pnf_list_mle <- Dat_pnf_mle(p, dat_mlelist, mleFit_list_new)
    # --- 将当前循环的所有结果存入 all_results_df ---
    all_results_df[k, "Iteration"] <- k
    
    # 新模型结果
    all_results_df[k, "AUC_New"] <- auc_list_mle_new$auc
    all_results_df[k, "AUC_New_se"] <- sqrt(auc_list_mle_new$var_auc)
    all_results_df[k, "PCF_New"] <- pcf_list_mle$pcf
    all_results_df[k, "PCF_New_se"] <- sqrt(pcf_list_mle$var_pcf)
    all_results_df[k, "PNF_New"] <- pnf_list_mle$pnf
    all_results_df[k, "PNF_New_se"] <- sqrt(pnf_list_mle$var_pnf)
    
    # 旧模型结果
    all_results_df[k, "AUC_Old"] <- auc_list_mle_old$auc
    all_results_df[k, "AUC_Old_se"] <- auc_list_mle_old$sq_var_auc
    all_results_df[k, "PCF_Old"] <- pcf_list_mle_old$pcf
    all_results_df[k, "PCF_Old_se"] <- pcf_list_mle_old$sq_var_pcf
    all_results_df[k, "PNF_Old"] <- pnf_list_mle_old$pnf
    all_results_df[k, "PNF_Old_se"] <- pnf_list_mle_old$sq_var_pnf
    
    # IDI 和 NRI 结果
    all_results_df[k, "IDI"] <- idi_list$idi
    all_results_df[k, "IDI_se"] <- sqrt(idi_list$var_idi)
    all_results_df[k, "NRI"] <- nri_list$nri
    all_results_df[k, "NRI_se"] <- sqrt(nri_list$var_nri)
    
    # Beta 参数估计值
    all_results_df[k, beta_names_new] <- mleFit_list_new$beta.mle
    all_results_df[k, beta_names_old] <- mleFit_list_old$beta.mle
    all_results_df[k, beta_new_se_names] <- sqrt(mleFit_list_new$beta.var.mle)
    all_results_df[k, beta_old_se_names] <- sqrt(mleFit_list_old$beta.var.full)
    # 打印进度
    
    if(k %% 10 == 0)
    {
      cat("Iteration:", k, "completed.\n")
    }
 
  }, error = function(e) {
    # 如果发生错误，打印错误信息和迭代次数
    cat("Error in iteration", k, ":", e$message, "\n")
  })
}

cat("Simulation finished.\n")

# --- 5. 汇总与分析结果 (简化版) ---
cat("\n\n--- Final Summary Table ---\n")

# 定义一个辅助函数来计算汇总统计量
get_summary_stats <- function(metric_name) {
  est_col <- all_results_df[[metric_name]]
  se_col <- all_results_df[[paste0(metric_name, "_se")]]
  
  mean_est <- mean(est_col, na.rm = TRUE)
  mean_se <- mean(se_col, na.rm = TRUE)
  empirical_sd <- sd(est_col, na.rm = TRUE)
  
  return(c(Mean_Estimate = mean_est, Mean_Theoretical_SD = mean_se, Empirical_SD = empirical_sd))
}

# 计算所有指标的统计量
summary_auc_new <- get_summary_stats("AUC_New")
summary_auc_old <- get_summary_stats("AUC_Old")
summary_pcf_new <- get_summary_stats("PCF_New")
summary_pcf_old <- get_summary_stats("PCF_Old")
summary_pnf_new <- get_summary_stats("PNF_New")
summary_pnf_old <- get_summary_stats("PNF_Old")
summary_idi <- get_summary_stats("IDI")
summary_nri <- get_summary_stats("NRI")

# 创建最终的汇总摘要表
summary_df <- data.frame(
  Metric = c("AUC (New Model)", "AUC (Old Model)",
             "PCF (New Model)", "PCF (Old Model)",
             "PNF (New Model)", "PNF (Old Model)",
             "IDI", "NRI"),
  rbind(summary_auc_new, summary_auc_old,
        summary_pcf_new, summary_pcf_old,
        summary_pnf_new, summary_pnf_old,
        summary_idi, summary_nri)
)

print(summary_df, row.names = FALSE, digits = 4)

# --- (可选) Beta 参数的总结 ---
cat("\n\n--- Beta Coefficient Estimation Summary ---\n")

# --- 新模型 Beta 总结 ---
beta_new_estimates <- all_results_df[, beta_names_new]
beta_new_theoretical_ses <- all_results_df[, beta_new_se_names]
summary_beta_new <- data.frame(
  Model               = "New Model",
  Parameter           = beta_names_new,
  True_Value          = beta_true,
  Mean_Estimate       = colMeans(beta_new_estimates, na.rm = TRUE),
  Mean_Theoretical_SD = colMeans(beta_new_theoretical_ses, na.rm = TRUE), # 新增
  Empirical_SD        = apply(beta_new_estimates, 2, sd, na.rm = TRUE)
)
print(summary_beta_new, row.names = FALSE, digits = 4)
cat("\n")

# --- 旧模型 Beta 总结 ---
beta_old_estimates <- all_results_df[, beta_names_old]
beta_old_theoretical_ses <- all_results_df[, beta_old_se_names]
summary_beta_old <- data.frame(
  Model               = "Old Model",
  Parameter           = beta_names_old,
  True_Value          = beta_true[1:numbeta_old],
  Mean_Estimate       = colMeans(beta_old_estimates, na.rm = TRUE),
  Mean_Theoretical_SD = colMeans(beta_old_theoretical_ses, na.rm = TRUE), # 新增
  Empirical_SD        = apply(beta_old_estimates, 2, sd, na.rm = TRUE)
)
print(summary_beta_old, row.names = FALSE, digits = 4)
# --- (可选) 保存详细结果到 Excel ---
#require(writexl)
#file_path <- "C:/Users/良/Desktop/文章/1000次模拟结果/case-control/cc_simulation_results_m-5-1.xlsx"
#write_xlsx(all_results_df, file_path)
#cat(paste("\nFull simulation results saved to:", file_path, "\n"))