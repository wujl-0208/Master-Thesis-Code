require(osDesign)
require(plyr)
require(MESS)
require(mvtnorm)
require(writexl)

# --- 1. 必需的函数 (来自您的代码) ---

Dat_gen <- function(n, beta, Rho){
  #n: total number of subjects in phase I
  #beta: log odds ratio
  
  # Phase I sample
  if (Rho == "Independent")
  {x1 <- rnorm(n)
  z <- rnorm(n)}
  if (Rho == "0.3")
  {sigma <- matrix(c(1, 0.3, 0.3, 1), nrow=2, byrow=TRUE)
  x1z <- rmvnorm(n, mean=c(0, 0), sigma=sigma)
  x1 <- x1z[,1]
  z <- x1z[,2]
  }
  if (Rho == "0.5")
  {sigma <- matrix(c(1, 0.5, 0.5, 1), nrow=2, byrow=TRUE) 
  x1z <- rmvnorm(n, mean=c(0, 0), sigma=sigma)
  x1 <- x1z[,1]
  z <- x1z[,2]
  }
  x2 <- runif(n)
  x3 <- rbinom(n, size=1, prob=0.2)
  X <- cbind(rep(1,n), x1, x2, x3, z)
  y <- rbinom(n, size=1, prob=exp(X %*% beta)/(1+exp(X %*% beta)))
  dat <- data.frame(cbind(y,x1,x2,x3,z))
  dat
}

Dat_fullFit <- function(data){
  # Directly fit the logistic model 
  m <- glm(y ~ x1 + x2 + x3 + z, family=binomial, data=data)
  yhat.full <- m$fitted.values
  beta.full <- m$coefficients
  beta.var.full <- summary(m)$cov.unscaled
  X <- cbind(rep(1,nrow(data)), data$x1, data$x2, data$x3, data$z) #Design matrix
  
  return(list(yhat.full=yhat.full, beta.full=beta.full, beta.var.full=diag(beta.var.full)))
}

# --- 2. 【修改后】计算 SNB, NB 及其理论方差的函数 ---

#' @param data 完整数据集 (来自 Dat_gen)
#' @param fullFit_list 完整模型拟合结果 (来自 Dat_fullFit)
#' @param threshold_t 用于计算 NB/SNB 的决策阈值 t
Dat_SNB_full <- function(data, fullFit_list, threshold_t){
  
  # --- 2.1 获取模型预测值 ---
  y_true <- data$y
  yhat_new <- fullFit_list$yhat.full # 新模型 (m1)
  
  # 拟合基准模型 m0 (不包含z)
  m0 <- glm(y ~ x1 + x2 + x3, family=binomial, data=data)
  yhat_old <- m0$fitted.values # 基准模型 (m0)
  
  # --- 2.2 基础统计量 ---
  N <- length(y_true)
  event_rate_y <- mean(y_true)
  n_cases <- sum(y_true == 1)
  n_controls <- sum(y_true == 0)
  t <- threshold_t
  
  # --- 2.3 内部辅助函数：计算点估计值 ---
  calculate_nb_snb_point <- function(yhat) {
    # 计算 Se 和 Sp
    true_positives_count <- sum(yhat > t & y_true == 1)
    false_positives_count <- sum(yhat > t & y_true == 0)
    true_negatives_count <- sum(yhat <= t & y_true == 0)
    
    Se <- true_positives_count / n_cases
    Sp <- true_negatives_count / n_controls
    
    # (1) NB(t)
    # 公式: Se * y - (1-Sp) * (1-y) * t/(1-t)
    # 或者是: (TP - FP * w) / N
    nb_model <- (true_positives_count / N) - (false_positives_count / N) * (t / (1 - t))
    
    # (2) SNB(t)
    nb_treat_all <- event_rate_y - (1 - event_rate_y) * (t / (1 - t))
    nb_strategy <- max(0, nb_treat_all)
    nb_perfect <- event_rate_y
    
    snb_denominator <- nb_perfect - nb_strategy
    if (snb_denominator < 1e-9) {
      snb_model <- NA 
    } else {
      snb_model <- (nb_model - nb_strategy) / snb_denominator
    }
    return(list(nb = nb_model, snb = snb_model, Se = Se, Sp = Sp))
  }
  
  metrics_new <- calculate_nb_snb_point(yhat_new)
  metrics_old <- calculate_nb_snb_point(yhat_old)
  
  # --- 2.4 计算 Se 和 Sp 的方差 (指示变量法) ---
  # Se 的指示变量 (仅在病例中): I(p > t)
  ind_se_new <- (yhat_new[y_true == 1] > t)
  ind_se_old <- (yhat_old[y_true == 1] > t)
  
  # Sp 的指示变量 (仅在对照中): I(p <= t)
  ind_sp_new <- (yhat_new[y_true == 0] <= t)
  ind_sp_old <- (yhat_old[y_true == 0] <= t)
  
  # 计算各个分量的方差 (Var(Mean) = Var(X) / n)
  var_se_new <- var(ind_se_new) / n_cases
  var_se_old <- var(ind_se_old) / n_cases
  
  var_sp_new <- var(ind_sp_new) / n_controls
  var_sp_old <- var(ind_sp_old) / n_controls
  
  # 计算差分的方差 (Delta)
  var_delta_se <- var(ind_se_new - ind_se_old) / n_cases
  var_delta_sp <- var(ind_sp_new - ind_sp_old) / n_controls
  
  # --- 2.5 应用方差公式计算 NB 和 SNB 的方差 ---
  
  # 公共系数 K (用于 NB)
  K_nb <- (1 - event_rate_y) * (t / (1 - t))
  
  # (A) NB 的方差: Var(NB) = y^2 * Var(Se) + K^2 * Var(Sp)
  var_nb_new <- (event_rate_y^2) * var_se_new + (K_nb^2) * var_sp_new
  var_nb_old <- (event_rate_y^2) * var_se_old + (K_nb^2) * var_sp_old
  
  # (B) SNB 和 Delta SNB 的方差
  # 根据 y <= t 或 y > t 选择权重
  var_snb_new <- NA
  var_snb_old <- NA
  var_delta_snb <- NA
  delta_snb_value <- NA
  
  delta_Se <- metrics_new$Se - metrics_old$Se
  delta_Sp <- metrics_new$Sp - metrics_old$Sp
  
  if (event_rate_y <= t) {
    # --- 情况 1: y <= t ---
    # 权重 W = t(1-y) / y(1-t)
    # Var(SNB) = Var(Se) + W^2 * Var(Sp)
    weight <- (t * (1 - event_rate_y)) / (event_rate_y * (1 - t))
    
    # Delta SNB 值
    delta_snb_value <- delta_Se + weight * delta_Sp
    
    # 方差计算
    var_snb_new <- var_se_new + (weight^2) * var_sp_new
    var_snb_old <- var_se_old + (weight^2) * var_sp_old
    var_delta_snb <- var_delta_se + (weight^2) * var_delta_sp
    
  } else {
    # --- 情况 2: y > t ---
    # 权重 W' = y(1-t) / t(1-y)
    # Var(SNB) = W'^2 * Var(Se) + Var(Sp)
    weight <- (event_rate_y * (1 - t)) / (t * (1 - event_rate_y))
    
    # Delta SNB 值
    delta_snb_value <- weight * delta_Se + delta_Sp
    
    # 方差计算
    var_snb_new <- (weight^2) * var_se_new + var_sp_new
    var_snb_old <- (weight^2) * var_se_old + var_sp_old
    var_delta_snb <- (weight^2) * var_delta_se + var_delta_sp
  }
  
  # --- 2.6 返回结果 ---
  return(list(
    NB_new = metrics_new$nb,
    NB_old = metrics_old$nb,
    SNB_new = metrics_new$snb,
    SNB_old = metrics_old$snb,
    Delta_SNB = delta_snb_value,
    
    # 返回标准差 (SE)
    NB_new_se = sqrt(var_nb_new),
    NB_old_se = sqrt(var_nb_old),
    SNB_new_se = sqrt(var_snb_new),
    SNB_old_se = sqrt(var_snb_old),
    Delta_SNB_se = sqrt(var_delta_snb)
  ))
}


# --- 3. 模拟设置 (Simulation Setup) ---
set.seed(0)
Prev <- "Moder"
Rho <- 0.5#"Independent"
n <- 3000
if (Prev == "Rare") {
  beta <- log(c(0.03, 0.6, 1.6, 0.6, 1.5))
}
if (Prev == "Moder") {
  beta <- log(c(0.13, 0.6, 1.6, 0.6, 1.5))
}
numbeta <- length(beta)

# 设定阈值
threshold_t_snb <- 0.3

K = 100 # 模拟次数

# --- 4. 初始化结果存储 (已扩充列名) ---
col_names <- c("Iteration", 
               "NB_new", "NB_new_se",
               "NB_old", "NB_old_se",
               "SNB_new", "SNB_new_se", 
               "SNB_old", "SNB_old_se",
               "Delta_SNB", "Delta_SNB_se")
all_results_df <- data.frame(matrix(NA, nrow = K, ncol = length(col_names)))
colnames(all_results_df) <- col_names

# --- 5. 主模拟循环 (Main Simulation Loop) ---
for (k in 1:K)
{
  tryCatch({
    data.full <- Dat_gen(n, beta, Rho)
    fullFit_list <- Dat_fullFit(data.full)
    
    # 调用修改后的函数
    snb_list_full <- Dat_SNB_full(data.full, fullFit_list, threshold_t = threshold_t_snb)
    
    # 存储结果
    all_results_df[k, "Iteration"] <- k
    
    all_results_df[k, "NB_new"] <- snb_list_full$NB_new
    all_results_df[k, "NB_new_se"] <- snb_list_full$NB_new_se
    
    all_results_df[k, "NB_old"] <- snb_list_full$NB_old
    all_results_df[k, "NB_old_se"] <- snb_list_full$NB_old_se
    
    all_results_df[k, "SNB_new"] <- snb_list_full$SNB_new
    all_results_df[k, "SNB_new_se"] <- snb_list_full$SNB_new_se
    
    all_results_df[k, "SNB_old"] <- snb_list_full$SNB_old
    all_results_df[k, "SNB_old_se"] <- snb_list_full$SNB_old_se
    
    all_results_df[k, "Delta_SNB"] <- snb_list_full$Delta_SNB
    all_results_df[k, "Delta_SNB_se"] <- snb_list_full$Delta_SNB_se
    
    if (round(k / 10) == k / 10)
    {
      print(k)
    }
  }, error = function(e) {
    cat("Error in iteration", k, ":", e$message, "\n")
  })
}

# --- 6. 模拟结果总结 ---
cat("\n--- Simulation Finished: Final Summary ---\n\n")

# 辅助函数
get_summary_stats <- function(metric_name) {
  est_col <- all_results_df[[metric_name]]
  
  # 自动查找对应的 se 列
  se_col_name <- paste0(metric_name, "_se")
  if (se_col_name %in% colnames(all_results_df)) {
    se_col <- all_results_df[[se_col_name]]
    mean_se <- mean(se_col, na.rm = TRUE)
  } else {
    mean_se <- NA 
  }
  
  mean_est <- mean(est_col, na.rm = TRUE)
  empirical_sd <- sd(est_col, na.rm = TRUE)
  
  return(c(Mean_Estimate = mean_est, Mean_Theoretical_SD = mean_se, Empirical_SD = empirical_sd))
}

# 计算所有指标的统计量
summary_nb_new <- get_summary_stats("NB_new")
summary_nb_old <- get_summary_stats("NB_old")
summary_snb_new <- get_summary_stats("SNB_new")
summary_snb_old <- get_summary_stats("SNB_old")
summary_delta_snb <- get_summary_stats("Delta_SNB")

# 创建指标汇总表
summary_table <- data.frame(
  Metric = c("NB_new", "NB_old", "SNB_new", "SNB_old", "Delta_SNB"),
  rbind(summary_nb_new, summary_nb_old, summary_snb_new, summary_snb_old, summary_delta_snb)
)
print(summary_table, row.names = FALSE, digits = 6)

# --- 7. 保存到 Excel (可选) ---
# file_path <- "SNB_simulation_results_complete.xlsx"
# write_xlsx(all_results_df, path = file_path)
# cat(paste("\n\n--- 详细结果已保存至 Excel ---\n"))