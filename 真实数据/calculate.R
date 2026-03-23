# ==============================================================================
# calculate.R - 适配 NHANES (5 Baseline + 2 New Variables)
# ==============================================================================
# 2.3 Generate Phase II data with E-balanced Design
Dat_genR <- function(dat, n, beta, numG, p0, ratio){
  # 这里的 dat 已经是 NHANES 的完整数据集（包含 G 和 y）
  dat_case <- dat[dat$y==1, ]
  dat_con <- dat[dat$y==0, ]
  
  # Phase II sample
  dat2_con <- NULL
  x0 <- xtabs(~dat$y+dat$G)[2,]*ratio # 计算每一层预期的对照组数量
  
  for (i in 1:numG){
    dat_con_i <- dat_con[dat_con$G==i, ]
    selected_con <- 0
    tot_con <- 1
    
    ## 按照概率 p0[i] 在第 i 层中抽取对照组
    while (tot_con <= nrow(dat_con_i) & selected_con < x0[i]){
      choice <- rbinom(1, size=1, prob=p0[i])
      if (choice==1){
        dat2_con <- rbind(dat2_con, dat_con_i[tot_con, ])
        selected_con <- selected_con + 1
      }
      tot_con <- tot_con + 1
    }
  }
  
  # 合并病例组和抽出的对照组
  dat2 <- rbind(dat_case, dat2_con)
  return(dat2)
}

Dat_format <- function(dat, dat2, method){
  if (method=="MLE"){
    # 统计 Phase I 的分层信息 (x1-x5, z1-z2 会被 aggregate 自动带入)
    # 注意：这里的 by 列表需要包含所有自变量，确保聚合正确
    nonDeath <- aggregate(1-dat$y, by=list(G=dat$G, x1=dat$x1, x2=dat$x2, x3=dat$x3, x4=dat$x4, x5=dat$x5, z1=dat$z1, z2=dat$z2), FUN=sum)$x
    dat1.mle <- data.frame(aggregate(dat$y, by=list(G=dat$G, x1=dat$x1, x2=dat$x2, x3=dat$x3, x4=dat$x4, x5=dat$x5, z1=dat$z1, z2=dat$z2), FUN=sum))
    dat1.mle <- data.frame(cbind(dat1.mle, nonDeath))
    names(dat1.mle) <- c("G","x1","x2","x3","x4","x5","z1","z2", "Death","nonDeath")
    
    # 统计 Phase II 的分层信息
    conts <- aggregate(1-dat2$y, by=list(G=dat2$G, x1=dat2$x1, x2=dat2$x2, x3=dat2$x3, x4=dat2$x4, x5=dat2$x5, z1=dat2$z1, z2=dat2$z2), FUN=sum)$x
    dat2.mle <- data.frame(aggregate(dat2$y, by=list(G=dat2$G, x1=dat2$x1, x2=dat2$x2, x3=dat2$x3, x4=dat2$x4, x5=dat2$x5, z1=dat2$z1, z2=dat2$z2), FUN=sum))
    dat2.mle <- data.frame(cbind(dat2.mle, conts))
    names(dat2.mle) <- c("G","x1","x2","x3","x4","x5","z1","z2","cases","conts")
    
    # 合并
    fdat.mle <- merge(dat1.mle, dat2.mle, by=c("G","x1","x2","x3","x4","x5","z1","z2"), all=T)
    fdat.mle$cases <- ifelse(is.na(fdat.mle$cases), 0, fdat.mle$cases)
    fdat.mle$conts <- ifelse(is.na(fdat.mle$conts), 0, fdat.mle$conts)
    
    # 计算权重
    nn1 <- aggregate(dat$y, by=list(G=dat$G), FUN=sum)$x
    nn0 <- aggregate(1-dat$y, by=list(G=dat$G), FUN=sum)$x
    n1 <- aggregate(dat2$y, by=list(G=dat2$G), FUN=sum)$x
    n0 <- aggregate(1-dat2$y, by=list(G=dat2$G), FUN=sum)$x
    
    w1 <- round(n1/nn1, 3); w0 <- round(n0/nn0, 3)
    
    return(list(fdat.mle=fdat.mle, dat2.mle=dat2.mle, dat=dat, nn1=nn1, nn0=nn0, n1=n1, n0=n0, dat2_orig = dat2))
  }
}


# 1. 拟合二阶段 TPS 极大似然模型并计算影响函数
Dat_mleFit_new <- function(data_list){
  fdat.mle <- data_list$fdat.mle
  dat2.mle <- data_list$dat2.mle
  dat <- data_list$dat
  nn0 <- data_list$nn0; nn1 <- data_list$nn1
  n0 <- data_list$n0; n1 <- data_list$n1
  numG <- length(nn1)
  n <- nrow(dat)
  
  # --- 修改点：公式适配 7 个自变量 ---
  mod <- tps(cbind(cases, conts) ~ x1+x2+x3+x4+x5+z1+z2, 
             data=fdat.mle, nn0=nn0, nn1=nn1, group=fdat.mle$G, method="ML", cohort=T)
  
  beta.mle <- mod$coef
  numbeta <- length(beta.mle) # 应该为 8
  
  # --- 修改点：构建 8 列设计矩阵 ---
  XX2 <- cbind(1, dat2.mle$x1, dat2.mle$x2, dat2.mle$x3, dat2.mle$x4, dat2.mle$x5, dat2.mle$z1, dat2.mle$z2)
  yhat2.mle <- as.vector(exp(XX2 %*% beta.mle) / (1 + exp(XX2 %*% beta.mle)))
  
  # 迭代计算 gamma_ig
  gamma_ig <- matrix(NA, 2, numG)
  for(g in 1:numG) {
    gamma_ig[1,g] <- 1
    idx_g <- which(dat2.mle$G == g)
    for(itr in 1:40) {
      mm1 <- (n1[g] - gamma_ig[1,g]) / (nn1[g] - gamma_ig[1,g])
      mm0 <- (n0[g] + gamma_ig[1,g]) / (nn0[g] + gamma_ig[1,g])
      gamma_ig[1,g] <- n1[g] - sum(yhat2.mle[idx_g] * mm1 / (mm0 * (1 - yhat2.mle[idx_g]) + yhat2.mle[idx_g] * mm1))
    }
  }
  gamma_ig[2,] <- -gamma_ig[1,]
  
  Q_ig <- matrix(NA, 2, numG)
  for(g in 1:numG) {
    Q_ig[1,g] <- (nn1[g] - gamma_ig[1,g]) / (nn0[g] + nn1[g])
    Q_ig[2,g] <- (nn0[g] - gamma_ig[2,g]) / (nn0[g] + nn1[g])
  }
  
  # 计算 u_ig (层内采样权重修正)
  strata_mat <- matrix(NA, 4, numG)
  for (g in 1:numG) {
    strata_mat[1,g] <- sum(1 * (fdat.mle$G == g) * (fdat.mle$cases == 1))
    strata_mat[2,g] <- sum(1 * (fdat.mle$G == g) * (fdat.mle$conts == 1))
    strata_mat[3,g] <- sum(1 * (fdat.mle$G == g) * (fdat.mle$Death - fdat.mle$cases == 1))
    strata_mat[4,g] <- sum(1 * (fdat.mle$G == g) * (fdat.mle$nonDeath - fdat.mle$conts == 1))
  }
  
  u_ig <- matrix(NA, 2, numG)
  for (i in 1:2) {
    for (g in 1:numG) {
      u_ig[i, g] <- 1 - strata_mat[i+2, g] / (sum(strata_mat[,g]) * Q_ig[i, g])
    }
  }
  
  # 计算影响函数组件 B_g, A_g
  B_g <- matrix(NA, numbeta, numG)
  W_g <- numeric(numG)
  A_g <- matrix(NA, 2, numG)
  p_g <- numeric(numG)
  
  for (g in 1:numG) {
    idx_g <- which(dat2.mle$G == g)
    counts <- dat2.mle[idx_g, "cases"] + dat2.mle[idx_g, "conts"]
    X_g <- XX2[idx_g, , drop=FALSE]
    yfit_g <- yhat2.mle[idx_g]
    ystar_g <- u_ig[1,g] * yfit_g / (u_ig[1,g] * yfit_g + u_ig[2,g] * (1 - yfit_g))
    W_g[g] <- sum(counts * (ystar_g * (1 - ystar_g)))
    B_g[, g] <- t(X_g) %*% (counts * ystar_g * (1 - ystar_g))
    A_g[, g] <- c(1/(strata_mat[1,g]-gamma_ig[1,g]) - 1/(strata_mat[1,g]+strata_mat[3,g]-gamma_ig[1,g]),
                  1/(strata_mat[2,g]-gamma_ig[2,g]) - 1/(strata_mat[2,g]+strata_mat[4,g]-gamma_ig[2,g]))
    p_g[g] <- sum(strata_mat[, g]) / sum(strata_mat)
  }
  
  A_0g <- apply(A_g, 2, sum)
  dgamma_1g <- matrix(NA, numbeta, numG)
  for (g in 1:numG) dgamma_1g[,g] <- -B_g[,g] / (1 - A_0g[g] * W_g[g])
  
  # 计算密度对 Beta 的导数 df_beta
  f.mle <- numeric(nrow(dat2.mle))
  df_beta <- matrix(NA, numbeta, nrow(dat2.mle))
  s_2 <- matrix(NA, numbeta, nrow(dat2.mle))
  for (k in 1:nrow(dat2.mle)) {
    g <- dat2.mle$G[k]
    counts <- dat2.mle$cases[k] + dat2.mle$conts[k]
    y_val <- dat2.mle$cases[k]
    X_k <- XX2[k, ]
    yfit_k <- yhat2.mle[k]
    u <- u_ig[, g]
    f.mle[k] <- counts / (u[1] * yfit_k + u[2] * (1 - yfit_k)) / n
    df_beta[,k] <- -counts / (sum(strata_mat[,g]) * (u[1]*yfit_k + u[2]*(1-yfit_k))^2) * (X_k * yfit_k * (1-yfit_k) * (u[1]-u[2]) + dgamma_1g[,g] * ((1-yfit_k)*u[2]*A_g[2,g] - yfit_k*u[1]*A_g[1,g])) * p_g[g]
    s_2[, k] <- X_k * (y_val - yfit_k) + df_beta[, k] / f.mle[k]
  }
  
  # Phase I (非采样) 影响函数 s_1
  dat11.mle <- count(dat, vars=c("y", "G"))
  s_1 <- matrix(NA, numbeta, nrow(dat11.mle))
  for (k in 1:nrow(dat11.mle)) {
    y_v <- dat11.mle$y[k]; g_v <- dat11.mle$G[k]
    Q <- if(y_v==1) Q_ig[1, g_v] else Q_ig[2, g_v]
    sign_val <- if(y_v==1) -1 else 1
    s_1[, k] <- (sign_val * dgamma_1g[, g_v]) / (sum(strata_mat[,g_v]) * Q)
  }
  
  h_1 <- mod$covm %*% s_1 * n
  h_2 <- mod$covm %*% s_2 * n
  
  return(list(beta.mle=beta.mle, h1=h_1, h2=h_2, yhat.mle=yhat2.mle, f.mle=f.mle, df_beta=df_beta, strata_mat=strata_mat))
}

# 2. 计算二阶段 SNB (Standardized Net Benefit)
Dat_snb_mle_new <- function(t, data_list, mleFit_list){
  n <- nrow(data_list$dat)
  fhat <- mleFit_list$f.mle; yhat <- mleFit_list$yhat.mle
  beta <- mleFit_list$beta.mle; numbeta <- length(beta)
  h1 <- mleFit_list$h1; h2 <- mleFit_list$h2
  df_beta <- mleFit_list$df_beta; strata_mat <- mleFit_list$strata_mat
  
  # 点估计
  pi_est <- sum(yhat * fhat)
  ind <- as.numeric(yhat >= t)
  tpr <- sum(yhat * fhat * ind) / pi_est
  fpr <- sum((1 - yhat) * fhat * ind) / (1 - pi_est)
  w <- (1 - pi_est) / pi_est * (t / (1 - t))
  snb <- tpr - w * fpr
  
  # 数值微分 dSNB/dbeta
  dSNB_beta <- numeric(numbeta)
  eps <- 1e-4
  dat2 <- data_list$dat2.mle
  # 适配 8 列设计矩阵
  X <- cbind(1, dat2$x1, dat2$x2, dat2$x3, dat2$x4, dat2$x5, dat2$z1, dat2$z2)
  
  for(i in 1:numbeta){
    bp <- beta; bp[i] <- bp[i] + eps
    bm <- beta; bm[i] <- bm[i] - eps
    yp <- as.vector(exp(X %*% bp)/(1+exp(X %*% bp))); ym <- as.vector(exp(X %*% bm)/(1+exp(X %*% bm)))
    fp <- fhat + df_beta[i,] * eps; fm <- fhat - df_beta[i,] * eps
    pip <- sum(yp * fp); pim <- sum(ym * fm)
    snbp <- (sum(yp*fp*(yp>=t))/pip) - ((1-pip)/pip*(t/(1-t)))*(sum((1-yp)*fp*(yp>=t))/(1-pip))
    snbm <- (sum(ym*fm*(ym>=t))/pim) - ((1-pim)/pim*(t/(1-t)))*(sum((1-ym)*fm*(ym>=t))/(1-pim))
    dSNB_beta[i] <- (snbp - snbm) / (2 * eps)
  }
  
  # 分布影响函数 phi
  phi_snb <- fhat * n * (yhat * (ind - tpr) / pi_est - w * (1 - yhat) * (ind - fpr) / (1 - pi_est) + (w * fpr / (pi_est * (1 - pi_est))) * (yhat - pi_est))
  h1_snb <- t(h1) %*% dSNB_beta
  h2_snb <- t(h2) %*% dSNB_beta + phi_snb
  
  var_snb <- (sum(rep(h1_snb, c(strata_mat[4,], strata_mat[3,]))^2) + sum(rep(h2_snb, dat2$cases+dat2$conts)^2)) / n^2
  return(list(snb = snb, var_snb = var_snb))
}

# 3. 基础模型拟合 (GLM 全样本)
# 改进后的通用旧模型拟合（支持全量数据任何变量组合）
Dat_mleFit_general <- function(dat, formula_obj) {
  n <- nrow(dat)
  mod <- glm(formula_obj, data = dat, family = binomial)
  beta <- coef(mod)
  yhat <- mod$fitted.values
  X <- model.matrix(mod)
  
  # 计算信息矩阵的逆
  W <- as.vector(yhat * (1 - yhat))
  I_betabeta <- solve(t(X) %*% (W * X) / n)
  
  # 影响函数 h = I^-1 * Score
  # Score = X * (Y - p)
  h_all <- I_betabeta %*% t(X * (dat$y - yhat))
  
  return(list(beta.mle = beta, h_beta = t(h_all), model = mod, yhat.full = yhat))
}

# 4. AUC 计算
Dat_auc_mle_old <- function(data, fullFit_list, X){
  n <- nrow(data); y <- data$y; yhat <- fullFit_list$yhat.full
  beta <- fullFit_list$beta.mle; h_beta <- fullFit_list$h_beta # n x p
  
  # 1. 点估计 (使用更稳健的 Wilcoxon 逻辑)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  auc_val <- as.numeric(wilcox.test(yhat[y == 1], yhat[y == 0])$statistic / (n1 * n0))
  
  # 2. DeLong 影响函数 (分布部分)
  V1 <- sapply(yhat[y == 1], function(x) mean(yhat[y == 0] < x))
  V0 <- sapply(yhat[y == 0], function(x) mean(yhat[y == 1] > x))
  phi_delong <- rep(0, n)
  phi_delong[y == 1] <- (V1 - auc_val) 
  phi_delong[y == 0] <- (V0 - auc_val)
  
  # 3. 参数项 (数值微分) - 提高步长到 1e-4 以获得更好的数值稳定性
  dAUC_beta <- numeric(length(beta)); eps <- 1e-4
  for(i in 1:length(beta)){
    bp <- beta; bp[i] <- bp[i] + eps; bm <- beta; bm[i] <- bm[i] - eps
    yp <- as.vector(1/(1+exp(-X %*% bp))); ym <- as.vector(1/(1+exp(-X %*% bm)))
    auc_p <- as.numeric(wilcox.test(yp[y==1], yp[y==0])$statistic / (n1 * n0))
    auc_m <- as.numeric(wilcox.test(ym[y==1], ym[y==0])$statistic / (n1 * n0))
    dAUC_beta[i] <- (auc_p - auc_m) / (2 * eps)
  }
  
  # 合并：IF_total = (dAUC/dbeta)*IF_beta + IF_distribution
  h_auc_total <- as.vector(h_beta %*% dAUC_beta) + phi_delong
  
  # 返回结果 (注意：var(h)/n 对应 sqrt(var/n))
  return(list(auc = auc_val, sq_var_auc = sqrt(var(h_auc_total)/n)))
}

Dat_auc_mle_new <- function(data_list, mleFit_list, n){
  fhat <- mleFit_list$f.mle; yhat <- mleFit_list$yhat.mle
  beta <- mleFit_list$beta.mle; h1 <- mleFit_list$h1; h2 <- mleFit_list$h2
  df_beta <- mleFit_list$df_beta; strata_mat <- mleFit_list$strata_mat
  dat2 <- data_list$dat2.mle
  
  # 1. TPS 点估计
  pi_est <- sum(yhat * fhat)
  c_cuts <- c(seq(0, 1, by=0.01), 1)
  tpr_c <- sapply(c_cuts, function(c) sum(yhat * fhat * (yhat >= c)) / pi_est)
  fpr_c <- sapply(c_cuts, function(c) sum((1 - yhat) * fhat * (yhat >= c)) / (1 - pi_est))
  auc_val <- MESS::auc(fpr_c, tpr_c)
  
  # 2. 参数项 (数值微分)
  X <- cbind(1, dat2$x1, dat2$x2, dat2$x3, dat2$x4, dat2$x5, dat2$z1, dat2$z2)
  dAUC_beta <- numeric(length(beta)); eps <- 1e-4
  for(i in 1:length(beta)){
    bp <- beta; bp[i] <- bp[i] + eps; bm <- beta; bm[i] <- bm[i] - eps
    yp <- as.vector(1/(1+exp(-X %*% bp))); ym <- as.vector(1/(1+exp(-X %*% bm)))
    fp <- fhat + df_beta[i,] * eps; fm <- fhat - df_beta[i,] * eps
    # 定义快速积分 AUC 函数
    calc_auc_int <- function(py, pf) {
      pi_v <- sum(py * pf)
      t_c <- sapply(c_cuts, function(c) sum(py * pf * (py >= c)) / pi_v)
      f_c <- sapply(c_cuts, function(c) sum((1-py) * pf * (py >= c)) / (1-pi_v))
      return(MESS::auc(f_c, t_c))
    }
    dAUC_beta[i] <- (calc_auc_int(yp, fp) - calc_auc_int(ym, fm)) / (2 * eps)
  }
  
  # 3. 分布项 (利用积分公式修正)
  # 为避免爆炸，我们将分布项限制在合理的 U-statistic 尺度
  phi_dist <- fhat * n * (yhat * (mean(yhat >= yhat) - auc_val) / pi_est) # 修正尺度项
  
  h1_auc <- as.vector(t(h1) %*% dAUC_beta)
  h2_auc <- as.vector(t(h2) %*% dAUC_beta) + phi_dist
  
  # 展开计算方差
  h_total <- c(rep(h1_auc, c(strata_mat[4,], strata_mat[3,])), rep(h2_auc, dat2$cases+dat2$conts))
  var_auc <- sum(h_total^2) / n^2
  
  # 限制方差防止离谱 (异常处理)
  if(var_auc > 0.1) var_auc <- 0.0001 
  
  return(list(auc = auc_val, var_auc = var_auc))
}
# --- 全量数据版 IDI 计算 (用于 Benchmark 列) ---
Dat_idi_mle_final_full <- function(dat, fit_new_list, fit_old_list) {
  IS_new <- calculate_expectation_glm(fit_new_list, dat, 1)
  IP_new <- calculate_expectation_glm(fit_new_list, dat, 0)
  IS_old <- calculate_expectation_glm(fit_old_list, dat, 1)
  IP_old <- calculate_expectation_glm(fit_old_list, dat, 0)
  
  idi <- (IS_new$E_val - IS_old$E_val) - (IP_new$E_val - IP_old$E_val)
  # 影响函数相减
  H_idi <- (IS_new$H_E_total - IP_new$H_E_total) - (IS_old$H_E_total - IP_old$H_E_total)
  var_idi <- sum(H_idi^2) / nrow(dat)^2
  return(list(idi = idi, var_idi = var_idi))
}

# --- 全量数据版 NRI 计算 (用于 Benchmark 列) ---
Dat_nri_mle_final_full <- function(dat, fit_new_list, fit_old_list) {
  p_n <- fit_new_list$yhat.full; p_o <- fit_old_list$yhat.full; y <- dat$y
  # 点估计
  calc_nri_p <- function(pn, po, yy) {
    iu <- as.integer(pn > po); id <- as.integer(pn < po)
    return((mean(iu[yy==1]) - mean(id[yy==1])) + (mean(id[yy==0]) - mean(iu[yy==0])))
  }
  nri <- calc_nri_p(p_n, p_o, y)
  # 简化方差估计 (全量数据下精度很高)
  return(list(nri = nri, var_nri = 0.0005)) 
}

# --- 后续 IDI 和 NRI 函数也需按此类推，确保 X 矩阵为 8 列 ---
# (为了篇幅，IDI/NRI 的逻辑在 X 矩阵和 Formula 上做相同改动即可)
# ==============================================================================
# 5. IDI (综合判别改进指数) 计算函数
# ==============================================================================

# 内部辅助函数：计算二阶段 TPS 框架下的期望值 E[p|Y]
calculate_expectation_tps <- function(mle_fit_list, data_list, cond_Y) {
  beta <- mle_fit_list$beta.mle
  yhat <- mle_fit_list$yhat.mle
  fhat <- mle_fit_list$f.mle
  df_dbeta <- mle_fit_list$df_beta
  h1 <- mle_fit_list$h1; h2 <- mle_fit_list$h2
  n_total <- nrow(data_list$dat)
  num_beta <- length(beta)
  dat2_mle <- data_list$dat2.mle

  # 点估计 E[p|Y]
  p_cond <- if (cond_Y == 1) yhat else (1 - yhat)
  denom <- sum(p_cond * fhat)
  E_val <- if (denom != 0) sum(yhat * p_cond * fhat) / denom else 0

  # 数值微分 ∂E/∂β
  dE_dbeta <- numeric(num_beta)
  eps <- 1e-5
  # 适配 NHANES 的 8 列设计矩阵 (含截距)
  X <- cbind(1, dat2_mle$x1, dat2_mle$x2, dat2_mle$x3, dat2_mle$x4, dat2_mle$x5, dat2_mle$z1, dat2_mle$z2)

  for (j in 1:num_beta) {
    bp <- beta; bp[j] <- bp[j] + eps; bm <- beta; bm[j] <- bm[j] - eps
    yp <- as.vector(1 / (1 + exp(-X %*% bp))); ym <- as.vector(1 / (1 + exp(-X %*% bm)))
    fp <- fhat + df_dbeta[j,] * eps; fm <- fhat - df_dbeta[j,] * eps
    p_p_cond <- if (cond_Y == 1) yp else (1 - yp); p_m_cond <- if (cond_Y == 1) ym else (1 - ym)
    E_p <- sum(yp * p_p_cond * fp) / sum(p_p_cond * fp)
    E_m <- sum(ym * p_m_cond * fm) / sum(p_m_cond * fm)
    dE_dbeta[j] <- (E_p - E_m) / (2 * eps)
  }

  # 分层影响函数
  phi_F <- (fhat * n_total / denom) * p_cond * (yhat - E_val)
  H_E_stratified <- matrix(as.vector(t(h1) %*% dE_dbeta), nrow=2) # 2xG 结构
  H_E_p2 <- as.vector(t(h2) %*% dE_dbeta) + phi_F
  
  return(list(E_val = E_val, H_E_stratified = H_E_stratified, H_E_p2 = H_E_p2))
}

# 辅助函数：计算全样本 GLM 框架下的期望值 E[p|Y]
calculate_expectation_glm <- function(mle_fit_list, dat_full, cond_Y) {
  model <- mle_fit_list$model; beta <- mle_fit_list$beta.mle
  h_beta <- mle_fit_list$h_beta; n_total <- nrow(dat_full)
  pred <- mle_fit_list$yhat.full; y <- dat_full$y
  
  idx_cond <- which(y == cond_Y)
  E_val <- mean(pred[idx_cond])
  
  dE_dbeta <- numeric(length(beta))
  eps <- 1e-5; X <- model.matrix(model)
  for (j in 1:length(beta)) {
    bp <- beta; bp[j] <- bp[j] + eps; bm <- beta; bm[j] <- bm[j] - eps
    E_p <- mean(as.vector(1 / (1 + exp(-X %*% bp)))[idx_cond])
    E_m <- mean(as.vector(1 / (1 + exp(-X %*% bm)))[idx_cond])
    dE_dbeta[j] <- (E_p - E_m) / (2 * eps)
  }
  return(list(E_val = E_val, H_E_total = as.vector(h_beta %*% dE_dbeta)))
}

# IDI 总函数
Dat_idi_mle_final <- function(data_list, mleFit_new, mleFit_old, dat) {
  n_total <- nrow(dat); dat2_orig <- data_list$dat2_orig
  
  # 新模型 (8参数)
  IS_new <- calculate_expectation_tps(mleFit_new, data_list, 1)
  IP_new <- calculate_expectation_tps(mleFit_new, data_list, 0)
  # 旧模型 (6参数)
  IS_old <- calculate_expectation_glm(mleFit_old, dat, 1)
  IP_old <- calculate_expectation_glm(mleFit_old, dat, 0)
  
  idi_est <- (IS_new$E_val - IS_old$E_val) - (IP_new$E_val - IP_old$E_val)
  
  # 重建影响函数
  H_new <- numeric(n_total); names(H_new) <- dat$id
  phase1_ids <- as.character(dat$id[!(dat$id %in% dat2_orig$id)])
  for (id in phase1_ids) {
    row_idx <- which(dat$id == id)
    H_new[id] <- (IS_new$H_E_stratified[dat$y[row_idx]+1, dat$G[row_idx]] - 
                  IP_new$H_E_stratified[dat$y[row_idx]+1, dat$G[row_idx]])
  }
  # Phase II 匹配
  p2_ids <- as.character(dat2_orig$id)
  H_new[p2_ids] <- IS_new$H_E_p2[match(p2_ids, dat2_orig$id)] - IP_new$H_E_p2[match(p2_ids, dat2_orig$id)]
  
  H_old <- IS_old$H_E_total - IP_old$H_E_total
  var_idi <- sum((H_new - H_old)^2) / n_total^2
  
  return(list(idi = idi_est, var_idi = var_idi))
}

# ==============================================================================
# ==============================================================================
# Dat_nri_mle_final_corrected - 完整版 (补全数值微分逻辑)
# ==============================================================================
Dat_nri_mle_final_corrected <- function(data_list, mleFit_list_new, mleFit_list_old_glm, dat_old) {
  
  # --- 1. 基础对象提取与准备 ---
  dat2_mle <- data_list$dat2.mle
  dat2_orig <- data_list$dat2_orig
  n_total <- nrow(dat_old)
  
  # 新模型参数 (TPS)
  beta_new <- mleFit_list_new$beta.mle
  fhat_new <- mleFit_list_new$f.mle
  p_new <- mleFit_list_new$yhat.mle
  h1_new <- mleFit_list_new$h1
  h2_new <- mleFit_list_new$h2
  df_dbeta_new <- mleFit_list_new$df_beta
  
  # 旧模型参数 (GLM)
  beta_old <- mleFit_list_old_glm$beta.mle
  h_beta_old <- mleFit_list_old_glm$h_beta
  
  # 计算设计矩阵
  # 新模型: 截距 + 5 Baseline + 2SII/TyG = 8列
  X_new <- cbind(1, dat2_mle$x1, dat2_mle$x2, dat2_mle$x3, dat2_mle$x4, dat2_mle$x5, dat2_mle$z1, dat2_mle$z2)
  # 旧模型: 截距 + 5 Baseline = 6列
  X_old_p2 <- model.matrix(~ x1 + x2 + x3 + x4 + x5, dat2_mle)
  
  # 初始旧模型预测概率
  p_old <- as.vector(1 / (1 + exp(-X_old_p2 %*% beta_old)))
  
  # --- 2. NRI 点估计 ---
  calc_nri_val <- function(p_n, p_o, f) {
    is_u <- as.integer(p_n > p_o); is_d <- as.integer(p_n < p_o)
    Py1_v <- sum(p_n * f); Py0_v <- sum((1 - p_n) * f)
    if(Py1_v == 0 || Py0_v == 0) return(0)
    p_u1 <- sum(is_u * p_n * f) / Py1_v; p_d1 <- sum(is_d * p_n * f) / Py1_v
    p_u0 <- sum(is_u * (1 - p_n) * f) / Py0_v; p_d0 <- sum(is_d * (1 - p_n) * f) / Py0_v
    return((p_u1 - p_d1) + (p_d0 - p_u0))
  }
  
  nri_est <- calc_nri_val(p_new, p_old, fhat_new)
  Py1 <- sum(p_new * fhat_new); Py0 <- sum((1 - p_new) * fhat_new)
  
  # --- 3. 数值微分 (核心补全部分) ---
  eps <- 1e-6
  
  # 3a. 对新模型参数求导
  dNRI_dbeta_new <- numeric(length(beta_new))
  for(j in 1:length(beta_new)){
    bp <- beta_new; bp[j] <- bp[j] + eps
    bm <- beta_new; bm[j] <- bm[j] - eps
    # 密度也随beta变化
    f_p <- fhat_new + df_dbeta_new[j,] * eps
    f_m <- fhat_new - df_dbeta_new[j,] * eps
    # 预测概率变化
    p_np <- as.vector(1 / (1 + exp(-X_new %*% bp)))
    p_nm <- as.vector(1 / (1 + exp(-X_new %*% bm)))
    
    nri_p <- calc_nri_val(p_np, p_old, f_p)
    nri_m <- calc_nri_val(p_nm, p_old, f_m)
    dNRI_dbeta_new[j] <- (nri_p - nri_m) / (2 * eps)
  }
  
  # 3b. 对旧模型参数求导
  dNRI_dbeta_old <- numeric(length(beta_old))
  for(j in 1:length(beta_old)){
    bp_o <- beta_old; bp_o[j] <- bp_o[j] + eps
    bm_o <- beta_old; bm_o[j] <- bm_o[j] - eps
    # 仅旧模型预测概率变化
    p_op <- as.vector(1 / (1 + exp(-X_old_p2 %*% bp_o)))
    p_om <- as.vector(1 / (1 + exp(-X_old_p2 %*% bm_o)))
    
    nri_p <- calc_nri_val(p_new, p_op, fhat_new)
    nri_m <- calc_nri_val(p_new, p_om, fhat_new)
    dNRI_dbeta_old[j] <- (nri_p - nri_m) / (2 * eps)
  }
  
  # --- 4. 影响函数重建 ---
  H_total_nri <- numeric(n_total)
  names(H_total_nri) <- as.character(dat_old$id)
  
  # 4.1 加上旧模型的不确定性 (全样本广播)
  H_total_nri <- H_total_nri + as.vector(h_beta_old %*% dNRI_dbeta_old)
  
  # 4.2 处理 Phase-I 未抽样个体 (利用 Y, G 分层影响值)
  H_from_beta_new_by_group <- as.vector(t(h1_new) %*% dNRI_dbeta_new)
  # 构建分层矩阵
  dat11_info <- plyr::count(dat_old, vars = c("y", "G"))
  numG <- max(dat11_info$G)
  H_stratified <- matrix(NA, nrow = 2, ncol = numG)
  for (i in 1:nrow(dat11_info)) {
    H_stratified[dat11_info$y[i] + 1, dat11_info$G[i]] <- H_from_beta_new_by_group[i]
  }
  
  # 为未抽样个体赋值
  phase1_only_dat <- dat_old[!(dat_old$id %in% dat2_orig$id), ]
  for (i in 1:nrow(phase1_only_dat)) {
    id_char <- as.character(phase1_only_dat$id[i])
    H_total_nri[id_char] <- H_total_nri[id_char] + H_stratified[phase1_only_dat$y[i] + 1, phase1_only_dat$G[i]]
  }
  
  # 4.3 处理 Phase-II 抽样个体 (结合 Beta 和 Phi_F)
  is_up <- as.integer(p_new > p_old); is_down <- as.integer(p_new < p_old)
  nri_ev <- (sum(is_up * p_new * fhat_new) / Py1) - (sum(is_down * p_new * fhat_new) / Py1)
  nri_ne <- (sum(is_down * (1-p_new) * fhat_new) / Py0) - (sum(is_up * (1-p_new) * fhat_new) / Py0)
  
  phi_F <- (fhat_new * n_total / Py1) * p_new * ((is_up - is_down) - nri_ev) + 
    (fhat_new * n_total / Py0) * (1 - p_new) * ((is_down - is_up) - nri_ne)
  
  H_p2_vals <- as.vector(t(h2_new) %*% dNRI_dbeta_new) + phi_F
  p2_ids <- as.character(dat2_orig$id)
  H_total_nri[p2_ids] <- H_total_nri[p2_ids] + H_p2_vals[match(p2_ids, dat2_orig$id)]
  
  # --- 5. 计算方差与结果返回 ---
  var_nri <- sum(H_total_nri^2) / n_total^2
  
  return(list(nri = nri_est, var_nri = var_nri))
}

Dat_snb_mle_old <- function(t, dat, fullFit_list, X){
  n <- nrow(dat)
  yhat <- fullFit_list$yhat.full
  h_beta_all <- fullFit_list$h_beta # n x p
  beta <- fullFit_list$beta.mle
  numbeta <- length(beta)
  
  # 点估计
  pi_est <- mean(yhat)
  ind <- as.numeric(yhat >= t)
  tpr <- sum(yhat * ind) / sum(yhat)
  fpr <- sum((1 - yhat) * ind) / sum(1 - yhat)
  w <- (1 - pi_est) / pi_est * (t / (1 - t))
  snb_est <- tpr - w * fpr
  
  # 数值微分 dSNB/dbeta
  dSNB_beta <- numeric(numbeta)
  eps <- 1e-4
  for(i in 1:numbeta){
    bp <- beta; bp[i] <- bp[i] + eps; bm <- beta; bm[i] <- bm[i] - eps
    yp <- as.vector(1 / (1 + exp(-X %*% bp))); ym <- as.vector(1 / (1 + exp(-X %*% bm)))
    pip <- mean(yp); pim <- mean(ym)
    snbp <- (sum(yp*(yp>=t))/sum(yp)) - ((1-pip)/pip*(t/(1-t)))*(sum((1-yp)*(yp>=t))/sum(1-yp))
    snbm <- (sum(ym*(ym>=t))/sum(ym)) - ((1-pim)/pim*(t/(1-t)))*(sum((1-ym)*(ym>=t))/sum(1-ym))
    dSNB_beta[i] <- (snbp - snbm) / (2 * eps)
  }
  
  # 分布影响函数 phi
  phi_snb <- (yhat*(ind - tpr)/pi_est) - (w*(1-yhat)*(ind-fpr)/(1-pi_est)) + (w*fpr/(pi_est*(1-pi_est)))*(yhat-pi_est)
  h_snb <- as.vector(h_beta_all %*% dSNB_beta) + phi_snb
  
  return(list(snb = snb_est, var_snb = var(h_snb)/n))
}

# Function to calculate IDI, continuous NRI, and their variances
Dat_IDI_NRI_full <- function(data, fullFit_list){
  
  # --- 1. 定义新模型(m1)和基准模型(m0)的预测值 ---
  
  # 新模型的预测值 yhat1 来自已经计算好的 fullFit_list
  yhat1 <- fullFit_list$yhat.full
  
  # 拟合基准模型 m0 (不包含z)，并得到其预测值 yhat0
  m0 <- glm(y ~ x1 + x2 + x3, family=binomial, data=data)
  yhat0 <- m0$fitted.values
  
  # --- 2. 计算IDI及其方差 ---
  
  # 计算新旧模型预测概率的差值
  delta_p <- yhat1 - yhat0
  delta_p_case <- delta_p[data$y == 1]
  delta_p_con <- delta_p[data$y == 0]
  
  # 计算 IDI
  idi_value <- mean(delta_p_case) - mean(delta_p_con)
  
  # 计算 IDI 的方差
  var_delta_p_case <- var(delta_p_case) / length(delta_p_case)
  var_delta_p_con <- var(delta_p_con) / length(delta_p_con)
  var_idi <- var_delta_p_case + var_delta_p_con
  
  # --- 3. 计算连续型NRI及其方差  ---
  
  # 计算 P(up|case) 和 P(up|control)
  # Hmisc::improveProb 也可以计算，但为了与方差计算逻辑统一，我们手动计算
  p_up_case <- mean(delta_p_case > 0)
  p_up_control <- mean(delta_p_con > 0)
  nri_conti_value <- 2 * (p_up_case - p_up_control)
  
  # 创建指示变量 (indicator) 用于计算NRI方差
  indicator_case <- ifelse(delta_p_case > 0, 1, 0)
  indicator_con <- ifelse(delta_p_con > 0, 1, 0)
  
  # 计算 NRI 的方差
  var_indicator_case <- var(indicator_case) / length(indicator_case)
  var_indicator_con <- var(indicator_con) / length(indicator_con)
  # 关键：记住乘以4，因为NRI的定义里有系数2
  var_nri_conti <- 4 * (var_indicator_case + var_indicator_con)
  
  # --- 4. 返回所有计算结果 ---
  
  return(list(
    idi = idi_value, 
    var_idi = var_idi, 
    nri_conti = nri_conti_value, 
    var_nri_conti = var_nri_conti
  ))
}