rm(list = ls())
# scenario B
ref <- "B002"
seed <- 2023
if(dir.exists(ref)) {
  cat(paste0("The directory ", ref, " already exists."))
} else {
  dir.create(ref)
}
file.copy("simB_code.R", paste0(ref, "/sim_code_", ref, ".R"))
rep <- 100 # replicas
set.seed(seed); it_s <- sample(seq_len(rep), 1)
est <- cvg <- list(g = matrix(NA, nrow = rep, ncol = 4), # estimates and coverage matrices
                   b = matrix(NA, nrow = rep, ncol = 4)) # [g]aussian and [b]eta models
colnames(est$g) <- colnames(est$b) <- c("Intercept", "Time", "Group", "Assoc")
colnames(cvg$g) <- colnames(cvg$b) <- c("Intercept", "Time", "Group", "Assoc")
summ <- list("prop_tevent" = rep(NA, rep), # proportion of terminal events
             "prop_tcensor" = rep(NA, rep), # proportion of censoring
             "N" = rep(NA, rep), # total number of longitudinal obs
             "med_ni" = rep(NA, rep), # median number of longitudinal obs per individual
             "afp" = rep(NA, rep), # aggregated follow-up
             "med_fpi" = rep(NA, rep), # median follow-up duration per individual
             "med_ttime" = rep(NA, rep), # median terminal event time
             "med_tcensor" = rep(NA, rep), # median censoring time
             "prop_group1" = rep(NA, rep)) # proportion of level 1 in categorical var
n_iter <- 5000
n_burnin <- 2500
n <- 1000 # number of individuals
n_i <- 20 # number of measurements per individual
t_max <- 10 # maximum follow-up time
betas <- c("Intercept" = 2, "Time" = -1)
D <- diag(c(0.25, 0.15)) # RE vcov matrix
gammas_t <- c("Group" = 0.25)
alpha_t <- -2
h_0 <- 0.1
true <- c(betas, gammas_t, alpha_t)
expit <- function(x) exp(x) / (1 + exp(x))
invS_t <- function(t, i, eta_t, u_t) {
  h <- function(s) {
    X_s <- cbind(1, s)
    Z_s <- cbind(1, s)
    eta_y <- as.vector(X_s %*% betas + rowSums(Z_s * b[rep(i, nrow(Z_s)), ]))
    exp(log(h_0) + eta_t + expit(eta_y) * alpha_t)
  }
  integrate(h, lower = 0, upper = t)$value + log(u_t)
}
library(JMbayes2)
tic <- Sys.time(); print(tic)
for (it in seq_len(rep)){
  set.seed(seed + it)
  cat(paste0("\r", it, "/", rep, "   "))
  ## longitudinal 1/2 ==========================================================
  b <- MASS::mvrnorm(n, rep(0, nrow(D)), D)
  ## terminal ==================================================================
  group <- rbinom(n, 1, 0.5)
  summ$prop_group1[it] <- mean(group)
  term <- data.frame(id = seq_len(n), group = group)
  W <- model.matrix(~ -1 + group, data = term)
  eta_t <- as.vector(W %*% gammas_t)
  u_t <- runif(n)
  tstop <- numeric(n)
  for(i in seq_len(n)) {
    root <- try(uniroot(invS_t, interval = c(1e-05, 100), i = i, eta_t[i], u_t[i])$root, TRUE)  
    tstop[i] <- if (!inherits(root, "try-error")) root else t_max
  }
  term$tstop <- pmin(t_max, tstop)
  term$status <- as.numeric(tstop < t_max)
  summ$prop_tevent[it] <- mean(term$status)
  summ$prop_tcensor[it] <- 1 - mean(term$status)
  summ$med_ttime[it] <- median(term$tstop[term$status == 1])
  summ$med_tcensor[it] <- median(term$tstop[term$status == 0])
  if(it != rep) remove(eta_t, group, i, root, u_t, W)
  ## longitudinal 2/2 ==========================================================
  time <- c(replicate(n, c(0, sort(runif(n_i - 1, 0, t_max)))))
  long <- data.frame(id   = rep(seq_len(n), each = n_i),
                     time = time)
  X <- model.matrix(~ 1 + time, data = long) # FE design matrix
  Z <- model.matrix(~ 1 + time, data = long) # RE design matrix
  eta_y <- as.vector(X %*% betas + rowSums(Z * b[long$id, ])) # linear predictor
  mu <- expit(eta_y) # expected value
  phi <- 10000 # Beta precision parameter
  shape1 <- phi * mu
  shape2 <- phi * (1 - mu) 
  y <- rbeta(length(shape1), shape1, shape2)
  y <- (y * (length(y) - 1) + 0.5) / length(y)
  long$y <- y
  long <- long[long$time <= term$tstop[long$id], ]
  summ$N[it] <- length(long$y)
  summ$med_ni[it] <- median(table(long$id))
  head_rows <- tapply(seq_along(long$id), long$id, head, 1)
  tail_rows <- tapply(seq_along(long$id), long$id, tail, 1)
  summ$afp[it] <- sum(long$time[tail_rows] - long$time[head_rows])
  summ$med_fpi[it] <- median(long$time[tail_rows] - long$time[head_rows])
  if(it != rep) remove(b, eta_y, head_rows, mu, shape1, shape2, tail_rows, time, 
                       X, y, Z)
  ## beta joint model ==========================================================
  m_b <- mixed_model(y ~ time, random = ~ time || id, data = long,
                     family = beta.fam())
  ph <- coxph(Surv(tstop, status) ~ group, data = term)
  jm_b <- jm(ph, m_b, time_var = "time", functional_forms = ~ vexpit(value(y)),
             n_iter = n_iter, n_burnin = n_burnin)
  if(it == it_s) {
    betas_b <- jm_b$statistics$Mean$betas1
    saveRDS(jm_b, paste0(ref, "/jm_b_", ref, ".rds"))
    pdf(paste0(ref, "/tcp_b_", ref, ".pdf"), width = 12, height = 16)
    {
      par(mfrow = c(4, 3))
      traceplot(jm_b, "all")
      mtext(paste0("Dataset ", it_s, "/", rep), side = 3 , line = -2, 
            outer = TRUE)
    }
    dev.off()
  }
  est$b[it, ] <- c(fixef(jm_b)[[1]], coef(jm_b)$gammas, coef(jm_b)$association)
  low_b <- unlist(jm_b$statistics$CI_low[c("betas1", "gammas", "alphas")])
  upp_b <- unlist(jm_b$statistics$CI_upp[c("betas1", "gammas", "alphas")])
  cvg$b[it, ] <- true <= upp_b & true >= low_b
  if(it != rep) remove(low_b, m_b, jm_b, upp_b)
  ## Gaussian joint model ======================================================
  m_g <- lme(y ~ time, random = list(id = pdDiag(~time)), data = long, 
             control = list(opt = "optim"))
  jm_g <- jm(ph, m_g, time_var = "time", n_iter = n_iter, n_burnin = n_burnin)
  if(it == it_s) {
    betas_g <- jm_g$statistics$Mean$betas1
    saveRDS(jm_g, paste0(ref, "/jm_g_", ref, ".rds"))
    pdf(paste0(ref, "/tcp_g_", ref, ".pdf"), width = 12, height = 16)
    {
      par(mfrow = c(4, 3))
      traceplot(jm_g, "all")
      mtext(paste0("Dataset ", it_s, "/", rep), side = 3 , line = -2, 
            outer = TRUE)
    }
    dev.off()
  }
  est$g[it, ] <- c(fixef(jm_g)[[1]], coef(jm_g)$gammas, coef(jm_g)$association)
  low_g <- unlist(jm_g$statistics$CI_low[c("betas1", "gammas", "alphas")])
  upp_g <- unlist(jm_g$statistics$CI_upp[c("betas1", "gammas", "alphas")])
  cvg$g[it, ] <- true <= upp_g & true >= low_g
  if(it != rep) remove(jm_g, low_g, m_g, ph, term, upp_g)
  ## plot ======================================================================
  if(it == it_s){
    rx <- 100 # resolution
    toy <- data.frame(time = seq(0, t_max, length.out = rx))
    X_toy <- model.matrix(~ 1 + time, data = toy)
    toy$ytrue <- expit(X_toy %*% betas)
    toy$yhat_b <- expit(X_toy %*% betas_b)
    toy$yhat_g <- X_toy %*% betas_g
    pdf(paste0(ref, "/long_prof_", ref, ".pdf"), width = 15, height = 5)
    {
      par(mfrow = c(1, 2))
      ### 1/2
      plot(NA, ylim = range(long$y), xlim = c(0, t_max), 
           xlab = "Time", ylab = "Y1", xaxs = "i", yaxs = "i")
      invisible(tapply(seq_along(long$id), long$id, function(rows) {
        lines(long$time[rows], long$y[rows], col = rgb(0, 0, 0, 0.05))
      }))
      lines(toy$time, toy$ytrue, lwd = 2, col = "green") # true model
      lines(toy$time, toy$yhat_g, lwd = 2, lty = 2, col = "yellow") # estimated
      lines(toy$time, toy$yhat_b, lwd = 2, lty = 2, col = "yellow") # estimated
      mtext(text = paste0("Dataset ", it_s, "/", rep), side = 3, line = 0)
      ### 2/2
      boxplot(tstop, horizontal = TRUE, ylab = "", xlab = "Survival time")
      abline(v = t_max, col = 2)
    }
    dev.off()
    remove(rx, toy, X_toy)
  }
  if(it != rep) remove(long, tstop)
}
toc <- Sys.time(); print(toc)
diff <- difftime(toc, tic); print(diff)
cvg <- lapply(cvg, colMeans)
bias <- lapply(est, function(x) colMeans(x) - true) 
mse <- lapply(est, function(x) colMeans(sweep(x, 2, true)^2))
res <- list(cvg = cvg, est = est, bias = bias, mse = mse, true = true, 
            time = diff, summ = summ)
saveRDS(res, file = paste0(ref, "/res_", ref,".rds"))

pdf(paste0(ref, "/res_", ref, ".pdf"), width = 10, height = 5) 
{
  par(mfrow = c(1, 4))
  for(i in seq_len(4)){
    par(font.main = 1)
    boxplot(c(est$b[, i], est$g[, i]) ~ rep(c("beta", "Gaussian"), each = rep), 
            main = colnames(est$b)[i], xlab = "", ylab = expression(hat(theta)))
    abline(h = true[i], col = 2, lwd = 2)
    label1 <- c("Bias:", round(c(bias$b[i], bias$g[i]), 2))
    mtext(label1, side = 1, at = 0:2, line = 3, cex = 0.7, adj = c(0, 0.5, 0.5))
    label2 <- c("CP (%):", round(c(cvg$b[i], cvg$g[i]) * 100, 2))
    mtext(label2, side = 1, at = 0:2, line = 4, cex = 0.7, adj = c(0, 0.5, 0.5))
  }
} 
dev.off()