rm(list = ls())
# scenario C
ref <- "A001"
seed <- 2023
if(dir.exists(ref)) {
  cat(paste0("The directory ", ref, " already exists."))
} else {
  dir.create(ref)
}
file.copy("simA_code.R", paste0(ref, "/sim_code_", ref, ".R"))
rep <- 100 # replicas
set.seed(seed); it_s <- sample(seq_len(rep), 1)
est <- cvg <- matrix(NA, nrow = rep, ncol = 15) # estimates and coverage matrices
colnames(est) <- colnames(est) <- c("Intercept1", "Time1", 
                                    "Intercept2", "Time2", 
                                    "Group1", "Group2", "Group3",
                                    "Assoc11", "Assoc12", # strata 1
                                    "Assoc21", "Assoc22", "Assoc23",  # strata 2
                                    "Assoc31", "Assoc32", "Assoc33")  # strata 3
summ <- list("prop_tevent1" = rep(NA, rep), # proportion of competing risk 1
             "prop_tevent2" = rep(NA, rep), # proportion of competing risk 2
             "prop_tcensor" = rep(NA, rep), # proportion of censoring
             "med_revents" = rep(NA, rep), # median number of recurrent events per individual
             "N1" = rep(NA, rep), # total number of long 1 obs
             "N2" = rep(NA, rep), # total number of long 2 obs
             "med_ni1" = rep(NA, rep), # median number of long 1 obs per individual
             "med_ni2" = rep(NA, rep), # median number of long 2 obs per individual
             "afp1" = rep(NA, rep), # aggregated follow-up long 1
             "afp2" = rep(NA, rep), # aggregated follow-up long 2
             "med_fpi1" = rep(NA, rep), # median follow-up duration per individual long 1
             "med_fpi2" = rep(NA, rep), # median follow-up duration per individual long 2
             "med_ttime1" = rep(NA, rep), # median competing risk 1 time
             "med_ttime2" = rep(NA, rep), # median competing risk 2 time
             "med_tcensor" = rep(NA, rep), # median censoring time
             "prop_group1" = rep(NA, rep)) # proportion of level 1 in categorical var
n_iter <- 10000
n_burnin <- 7500
n <- 1000 # number of individuals
n_i <- 20 # number of measurements per individual
t_max <- 10 # maximum follow-up time
betas1 <- c("Intercept" = 2,   "Time" = -1.5) # beta model
betas2 <- c("Intercept" = 0.8, "Time" = -0.05) # Gauss model
D1 <- diag(c(0.25, 0.15)) # RE vcov matrix
D2 <- diag(c(0.01, 0.01)) # RE vcov matrix #?? c(0.001, 0.00001)
gammas_t1 <- c("Group" = 0.25)
gammas_t2 <- c("Group" = 0.25)
gammas_r  <- c("Group" = 0.25)
alpha_t1 <- c("Long1" = -2, "Long2" = -1, "Frailty" = 1)
alpha_t2 <- c("Long1" = -2, "Long2" = -1, "Frailty" = 1)
alpha_r  <- c("Long1" = -2, "Long2" = -1)
phi1 <- 10000 # Beta precision parameter
sigma2 <- 0.005 # error term SD
sigmaF <- 0.5 # frailty SD
h1_0 <- 0.2
h2_0 <- 0.2
hr_0 <- 2.0
true <- c(betas1, betas2,
          gammas_r, gammas_t1, gammas_t2,
          alpha_r, alpha_t1, alpha_t2)
expit <- function(x) exp(x) / (1 + exp(x))
invS <- function(t, i, h0, tstart, alpha, eta_i, u_i, f_i) {
  h <- function(s) {
    # for the competing risks tstart is allways 0
    X_s <- cbind(1, s + tstart)
    Z_s <- cbind(1, s + tstart)
    eta_y1 <- as.vector(X_s %*% betas1 + rowSums(Z_s * b1[rep(i, nrow(Z_s)), ]))
    eta_y2 <- as.vector(X_s %*% betas2 + rowSums(Z_s * b2[rep(i, nrow(Z_s)), ]))
    # the difference between calendar and gap timescales is in the time var 
    # in the baseline hazard: 's' in gap becomes 's + tstart' in calendar. 
    # Here, because we use a constant baseline hazard, there is no s below, 
    # they are the same 
    exp(log(h0) + eta_i + expit(eta_y1) * alpha[1] + eta_y2 * alpha[2] + f_i * alpha[3])
  }
  integrate(h, lower = 0, upper = t)$value + log(u_i)
}
library(JMbayes2)
tic <- Sys.time(); print(tic)
for (it in seq_len(rep)){
  set.seed(seed + it)
  cat(paste0("\r", it, "/", rep, "   "))
  ## longitudinal 1/2 ==========================================================
  b1 <- MASS::mvrnorm(n, numeric(nrow(D1)), D1)
  b2 <- MASS::mvrnorm(n, numeric(nrow(D2)), D2)
  ## recurrent 1/2 =============================================================
  f <- rnorm(n, mean = 0, sd = sigmaF) # frailty
  ## terminal 1/2 ==============================================================
  group <- rbinom(n, 1, 0.5)
  summ$prop_group1[it] <- mean(group)
  term <- data.frame(id = seq_len(n), group = group, tstart = 0)
  W_t <- model.matrix(~ -1 + group, data = term)
  eta_t1 <- as.vector(W_t %*% gammas_t1)
  eta_t2 <- as.vector(W_t %*% gammas_t2)
  u_t1 <- runif(n)
  u_t2 <- runif(n)
  tstop1 <- tstop2 <- numeric(n)
  for(i in seq_len(n)) {
    root1 <- try(uniroot(invS, interval = c(1e-05, 100), i = i, h0 = h1_0, 
                         tstart = 0, alpha = alpha_t1, eta_i = eta_t1[i], 
                         u_i = u_t1[i], f_i = f[i])$root, TRUE)
    tstop1[i] <- if (!inherits(root1, "try-error")) root1 else t_max
    root2 <- try(uniroot(invS, interval = c(1e-05, 100), i = i, h0 = h2_0,
                         tstart = 0, alpha = alpha_t2, eta_i = eta_t2[i], 
                         u_i = u_t2[i], f_i = f[i])$root, TRUE)  
    tstop2[i] <- if (!inherits(root2, "try-error")) root2 else t_max
  }
  term$tstop <- pmin(tstop1, tstop2)
  term$strata <- as.numeric(tstop1 > tstop2) + 1
  term$status <- term$tstop < t_max
  term$tstop <- pmin(t_max, term$tstop)
  summ$prop_tcensor[it] <- mean(!term$status)
  summ$prop_tevent1[it] <- mean(term$strata == 1 & term$status)
  summ$prop_tevent2[it] <- mean(term$strata == 2 & term$status)
  summ$med_ttime1[it] <- median(term$tstop[term$strata == 1 & term$status])
  summ$med_ttime2[it] <- median(term$tstop[term$strata == 2 & term$status])
  summ$med_tcensor[it] <- median(term$tstop[!term$status])
  if(it != rep) remove(eta_t1, eta_t2, i, root1, root2, u_t1, u_t2, W_t)
  ## recurrent 2/2 =============================================================
  rec <- data.frame(id = seq_len(n), group = group)
  W_r <- model.matrix(~ -1 + group, data = rec)
  eta_r <- as.vector(W_r %*% gammas_r) # we include the frailty in invS()
  tstarts <- tstops <- ids <- list()
  j <- 1
  for(i in seq_len(n)) {
    tstart <- 0
    while(tstart < term$tstop[i]){
      u_r <- runif(1)
      root <- try(uniroot(invS, interval = c(1e-05, 100), i = i, h0 = hr_0, 
                          tstart = tstart, alpha = c(alpha_r, 1), 
                          eta_i = eta_r[i], u_i = u_r, f_i = f[i])$root, TRUE)  
      tstop <- if (!inherits(root, "try-error")) root else t_max 
      tstarts[[j]] <- tstart
      tstops[[j]] <- tstart + tstop
      dur <- 0 # we can add (random or deterministic) event duration
      tstart <- tstart + tstop + dur
      ids[[j]] <- i
      j <- j + 1
    }
  }
  ids <- unlist(ids)
  summ$med_revents[it] <- median(table(ids) - 1) # summary(c(table(ids) - 1))
  rec <- rec[ids, ]
  rec$tstart <- unlist(tstarts)
  rec$tstop <- unlist(tstops)
  rec$status <- rec$tstop < term$tstop[ids]
  rec$tstop <- pmin(rec$tstop, term$tstop[ids])
  rec$strata <- as.factor(1)
  if(it != rep) remove(eta_r, f, group, i, j, root, tstart, tstop, tstarts, 
                       tstops, dur, u_r, W_r)
  ## longitudinal 2/2 ==========================================================
  time <- c(replicate(n, c(0, sort(runif(n_i - 1, 0, t_max)))))
  long <- data.frame(id   = rep(seq_len(n), each = n_i),
                     time = time)
  X <- model.matrix(~ 1 + time, data = long) # FE design matrix 1/2
  Z <- model.matrix(~ 1 + time, data = long) # RE design matrix 1/2
  eta_y1 <- as.vector(X %*% betas1 + rowSums(Z * b1[long$id, ])) # linear predictor 1
  eta_y2 <- as.vector(X %*% betas2 + rowSums(Z * b2[long$id, ])) # linear predictor 2
  mu1 <- expit(eta_y1) # expected value 1
  mu2 <- eta_y2 # expected value 2
  shape1 <- phi1 * mu1
  shape2 <- phi1 * (1 - mu1) 
  y1 <- rbeta(length(shape1), shape1, shape2)
  y1 <- (y1 * (length(y1) - 1) + 0.5) / length(y1)
  long$y1 <- y1
  long$y2 <- rnorm(length(mu2), mu2, sigma2)
  long <- long[long$time <= term$tstop[long$id], ]
  summ$N1[it] <- summ$N2[it] <- length(long$y1) # y1 and y2 share the obs times
  summ$med_ni1[it] <- summ$med_ni2[it] <- median(table(long$id)) # y1 and y2 share the obs times
  head_rows <- tapply(seq_along(long$id), long$id, head, 1)
  tail_rows <- tapply(seq_along(long$id), long$id, tail, 1)
  time_diff <- long$time[tail_rows] - long$time[head_rows]
  summ$afp1[it] <- summ$afp2[it] <- sum(time_diff)  # y1 and y2 share the obs times
  summ$med_fpi1[it] <- summ$med_fpi2[it] <- median(time_diff) # y1 and y2 share the obs times
  if(it != rep) remove(b1, b2, head_rows, eta_y1, eta_y2, mu1, mu2, shape1, 
                       shape2, time, time_diff, X, y1, Z)
  ## terminal 2/2 ==============================================================
  term <- term[rep(term$id, each = 2), ]
  rank <- rep(seq_len(2), times = n)
  term$status <- term$strata == rank & term$status
  term$strata <- as.factor(rank + 1) # + 1, strata 1 is the recurrent events
  if(it != rep) remove(rank)
  ## data ======================================================================
  surv <- rbind(rec, term)
  surv <- surv[order(surv$id, surv$strata, surv$tstart), ]
  surv$status <- as.numeric(surv$status)
  rownames(surv) <- seq_along(surv$id)
  if(it != rep) remove(term, rec)
  ## joint model ===============================================================
  m1 <- mixed_model(y1 ~ time, random = ~ 1 + time || id, data = long,
                    family = beta.fam())
  m2 <- lme(y2 ~ time, random = list(id = pdDiag(~ 1 + time)), data = long,
            control = list(opt = "optim"))
  ph <- coxph(Surv(tstart, tstop, status) ~ group : strata(strata), data = surv)
  jm <- jm(ph, list(m1, m2), time_var = "time", recurrent = "gap",
           functional_forms = ~ (vexpit(value(y1)) + value(y2)):strata,
           n_iter = n_iter, n_burnin = n_burnin)
  if(it == it_s) {
    saveRDS(jm, paste0(ref, "/jm", ref, ".rds"))
    pdf(paste0(ref, "/tcp_", ref, ".pdf"), width = 12, height = 16)
    {
      par(mfrow = c(4, 3))
      traceplot(jm, "all")
      mtext(paste0("Dataset ", it_s, "/", rep), side = 3 , line = -2, 
            outer = TRUE)
    }
    dev.off()
  }
  if(it != rep) remove(m1, m2, ph)
  nms <- c("betas1", "betas2", "gammas", "alphas", "alphaF")
  ord <- c(1:8, 11,  9, 12, 14, 10, 13, 15)
  est[it, ] <- unlist(jm$statistics$Mean[nms])[ord]
  low <- unlist(jm$statistics$CI_low[nms])[ord]
  upp <- unlist(jm$statistics$CI_upp[nms])[ord]
  cvg[it, ] <- true <= upp & true >= low
  if(it != rep) remove(low, upp, nms, ord)
  ## plot ======================================================================
  if(it == it_s){
    rx <- 100 # resolution
    toy <- data.frame(time = seq(0, t_max, length.out = rx))
    X_toy <- model.matrix(~ 1 + time, data = toy)
    toy$ytrue1 <- expit(X_toy %*% betas1)
    toy$ytrue2 <- X_toy %*% betas2
    toy$yhat1 <- expit(X_toy %*% jm$statistics$Mean$betas1)
    toy$yhat2 <- X_toy %*% jm$statistics$Mean$betas2
    pdf(paste0(ref, "/long_prof_", ref, ".pdf"), width = 15, height = 5)
    {
      par(mfrow = c(1, 4))
      ### 1/4
      plot(NA, ylim = range(long$y1), xlim = c(0, t_max), 
           xlab = "Time", ylab = "Y1", xaxs = "i", yaxs = "i")
      invisible(tapply(seq_along(long$id), long$id, function(rows) {
        lines(long$time[rows], long$y1[rows], col = rgb(0, 0, 0, 0.05))
      }))
      lines(toy$time, toy$ytrue1, lwd = 2, col = "green") # true model
      lines(toy$time, toy$yhat1, lwd = 2, lty = 2, col = "yellow") # estimated
      mtext(text = paste0("Dataset ", it_s, "/", rep), side = 3, line = 0)
      
      ### 2/4
      plot(NA, ylim = range(long$y2), xlim = c(0, t_max), 
           xlab = "Time", ylab = "Y2", xaxs = "i", yaxs = "i")
      invisible(tapply(seq_along(long$id), long$id, function(rows) {
        lines(long$time[rows], long$y2[rows], col = rgb(0, 0, 0, 0.05))
      }))
      lines(toy$time, toy$ytrue2, lwd = 2, col = "green") # true model
      lines(toy$time, toy$yhat2, lwd = 2, lty = 2, col = "yellow") # estimated
      ### 3/4
      boxplot(c(tstop1, tstop2) ~ rep(1:2, each = n), horizontal = TRUE, 
              ylab = "", xlab = "Survival time")
      abline(v = t_max, col = 2)
      ### 4/4
      bp <- boxplot(table(ids), ylab = "# recurrent events/id", xlab = "Time")
      text(1.35, bp$stats, labels = round(bp$stats, 2))
    }
    dev.off()
    remove(rx, toy, X_toy)
  }
  if(it != rep) remove(ids, long, jm, tstop1, tstop2)
}
toc <- Sys.time(); print(toc)
diff <- difftime(toc, tic); print(diff)
cvg <- colMeans(cvg)
bias <- colMeans(est) - true 
mse <- colMeans(sweep(est, 2, true)^2)
res <- list(cvg = cvg, est = est, bias = bias, mse = mse, true = true, 
            time = diff, summ = summ)
saveRDS(res, file = paste0(ref, "/res_", ref,".rds"))

pdf(paste0(ref, "/res_", ref, ".pdf"), width = 10, height = 10) 
{
  par(mfrow = c(3, 5))
  for(i in seq_len(15)){
    par(font.main = 1)
    boxplot(est[, i], main = colnames(est)[i], xlab = "", 
            ylab = expression(hat(theta)))
    abline(h = true[i], col = 2, lwd = 2)
    label1 <- c("Bias:", round(bias[i], 2))
    mtext(label1, side = 1, at = c(0.5, 1), line = 3, cex = 0.7)
    label2 <- c("CP (%):", round(cvg[i] * 100, 2))
    mtext(label2, side = 1, at = c(0.5, 1), line = 4, cex = 0.7)
  }
} 
dev.off()