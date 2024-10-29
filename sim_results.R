rm(list = ls())
# data =========================================================================
## scenario A ==================================================================
ref_A <- "A2001"
res_A  <- readRDS(paste0(ref_A, "/res_", ref_A,".rds"))
jm_A <- readRDS(paste0(ref_A, "/jm", ref_A, ".rds"))
## scenario B ==================================================================
ref_B <- "B2001"
res_B  <- readRDS(paste0(ref_B, "/res_", ref_B,".rds"))
jm_b_B <- readRDS(paste0(ref_B, "/jm_b_", ref_B, ".rds"))
jm_g_B <- readRDS(paste0(ref_B, "/jm_g_", ref_B, ".rds"))
# functions ====================================================================
iqr <- function(x, dig = 2) {
  q <- round(quantile(x, prob = c(0.5, 0.25, 0.75)), dig)
  paste0(q[1], " (", q[2], "--", q[3], ")")
}
hist2 <- function(v1, v2 = v1, lwd1 = 1.5, lwd2 = 1.5, col1 = 1, col2 = 2) {
  h1 <- hist(v1, plot = FALSE)
  x1 <- h1$breaks
  y1 <- h1$density
  h2 <- hist(v2, plot = FALSE)
  x2 <- h2$breaks
  y2 <- h2$density
  xlim <- range(x1, x2)
  ylim <- range(y1, y2)
  plot(NA, xlim = xlim, ylim = ylim, xlab = "", ylab = "", bty = "n", xaxt = "n",
       yaxt = "n")
  box(lwd = 0.5)
  axis(1, col.axis = 1, col = NA, col.ticks = 1)
  axis(2, col.axis = 1, col = NA, col.ticks = 1, las = 2)
  segments(x0 = head(x1, -1), x1 = tail(x1, -1), y0 = y1, lwd = lwd1, col = col1) # horizontal lines 1
  segments(y0 = c(par("usr")[3], y1), y1 = c(y1, par("usr")[3]), x0 = x1, 
           lwd = lwd1, col = col1) # vertical lines 1
  if(!identical(v1, v2)) {
    segments(x0 = head(x1, -1), x1 = tail(x1, -1), y0 = y1, lwd = lwd2, col = col2) # horizontal lines 2
    segments(y0 = c(par("usr")[3], y1), y1 = c(y1, par("usr")[3]), x0 = x1, 
             lwd = lwd2, col = col2) # vertical lines 2
  }
}
# tab main document ============================================================
## scenario A ==================================================================
round(colMeans(res_A$est), 3)[c(1:5, 8:9, 6, 10:12, 7, 13:15)] # estimate
round(res_A$bias, 3)[c(1:5, 8:9, 6, 10:12, 7, 13:15)] # bias
round(res_A$mse, 3)[c(1:5, 8:9, 6, 10:12, 7, 13:15)] # mse
round(res_A$ese, 3)[c(1:5, 8:9, 6, 10:12, 7, 13:15)] # ese
round(res_A$psd, 3)[c(1:5, 8:9, 6, 10:12, 7, 13:15)] # sd
round(res_A$cvg, 3)[c(1:5, 8:9, 6, 10:12, 7, 13:15)] # cvg
## scenario B ==================================================================
round(colMeans(res_B$est$b), 3) # estimates beta
round(res_B$bias$b, 3) # bias beta
round(res_B$mse$b, 3) # mse beta
round(res_B$ese$b, 3) # ese beta
round(res_B$psd$b, 3) # sd beta
round(res_B$cvg$b, 3) # cvg beta
round(res_B$bias$g, 3) # bias Gaussian
round(colMeans(res_B$est$g), 3) # estimates Gaussian
round(res_B$mse$g, 3) # mse Gaussian
round(res_B$ese$g, 3) # ese Gaussian
round(res_B$psd$g, 3) # sd Gaussian
round(res_B$cvg$g, 3) # cvg Gaussian
# tab web material =============================================================
## scenario A ==================================================================
iqr(res_A$summ$N1) # number of observations 1
iqr(res_A$summ$N2) # number of observations 2
iqr(res_A$summ$med_ni1) # number of observations/individual 1
iqr(res_A$summ$med_ni2) # number of observations/individual 2
iqr(res_A$summ$afp1) # aggregated follow-up duration 1
iqr(res_A$summ$afp2) # aggregated follow-up duration 2
iqr(res_A$summ$med_fpi1) # follow-up duration/individual 1
iqr(res_A$summ$med_fpi2) # follow-up duration/individual 2
iqr(res_A$summ$med_ttime1) # terminal event time 1
iqr(res_A$summ$med_ttime2) # terminal event time 2
iqr(res_A$summ$med_tcensor) # censoring event time
iqr(res_A$summ$prop_tevent1) # proportion of event 1
iqr(res_A$summ$prop_tevent2) # proportion of event 2
iqr(res_A$summ$prop_tcensor) # proportion of censoring
iqr(res_A$summ$med_revents) # number of recurrent events
iqr(res_A$summ$prop_group1) # proportion of group 1
cairo_ps(paste0("hist_ni_A.eps"), # number of observations/individual
         height = 3 * 1, width = 5 * 1, family = "helvetica")
{
  par(mar = c(5.1 - 1.5, 4.1 - 0.5, 4.1 - 3.5, 2.1 - 2.1))
  hist2(res_A$summ$med_ni1, res_A$summ$med_ni2)
  mtext("Number of observations/individual", 1, line = 2.5)
  mtext("Density", 2, line = 2.5)
}
dev.off()
cairo_ps(paste0("hist_fpi_A.eps"), # follow-up duration/individual
         height = 3 * 1, width = 5 * 1, family = "helvetica")
{
  par(mar = c(5.1 - 1.5, 4.1 - 0.5, 4.1 - 3.5, 2.1 - 2.1))
  hist2(res_A$summ$med_fpi1, res_A$summ$med_fpi2)
  mtext("Median follow-up duration/individual", 1, line = 2.5)
  mtext("Density", 2, line = 2.5)
}
dev.off()
cairo_ps(paste0("hist_ttime_A.eps"), # event time
         height = 3 * 1, width = 5 * 1, family = "helvetica")
{
  par(mar = c(5.1 - 1.5, 4.1 - 0.5, 4.1 - 3.5, 2.1 - 2.1))
  hist2(res_A$summ$med_ttime1)
  mtext("Event time", 1, line = 2.5)
  mtext("Density", 2, line = 2.5)
}
dev.off()
## scenario B ==================================================================
iqr(res_B$summ$N) # number of observations
iqr(res_B$summ$med_ni) # number of observations/individual
iqr(res_B$summ$afp) # aggregated follow-up duration
iqr(res_B$summ$med_fpi) # follow-up duration/individual
iqr(res_B$summ$med_ttime) # terminal event time
iqr(res_B$summ$med_tcensor) # censoring event time
iqr(res_B$summ$prop_tevent) # proportion of event
iqr(res_B$summ$prop_tcensor) # proportion of censoring
iqr(res_B$summ$prop_group1) # proportion of group 1
cairo_ps(paste0("hist_ni_B.eps"), # number of observations/individual
         height = 3 * 1, width = 5 * 1, family = "helvetica")
{
  par(mar = c(5.1 - 1.5, 4.1 - 0.5, 4.1 - 3.5, 2.1 - 2.1))
  hist2(res_B$summ$med_ni)
  mtext("Number of observations/individual", 1, line = 2.5)
  mtext("Density", 2, line = 2.5)
}
dev.off()
cairo_ps(paste0("hist_fpi_B.eps"), # follow-up duration/individual
         height = 3 * 1, width = 5 * 1, family = "helvetica")
{
  par(mar = c(5.1 - 1.5, 4.1 - 0.5, 4.1 - 3.5, 2.1 - 2.1))
  hist2(res_B$summ$med_fpi)
  mtext("Median follow-up duration/individual", 1, line = 2.5)
  mtext("Density", 2, line = 2.5)
}
dev.off()
cairo_ps(paste0("hist_ttime_B.eps"), # event time
         height = 3 * 1, width = 5 * 1, family = "helvetica")
{
  par(mar = c(5.1 - 1.5, 4.1 - 0.5, 4.1 - 3.5, 2.1 - 2.1))
  hist2(res_B$summ$med_ttime)
  mtext("Event time", 1, line = 2.5)
  mtext("Density", 2, line = 2.5)
}
dev.off()
# traceplot ====================================================================
## scenario A ==================================================================
names <- c(expression(beta["1,0"]), expression(beta["1,t"]), 
           expression(beta["2,0"]), expression(beta["2,t"]), 
           expression(gamma^"R"), expression(gamma["1"]^"T"), 
           expression(gamma["2"]^"T"),
           expression(alpha["1"]^"R"), expression(alpha["2"]^"R"), 
           expression(alpha["1,1"]^"T"), expression(alpha["1,2"]^"T"),
           expression(alpha["2,1"]^"T"), expression(alpha["2,2"]^"T"),
           expression(alpha["1"]^upsilon), expression(alpha["2"]^upsilon))
mcmc <- list()
mcmc <- jm_A$mcmc[c("betas1", "betas2", "gammas", "alphas", "alphaF")]
Rhat <- do.call(rbind, jm_A$statistics$Rhat[c("betas1", "betas2", "gammas", "alphas", "alphaF")])[, 1] 
cols <- rgb(0, 0, 0, alpha = c(0.9, 0.6, 0.3))
cairo_ps(paste0("tcp_A.eps"), height = 12 * 0.8, width = 12 * 0.8, family = "helvetica")
{
  par(mfrow = c(4, 4), mar = c(5.1, 4.1 + 0.5, 4.1 + 1, 2.1))
  k <- 1
  for(i in seq_len(5)) {
    for(j in seq_len(ncol(mcmc[[i]][[1]]))) {
      ylim <- range(do.call(rbind, mcmc[[i]])[, j])
      ylim[2] <- ylim[2] + 0.2 * abs(diff(ylim))
      xlim <- c(1, nrow(mcmc[[1]][[1]]))
      plot(NA, xlim = xlim, ylim = ylim, bty = "n", xlab = "", ylab = "", 
           xaxt = "n", yaxt = "n", main = "")
      box(lwd = 0.5)
      axis(1, at = c(1, 625, 1250, 1875, 2500), labels = c(1, NA, 1250, NA, 2500),
           col.axis = 1, col = NA, col.ticks = 1)
      axis(2, col.axis = 1, col = NA, col.ticks = 1, las = 2)
      mtext("Iterations", side = 1, line = 2, cex = 0.65)
      mtext("Sample value", side = 2, line = 3.5, cex = 0.65)
      mtext(names[k], side = 3, line = 0.5, cex = 0.8)
      lines(x = seq_len(xlim[2]), y = mcmc[[i]][[1]][, j], col = cols[1])
      lines(x = seq_len(xlim[2]), y = mcmc[[i]][[2]][, j], col = cols[2])
      lines(x = seq_len(xlim[2]), y = mcmc[[i]][[3]][, j], col = cols[3])
      text(x = xlim[1], y = ylim[2], adj = c(0, 1), 
           labels = bquote(hat(R) == .(sprintf("%.2f", Rhat[k]))))
      k <- k + 1
    }
  }
  ## legend
  par(fig = c(0, 1, 0, 0.5), mar = c(0, 0, 0, 0), new = TRUE)
  plot(NA, xlim = 0:1, ylim = 0:1, bty = "n", xaxt = "n", yaxt = "n",
       xaxs = "i", yaxs = "i")
  text(x = 0.999, y = 0, adj = c(0, 0), cex = 0.9, srt = 90,
       labels = paste0("Dataset: 66/500"))
}
dev.off()
## scenario B ==================================================================
names <- c(expression(beta["1,0"]), expression(beta["1,t"]), 
           expression(gamma["1"]^"T"), expression(alpha["1,1"]^"T"))
mcmc <- Rhat <- list()
mcmc$g <- jm_g_B$mcmc[c("betas1", "gammas", "alphas")]
mcmc$b <- jm_b_B$mcmc[c("betas1", "gammas", "alphas")]
Rhat$g <- do.call(rbind, jm_g_B$statistics$Rhat[c("betas1", "gammas", "alphas")])[, 1] 
Rhat$b <- do.call(rbind, jm_b_B$statistics$Rhat[c("betas1", "gammas", "alphas")])[, 1] 
cols <- rgb(0, 0, 0, alpha = c(0.9, 0.6, 0.3))
cairo_ps(paste0("tcp_B.eps"), height = 6 * 0.8, width = 12 * 0.8, family = "helvetica")
{
  par(mfrow = c(2, 4), mar = c(5.1, 4.1 + 0.5, 4.1 + 1, 2.1))
  ## top - beta
  k <- 1
  for(i in seq_len(3)) {
    for(j in seq_len(ncol(mcmc$b[[i]][[1]]))) {
      ylim <- range(do.call(rbind, mcmc$b[[i]])[, j])
      ylim[2] <- ylim[2] + 0.2 * abs(diff(ylim))
      xlim <- c(1, nrow(mcmc$b[[1]][[1]]))
      plot(NA, xlim = xlim, ylim = ylim, bty = "n", xlab = "", ylab = "", 
           xaxt = "n", yaxt = "n", main = "")
      box(lwd = 0.5)
      axis(1, at = c(1, 625, 1250, 1875, 2500), labels = c(1, NA, 1250, NA, 2500),
           col.axis = 1, col = NA, col.ticks = 1)
      axis(2, col.axis = 1, col = NA, col.ticks = 1, las = 2)
      mtext("Iterations", side = 1, line = 2, cex = 0.65)
      mtext("Sample value", side = 2, line = 3.5, cex = 0.65)
      mtext(names[k], side = 3, line = 0.5, cex = 0.8)
      lines(x = seq_len(xlim[2]), y = mcmc$b[[i]][[1]][, j], col = cols[1])
      lines(x = seq_len(xlim[2]), y = mcmc$b[[i]][[2]][, j], col = cols[2])
      lines(x = seq_len(xlim[2]), y = mcmc$b[[i]][[3]][, j], col = cols[3])
      text(x = xlim[1], y = ylim[2], adj = c(0, 1), 
           labels = bquote(hat(R) == .(sprintf("%.2f", Rhat$b[k]))))
      k <- k + 1
    }
  }
  ## bottom - Gaussian
  k <- 1
  for(i in seq_len(3)) {
    for(j in seq_len(ncol(mcmc$g[[i]][[1]]))) {
      ylim <- range(do.call(rbind, mcmc$g[[i]])[, j])
      ylim[2] <- ylim[2] + 0.2 * abs(diff(ylim))
      xlim <- c(1, nrow(mcmc$g[[1]][[1]]))
      plot(NA, xlim = xlim, ylim = ylim, bty = "n", xlab = "", ylab = "", 
           xaxt = "n", yaxt = "n", main = "")
      box(lwd = 0.5)
      axis(1, at = c(1, 625, 1250, 1875, 2500), labels = c(1, NA, 1250, NA, 2500),
           col.axis = 1, col = NA, col.ticks = 1)
      axis(2, col.axis = 1, col = NA, col.ticks = 1, las = 2)
      mtext("Iterations", side = 1, line = 2, cex = 0.65)
      mtext("Sample value", side = 2, line = 3.5, cex = 0.65)
      mtext(names[k], side = 3, line = 0.5, cex = 0.8)
      lines(x = seq_len(xlim[2]), y = mcmc$g[[i]][[1]][, j], col = cols[1])
      lines(x = seq_len(xlim[2]), y = mcmc$g[[i]][[2]][, j], col = cols[2])
      lines(x = seq_len(xlim[2]), y = mcmc$g[[i]][[3]][, j], col = cols[3])
      text(x = xlim[1], y = ylim[2], adj = c(0, 1), 
           labels = bquote(hat(R) == .(sprintf("%.2f", Rhat$g[k]))))
      k <- k + 1
    }
  }
  ## legend
  par(fig = c(0, 1, 0, 0.5), mar = c(0, 0, 0, 0), new = TRUE)
  plot(NA, xlim = 0:1, ylim = 0:1, bty = "n", xaxt = "n", yaxt = "n",
       xaxs = "i", yaxs = "i")
  mtext("Gaussian", side = 3, at = 0.5, line = -1.5, cex = 0.9)
  text(x = 0.999, y = 0, adj = c(0, 0), cex = 0.9, srt = 90,
       labels = paste0("Dataset: 66/500"))
  par(fig = c(0, 1, 0.5, 1), mar = c(0, 0, 0, 0), new = TRUE)
  plot(NA, xlim = 0:1, ylim = 0:1, bty = "n", xaxt = "n", yaxt = "n",
       xaxs = "i", yaxs = "i")
  mtext("Beta", side = 3, at = 0.5, line = -1.5, cex = 0.9)
}
dev.off()
# post mean ====================================================================
## scenario A ==================================================================
names <- c(expression(beta["1,0"]), expression(beta["1,t"]), 
           expression(beta["2,0"]), expression(beta["2,t"]), 
           expression(gamma^"R"), expression(gamma["1"]^"T"), 
           expression(gamma["2"]^"T"),
           expression(alpha["1"]^"R"), expression(alpha["2"]^"R"), 
           expression(alpha["1,1"]^"T"), expression(alpha["1,2"]^"T"),
           expression(alpha["1"]^upsilon),
           expression(alpha["2,1"]^"T"), expression(alpha["2,2"]^"T"),
           expression(alpha["2"]^upsilon))
rep <- 500
true <- c(2, -1.5, 0.8, -0.05, 0.25, 0.25, 0.25, -2, -1, -2, -1, 1, -2, -1, 1)
cairo_ps(paste0("pmean_A.eps"), height = 12 * 0.8, width = 12 * 0.8, family = "helvetica")
{
  par(mfrow = c(3, 5))
  par(mar = c(5.1, 4.1 + 0.5, 4.1, 2.1 + 1.75))
  for(i in seq_len(15)){
    boxplot(res_A$est[, i], main = names[i], xlab = "", ylab = "",
            frame = FALSE, xaxt = "n", yaxt = "n")
    box(lwd = 0.5)
    axis(2, col.axis = 1, col = NA, col.ticks = 1, las = 2,
         at = axTicks(2), labels = round(axTicks(2), 2))
    axis(4, col.axis = 1, col = NA, col.ticks = 1, las = 2,
         at = axTicks(4), labels = round(axTicks(4) - true[i], 2))
    mtext("Posterior mean", side = 2, cex = 0.65, line = 3)
    mtext("Error", side = 4, cex = 0.65, line = 3)
    abline(h = true[i], col = 1, lwd = 1.5)
  }
}
dev.off()
## scenario B ==================================================================
names <- c(expression(beta["1,0"]), expression(beta["1,t"]), 
           expression(gamma["1"]^"T"), expression(alpha["1,1"]^"T"))
rep <- 500
true <- c(2, -1, 0.25, -2)
cairo_ps(paste0("pmean_B.eps"), height = 4 * 0.8, width = 12 * 0.8, family = "helvetica")
{
  par(mfrow = c(1, 4))
  par(mar = c(5.1, 4.1 + 0.5, 4.1, 2.1 + 1.75))
  for(i in seq_len(4)){
    boxplot(c(res_B$est$b[, i], res_B$est$g[, i]) ~ rep(1:2, each = rep), 
            main = names[i], xlab = "", ylab = "",
            frame = FALSE, xaxt = "n", yaxt = "n")
    box(lwd = 0.5)
    axis(2, col.axis = 1, col = NA, col.ticks = 1, las = 2,
         at = axTicks(2), labels = round(axTicks(2), 2))
    axis(4, col.axis = 1, col = NA, col.ticks = 1, las = 2,
         at = axTicks(4), labels = round(axTicks(4) - true[i], 2))
    axis(1, at = 1:2, labels = c("Beta\nmodel", "Gaussian\nmodel"), col.axis = 1, 
         col = NA, col.ticks = 1, padj = 0.5)
    mtext("Posterior mean", side = 2, cex = 0.65, line = 3)
    mtext("Error", side = 4, cex = 0.65, line = 3)
    abline(h = true[i], col = 1, lwd = 1.5)
  }
}
dev.off()
# computational time ===========================================================
round(quantile(res_A$rtime / 60, prob = c(0.25, 0.5, 0.75)), 2)