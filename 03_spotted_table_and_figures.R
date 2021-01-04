# Spotted Results processing 
library(moveds)
library(Distance)
library(ggplot2)
theme_set(theme_bw())

# Make Table 1 --------------------------------------------------------------
by_year <- readRDS("byyear.Rds")
tab <- matrix(0, nr = 5, nc = 6)

# titles
rownames(tab) <- c("Year", 1999, 2000, 2003, 2006)
colnames(tab) <- c("CDS1D", " ",  "CDS2D", " ", "MDS2D", " ")
tab[1, ] <- rep(c("Est.", "Var."), 3)

# multipliers
A <- 1
S <- 1

# rounding set 
rnd <- 0

## cds1d
cds1d <- by_year[[1]]
cds1d <- cds1d[-5,]
# rounding
cds1d[,2] <- round(cds1d[,2], rnd)
cds1d[,4] <- round(cds1d[,4]*100*S/A, rnd) 
cds1d[,5:6] <- round(cds1d[,5:6], rnd)
# means
tab[2:5, 1] <- cds1d[,2] 
# variances
tab[2:5, 2] <- paste0(cds1d[,4], "\\% (", cds1d[,5], ", ", cds1d[,6], ")")

## cds2d 
cds2d <- by_year[[2]]
restab <- t(sapply(cds2d, FUN = function(x){as.numeric(x$result[3,])}))
# scale
restab <- restab * A / S
# CV
restab[,2] <- restab[,2] / restab[,1]
# rounding 
restab[,1] <- round(restab[,1], rnd)
restab[,2] <- round(restab[,2]*100, rnd)
restab[,3:4] <- round(restab[,3:4], rnd)
# means
tab[2:5, 3] <- restab[, 1]
# variances
tab[2:5, 4] <- paste0(restab[,2], "\\% (", restab[,3], ", ", restab[,4], ")")

## mds2d 
mds2d <- by_year[[3]]
restab <- t(sapply(mds2d, FUN = function(x){as.numeric(x$result[4,])}))
# scale
restab <- restab *  A / S
# CV
restab[,2] <- restab[,2] / restab[,1]
# rounding 
restab[,1] <- round(restab[,1], rnd)
restab[,2] <- round(restab[,2]*100, rnd)
restab[,3:4] <- round(restab[,3:4], rnd)
# means
tab[2:5, 5] <- restab[, 1]
# variances
tab[2:5, 6] <- paste0(restab[,2], "\\% (", restab[,3], ", ", restab[,4], ")")

saveRDS(tab, "Table1.Rds")


# Make Figure 3 --------------------------------------------------------
mods <- readRDS("final_mods.Rds")

pred.cds1d <- ds.gof(mods$cds1d, breaks = seq(0, 5.5, 0.5))$chisquare$chi1$expected

xseq <- seq(-5.5, 5.0, 0.5)
cds2d <- mods$cds2d
pred.cds2d.full <- rowSums(plot(cds2d, extract = TRUE))
mat <- matrix(pred.cds2d.full, nr = 11, nc = 2)
mat[,1] <- rev(mat[,1])
pred.cds2d <- rowSums(mat) * cds2d$result[3,1]

mds2d <- mods$mds2d
pred.mds2d.full <- rowSums(plot(mds2d, delta = c(cds2d$ds$delta[1], mds2d$ds$delta[2]), extract = TRUE))
pred.mds2d.sub <- pred.mds2d.full[-1]
mat <- matrix(pred.mds2d.sub, nr = 11, nc = 2)
mat[,1] <- rev(mat[,1])
pred.mds2d <- rowSums(mat) * mds2d$result[4,1]

xs<- seq(0, 5.0, 0.5)
fit0 <- data.frame(x = xs, xend = xs, yend = pred.cds1d, y = 0)
fit1 <- data.frame(x = xs, xend = xs, yend = pred.cds2d, y = 0)
fit2 <- data.frame(x = xs, xend = xs, yend = pred.mds2d, y = 0)
dat <- cds2d$ds$data 
dat$x <- abs(dat$x)
dat <- dat[dat$x<=5.5,]
gr_hist <- ggplot(dat, aes(x = abs(x))) + 
  geom_histogram(binwidth = 0.5, center = 0.25, color = "white", fill = "grey90")  

pdf("Figure_3.pdf")
gr_hist + 
  geom_segment(data = fit0, aes(x = x + 0.5/4, y = y, xend = xend + 0.5/4, yend= yend), col = "blue", lwd = 1.0) + 
  geom_segment(data = fit1, aes(x = x + 2 * 0.5/4, y = y, xend = xend + 2*0.5/4, yend= yend), col = "black", 
               linetype = "dashed", lwd = 1.0) + 
  geom_segment(data = fit2, aes(x = x + 3 * 0.5/4, y = y, xend = xend + 3*0.5/4, yend= yend), col = "red", 
               linetype = "dotted", lwd = 1.0) + 
  geom_point(data = fit0, aes(x = x + 0.5/4, y = yend), col = "blue", size = 2.5) + 
  geom_point(data = fit1, aes(x = x + 2 * 0.5 / 4, y = yend), col = "black", size = 2.5) + 
  geom_point(data = fit2, aes(x = x + 3 * 0.5 / 4, y = yend), col = "red", size = 2.5) + 
  scale_x_continuous("Perpendicular distance (km)", breaks = seq(0, 5.5, 0.5), label = seq(0, 5.5, 0.5)) + 
  scale_y_continuous("Number of sightings", breaks = seq(0, 50, 5), label = c("0", "", "10", "", "20", "", "30", "", "40", "", "50")) + 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        strip.text.x = element_blank(), axis.text=element_text(size =20), axis.title=element_text(size = 20)) 
dev.off()


# Make Figure 4 ---------------------------------------------------------

xseq <- seq(0, 5.5, 0.01)
det.est <- as.numeric(exp(mods$cds1d$ddf$ds$par))
cds2d.est <- s2sigmab(cds2d$result[1,1], cds2d$result[2,1])
mds2d.est <- s2sigmab(mds2d$result[1,1], mds2d$result[2,1])
g0 <- data.frame(x = xseq, y = 1-exp(-(xseq/det.est[2])^(-det.est[1])))
g1 <- data.frame(x = xseq, y = 1-exp(-(xseq/(cds2d.est[2]))^(-cds2d.est[1])))
g2 <- data.frame(x = xseq, y = 1-exp(-(xseq/mds2d.est[2])^(-mds2d.est[1])))
pdf("Figure_4.pdf")
ggplot(g0) + 
  geom_line(data = g0, aes(x = x, y = y), lwd = 1.0, col = "blue") + 
  geom_line(data = g1, aes(x = x, y = y), lwd = 1.0, col = "black", linetype = "dashed") + 
  geom_line(data = g2, aes(x = x, y = y), lwd = 1.0, col = "red", linetype = "dotted") + 
  scale_x_continuous("Perpendicular distance (km)", breaks = seq(0, 5.5, 0.5), label = seq(0, 5.5, 0.5)) + 
  scale_y_continuous("Detection probability", breaks = seq(0, 1, 0.1), label = seq(0, 1, 0.1)) + 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        strip.text.x = element_blank(), axis.text=element_text(size =20), axis.title=element_text(size = 20)) 
dev.off()


