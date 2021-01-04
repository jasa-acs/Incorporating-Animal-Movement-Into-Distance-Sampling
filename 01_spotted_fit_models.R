## moveds example
library(moveds)
library(Distance)

# Read in Data ------------------------------------------------------------

spotted <- readRDS("spotted_data.rds")

aux <- c(1000,1000, 5.5, mean(spotted$transect[,6]), 0)
obs <- as.data.frame(spotted$enc)[,-1]
transdat <- spotted$transect
n.transects <- nrow(transdat)
distances <- abs(spotted$enc[,3])
movedat <- spotted$tags


# 1D CDS model ------------------------------------------------
names(obs) <- c("transect", "x", "y", "t")
region.table <- data.frame(Region.Label = 1, Area = prod(aux[1:2]))
sample.table <- data.frame(Sample.Label = transdat[,1], Region.Label = 1,
                           Effort = transdat[, 2])
obs.table <- data.frame(object = 1:nrow(obs), Region.Label = 1, Sample.Label = obs[, 1])
distances <- abs(obs$x)

# fit 1D CDS model
cds1d <- ds(distances,
            truncation = aux[3],
            key = "hr",
            adjustment = NULL,
            region.table = region.table,
            sample.table = sample.table,
            obs.table = obs.table)

# get estimated detection parameters
det.est <- as.numeric(exp(cds1d$ddf$ds$par))
# plot estimated detection function
plot(cds1d)
# check gof
cds1d.gof <- ds.gof(cds1d, breaks = seq(0, aux[3], 0.5))
# get estimated density
N.cds1d <- as.numeric(cds1d$dht$individuals$N$Estimate)


# 2D CDS model ------------------------------------------------------------
# fit base 2d CDS model 
ds <- list(data = obs,
           transect = transdat[,1:2],
           aux = aux,
           delta = c(0.25, 9),
           buffer = 0,
           hazardfn = 1,
           move = 0)
move <- list(data = movedat)
cds2d <- mds(ds, move, start = c(s = 1, d = 1), print = TRUE)
summary(cds2d)

plot(cds2d)
cds2d.gof <- mds.gof(cds2d)

# 2D MDS model ------------------------------------------------------------
ds$move <- 1
ds$delta <- c(0.5, 0.5)
ds$buffer <- 1
mds2d <- mds(ds, move, start = c(s = 1.0, d = 2.0, sd = 0.6), print = TRUE)
summary(mds2d)

# check dx  
dx_conv <- check.dx(ds, move, par = c(s = mds2d$result[1,1], d = mds2d$result[2,1], sd = mds2d$result[3,1]), dx = c(0.1, 1), 0.05, print = TRUE)

# converged < 0.2 for sure 
ds$delta <- c(0.2, 0.5)
# check dt, <1%
dt_conv <- check.dt(ds, move, par = c(s = mds2d$result[1,1], d = mds2d$result[2,1], sd = mds2d$result[3,1]), dt = c(0.1, 1), by = 0.1)

# fit model with stable times-step
ds$delta <- c(0.2, 0.2)
mds2d <- mds(ds, move, start = c(s = mds2d$result[1,1], d = mds2d$result[2,1], sd = mds2d$result[3,1]), print = TRUE)

mods <- list(cds1d = cds1d, cds2d = cds2d, mds2d = mds2d)

saveRDS(mods, "final_mods.Rds") 

######

summary(mds2d)
plot(mds2d, delta = c(0.25, 0.2))

mds2d.gof <- mds.gof(mds2d, delta = c(0.5, 0.5))
xpdf <- rowSums(mds2d.gof)
xpdf <- matrix(xpdf, nc = 2)
xpdf[,2] <- rev(xpdf[,2])
xpdf <- rev(rowMeans(xpdf))
xpdf <- xpdf/sum(xpdf)
obsn <- hist(distances[distances<=5.5])$count
npdf2 <- xpdf * 1256.8024 * 5.5 * sum(transdat[,2]) / 1000^2
chisq.test(npdf[-12], obsn)

# Plot estimated detection functions --------------------------------------
xseq <- seq(0, 5.5, 0.01)
g.cds1d <- 1-exp(-(xseq/det.est[2])^(-det.est[1]))
cds2d.est <- s2sigmab(cds2d$result[1,1], cds2d$result[2,1], v = 1)
g.cds2d <- 1-exp(-(xseq/(cds2d.est[2]))^(-cds2d.est[1]))
mds2d.est <- s2sigmab(mds2d$result[1,1], mds2d$result[2,1], v = 1)
g.mds2d <- 1-exp(-(xseq/mds2d.est[2])^(-mds2d.est[1]))
plot(xseq, g.cds1d, type = "l", lwd = 1.5, col = "blue", xlab = "Perpendicular Distance (km)", ylab = "Detection Probability")
lines(xseq, g.cds2d, col = "black", lty = "dashed", lwd = 1.5)
lines(xseq, g.mds2d, col = "red", lty = "dotted", lwd = 1.5)

int1 <- sum(g.cds1d) * 0.01
int2 <- sum(g.cds2d) * 0.01
plot(xseq, g.cds2d/int1, type = "l", lwd = 1.5, col = "blue", xlab = "Perpendicular Distance", ylab = "Detection Probability", ylim = c(0.1, 0.35))
lines(xseq, g.cds1d/int2, col = "red", lwd = 1.5)

