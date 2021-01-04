## moveds example
library(moveds)

# Read in Data ------------------------------------------------------------

spotted <- readRDS("data/spotted_dolphin/moveds_dat/spotted.rds")
aux <- c(1000,1000, 5.5, mean(spotted$transect[,6]), 0)
obs <- as.data.frame(spotted$enc)[,-1]
transdat <- spotted$transect
n.transects <- nrow(transdat)
distances <- abs(spotted$enc[,3])
movedat <- spotted$tags


# read in models ----------------------------------------------------------

mods <- readRDS("results/spotted_dolphin/final_mods.Rds")

# 1D CDS model ------------------------------------------------
library(Distance)
names(obs) <- c("transect", "x", "y", "t")
yrs <- unique(transdat[,8])
obs$year <- 1 
for (i in seq(nrow(obs))) obs$year[i] <- transdat[obs$transect[i]==transdat[,1], 8]
obs$len <- obs$sp <- rep(0, nrow(obs))
for (i in 1:nrow(obs)) {
  len <- transdat[which(transdat[,1] == obs$transect[i]), 2]
  sp <- transdat[which(transdat[,1]==obs$transect[i]), 6]
  obs$len[i] <- len
  obs$sp[i] <- sp
}
obs <- obs[obs$y + obs$t*obs$sp <= obs$len,]
region.table <- data.frame(Region.Label = yrs, Area = prod(aux[1:2]))
sample.table <- data.frame(Sample.Label = transdat[,1], Region.Label = transdat[,8],
                           Effort = transdat[, 2])
obs.table <- data.frame(object = 1:nrow(obs), Region.Label = obs$year, Sample.Label = obs[, 1])
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
# get estimated density
N.cds1d <- as.numeric(cds1d$dht$individuals$N$Estimate)[-5]
cds1d.byyr <- cds1d$dht$individuals$N

# 2D CDS model ------------------------------------------------------------
cds2d <- mods$cds2d
cds2d.byyr <- vector(mode = "list", length = 4)
trans <- cds2d$ds$transect
for (i in 1:4) {
  subtrans <- trans[transdat[,8] == yrs[i],]
  cds2d.byyr[[i]] <- predict(cds2d, subtrans, yrs[i])
  cat("i = ", i, " N = ", cds2d.byyr[[i]]$result[3,1], "\n")
}


# 2D MDS model ------------------------------------------------------------
mds2d <- mods$mds2d
mds2d.byyr <- vector(mode = "list", length = 4)
trans <- mds2d$ds$transect
for (i in 1:4) {
  subtrans <- trans[transdat[,8] == yrs[i],]
  mds2d.byyr[[i]] <- predict(mds2d, subtrans, yrs[i])
}

saveRDS(list(cds1d = cds1d.byyr, cds2d = cds2d.byyr, mds2d = mds2d.byyr), "byyear.Rds")
