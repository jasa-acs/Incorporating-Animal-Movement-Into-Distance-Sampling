## line transect simulation study
library(moveds)
library(Distance)
set.seed(15032)

# Set truth ---------------------------------------------------------------
# scale, shape, diffusion rates
true.par <- c(5, 3, 0)
diffusions <- seq(0.5, 4.0, 0.5)
ndiff <- length(diffusions)
hzfn <- 1
hb.par <- s2sigmab(true.par[1], true.par[2])
# auxiliary information: region width, region length, half width, observer_speed, transect type
aux <- c(1000, 1000, 30, 1.0, 0)
# simulation time step
sim.dt <- 1
# true abundance
true.N <- 100
# simulate movement or not?
sim.move <- 1
# movement observations
num.tagged <- 5
observation.times <- seq(0, 1000, 10)
# transects
n.transects <- 50
transdat <- matrix(0, nr = n.transects, nc = 2)
transdat[,1] <- 1:n.transects
transdat[,2] <- rep(1000, n.transects)
# number of simulations 
nsims <- 100

# setup results matrix ----------------------------------------------------

res.cds <- res.mds <- vector(mode = "list", length = ndiff)

# run simulation ----------------------------------------------------------

for (diff in 1:ndiff) {
  true.par[3] <- diffusions[diff]
  # save estimated N, SE, LCL, and UCL 
  res.cds[[diff]] <- matrix(0, nr = nsims, ncol = 4)
  res.mds[[diff]] <- matrix(0, nr = nsims, ncol = 4)
  for (sim in 1:nsims) {
    cat("Simulation ", sim, " / ", nsims, "\n")
    # Run simulation function to output data file
    SimulateDsData(c(rev(hb.par), true.par[3]),
                 true.N,
                 c(1000, 1000, 60, 1000, 1000, 1, n.transects, 0),
                 sim.dt,
                 sim.move)
    
    # Read in observations data  ----------------------------------------------
    obs <- read.csv("simulated_dsdata.csv", header = FALSE)
    names(obs) <- c("transect", "x", "y", "t")
    
    # Simulate movement data  -------------------------------------------------
    movedat <- SimulateMovementData(num.tagged, observation.times, true.par[3])
  
    # CDS model ------------------------------------------------------------
    cat("CDS....")
    # format for mds function
    ds <- list(data = obs,
             transect = transdat,
             aux = aux,
             delta = c(2.5, 1000),
             buffer = 0,
             hazardfn = 1,
             move = 0)
    move <- list(data = movedat)
  
    cds2d <- mds(ds, move, start = c(s = 5, d = 3), print = TRUE)
    res.cds[[diff]][sim, ] <- as.numeric(cds2d$result[3,])
  
    # MDS model ------------------------------------------------------------
    cat("MDS....")
    #ds$move <- 1
    #ds$delta <- c(2.5, 1)
    #mds2d <- mds(ds, move, start = c(s = 5, d = 3, sd = 2.5), print = TRUE)
    #res.mds[[diff]][sim,] <- as.numeric(mds2d$result[4,])
    
    saveRDS(list(cds2d=res.cds, mds2d=res.mds), "line_sim_res.Rds")
    cat("\n")
  }
}

# format results 
cds.bias <- sapply(res.cds, FUN = function(x) {(mean(x[,1]) - 100) / 100})
mds.bias <- sapply(res.mds, FUN = function(x) {(mean(x[,1]) - 100) / 100})

bias <- data.frame(diff = diffusions*sqrt(2/pi), 
                   mds = mds.bias, 
                   cds = cds.bias)

saveRDS(bias, "line_bias.Rds") 

# create Figure 1  
library(ggplot2)
theme_set(theme_bw())
bias.grp <- ggplot(bias, aes(x = diff)) + 
  geom_ribbon(data = data.frame(x = seq(0,3.5,0.1)), aes(x = x, ymin = -5), ymax = 5, fill = "grey80") + 
  geom_segment(x = 0, xend = 0.8, y = 5, yend = 5, col = "grey70", linetype = "dashed", lwd = 1) +
  geom_segment(x = 0.8, xend = 0.8, y = -10, yend = 5, col = "grey70", linetype = "dashed", lwd = 1) +   
  geom_segment(x = 0, xend = 1.55, y = 20, yend = 20, col = "grey70", linetype = "dashed", lwd = 1) +
  geom_segment(x = 1.55, xend = 1.55, y = -10, yend = 20, col = "grey70", linetype = "dashed", lwd = 1) +   
  geom_segment(x = 0, xend = 2.11, y = 50, yend = 50, col = "grey70", linetype = "dashed", lwd = 1) +
  geom_segment(x = 2.11, xend = 2.11, y = -10, yend = 50, col = "grey70", linetype = "dashed", lwd = 1) +  
  geom_segment(x = 0, xend = 3, y = 100, yend = 100, col = "grey70", linetype = "dashed", lwd = 1) +
  geom_segment(x = 3, xend = 3, y = -10, yend = 100, col = "grey70", linetype = "dashed", lwd = 1) + 
  geom_line(aes(y = cds, linetype = "dashed", col = "blue"), lwd = 1.5) + 
  geom_point(aes(y =  cds), cex = 3) + 
  geom_line(aes(y = mds, linetype = "solid", col = "black"), lwd = 1.5) + 
  geom_point(aes(y = mds), col = "blue", cex = 3) + 
  scale_x_continuous("Animal Speed (as % of Observer Speed)", breaks = seq(0.5, 3.5, 0.5), labels = seq(0.5, 3.5, 0.5)*100, limits = c(0.35,3.0)) + 
  scale_y_continuous("Relative Bias (%)", breaks = seq(-10, 100, 10), labels = seq(-10, 100, 10), limits = c(-5, 100)) + 
  scale_linetype_discrete(guide = FALSE, "parameter") + 
  scale_color_manual(guide = FALSE, "model", values = c("blue", "black"), labels = c("moveds", "cds")) +
  geom_text(x = 2.15, y = 65, label = "DS", size = 12, family = "sans") + 
  geom_text(x = 2.5, y = 10, label = "MDS", size = 12, family = "sans", col = "blue") +  
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        strip.text.x = element_blank(), axis.text=element_text(size =20), axis.title=element_text(size = 20)) 
bias.grp



