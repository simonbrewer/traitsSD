###############################################################################	
#
# Script to calculate hypervolumes from fossil pollen samples
#
# This ONLY calculates the HV and does a set of random resamples
# HVs are then stored for alpha and beta calcs
###############################################################################	

start.time = Sys.time()

library(hypervolume)
library(dplyr)
library(foreach)
library(parallel)

# cl <- parallel::makeCluster(4, setup_timeout = 0.5)
# doParallel::registerDoParallel(cl)

randomTraits <- function(n = 1, means, sds) {
  if (length(means) != length(sds)) {stop("Means and SDs not equal")}
  ntraits = length(means)
  trait.out = rep(NA, ntraits)
  
  for (i in 1:ntraits) {
    if (sds[i] > 0 & !is.na(sds[i])) {
      trait.out[i] <- rnorm(n, means[i], sds[i])
    } else {
      trait.out[i] <- means[i]
    }
  }
  return(trait.out)
}

## First get the trait values
load("../../pollen_traits/epdPollTraitVals.RData")

## Threshold for taxa selection (5%)
# tth = 0.005
tth = read.csv("../../pollen_traits/pvarCutOff.csv")

## Means and standard deviations
trait.mn = c(mean(dat$lsla.mn, na.rm = TRUE), 
             mean(dat$lhgt.mn, na.rm = TRUE), 
             mean(dat$lsdm.mn, na.rm = TRUE))

trait.sd = c(sd(dat$lsla.mn, na.rm = TRUE), 
             sd(dat$lhgt.mn, na.rm = TRUE), 
             sd(dat$lsdm.mn, na.rm = TRUE))

## Keep only 15000 years
dat = dat %>% 
  filter(agebp >= -50, agebp <= 15250)

## Missing values
dat = dat %>%
  filter(!is.na(lsla.mn) & !is.na(lhgt.mn) & !is.na(lsdm.mn))

## All terrestrial taxa
dat = dat %>% 
  filter(grp=="TRSH" | grp=="LIAN" | grp=="DWAR" | grp=="HERB")

## Binages (crude approach)
tslice = seq(-250,15250,by=500)
tmid = tslice[-length(tslice)] + 250
ntslice = length(tmid)

dat$age.bin = as.numeric(cut(dat$agebp, tslice) )
dat$age.bin = (dat$age.bin * 500) - 500

sites = unique(sort(dat$ent))
nsites = length(sites)

## Resampling iterations
nbit = 100

## Now loop through by site
for (i in 1:2) {
  # for (i in 1:2) {
  
  ## Output
  out.hv = list()
  out.rnd = list()
  out.crds = data.frame(lon=dat$lon[i],
                        lat=dat$lat[i],
                        alt=dat$alt[i])
  trait.mean = array(NA, dim = c(ntslice, 3))
  trait.rnd = array(NA, dim = c(ntslice, 3, nbit))
  
  site.dat = dat %>% 
    filter(ent == sites[i])
  
  ## Loop to remove low values
  nsd = dim(site.dat)[1]
  keepID = NULL
  
  for (j in 1:nsd) {
    if (site.dat$val[j] > tth[which(tth$pvar == site.dat$pvar[j]),]$cutoff) {
      keepID = c(keepID,j)
    }
  }
  # site.dat = site.dat[which(site.dat$val >= tth),]
  site.dat = site.dat[keepID,]
  
  for (j in 1:ntslice) {
  # for (j in 1:2) {
    print(paste("Doing",i,j))
    
    time.tmp.hv = list()
    time.dat = site.dat %>% 
      filter(age.bin == tmid[j]) 
    
    trait.dat = data.frame(lsla = time.dat$lsla.mn,
                           lhgt = time.dat$lhgt.mn,
                           lsdm = time.dat$lsdm.mn)
    trait.dat = unique(trait.dat)
    trait.mean[j,] <- apply(trait.dat, 2, mean, na.rm = TRUE)
    if (dim(trait.dat)[1] >= 5) {
      
      # samp.hv = hypervolume_svm(trait.dat,
      #                            samples.per.point = 1000, verbose = FALSE)
      # samp.hv = hypervolume_gaussian(trait.dat,
      #                                 kde.bandwidth = 0.5,
      #                                 quantile.requested = 0.95,
      #                                 quantile.requested.type = "probability",
      #                                 samples.per.point = 1000, verbose = FALSE)
      out.hv[[j]] <- hypervolume_box(trait.dat,
                                     samples.per.point = 1000, verbose = FALSE)
      
      ## Resampling loop - runs in parallel using foreach 
      ## Set up cores above ^
      # out.rnd[i,j,] <- foreach(k=1:nbit, .combine = 'c') %do% {
      oper <- foreach (k=1:nbit, .combine = 'c', .packages='hypervolume') %do% {
        
        # for (k in 1:nbit) {
        ## Get random traits
        trait.dat = data.frame(lsla = randomTraits(means = time.dat$lsla.mn,
                                                   sds = time.dat$lsla.sd),
                               lhgt = randomTraits(means = time.dat$lhgt.mn,
                                                   sds = time.dat$lhgt.sd),
                               lsdm = randomTraits(means = time.dat$lsdm.mn,
                                                   sds = time.dat$lsdm.sd))
        trait.rnd.tmp <- apply(trait.dat, 2, mean, na.rm = TRUE)
        #   # samp.hv = hypervolume_svm(trait.dat,
        #   #                           samples.per.point = 1000, verbose = FALSE)
        #   
        #   # samp.hv = hypervolume_gaussian(trait.dat,
        #   #                                kde.bandwidth = 0.5,
        #   #                                quantile.requested = 0.95,
        #   #                                quantile.requested.type = "probability",
        #   #                                samples.per.point = 1000, verbose = FALSE)
        #   # 
        tmp.hv <- hypervolume_box(trait.dat,
                                       samples.per.point = 1000, verbose = FALSE)
        #   
        #   # out.rnd[i,j,k] = samp.hv@Volume
        #   # samp.hv@Volume
        list(tmp.hv, trait.rnd.tmp)
        #   
      } # Resample loop
      
      ## Reformat output
      hvpos = seq(1, nbit*2, by = 2)
      trpos = seq(2, nbit*2, by = 2)
      for (k in 1:nbit) {
        time.tmp.hv[[k]] <- oper[[hvpos[k]]] 
        trait.rnd[j,,k] <- oper[[trpos[k]]] 
      }
      out.rnd[[j]] <- time.tmp.hv
    } else {
      out.rnd[[j]] <- 0
    } # Sample check

  } # Time loop
  save(out.crds, out.hv, out.rnd, trait.mean, trait.rnd, 
       file = paste0("./hv_out/site_",sites[i],".RData"))

} # Site loop

end.time = Sys.time()

# parallel::stopCluster(cl)

library(RPushbullet)
pbPost("note","Desktop", paste("calc_hvs: Finished in\n",
                               format(end.time - start.time)))
