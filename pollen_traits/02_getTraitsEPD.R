library(dplyr)

## Trait values by pvar
ptraitvals = read.csv("pvarTraitVals.csv")
## Taxonomic groups assignments
pgrp = read.csv("p_group.csv")
## Sites list to assign traits
sites = read.csv("entityDetails_v3.csv")

## Output data frame
out.df = NULL

## Loop through sites
filedir = "~/Dropbox/Data/epd/mapping/maps2010/csvfiles/"
filelist = list.files(filedir, pattern = "csv")
nsites = length(filelist)

## Log file
if (file.exists("taxontrait.log")) {
  file.remove("taxontrait.log")
}

log_con <- file("taxontrait.log", open="a")

missing_taxa = NULL

for (i in 1:nsites) {
# for (i in 1:10) {
    
  myent = strsplit(filelist[i], "_")[[1]][1]
  print(paste("Doing", myent, filelist[i]))
  cat(paste("Doing", myent, filelist[i]), file = log_con, sep="\n")
  siteID = which(sites$E. == myent)
  
  dat = read.csv(paste0(filedir,filelist[i]))
  poll = dat %>%
    select(-seq(1,7), -sum)
  pvars = as.numeric(substr(names(poll), 2, 10))
  nvars = length(pvars)
  
  ## Make up temporary ID list
  traitID = rep(NA, nvars)
  grps = rep(0, nvars)
  for (j in 1:nvars) {
    tmpID = which(ptraitvals$pvar == pvars[j])
    grpID = which(pgrp$Var == pvars[j])
    if (length(tmpID) != 1) { 
      # stop("HELP!!!")
      print(paste("Missing taxon:", pvars[j]))
      cat(paste("Missing taxon:", pvars[j]), file = log_con, sep="\n")
      missing_taxa = c(missing_taxa, pvars[j])
    } else {
      traitID[j] = tmpID
      grps[j] = pgrp$Group[grpID]
      # print(paste(pgrp$Group[grpID], pgrp$VarName[grpID]))
    }
  }

  ## Vectorize the pollen
  nsamp = dim(dat)[1]
  poll.df = data.frame(ent = rep(myent, nsamp * nvars),
                       sample = rep(seq(1:nsamp), each = nvars),
                       agebp = rep(dat$agebp, each = nvars),
                       lon = rep(sites$LonDD[siteID], nsamp * nvars),
                       lat = rep(sites$LatDD[siteID], nsamp * nvars),
                       alt = rep(sites$Elevation[siteID], nsamp * nvars),
                       pvar = rep(pvars, nsamp), 
                       val = c(t(as.matrix(poll))),
                       grp = rep(grps, nsamp),
                       ptraitvals[traitID,-c(1,2,3)], row.names = NULL)
  
  # write.csv(poll.df, "test.csv", row.names = FALSE)
  ## Cut samples with 0 abundance
  poll.df = poll.df %>%
    filter(val > 0)
  
  ## Write out
  if (i == 1) {
    write.table(poll.df, "epdPollTraitVals.csv", sep = ",", row.names = FALSE)
  } else {
    write.table(poll.df, "epdPollTraitVals.csv", sep = ",", 
                col.names = FALSE, row.names = FALSE, append = TRUE)
  }
  # ## Now concatentate
  # out.df = rbind(out.df, poll.df)
  
}
close(log_con)
