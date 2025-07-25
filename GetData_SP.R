library(tidyverse)
library(parallel)

subPath <- "./Subtyping_Results"
datPath <- "./Data/DSCC_Main"
savePath <- "./SubSurvClin"

datasets <- list.files(datPath)
datasets <- strsplit(datasets, ".rds")
datasets <- lapply(datasets, function(elm) { elm[1] }) %>% unlist()
datasets <- setdiff(datasets, c("P23918603", "P38007532"))

methods <- c("nosubtype", "DSCC", "CC", "CIMLR", "SNF", "LRACluster", "ANF", "IntNMF")

lapply(methods, function(method){
  print(method)
  if (!file.exists(file.path(savePath, method))){
    dir.create(file.path(savePath, method))
  }

  mclapply(datasets, mc.cores = 4, function(dataset){
    print(dataset)
    clindat <- readRDS(file.path(datPath, paste0(dataset, ".rds")))
    survival <- clindat$survival
    clindat <- clindat$clinicalImputedV2

    if(method == "nosubtype"){
      subtypes <- NULL
    }else{
      subtypes <- readRDS(file.path(subPath, paste0(method, "-", dataset, ".rds")))
      subtypes <- subtypes$cluster

      if(length(subtypes) <= 1){
        subtypes <- NA
      }
    }

    rdsSaved <- list(cluster = subtypes, clinical = clindat, survival = survival)
    saveRDS(rdsSaved, file.path(savePath, method, paste0(dataset, ".rds")))
    return(NULL)
  })
  return(NULL)
})


