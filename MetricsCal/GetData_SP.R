library(tidyverse)
library(parallel)

# subPath <- "/data/daotran/Cancer_RP/Subtyping/BiB_Subtyping/Results/allbutClin"
# datPath <- "/data/daotran/Cancer_RP/Subtyping/Data/TCGA-mapped-2"
# savePath <- "/data/daotran/Cancer_RP/Subtyping/BiB_Subtyping/Data/SubSurvClin"

subPath <- "/data/daotran/Cancer_RP/Subtyping/BiB_Subtyping/Results/Mtb_results"
datPath <- "/data/daotran/Cancer_RP/Subtyping/Metabolite_Analysis/Data/ProcessedData_Map"
savePath <- "/data/daotran/Cancer_RP/Subtyping/BiB_Subtyping/Data/SubSurvClin"

datasets <- list.files(datPath)
datasets <- strsplit(datasets, ".rds")
datasets <- lapply(datasets, function(elm) { elm[1] }) %>% unlist()
datasets <- setdiff(datasets, c("P23918603", "P38007532", "P25091696"))
# datasets <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-OV")

# allstat <- lapply(datasets, function(file){
#   rds <- readRDS(file.path(datPath, paste0(file, ".rds")))
#   names(rds)
# })

# methods <- c("DSCCM.sig05", "CC", "CIMLR", "SNF", "NEMO", "LRACluster", "IntNMF", "ANF")
methods <- c("nosubtype", "DSCC", "CC", "CIMLR", "SNF", "NEMO", "LRACluster", "ANF", "IntNMF")
# methods <- c("DSCC")


lapply(methods, function(method){
  print(method)
  if (!file.exists(file.path(savePath, method))){
    dir.create(file.path(savePath, method))
  }

  mclapply(datasets, mc.cores = 16, function(dataset){
    print(dataset)
    clindat <- readRDS(file.path(datPath, paste0(dataset, ".rds")))
    survival <- clindat$survival
    clindat <- clindat$clinicalImputedV2

     # colkeeps <- c("age_at_initial_pathologic_diagnosis", "race", "gender")
    # keeps <- lapply(colnames(clindat), function(name){
    #   countmatch <- sapply(colkeeps, function(col) {
    #     length(grep(col, name))
    #   }) %>% sum()
    #   if(countmatch == 0){
    #     return(NULL)
    #   }else{
    #     return(name)
    #   }
    # }) %>% do.call(what = c) %>% unique()
    # keeps <- keeps[!is.null(keeps)]
    #
    # clindat <- as.data.frame(as.matrix(clindat[, keeps])) %>% `rownames<-` (rownames(clindat))

    # cSamples <- intersect(rownames(survival), rownames(clindat))

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


