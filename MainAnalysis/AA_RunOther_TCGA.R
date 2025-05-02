source("/data/daotran/Cancer_RP/Subtyping/BiB_Subtyping/Code/AA_Other_helper.R")
library(data.table)
library(matrixStats)
library(mltools)
library(survminer)
library(metap)

dataPath <- "/data/daotran/Cancer_RP/Subtyping/Data/TCGA_original/"
ptoPath <- "/data/daotran/Cancer_RP/Subtyping/Data/TCGA_Proteomics"
resPath <- "/data/daotran/Cancer_RP/Subtyping/BiB_Subtyping/Results/"
# dir.create(resPath, recursive = T)
datasets <- list.files(dataPath)
datasets <- strsplit(datasets, ".rds")
datasets <- lapply(datasets, function(elm) { elm[1] }) %>% unlist()
# datasets <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-OV")
# datasets <- datasets[33]
# dataTypes <- c("mRNATPM", "mRNAFPKM", "mRNAFPKMuq", "mRNAUnstranded", "miRNA", "miRNAiso",
#                "meth450", "cnv", "SNVMuseSignature", "SNVVarscanSignature",
#                "SNVSniperSignature", "SNVMutectSignature")
dataTypes <- c("mRNATPM", "mRNAFPKM", "mRNAFPKMuq", "mRNAUnstranded", "miRNA", "miRNAiso",
               "meth450", "cnv", "SNVMuseSignature", "SNVVarscanSignature",
               "SNVSniperSignature", "SNVMutectSignature", "Metabo")

ptoTypes <- c("PTO_Sum_SpectralCounts", "PTO_Sum_DistinctPeptides", "PTO_Sum_UnsharedPeptides", "PTO_Peptides_SpectralCount",
               "PTO_Sum_SpectralCounts_JHU", "PTO_Sum_DistinctPeptides_JHU", "PTO_Sum_UnsharedPeptides_JHU", "PTO_Peptides_SpectralCount_JHU")

# fold <- paste(dataTypes, collapse = "_")
# fold <- paste0(fold, "_all")
# resPath <- paste0(resPath, fold,  "/")
resPath <- paste0(resPath, "allbutClin",  "/")

if (!file.exists(resPath)){
  dir.create(resPath)
}

# methods <- c("CC", "SNF", "CIMLR", "iCB", "NEMO", "moCLuster", "PINSPlus", "LRACluster", "MC", "MCCA", "IntNMF", "ANF")
# methods <- c("CC", "SNF", "CIMLR", "NEMO", "LRACluster", "ANF", "IntNMF")
# methods <- c("CC", "SNF", "CIMLR", "IntNMF")
methods <- c("CC")
union_methods <- c("NEMO", "PINSPlus")
intersect_methods <- c("CC", "SNF", "CIMLR", "iCB", "moCLuster", "LRACluster", "MC", "MCCA", "IntNMF", "ANF")
# datasets <- c("P23236214")
# datasets <- setdiff(datasets, c("ST001236", "P38007532"))

for (dataset in datasets) {
  (function() {
    message("Subtyping on: ", dataset)
    # toUseDataList <- readRDS(paste0(dataPath, dataset, ".rds"))
    rds <- readRDS(file.path(dataPath, paste0(dataset, ".rds")))

    # Calculate survival information
    survival <- rds$survival
    survival <- survival[rowSums(is.na(survival)) == 0,]
    survival <- survival[survival$os >= 0,]
    survival$os[survival$os == 0] <- 1

    toUseDataList <- rds[setdiff(names(rds), c("survival"))]

    toUseDataList <- toUseDataList[names(toUseDataList) %in% dataTypes]

    # keeps <- lapply(names(toUseDataList), function(name) {
    #   countmatch <- sapply(dataTypes, function(dataType) {
    #     length(grep(dataType, name))
    #   }) %>% sum()
    #   if (countmatch == 0) {
    #     return(NULL)
    #   }else {
    #     return(name)
    #   }
    # }) %>%
    #   do.call(what = c) %>%
    #   unique()
    # keeps <- keeps[!is.null(keeps)]
    # keeps <- setdiff(keeps, "miRNAiso")
    # toUseDataList <- toUseDataList[names(toUseDataList) %in% dataTypes]
    # toUseDataList <- toUseDataList[keeps]

    # toUseDataList <- toUseDataList[setdiff(names(toUseDataList), "GeneExp_1")]
    # removeid <- grep("Metabo", names(toUseDataList))
    # toUseDataList <- toUseDataList[setdiff(names(toUseDataList), names(toUseDataList)[removeid])]
    toUseDataList <- toUseDataList[!sapply(toUseDataList, is.null)]

    toUseDataList <- lapply(names(toUseDataList), function(dataType) {
      data <- toUseDataList[[dataType]]
      # if (dataType != "meth450"){
      #   colnames(data) <- paste0("Var", c(1:ncol(data)))
      # }
      # if(dataType %in% c("clinical", "clinicalImputedV2")){
      #   colkeeps <- c("age_at_initial_pathologic_diagnosis", "race", "gender", "ethnicity", "height", "weight")
      #   keeps <- lapply(colnames(data), function(name){
      #     countmatch <- sapply(colkeeps, function(col) {
      #       length(grep(col, name))
      #     }) %>% sum()
      #     if(countmatch == 0){
      #       return(NULL)
      #     }else{
      #       return(name)
      #     }
      #   }) %>% do.call(what = c) %>% unique()
      #   keeps <- keeps[!is.null(keeps)]
      #   data <- data[, keeps, drop=FALSE]
      # }
      data
    }) %>% `names<-`(names(toUseDataList))

    if (file.exists(file.path(ptoPath, paste0(dataset, ".rds")))){
      ptodata <- readRDS(file.path(ptoPath, paste0(dataset, ".rds")))
      ptodata <- ptodata[names(ptodata) %in% ptoTypes]
      toUseDataList <- c(toUseDataList, ptodata)
    }

    rm(rds)
    gc()

    # toUseDataList <- lapply(toUseDataList, function(df) {
    #   df <- as.matrix(df)
    #   df[is.na(df)] <- 0
    #   # log2(df + 1)
    #   if (max(df) > 100) {
    #     df <- log2(df + 1)
    #   }
    #   df
    # })

    for (method in methods) {
      (function() {
        message(method, "\t", dataset)

        if (method %in% union_methods) {
          unionSamples <- Reduce(c, lapply(toUseDataList, rownames)) %>% unique()
          toUseSamples <- intersect(unionSamples, rownames(survival))

          toUseDataList <- lapply(toUseDataList, function(data) {
            sampkeep <- rownames(data)[rownames(data) %in% toUseSamples]
            if (length(sampkeep) < 2) {
              data <- NULL
            }else {
              data <- data[sampkeep,]
            }
            data
          }) %>% `names<-`(names(toUseDataList))

          toUseDataList <- toUseDataList[!sapply(toUseDataList, is.null)]
          unionSamples <- Reduce(c, lapply(toUseDataList, rownames)) %>% unique()
          toUseSamples <- intersect(unionSamples, rownames(survival))


        }else {
          # if (dataset == "TCGA-OV") {
          #   # keepid <- grep("Metabo", names(toUseDataList))
          #   toUseDataList <- toUseDataList[setdiff(names(toUseDataList), "meth450")]  ## low number of samples
          # }
          intersectSamples <- Reduce(intersect, lapply(toUseDataList, rownames)) %>% unique()
          # if (length(intersectSamples) < 10){
          #   toUseDataList <- toUseDataList[!names(toUseDataList) %in% ptoTypes]
          #   intersectSamples <- Reduce(intersect, lapply(toUseDataList, rownames)) %>% unique()
          # }
          toUseSamples <- intersect(intersectSamples, rownames(survival))

          toUseDataList <- lapply(toUseDataList, function(data) {
            sampkeep <- rownames(data)[rownames(data) %in% toUseSamples]
            sampkeep <- sort(sampkeep, decreasing = T)
            if (length(sampkeep) < 2) {
              data <- NULL
            }else {
              data <- data[sampkeep,]
            }
            data
          }) %>% `names<-`(names(toUseDataList))

          toUseDataList <- toUseDataList[!sapply(toUseDataList, is.null)]
          intersectSamples <- Reduce(intersect, lapply(toUseDataList, rownames)) %>% unique()
          toUseSamples <- intersect(intersectSamples, rownames(survival))
        }

        ###
        toUseDataList <- lapply(names(toUseDataList), function(dataType) {
          df <- as.matrix(toUseDataList[[dataType]])
          df[is.na(df)] <- 0
          # log2(df + 1)
          if(!dataType %in% c("clinical", "clinicalImputedV2", "clinicalMCA")){
            if (min(df) >= 0){
              if (max(df) > 100) {
                df <- log2(df + 1)
              }
            }
            
            df <- df[, colVars(df) > 0, drop = FALSE]
            df <- df[, order(colVars(df), decreasing = TRUE)[1:min(8000, ncol(df))]]
          }
          df
        }) %>% `names<-`(names(toUseDataList))

        resFile <- paste0(resPath, method, "-", dataset, ".rds")
        runningTime <- Sys.time()
        cluster <- try({ get(paste0("run", method))(toUseDataList) })
        ncluster <- max(cluster[toUseSamples])
        runningTime <- Sys.time() - runningTime
        if (!is(cluster, "try-error")) {
          coxFit <- try({ summary(coxph(Surv(time = os, event = isDead) ~ as.factor(cluster[toUseSamples]), data = survival[toUseSamples,], ties = "exact")) })
          if (!is(coxFit, "try-error")) {
            pval <- round(coxFit$sctest[3], digits = 40)
            # pval <- coxFit$sctest[3]
            concordance <- coxFit$concordance[1]
            coxFit = summary(coxFit)
          } else {
            pval <- NA
            concordance <- NA
            coxFit <- NA
          }
          data_tmp <- as.data.frame(cbind(survival[toUseSamples, c("os", "isDead")], as.factor(cluster[toUseSamples])))
          colnames(data_tmp) <- c("os", "isDead", "cluster")

          pairwise_results <- try({ pairwise_survdiff(Surv(time = os, event = isDead) ~ cluster,
                                                      data = data_tmp,
                                                      p.adjust.method = "BH")  # Benjamini-Hochberg
          })
          if (!is(pairwise_results, "try-error")) {
            pval_lr <- pairwise_results$p.value
            pval_lr <- pval_lr[!is.na(pval_lr)]
            if (length(pval_lr) == 1) {
              fisher_pv <- stouffer_pv <- pval_lr
            }else {
              fisher_pv <- sumlog(pval_lr)$p
              stouffer_pv <- sumz(pval_lr)$p
            }
          }else {
            fisher_pv <- NA
            stouffer_pv <- NA
          }

          message("Result for", "\t", dataset, "\t", pval, "\t", fisher_pv, "\t", stouffer_pv, "\t", concordance, "\t", ncluster)


          saveRDS(
            list(
              cluster = cluster[toUseSamples],
              nclusters = ncluster,
              pval_cox = pval,
              pval_lr_fisher = fisher_pv,
              pval_lr_stf = stouffer_pv,
              noSamples = length(toUseSamples),
              nDeaths = sum(survival[toUseSamples, "isDead"]),
              noDataTypes = length(toUseDataList),
              concordance = concordance,
              coxFit = coxFit,
              dataset = dataset,
              method = method,
              runningTime = runningTime
            ), file = resFile)
        } else {
          saveRDS(
            list(
              cluster = NA,
              nclusters = NA,
              pval_cox = NA,
              pval_lr_fisher = NA,
              pval_lr_stf = NA,
              noSamples = length(toUseSamples),
              nDeaths = sum(survival[toUseSamples, "isDead"]),
              noDataTypes = length(toUseDataList),
              concordance = NA,
              coxFit = NA,
              dataset = dataset,
              method = method,
              runningTime = runningTime
            ), file = resFile)
        }
      })()
      gc()
    }

  })()
  gc()
}
