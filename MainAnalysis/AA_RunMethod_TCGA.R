### Use only the one distance among: angular, euclidean, correlation
source("/data/daotran/Cancer_RP/Subtyping/BiB_Subtyping/Code_Test_Mtb/AA_Methodhelper_test.R")

library(parallel)
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)
Sys.setenv(OMP_NUM_THREADS = 1,
           OPENBLAS_NUM_THREADS = 1,
           MKL_NUM_THREADS = 1,
           VECLIB_MAXIMUM_THREADS = 1,
           NUMEXPR_NUM_THREADS = 1)
library(tidyverse)
# data_root <- "/nfs/blanche/share/daotran/Subtyping/data-analysis/TCGA/"
# data_root <- processedDataPath <- "/nfs/blanche/share/daotran/Subtyping/data-analysis/TCGA/"
# data_root <- processedDataPath <- "/nfs/blanche/share/daotran/SurvivalPrediction/RP_Data_rds/"
data_root <- processedDataPath <- "/data/daotran/Cancer_RP/Subtyping/Data/TCGA-mapped-2/"
pto_path <- "/data/daotran/Cancer_RP/Subtyping/Data/TCGA_Proteomics_map"

datasets <- list.files(data_root)
datasets <- strsplit(datasets, ".rds") %>% do.call(what = c)
# datasets <- setdiff(datasets, c("TCGA-BRCA"))
# datasets <- c("TCGA-BRCA")

# datasets <- c("TCGA-BRCA", "TCGA-UCEC", "TCGA-KIRC", "TCGA-HNSC", "TCGA-LUAD", "TCGA-THCA", "TCGA-PRAD", "TCGA-LGG",
#               "TCGA-OV", "TCGA-LUSC", "TCGA-COAD", "TCGA-STAD", "TCGA-BLCA", "TCGA-LIHC")

# datasets <- c("METABRIC")

numClusters <- length(datasets)
# numClusters <- 8

## Name of the folder to save the result
config <- "NEMO-DIST1-GENE-NETWORK"
# config <- "NEMO-KEGG-PWS-RERUN-v2"

jobs_root <- "/data/daotran/Cancer_RP/Subtyping/Results/NEMO_GN_jobs2/"
jobFilePath <- paste0(jobs_root, config, ".rds")
file.remove(jobFilePath)

## Root folder to save the result folder
tmpR <- "/data/daotran/Cancer_RP/Subtyping/Results/NEMO_GN2/"
tmpRoot <- paste0(tmpR, config, "/")
unlink(tmpRoot, recursive = T)

dir.create(tmpRoot, recursive = TRUE, showWarnings = F)

# dataTypes <- c("meth450", "miRNA", "miRNAiso", "mRNAFPKM", "mRNAFPKMuq", "cnv",
#                "mRNATPM", "mRNAUnstranded")
# dataTypes <- c("mRNA", "miRNA", "meth", "cnv", "SNV", "cna", "mutation", "snv")
# dataTypes <- c("mRNATPM_map", "miRNA_map", "meth450_map")
# dataTypes <- c("mRNATPM", "mRNAFPKM", "mRNAFPKMuq", "mRNAUnstranded", "miRNA", "miRNAiso",
#                "meth450", "cnv", "clinicalMCA", "SNVMuseSignature", "SNVVarscanSignature",
#                "SNVSniperSignature", "SNVMutectSignature")


# dataTypes <- c("mRNATPM", "mRNAFPKM", "mRNAFPKMuq", "mRNAUnstranded", "miRNA", "miRNAiso",
#                "meth450", "cnv", "SNVMuseSignature", "SNVVarscanSignature",
#                "SNVSniperSignature", "SNVMutectSignature", 
#                "PPPTO_Peptides_SpectralCount", "PPPTO_Peptides_AmbigSpectralCount", "PPPTO_PPS_ITRAQ", "PPPTO_PPPeptides_ITRAQ", "PPPTO_Sum_SpectralCounts", "PPPTO_Sum_DistinctPeptides", 
#                "PPPTO_Sum_UnsharedPeptides", 
#                "PTO_Sum_SpectralCounts", "PTO_Sum_DistinctPeptides", "PTO_Sum_UnsharedPeptides", "PTO_Peptides_SpectralCount", "PTO_Peptides_AmbigSpectralCount", "PTO_ITRAQ_UnsharedLogRatio", 
#                "PTO_ITRAQ_LogRatio")

dataTypes <- c("mRNATPM", "mRNAFPKM", "mRNAFPKMuq", "mRNAUnstranded", "miRNA", "miRNAiso",
               "meth450", "cnv", "SNVMuseSignature", "SNVVarscanSignature",
               "SNVSniperSignature", "SNVMutectSignature", "Metabo")

ptoTypes <- c(
              #  "PPPTO_Peptides_SpectralCount", "PPPTO_Sum_SpectralCounts", "PPPTO_Sum_DistinctPeptides", "PPPTO_Sum_UnsharedPeptides", "PPPTO_Peptides_AmbigSpectralCount", 
               "PTO_Sum_SpectralCounts", "PTO_Sum_DistinctPeptides", "PTO_Sum_UnsharedPeptides", "PTO_Peptides_SpectralCount",
               "PTO_Sum_SpectralCounts_JHU", "PTO_Sum_DistinctPeptides_JHU", "PTO_Sum_UnsharedPeptides_JHU", "PTO_Peptides_SpectralCount_JHU"
               )

# dataTypes <- c("mRNATPM", "miRNA", "cnv")

# dataTypes <- c("mRNAFPKM", "miRNA", "meth450")

# aging_Genes <- read.csv("/data/daotran/Cancer_RP/Subtyping/Tmp_Data/Aging-relevant-genes.csv")
# aging_Genes <- aging_Genes$Symbol
# allPathways <- lapply(names(allPathways), function(pw){
#   GS <- allPathways[[pw]]
#   GS <- GS[!GS %in% aging_Genes]
#   GS
# }) %>% `names<-` (names(allPathways))

# allPathways$Agingpw <- aging_Genes
# allPathways <- allPathways[c("path:hsa04110", "path:hsa04114", "path:hsa04210", "path:hsa04216", "path:hsa04217", "path:hsa04115", "path:hsa04218", "path:hsa04211")]

## Final genesets variable must be a list of gene sets/pathways, each element is a list of genes

# keggGeneSet <- readRDS("/data/daotran/Cancer_RP/Test23.08.18/data/KEGGPathways.rds")
# allPathways <- keggGeneSet

genesets <- list("This is not important")

### Functions
if (file.exists(jobFilePath)) {
  jobs <- readRDS(jobFilePath)
  fs <- list.files(tmpRoot, full.names = T)
  jobs <- jobs[!(jobs$savePath %in% fs),]

} else {

  jobs <- expand.grid(datasets = datasets, pwID = genesets, method = "DSCC", stringsAsFactors = FALSE)
  jobs$savePath <- lapply(1:nrow(jobs), function(i) {

    tmpDir <- file.path(tmpRoot, paste0(sample(c(LETTERS, letters), 10, TRUE), collapse = ""))
    paste0(tmpDir, ".rds")
  }) %>% unlist()

  jobs$pwID2 <- NA
  count <- 0
  for (i in 1:nrow(jobs)) {
    if (jobs$datasets[i] == datasets[1]) {
      count <- count + 1
    }
    jobs$pwID2[i] <- count
  }

  saveRDS(jobs, jobFilePath)
}

print(nrow(jobs))


# cl <- makeCluster(numClusters)
# clusterExport(cl, "jobs")
# clusterEvalQ(cl, {

# Wrapper to run each method
suppressMessages({ suppressWarnings({
  library(wordspace)
  library(tidyverse)
  library(survival)
  library(matrixStats)
  library(SNFtool)
  library(igraph)
  library(NEMO)
  library(psych)
  library(parallel)
  library(cluster)
}) })

# processedDataPath <- "/nfs/blanche/share/daotran/Subtyping/data-analysis/TCGA/"

# })

cl <- makeCluster(numClusters, type = 'FORK', outfile = sprintf("%s.txt", config))

runParallel <- function(j) {

  dataset <- jobs[j, 1]
  cat("Loading processed data\n")
  toUseDataList <- readRDS(paste0(processedDataPath, dataset, ".rds"))
  survival <- toUseDataList$survival

  # survival <- survival[, c("time", "status")]
  survival <- survival[, c("os", "isDead")]

  colnames(survival) <- c("OStime", "OSstatus")

  # Get union of sample
  dataType.used <- setdiff(names(toUseDataList), c("survival"))
  toUseDataList <- toUseDataList[dataType.used]
  # toUseDataList <- toUseDataList[grep("mRNA", names(toUseDataList))]

  # keeps <- lapply(names(toUseDataList), function(name){
  #   countmatch <- sapply(dataTypes, function(dataType) {
  #     length(grep(dataType, name))
  #   }) %>% sum()
  #   if(countmatch == 0){
  #     return(NULL)
  #   }else{
  #     return(name)
  #   }
  # }) %>% do.call(what = c) %>% unique()
  # keeps <- keeps[!is.null(keeps)]
  # keeps <- setdiff(keeps, "miRNAiso")
  toUseDataList <- toUseDataList[names(toUseDataList) %in% dataTypes]

  if (file.exists(file.path(pto_path, paste0(dataset, ".rds")))){
    ptodata <- readRDS(file.path(pto_path, paste0(dataset, ".rds")))
    ptodata <- ptodata[names(ptodata) %in% ptoTypes]
    toUseDataList <- c(toUseDataList, ptodata)
  }

  # toUseDataList <- toUseDataList[keeps]

  # toUseDataList <- lapply(names(toUseDataList), function(dataType) {
  #     df <- as.matrix(toUseDataList[[dataType]])
  #     df[is.na(df)] <- 0
  #     # log2(df + 1)
  #     if (max(df) > 100) {
  #       df <- log2(df + 1)
  #     }
  #     df
  #   }) %>% `names<-` (names(toUseDataList))
  #
  # if ("mRNATPM" %in% names(toUseDataList)) {
  #     tmp <- as.matrix(toUseDataList$mRNATPM)
  #     tmp <- tmp[, order(colSds(tmp), decreasing = TRUE)[1:10000]]
  #     toUseDataList$mRNATPM <- tmp
  #     rm(tmp)
  #   }
  #
  # if ("meth450" %in% names(toUseDataList)) {
  #     tmp <- as.matrix(toUseDataList$meth450)
  #     tmp <- tmp[, order(colSds(tmp), decreasing = TRUE)[1:10000]]
  #     toUseDataList$meth450 <- tmp
  #     rm(tmp)
  #   }

  # if ("meth450" %in% names(toUseDataList)){
  #   tmp <- as.matrix(toUseDataList$meth450)
  #   tmp <- tmp[, order(apply(tmp, 2, var), decreasing=TRUE)[1:5000]]
  #   toUseDataList$meth450 <- tmp
  #   rm(tmp)
  # }

  toUseDataList <- toUseDataList[!sapply(toUseDataList, is.null)]

  # toUseDataList <- lapply(names(toUseDataList), function(dataType) {
  #   data <- toUseDataList[[dataType]]
    # if (dataType != "meth450"){
    #   colnames(data) <- paste0("Var", c(1:ncol(data)))
    # }
    # if (dataType %in% c("clinical", "clinicalImputedV2")) {
    #   colkeeps <- c("age_at_initial_pathologic_diagnosis", "race", "gender", "ethnicity", "height", "weight")
    #   keeps <- lapply(colnames(data), function(name) {
    #     countmatch <- sapply(colkeeps, function(col) {
    #       length(grep(col, name))
    #     }) %>% sum()
    #     if (countmatch == 0) {
    #       return(NULL)
    #     }else {
    #       return(name)
    #     }
    #   }) %>%
    #     do.call(what = c) %>%
    #     unique()
    #   keeps <- keeps[!is.null(keeps)]
    #   data <- data[, keeps, drop = FALSE]
    # }
  #   data
  # }) %>% `names<-`(names(toUseDataList))

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

  ###
  toUseDataList <- lapply(names(toUseDataList), function(dataType) {
    data <- toUseDataList[[dataType]]
    data <- as.matrix(data)
    data[is.na(data)] <- 0
  
    if (min (data) >= 0){
      if (max(data) > 100) {
        data <- log2(data + 1)
      }
    }
  
    # if (dataType %in% c("mRNATPM", "miRNA")){
    # if (!dataType %in% c("clinical", "clinicalImputedV2")) {
    #   if (max(data) > 100) {
    #     data <- log2(data + 1)
    #   }
    # }else {
    #   data <- data[, colSds(data) > 0, drop=FALSE]
    #   data <- apply(data, 2, min_max_norm)
    # }
    #   num_cols <- c("age_at_initial_pathologic_diagnosis", "height", "weight")
    #   keeps <- lapply(colnames(data), function(name) {
    #     countmatch <- sapply(num_cols, function(col) {
    #       length(grep(col, name))
    #     }) %>% sum()
    #     if (countmatch == 0) {
    #       return(NULL)
    #     }else {
    #       return(name)
    #     }
    #   }) %>% do.call(what = c) %>% unique()
    #   keeps <- keeps[!is.null(keeps)]
    #   if (length(keeps) > 0) {
    #     for (i in keeps) {
    #       if(max(data[, i]) > 1){
    #         data[, i] <- log2(data[, i] + 1)
    #       }
    #     }
    #   }
    # }
    # data <- data[, colSds(data) > 0]
    # data <- scale(data)
    # }
    data
  }) %>% `names<-`(names(toUseDataList))

  noDeath <- survival[toUseSamples, "OSstatus"]
  noDeath <- sum(noDeath, na.rm = T)
  cat("Finished loading data\n")

  resFile <- jobs[j, 4]
  gs <- unlist(jobs[j, 2])
  cat(gs)
  method <- jobs[j, 3]
  pwID <- jobs[j, 5]
  runningTime <- Sys.time()

  cluster <- try({ runDSCC(dataList = toUseDataList, nSamples = length(toUseSamples)) })
  runningTime <- Sys.time() - runningTime
  if (!is(cluster, "try-error")) {
    nCluster <- max(cluster[toUseSamples])
    coxFit <- try({ summary(coxph(Surv(time = OStime, event = OSstatus) ~ as.factor(cluster[toUseSamples]), data = survival[toUseSamples,], ties = "exact")) })
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
    message(method, "\t", dataset, "\t", pval, "\t", concordance)
    saveRDS(
      list(
        pwID = pwID,
        pwList = gs,
        cluster = cluster,
        ncluster = nCluster,
        pval = pval,
        noSamples = length(toUseSamples),
        noDeath = noDeath,
        noDataTypes = length(toUseDataList),
        concordance = concordance,
        coxFit = coxFit,
        dataset = dataset,
        method = method,
        runningTime = runningTime
      ), file = resFile)
  }else {
    saveRDS(
      list(
        pwID = pwID,
        pwList = gs,
        cluster = NA,
        ncluster = NA,
        pval = NA,
        noSamples = length(unionSamples),
        noDeath = noDeath,
        noDataTypes = length(toUseDataList),
        concordance = NA,
        coxFit = NA,
        dataset = dataset,
        method = method,
        runningTime = runningTime
      ), file = resFile)
  }
  rm(toUseDataList)
  gc()
  return(NULL)
}

tmp <- parallel::parLapply(cl, 1:nrow(jobs), runParallel)

stopCluster(cl)



