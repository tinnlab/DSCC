source("/data/daotran/Cancer_RP/Subtyping/BiB_Subtyping/Code_Test_Mtb/AA_Methodhelper_test.R")
library(data.table)
library(matrixStats)
library(mltools)
library(survminer)
library(metap)
library(KEGGREST)

dataPath <- "/data/daotran/Cancer_RP/Subtyping/Metabolite_Analysis/Data/ProcessedData_Map_2/"
# dataPath <- "/nfs/blanche/share/daotran/Subtyping/data-analysis/TCGA/"

resPath <- "/data/daotran/Cancer_RP/Subtyping/BiB_Subtyping/Results/Mtb_results/"
# resPath <- "/data/daotran/Cancer_RP/Subtyping/Results/QuickTest_TCGA"
# dir.create(resPath, recursive = T)

keggGeneSet <- readRDS("/data/daotran/Cancer_RP/Test23.08.18/data/KEGGPathways.rds")

# aging_Genes <- read.csv("/data/daotran/Cancer_RP/Subtyping/Tmp_Data/Aging-relevant-genes.csv")
# aging_Genes <- aging_Genes$Symbol
# keggGeneSet$Agingpw <- aging_Genes
# cellGDPathways <- c("path:hsa04110", "path:hsa04210", "path:hsa04216", "path:hsa04115", "path:hsa04114", "path:hsa04217",
#                     "path:hsa04218", "path:hsa04211")
# cellGDPathways <- c("path:hsa04110", "path:hsa04210", "path:hsa04216", "path:hsa04115")
# allcellPWs <- c("4110", "4114", "4210", "4216", "4217", "4115", "4218", "4144", "4145", "4142", "4146", "4140", "4137",
#                 "4139", "4148", "4510", "4520", "4530", "4540", "4550", "2024", "5111", "2025", "2026", "2030", "2040",
#                 "4814", "4820", "4810")
# allcellPWs <- paste0("path:hsa0", allcellPWs)
# genePws <- keggGeneSet[names(keggGeneSet) %in% cellGDPathways]
genePws <- keggGeneSet

# metaData <- readRDS("/data/dungp/projects/subtyping/DSCC-subtyping/DSCC-Subtyping/tmp/allMetaData.rds")  %>%
#   dplyr::filter(!(dataType %in% c("survival", "clinical"))) %>%
#   dplyr::select(dataSet, nSample) %>%
#   group_by(dataSet) %>%
#   summarize(nSample = max(nSample)) %>%
#   arrange(nSample)
datasets <- list.files(dataPath)
datasets <- strsplit(datasets, ".rds")
datasets <- lapply(datasets, function(elm) { elm[1] }) %>% unlist()
datasets <- setdiff(datasets, c("P23236214", "P25091696"))
# datasets <- c("ST001237")


# methods <- c("DSCC")

for (dataset in datasets) {
  message("Subtyping on: ", dataset)
  # toUseDataList <- readRDS(paste0(dataPath, dataset, ".rds"))
  rds <- readRDS(file.path(dataPath, paste0(dataset, ".rds")))

  # Get union of sample
  survival <- rds$survival

  # print(length(metaboPws))
  # # print(metaboPws[[1]])
  # allmetabos <- colnames(rds[['Metabo']])
  # # print(allmetabos)
  # mtbkeep <- metaboPws %>% unlist() %>% unique()
  # # print(mtbkeep)
  # percentkeep <- length(intersect(mtbkeep, allmetabos)) / length(allmetabos)
  # # print(intersect())
  # print(percentkeep)

  # nGenes <- lapply(metaboPws, length) %>% do.call(what = c)
  # metaboPws <- metaboPws[which.max(nGenes)[1]]

  # toUseDataList <- rds[setdiff(names(rds), c("survival", "clinical"))]
  toUseDataList <- rds[setdiff(names(rds), c("survival", "clinical", "clinicalImputedV2", "clinicalMCA",
                                              "clinicalImputedV3", "clinicalMCAV2"))]
  toUseDataList <- toUseDataList[lapply(toUseDataList, ncol) %>% do.call(what = c) >= 10]

  # if (dataset == "P23236214"){
  #   toUseDataList <- toUseDataList[c("GeneExp_2", "Metabo")]
  # }

  toUseDataList <- toUseDataList[!sapply(toUseDataList, is.null)]
  unionSamples <- Reduce(c, lapply(toUseDataList, rownames)) %>% unique()
  toUseSamples <- intersect(unionSamples, rownames(survival))

  toUseDataList <- lapply(toUseDataList, function(data){
    sampkeep <- rownames(data)[rownames(data) %in% toUseSamples]
    if (length(sampkeep) < 2){
      data <- NULL
    }else{
      data <- data[sampkeep, ]
    }
    data
  }) %>% `names<-` (names(toUseDataList))

  toUseDataList <- toUseDataList[!sapply(toUseDataList, is.null)]
  unionSamples <- Reduce(c, lapply(toUseDataList, rownames)) %>% unique()
  toUseSamples <- intersect(unionSamples, rownames(survival))

  print(names(toUseDataList))

  rm(rds)
  gc()

  toUseDataList <- lapply(names(toUseDataList), function(dataType) {
    df <- as.matrix(toUseDataList[[dataType]])
    df[is.na(df)] <- 0
    if (min(df) >= 0){
      if (max(df) > 100) {
        df <- log2(df + 1)
      }

      # if (dataType == "clinicalImputedV2"){
      #   df <- log2(df + 1)
      # }

      # if (length(grep("GeneExp", dataType)) > 0){
      # df <- df[, colSds(df) > 0]
      # df <- scale(df)
      # }
      # df <- scale(df)
    }
    df
  }) %>% `names<-` (names(toUseDataList))

  resFile <- paste0(resPath, "DSCC", "-", dataset, ".rds")
  runningTime <- Sys.time()
  cluster <- try({ runDSCC(dataList = toUseDataList, nSamples = length(toUseSamples)) })
  runningTime <- Sys.time() - runningTime
  if (!is(cluster, "try-error")) {
    nCluster <- max(cluster[toUseSamples])
    coxFit <- try({ summary(coxph(Surv(time = os, event = isDead) ~ as.factor(cluster[toUseSamples]), data = survival[toUseSamples,], ties = "exact")) })
    # coxFit <- try({ summary(coxph(Surv(time = OStime, event = OSstatus) ~ as.factor(cluster[toUseSamples]), data = survival[toUseSamples,], ties = "exact")) })
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
      if (length(pval_lr) == 1){
        fisher_pv <- stouffer_pv <- pval_lr
      }else{
        fisher_pv <- sumlog(pval_lr)$p
        stouffer_pv <- sumz(pval_lr)$p
      }
    }else{
      fisher_pv <- NA
      stouffer_pv <- NA
    }

    message("Result for", "\t", dataset, "\t", pval, "\t", fisher_pv, "\t", stouffer_pv, "\t", concordance, "\t", nCluster)
    saveRDS(
      list(
        cluster = cluster,
        nclusters = nCluster,
        pval_cox = pval,
        pval_lr_fisher = fisher_pv,
        pval_lr_stf = stouffer_pv,
        noSamples = length(toUseSamples),
        nDeaths = sum(survival[toUseSamples, "isDead"]),
        noDataTypes = length(toUseDataList),
        concordance = concordance,
        coxFit = coxFit,
        dataset = dataset,
        method = "DSCC",
        runningTime = runningTime
      ), file = resFile)
  } else {
    print("ERROR OCCURED")
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
        method = "DSCC",
        runningTime = runningTime
      ), file = resFile)
  }
  gc()
}
