library(tidyverse)
library(parallel)
library(scales)
library(ggplot2)
library(tidyr)
library(matrixStats)
library(grid)
library(gridExtra)

resPath1 <- "/data/daotran/Cancer_RP/Subtyping/BiB_Subtyping/Results/allbutClin"
datasets1 <- list.files(resPath1)
datasets1 <- datasets1[grep("DSCC-", datasets1)]
datasets1 <- gsub("DSCC-", "", datasets1)
datasets1 <- strsplit(datasets1, ".rds")
datasets1 <- lapply(datasets1, function(elm) { elm[1] }) %>% unlist()

resPath2 <- "/data/daotran/Cancer_RP/Subtyping/BiB_Subtyping/Results/Mtb_results"
datasets2 <- list.files(resPath2)
datasets2 <- datasets2[grep("DSCC-", datasets2)]
datasets2 <- gsub("DSCC-", "", datasets2)
datasets2 <- strsplit(datasets2, ".rds")
datasets2 <- lapply(datasets2, function(elm) { elm[1] }) %>% unlist()
datasets2 <- setdiff(datasets2, c("P23236214", "P25091696"))

alldatasets <- union(datasets1, datasets2)
# alldatasets <- datasets2

# datasets <- datasets[1:32]
# datasets <- setdiff(datasets, c("P38007532"))
# datasets <- datasets[2:6]

# methods <- c("NEMO", "PINSPlus", "CC", "SNF", "CIMLR", "iCB")
# methods <- c("NEMO", "CC", "CIMLR")
# methods <- c("DSCC", "CC", "CIMLR", "NEMO", "SNF", "LRACluster", "MC", "MCCA", "IntNMF", "ANF")

# DSCC_res <- readRDS("/data/daotran/Cancer_RP/Subtyping/Metabolite_Analysis/Tmp/DSCC_TCGA.rds")
# colnames(DSCC_res) <- c("pvalue_cox", "c_index", "ncluster", "nSamples", "nDeath")

# methods <- c("DSCC", "CC", "CIMLR", "SNF", "LRACluster", "IntNMF", "ANF", "NEMO")
methods <- c("DSCC", "CC", "CIMLR", "SNF", "LRACluster", "IntNMF", "ANF")
# methods <- "DSCC"
# methods <- c("DSCC.sig02", "DSCC.sig05", "DSCC.sigmed")


# methods <- c("DSCC", "NEMO")
allRes_methods <- lapply(methods, function(method) {

  allRes <- lapply(alldatasets, function(dataset) {
    if (dataset %in% datasets1) {
      ResList <- readRDS(file.path(resPath1, paste0(method, "-", dataset, ".rds")))
    }else {
      ResList <- readRDS(file.path(resPath2, paste0(method, "-", dataset, ".rds")))
    }
    ResList
  }) %>% `names<-`(alldatasets)

  ResTable <- lapply(allRes, function(ResList) {
    data.frame(pvalue_cox = ResList$pval_cox,
               # pvalue_fisher = ResList$pval_lr_fisher,
               # pvalue_stf = ResList$pval_lr_stf,
               c_index = ResList$concordance,
               ncluster = ResList$ncluster,
               nSamples = ResList$noSamples,
               nDeath = ResList$nDeaths)
  }) %>%
    do.call(what = rbind) %>%
    `rownames<-`(names(allRes)) %>%
    as.data.frame()

  ResTable
}) %>% `names<-`(methods)

pva_cox_tab <- lapply(allRes_methods, function(table) {
  table$pvalue_cox
}) %>%
  do.call(what = cbind) %>%
  `colnames<-`(names(allRes_methods)) %>%
  `rownames<-`(rownames(allRes_methods[[1]]))
# pva_cox_tab[is.na(pva_cox_tab)] <- 0.5
nSigs <- colSums(pva_cox_tab < 0.05, na.rm = T)
# colSums(pva_cox_tab < 0.05)

minus_pvalue <- -log10(pva_cox_tab)
# minus_pvalue[!is.finite(minus_pvalue)] <- NA
# minus_pvalue[!is.finite(minus_pvalue)] <- max(minus_pvalue, na.rm = T) + 9

# pva_cox_tab <- -log10(pva_cox_tab)

cindex_tab <- lapply(allRes_methods, function(table) {
  table$c_index
}) %>%
  do.call(what = cbind) %>%
  `colnames<-`(names(allRes_methods)) %>%
  `rownames<-`(rownames(allRes_methods[[1]]))
# cindex_tab[is.na(cindex_tab)] <- 0.5
# cindex_tab <- cindex_tab[c(1:6, 9:12), ]

## number of clusters
nclus_tab <- lapply(allRes_methods, function(table) {
  table$ncluster
}) %>%
  do.call(what = cbind) %>%
  `colnames<-`(names(allRes_methods)) %>%
  `rownames<-`(rownames(allRes_methods[[1]]))

# for (i in 1:ncol(nclus_tab)){
#   tmp <- nclus_tab[, i]
#   tmp[is.na(tmp)] <- mean(tmp, na.rm = T)
#   nclus_tab[, i] <- tmp
# }

format_scientific_vector <- function(nums, digits = 2) {
  formatted <- sapply(nums, function(num) {
    fmt <- format(num, scientific = TRUE, digits = digits)
    fmt <- gsub("e", "E", fmt)
    return(fmt)
  })
  return(formatted)
}

pval_nclus_tab <- lapply(1:ncol(nclus_tab), function(i) {
  pval_col <- pva_cox_tab[, i]
  pval_col[!is.na(pval_col)] <- format_scientific_vector(pval_col[!is.na(pval_col)])
  pval_col[is.na(pval_col)] <- "NA"
  nclus_col <- nclus_tab[, i]
  nclus_col[is.na(nclus_col)] <- "NA"
  cb <- paste0(pval_col, " (", nclus_col, ")")
}) %>%
  do.call(what = cbind) %>%
  `colnames<-`(colnames(nclus_tab)) %>%
  `rownames<-`(rownames(nclus_tab))

####

# minus log-pval plot
# {
#   # Assuming nSigs is a vector with values for each method
#   # Make sure nSigs has names that match the method names
#
#   allValue <- pivot_longer(as.data.frame(minus_pvalue),
#                            cols = everything(),
#                            names_to = "Method",
#                            values_to = "Value")
#   allValue <- hablar::retype(allValue)
#   allValue$Value <- as.numeric(allValue$Value)
#   allValue$Method <- factor(allValue$Method, levels = methods)
#
#   # Calculate median values for each method to position the text
#   medians <- aggregate(Value ~ Method, data = allValue, FUN = median)
#
#   # Create a data frame for the nSigs values
#   # Ensure nSigs has the same order as methods
#   nSigs_df <- data.frame(
#     Method = factor(names(nSigs), levels = methods),
#     nSigs = nSigs
#   )
#
#   # Merge the medians and nSigs data frames
#   annotations <- merge(medians, nSigs_df, by = "Method")
#
#   # out_name <- "minuslog_pval_TCGA"
#   # pdf(sprintf("/data/daotran/Cancer_RP/Subtyping/Metabolite_Analysis/Plots/%s.pdf", out_name), width = 5, height = 4)
#   plt_pv <- ggplot(allValue, aes(x = Method, y = Value, fill = Method)) +
#     theme_classic() +
#     # scale_fill_npg() +
#     scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF",
#                                  "#E64B357F", "#4DBBD57F", "#00A0877F", "#3C54887F", "#F39B7F7F", "#8491B47F", "#91D1C27F", "#DC00007F", "#7E61487F", "#B09C857F")) +
#     labs(x = "", y = "-log10(p-value)", color = "black", size = 14) +
#     geom_boxplot(outlier.shape = NA, width = 0.4, position = position_dodge(width = 0.1)) +
#     # Add the horizontal line at y=1.3
#     geom_hline(yintercept = 1.3, linetype = "dashed") +
#     # Add text labels for nSigs values at the median for each method
#     geom_text(data = annotations,
#               aes(x = Method, y = Value, label = sprintf("%d", medians)),
#               vjust = -0.5, size = 3, fontface = "bold") +
#     theme(axis.text.x = element_text(angle = 35, vjust = 0.5, hjust = 0.4, size = 9, color = "black"),
#           axis.text.y = element_text(size = 9, color = "black"),
#           axis.title.y = element_text(size = 10, color = "black", face = "bold")) +
#     theme(legend.position = "none") +
#     scale_y_continuous(limits = c(0,7), oob = rescale_none, breaks = seq(0, 7, by=1))
#
#   print(plt_pv)
#   # dev.off()
#   # print(plt)
# }

{
  # Assuming nSigs is a vector with values for each method
  # Make sure nSigs has names that match the method names
  allValue <- pivot_longer(as.data.frame(minus_pvalue),
                           cols = everything(),
                           names_to = "Method",
                           values_to = "Value")
  allValue <- hablar::retype(allValue)
  allValue$Value <- as.numeric(allValue$Value)
  allValue$Method <- factor(allValue$Method, levels = methods)

  # Calculate median values for each method to position the text
  medians <- aggregate(Value ~ Method, data = allValue, FUN = median)


  # Create a data frame for the median values
  annotations <- medians  # Now using just the medians dataframe

  # out_name <- "minuslog_pval_TCGA"
  # pdf(sprintf("/data/daotran/Cancer_RP/Subtyping/Metabolite_Analysis/Plots/%s.pdf", out_name), width = 5, height = 4)
  plt_pv <- ggplot(allValue, aes(x = Method, y = Value, fill = Method)) +
    theme_classic() +
    # scale_fill_npg() +
    scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF",
                                 "#E64B357F", "#4DBBD57F", "#00A0877F", "#3C54887F", "#F39B7F7F", "#8491B47F", "#91D1C27F", "#DC00007F", "#7E61487F", "#B09C857F")) +
    labs(x = "", y = "-log10(p-value)", color = "black", size = 14) +
    geom_boxplot(outlier.shape = NA, width = 0.4, position = position_dodge(width = 0.1)) +
    # Add the horizontal line at y=1.3
    geom_hline(yintercept = 1.3, linetype = "dashed") +
    # Add text labels for median values above the median for each method
    geom_text(data = annotations,
              aes(x = Method, y = Value, label = sprintf("%.2f", Value)),  # Display the median value with 2 decimal places
              vjust = -0.5, size = 3, fontface = "bold") +
    theme(axis.text.x = element_text(angle = 35, vjust = 0.5, hjust = 0.4, size = 9, color = "black"),
          axis.text.y = element_text(size = 9, color = "black"),
          axis.title.y = element_text(size = 10, color = "black", face = "bold")) +
    theme(legend.position = "none") +
    scale_y_continuous(limits = c(0, 4), oob = rescale_none, breaks = seq(0, 5, by = 1))
  print(plt_pv)
  # dev.off()
  # print(plt)
}
ggsave("/data/daotran/Cancer_RP/Subtyping/BiB_Subtyping/Plot/allbutClin_minpval.pdf", plot = plt_pv, width = 6, height = 3)


colMeans(minus_pvalue)
colMedians(minus_pvalue)

#cindex plot
{
  allValue <- pivot_longer(as.data.frame(cindex_tab),
                           cols = everything(),
                           names_to = "Method",
                           values_to = "Value")
  allValue <- hablar::retype(allValue)
  allValue$Value <- as.numeric(allValue$Value)
  allValue$Method <- factor(allValue$Method, levels = methods)

  # Calculate statistics for each method
  mean_values <- aggregate(Value ~ Method, data = allValue, FUN = mean)
  median_values <- aggregate(Value ~ Method, data = allValue, FUN = median)
  sd_values <- aggregate(Value ~ Method, data = allValue, FUN = sd)

  # Combine statistics
  stats_values <- mean_values
  stats_values$median <- median_values$Value
  stats_values$sd <- sd_values$Value

  # out_name <- "cindex_TCGA_tmp"
  # pdf(sprintf("/data/daotran/Cancer_RP/Subtyping/Metabolite_Analysis/Plots/%s.pdf", out_name), width = 5, height = 4)
  plt_ci <- ggplot(allValue, aes(x = Method, y = Value, fill = Method)) +
    theme_classic() +
    # scale_fill_npg() +
    scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF",
                                 "#E64B357F", "#4DBBD57F", "#00A0877F", "#3C54887F", "#F39B7F7F", "#8491B47F", "#91D1C27F", "#DC00007F", "#7E61487F", "#B09C857F")) +
    labs(x = "", y = "C-Index", color = "black", size = 14) +
    geom_boxplot(outlier.shape = NA, width = 0.4, position = position_dodge(width = 0.1)) +
    # Add mean values above the median
    # geom_text(data = stats_values,
    #           aes(x = Method, y = median, label = sprintf("%.3f", Value)),
    #           position = position_dodge(width = 0.1),
    #           vjust = -1.0, size = 3) +
    # Add standard deviation below the median
    # geom_text(data = stats_values,
    #           aes(x = Method, y = median, label = sprintf("%.3f", sd)),
    #           position = position_dodge(width = 0.1),
    #           vjust = 2.0, size = 3) +
    theme(axis.text.x = element_text(angle = 35, vjust = 0.5, hjust = 0.4, size = 09, color = "black"),
          axis.text.y = element_text(size = 09, color = "black"),
          axis.title.y = element_text(size = 10, color = "black", face = "bold")) +
    geom_hline(yintercept = 1.3, linetype = "dashed") +
    theme(legend.position = "none") +
    scale_y_continuous(limits = c(0.4, 0.8), oob = rescale_none, breaks = c(0.4, 0.6, 0.8))
  print(plt_ci)
  # dev.off()
  # print(plt)
}

colMeans(cindex_tab, na.rm = T)
colMedians(cindex_tab, na.rm = T)

# nclus plot
{
  allValue <- pivot_longer(as.data.frame(nclus_tab),
                           cols = everything(),
                           names_to = "Method",
                           values_to = "Value")

  allValue <- hablar::retype(allValue)
  allValue$Value <- as.numeric(allValue$Value)
  allValue$Method <- factor(allValue$Method, levels = methods)

  # out_name <- "minuslog_pval_TCGA"

  # pdf(sprintf("/data/daotran/Cancer_RP/Subtyping/Metabolite_Analysis/Plots/%s.pdf", out_name), width = 5, height = 4)

  plt_nc <- ggplot(allValue, aes(x = Method, y = Value, fill = Method)) +
    theme_classic() +
    # scale_fill_npg() +
    scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF",
                                 "#E64B357F", "#4DBBD57F", "#00A0877F", "#3C54887F", "#F39B7F7F", "#8491B47F", "#91D1C27F", "#DC00007F", "#7E61487F", "#B09C857F")) +
    labs(x = "", y = "cluster number", color = "black", size = 14) +
    geom_boxplot(outlier.shape = NA, width = 0.4, position = position_dodge(width = 0.1)) +
    theme(axis.text.x = element_text(angle = 35, vjust = 0.5, hjust = 0.4, size = 09, color = "black"),
          axis.text.y = element_text(size = 09, color = "black"),
          axis.title.y = element_text(size = 10, color = "black", face = "bold")) +
    # geom_hline(yintercept = 1.3, linetype = "dashed") +
    theme(legend.position = "none") +
    scale_y_continuous(limits = c(1, 10), oob = rescale_none, breaks = seq(1, 10, by = 1))

  print(plt_nc)
  # dev.off()
  # print(plt)
}
colMeans(nclus_tab, na.rm = T)
colMedians(nclus_tab, na.rm = T)

combined_plot <- grid.arrange(plt_pv, plt_ci, plt_nc, ncol = 1)

ggsave("/data/daotran/Cancer_RP/Subtyping/BiB_Subtyping/Plot/allbutClin_alldts.pdf", plot = combined_plot, width = 6, height = 11)



