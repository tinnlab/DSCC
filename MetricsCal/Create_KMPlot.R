library(gridExtra)
library(tidyverse)
library(survival)
library(survminer)
library(grid)

tmpRoot <- "/data/daotran/Cancer_RP/Subtyping/BiB_Subtyping/Results/allbutClin"
processedDataPath <- "/data/daotran/Cancer_RP/Subtyping/Data/TCGA-mapped-2"

tmpRoot2 <- "/data/daotran/Cancer_RP/Subtyping/BiB_Subtyping/Results/Mtb_results"
processedDataPath2 <- "/data/daotran/Cancer_RP/Subtyping/Metabolite_Analysis/Data/ProcessedData_Map"

allFiles <- list.files(tmpRoot)
method <- "DSCC-"
allFiles <- allFiles[grep(method, allFiles)]
# allFiles <- setdiff(allFiles, c("DSCCM-ST001236.rds", "DSCCM-P38007532.rds", "DSCCM-P33577785.rds",
#                                 "DSCCM-P38395893.rds"))
allFiles <- setdiff(allFiles, c("DSCC-P23236214.rds", "DSCC-P25091696.rds"))

# allFiles <- allFiles[25:33]

allFiles2 <- list.files(tmpRoot2)
method <- "DSCC-"
allFiles2 <- allFiles2[grep(method, allFiles2)]
# allFiles <- setdiff(allFiles, c("DSCCM-ST001236.rds", "DSCCM-P38007532.rds", "DSCCM-P33577785.rds",
#                                 "DSCCM-P38395893.rds"))
allFiles2 <- setdiff(allFiles2, c("DSCC-P23236214.rds", "DSCC-P25091696.rds"))

allFiles3 <- union(allFiles, allFiles2)
allFiles3 <- allFiles3[37:43]

allFiles3 <- c("DSCC-TCGA-ACC.rds", "DSCC-P24316975.rds")


# allFiles <- allFiles[3] #  Adjust ncol for the number of columns
{
  allRes <- lapply(allFiles3, function(file) {
    if (file %in% allFiles){
      readRDS(file.path(tmpRoot, file))
    }else{
      readRDS(file.path(tmpRoot2, file))
    }
  })
  # datasets <- lapply(1:length(allRes), function(i) {
  #   allRes[[i]]$dataset
  # }) %>% unlist()

  datasets <- strsplit(allFiles3, "DSCC-")
  datasets <- lapply(datasets, function(elm){elm[2]}) %>% unlist()
  datasets <- gsub(".rds", "",  datasets)

  names(allRes) <- datasets
  # datasets <- sort(datasets)
  # allRes <- allRes[datasets]

  colors <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF", "#E64B357F", "#4DBBD57F", "#00A0877F", "#3C54887F", "#F39B7F7F", "#8491B47F", "#91D1C27F", "#DC00007F", "#7E61487F", "#B09C857F")
  plots <- list()
  for (dataset in datasets) {
    print(dataset)
    # title <- strsplit(dataset, "-BC-GENE-MAP")
    title <- dataset
    # survival <- readRDS(paste0(processedDataPath, dataset, "-survival.rds"))
    if (length(grep(dataset, allFiles)) > 0){
      survival <- readRDS(file.path(processedDataPath, paste0(dataset, ".rds")))$survival # read survival information
    }else{
      survival <- readRDS(file.path(processedDataPath2, paste0(dataset, ".rds")))$survival # read survival information
    }
    # cluster <- readRDS(file.path(tmpRoot, paste0("DSCC-", dataset, ".rds")))$cluster # read cluster result
    cluster <- allRes[[dataset]]$cluster

    if (length(cluster) == 1) {
      next
    }
    toUseSamples <- rownames(survival)
    # toUseSamples <- names(cluster)

    # creating kaplan meier plot
    survival_sub <- survival[intersect(toUseSamples, names(cluster)),]
    cluster_sub <- cluster[intersect(toUseSamples, names(cluster))]
    fit1 <- survfit(Surv(time = os, event = isDead) ~ as.factor(cluster_sub), data = survival_sub)
    pval <- surv_pvalue(fit1)$pval
    p <- ggsurvplot(fit1,
                    data = survival_sub,
                    palette = colors,
                    conf.int = FALSE,
                    xlab = "Days",
                    pval = signif(pval, digits = 3),
                    pval.size = 3,
                    # pval.size = 2.5,
                    # size = 0.5,
                    size = 1,
                    legend = 'right',
                    legend.title = "Subtype",  # Change legend title to "Subtype"
                    legend.labs = paste0("", 1:length(unique(cluster_sub))),
                    censor = FALSE)$plot

    consistent_size <- 8  # Choose the size you want for all elements
    p <- p +
      ggtitle(dataset) +
      theme(
        # axis.text = element_text(size = 2)
        plot.title = element_text(size = consistent_size),  # Smaller title
        axis.title.x = element_text(size = consistent_size),  # Smaller x-axis title
        axis.title.y = element_text(size = consistent_size),  # Smaller y-axis title
        axis.text.x = element_text(size = consistent_size),   # Smaller x-axis values
        axis.text.y = element_text(size = consistent_size),   # Smaller y-axis values
        # legend.position = "bottom",
        legend.title = element_text(size = consistent_size),  # Smaller legend title
        legend.text = element_text(size = consistent_size), # Smaller legend text
        legend.key.height = unit(0.6, "lines"),
        legend.key.width = unit(0.5, "lines"),
        legend.spacing.y = unit(0.1, "cm"),
        # legend.margin = margin(0, 0, 0, 0, "cm"),
        # plot.margin = margin(5, 0, 5, 5, "pt"),
        # legend.box.margin = margin(-10, 0, 0, -10, "pt")


        # legend.key.size = unit(0.4, "lines"),   # Smaller legend keys
        # legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
        # legend.box.margin = margin(t = -10),
        # panel.grid.minor = element_blank(),
        # panel.grid.major = element_line(color = "gray90", size = 0.1)
      )
    # ggsave(paste0("/data/dungp/projects/subtyping/DSCC-subtyping/DSCC-Subtyping/src/new-code/KEGG-Combinations/new-data-processing/figs/kaplan-meier/", dataset, "-cluster-survival.pdf"), p, width = 6, height = 6)
    plots[[dataset]] = p
  }

  # Initialize min and max limits for x and y axes
  global_xlim <- c(Inf, -Inf)
  global_ylim <- c(Inf, -Inf)
  # Find number of plots
  n_plots <- length(plots)

  p <- grid.arrange(grobs = plots, ncol = 2
                    # top = textGrob("DSCCM", gp = gpar(fontsize = 20, fontface = "bold"))
  ) #
}

# method <- "NEMO"
# ggsave(paste0("/data/daotran/Cancer_RP/Subtyping/BiB_Subtyping/Plot/KMPlots/DSCCM_km_mtb2.pdf"), p, width = 7, height = 5)
# ggsave(paste0("/data/daotran/Cancer_RP/Subtyping/BiB_Subtyping/Plot/KMPlots/DSCCM_km_tcga3_mtb1.pdf"), p, width = 7, height = 7)
ggsave(paste0("/data/daotran/Cancer_RP/Subtyping/BiB_Subtyping/Plot/KMPlots/DSCCM_km_2dts.pdf"), p, width = 6.5, height = 3)

