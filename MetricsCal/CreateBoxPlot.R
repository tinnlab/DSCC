library(scales)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(matrixStats)
library(grid)
library(gridExtra)

resPath <- "/data/daotran/Cancer_RP/Subtyping/BiB_Subtyping/Results/AA_SP"
ResTable_CI <- read.csv(file.path(resPath, 'CIndex_Table.csv'), row.names=1)
ResTable_tdCI <- read.csv(file.path(resPath, 'tdCIndex_Table.csv'), row.names=1)
ResTable_IBS <- read.csv(file.path(resPath, 'IBS_Table.csv'), row.names=1)
ResTable_Cal <- read.csv(file.path(resPath, 'Calibration_Table.csv'), row.names=1)

ResTable_IBS <- 1-ResTable_IBS
ResTable_Cal <- ResTable_Cal/25

alldatasets <- rownames(ResTable_CI)
alldatasets <- alldatasets[grep("TCGA", alldatasets)]
alldatasets <- union(alldatasets, rownames(ResTable_CI))

# rownames(ResTable) <- paste0(rownames(ResTable), "_No Standardization")

# resPath <- "/data/daotran/Cancer_RP/Subtyping/Metabolite_Analysis/Results_SP_scale"
# ResTable_sc <- read.csv(file.path(resPath, 'CIndex_Table.csv'), row.names=1)
# rownames(ResTable_sc) <- paste0(rownames(ResTable_sc), "_With Standardization")

# cindex_cb <- rbind(ResTable, ResTable_sc)
#
# colMedians(as.matrix(ResTable))
# colMeans(as.matrix(ResTable))
# colSds(as.matrix(ResTable))

# methods <- c("DSCCM", "CC", "CIMLR", "SNF", "NEMO", "LRACluster", "IntNMF", "ANF")
methods <- c("nosubtype", "DSCC", "CC", "CIMLR", "SNF", "LRACluster", "IntNMF", "ANF")
ResTable_CI <- ResTable_CI[alldatasets, methods]

colMeans(ResTable_CI, na.rm = T)
colMedians(as.matrix(ResTable_CI), na.rm = T)

# {
#   allValue <- pivot_longer(as.data.frame(ResTable_CI[, methods]),
#                            cols = everything(),
#                            names_to = "Method",
#                            values_to = "Value")
#   allValue <- hablar::retype(allValue)
#
#   # Ensure Value is numeric - force conversion and handle any issues
#   allValue$Value <- as.numeric(as.character(allValue$Value))
#
#   # Check for conversion issues
#   # print(paste("Number of NA values after conversion:", sum(is.na(allValue$Value))))
#
#   # Drop NA values if needed
#   # allValue <- allValue[!is.na(allValue$Value), ]
#
#   allValue$Method <- factor(allValue$Method, levels = methods)
#
#   # Calculate mean values for each method
#   mean_values <- aggregate(Value ~ Method, data = allValue, FUN = mean)
#
#   # Define a function to calculate median and add mean label
#   median_plus_mean_label <- function(x) {
#     m <- median(x)
#     mean_val <- mean(x)
#     return(data.frame(y = m + 0.02, label = sprintf("%.3f", mean_val)))
#   }
#
#   # Create plot
#   plt_ci <- ggplot(allValue, aes(x = Method, y = Value, fill = Method)) +
#     theme_classic() +
#     # scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF",
#     #                              "#E64B357F", "#4DBBD57F", "#00A0877F", "#3C54887F", "#F39B7F7F", "#8491B47F", "#91D1C27F", "#DC00007F", "#7E61487F", "#B09C857F")) +
#     scale_fill_manual(values = c("#7E6148FF", "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#B09C85FF",
#                                  "#7E6148FF", "#E64B357F", "#4DBBD57F", "#00A0877F", "#3C54887F", "#F39B7F7F", "#8491B47F", "#91D1C27F", "#B09C857F")) +
#     labs(x = "", y = "Concordance Index", color = "black", size = 14) +
#     geom_boxplot(outlier.shape = NA, width = 0.4, position = position_dodge(width = 0.1)) +
#     theme(axis.text.x = element_text(angle = 35, vjust = 0.5, hjust = 0.4, size = 09, color = "black"),
#           axis.text.y = element_text(size = 09, color = "black"),
#           axis.title.y = element_text(size = 10, color = "black", face = "bold")) +
#     # geom_hline(yintercept = 0.9, linetype = "dashed") +
#     theme(legend.position = "none") +
#     scale_y_continuous(limits = c(0.4, 1), oob = rescale_none, breaks = c(0.4, 0.6, 0.8, 1)) +
#     # Add mean points
#     # stat_summary(fun = "mean", geom = "point", shape = 18, size = 2, color = "black") +
#     # Add mean labels above median line
#     stat_summary(fun.data = median_plus_mean_label, geom = "text", size = 2.5)
#
#   print(plt_ci)
# }

{
  allValue <- pivot_longer(as.data.frame(ResTable_CI[, methods]),
                           cols = everything(),
                           names_to = "Method",
                           values_to = "Value")
  allValue <- hablar::retype(allValue)

  # Ensure Value is numeric - force conversion and handle any issues
  allValue$Value <- as.numeric(as.character(allValue$Value))

  # Check for conversion issues
  # print(paste("Number of NA values after conversion:", sum(is.na(allValue$Value))))

  # Drop NA values if needed
  # allValue <- allValue[!is.na(allValue$Value), ]

  allValue$Method <- factor(allValue$Method, levels = methods)

  # Calculate median values for each method to position the text
  medians <- aggregate(Value ~ Method, data = allValue, FUN = median)

  # Create a data frame for the median values
  annotations <- medians  # Now using just the medians dataframe

  # out_name <- "minuslog_pval_TCGA"
  # pdf(sprintf("/data/daotran/Cancer_RP/Subtyping/Metabolite_Analysis/Plots/%s.pdf", out_name), width = 5, height = 4)
  plt_ci <- ggplot(allValue, aes(x = Method, y = Value, fill = Method)) +
    theme_classic() +
    # scale_fill_npg() +
    scale_fill_manual(values = c("#7E6148FF", "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#B09C85FF",
                                 "#7E6148FF", "#E64B357F", "#4DBBD57F", "#00A0877F", "#3C54887F", "#F39B7F7F", "#8491B47F", "#91D1C27F", "#B09C857F")) +
    labs(x = "", y = "Concordance Index", color = "black", size = 14) +
    geom_boxplot(outlier.shape = NA, width = 0.4, position = position_dodge(width = 0.1)) +
    # Add the horizontal line at y=1.3
    geom_hline(yintercept = 1.3, linetype = "dashed") +
    # Add text labels for median values above the median for each method
    geom_text(data = annotations,
              aes(x = Method, y = Value, label = sprintf("%.3f", Value)),  # Display the median value with 2 decimal places
              vjust = -0.5, size = 2.5, fontface = "bold") +
    theme(axis.text.x = element_text(angle = 35, vjust = 0.5, hjust = 0.4, size = 9, color = "black"),
          axis.text.y = element_text(size = 9, color = "black"),
          axis.title.y = element_text(size = 10, color = "black", face = "bold")) +
    theme(legend.position = "none") +
    scale_y_continuous(limits = c(0.4, 1), oob = rescale_none, breaks = seq(0.4, 1, by = 0.2))
  print(plt_ci)
}
ggsave("/data/daotran/Cancer_RP/Subtyping/BiB_Subtyping/Plot/allbutClin_cindex.pdf", plot = plt_ci, width = 6, height = 3)


{
  allValue <- pivot_longer(as.data.frame(ResTable_tdCI),
                           cols = everything(),
                           names_to = "Method",
                           values_to = "Value")
  allValue <- hablar::retype(allValue)

  # Ensure Value is numeric - force conversion and handle any issues
  allValue$Value <- as.numeric(as.character(allValue$Value))

  # Check for conversion issues
  # print(paste("Number of NA values after conversion:", sum(is.na(allValue$Value))))

  # Drop NA values if needed
  # allValue <- allValue[!is.na(allValue$Value), ]

  allValue$Method <- factor(allValue$Method, levels = methods)

  # Calculate mean values for each method
  mean_values <- aggregate(Value ~ Method, data = allValue, FUN = mean)

  # Define a function to calculate median and add mean label
  median_plus_mean_label <- function(x) {
    m <- median(x)
    mean_val <- mean(x)
    return(data.frame(y = m + 0.02, label = sprintf("%.3f", mean_val)))
  }

  # Create plot
  plt_tdci <- ggplot(allValue, aes(x = Method, y = Value, fill = Method)) +
    theme_classic() +
    scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF",
                                 "#E64B357F", "#4DBBD57F", "#00A0877F", "#3C54887F", "#F39B7F7F", "#8491B47F", "#91D1C27F", "#DC00007F", "#7E61487F", "#B09C857F")) +
    labs(x = "", y = "td C-Index", color = "black", size = 14) +
    geom_boxplot(outlier.shape = NA, width = 0.4, position = position_dodge(width = 0.1)) +
    theme(axis.text.x = element_text(angle = 35, vjust = 0.5, hjust = 0.4, size = 09, color = "black"),
          axis.text.y = element_text(size = 09, color = "black"),
          axis.title.y = element_text(size = 10, color = "black", face = "bold")) +
    # geom_hline(yintercept = 0.9, linetype = "dashed") +
    theme(legend.position = "none") +
    scale_y_continuous(limits = c(0.4, 0.8), oob = rescale_none, breaks = c(0.4, 0.6, 0.8)) +
    # Add mean points
    # stat_summary(fun = "mean", geom = "point", shape = 18, size = 2, color = "black") +
    # Add mean labels above median line
    stat_summary(fun.data = median_plus_mean_label, geom = "text", size = 2.5)

  print(plt_tdci)
}


{
  allValue <- pivot_longer(as.data.frame(ResTable_IBS),
                           cols = everything(),
                           names_to = "Method",
                           values_to = "Value")
  allValue <- hablar::retype(allValue)

  # Ensure Value is numeric - force conversion and handle any issues
  allValue$Value <- as.numeric(as.character(allValue$Value))

  # Check for conversion issues
  # print(paste("Number of NA values after conversion:", sum(is.na(allValue$Value))))

  # Drop NA values if needed
  # allValue <- allValue[!is.na(allValue$Value), ]

  allValue$Method <- factor(allValue$Method, levels = methods)

  # Calculate mean values for each method
  mean_values <- aggregate(Value ~ Method, data = allValue, FUN = mean)

  # Define a function to calculate median and add mean label
  median_plus_mean_label <- function(x) {
    m <- median(x)
    mean_val <- mean(x)
    return(data.frame(y = m + 0.02, label = sprintf("%.3f", mean_val)))
  }

  # Create plot
  plt_ibs <- ggplot(allValue, aes(x = Method, y = Value, fill = Method)) +
    theme_classic() +
    scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF",
                                 "#E64B357F", "#4DBBD57F", "#00A0877F", "#3C54887F", "#F39B7F7F", "#8491B47F", "#91D1C27F", "#DC00007F", "#7E61487F", "#B09C857F")) +
    labs(x = "", y = "1-IBS", color = "black", size = 14) +
    geom_boxplot(outlier.shape = NA, width = 0.4, position = position_dodge(width = 0.1)) +
    theme(axis.text.x = element_text(angle = 35, vjust = 0.5, hjust = 0.4, size = 09, color = "black"),
          axis.text.y = element_text(size = 09, color = "black"),
          axis.title.y = element_text(size = 10, color = "black", face = "bold")) +
    # geom_hline(yintercept = 0.9, linetype = "dashed") +
    theme(legend.position = "none") +
    scale_y_continuous(limits = c(0.7, 1), oob = rescale_none, breaks = c(0.7, 0.8,0.9, 1)) +
    # Add mean points
    # stat_summary(fun = "mean", geom = "point", shape = 18, size = 2, color = "black") +
    # Add mean labels above median line
    stat_summary(fun.data = median_plus_mean_label, geom = "text", size = 2.5)

  print(plt_ibs)
}


{
  allValue <- pivot_longer(as.data.frame(ResTable_Cal),
                           cols = everything(),
                           names_to = "Method",
                           values_to = "Value")
  allValue <- hablar::retype(allValue)

  # Ensure Value is numeric - force conversion and handle any issues
  allValue$Value <- as.numeric(as.character(allValue$Value))

  # Check for conversion issues
  # print(paste("Number of NA values after conversion:", sum(is.na(allValue$Value))))

  # Drop NA values if needed
  # allValue <- allValue[!is.na(allValue$Value), ]

  allValue$Method <- factor(allValue$Method, levels = methods)

  # Calculate mean values for each method
  mean_values <- aggregate(Value ~ Method, data = allValue, FUN = mean)

  # Define a function to calculate median and add mean label
  median_plus_mean_label <- function(x) {
    m <- median(x)
    mean_val <- mean(x)
    return(data.frame(y = m + 0.02, label = sprintf("%.3f", mean_val)))
  }

  # Create plot
  plt_cal <- ggplot(allValue, aes(x = Method, y = Value, fill = Method)) +
    theme_classic() +
    scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF",
                                 "#E64B357F", "#4DBBD57F", "#00A0877F", "#3C54887F", "#F39B7F7F", "#8491B47F", "#91D1C27F", "#DC00007F", "#7E61487F", "#B09C857F")) +
    labs(x = "", y = "Calibration", color = "black", size = 14) +
    geom_boxplot(outlier.shape = NA, width = 0.4, position = position_dodge(width = 0.1)) +
    theme(axis.text.x = element_text(angle = 35, vjust = 0.5, hjust = 0.4, size = 09, color = "black"),
          axis.text.y = element_text(size = 09, color = "black"),
          axis.title.y = element_text(size = 10, color = "black", face = "bold")) +
    # geom_hline(yintercept = 0.9, linetype = "dashed") +
    theme(legend.position = "none") +
    scale_y_continuous(limits = c(0.5,1), oob = rescale_none, breaks = seq(0.5,1,0.1)) +
    # Add mean points
    # stat_summary(fun = "mean", geom = "point", shape = 18, size = 2, color = "black") +
    # Add mean labels above median line
    stat_summary(fun.data = median_plus_mean_label, geom = "text", size = 2.5)

  print(plt_cal)
}


combined_plot <- grid.arrange(plt_ci, plt_tdci, plt_ibs, plt_cal, ncol = 2)


# {
#   allValue <- pivot_longer(as.data.frame(ResTable),
#                            cols = everything(),
#                            names_to = "Method",
#                            values_to = "Value")
#   allValue <- hablar::retype(allValue)

#   # Ensure Value is numeric - force conversion and handle any issues
#   allValue$Value <- as.numeric(as.character(allValue$Value))

#   # Check for conversion issues
#   # print(paste("Number of NA values after conversion:", sum(is.na(allValue$Value))))

#   # Drop NA values if needed
#   # allValue <- allValue[!is.na(allValue$Value), ]

#   allValue$Method <- factor(allValue$Method, levels = methods)

#   # Calculate mean values for each method
#   mean_values <- aggregate(Value ~ Method, data = allValue, FUN = mean)

#   # Define a function to calculate median and add mean label
#   median_plus_mean_label <- function(x) {
#     m <- median(x)
#     mean_val <- mean(x)
#     return(data.frame(y = m + 0.02, label = sprintf("%.3f", mean_val)))
#   }

#   # Create plot
#   plt <- ggplot(allValue, aes(x = Method, y = Value, fill = Method)) +
#     theme_classic() +
#     scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF",
#                                  "#E64B357F", "#4DBBD57F", "#00A0877F", "#3C54887F", "#F39B7F7F", "#8491B47F", "#91D1C27F", "#DC00007F", "#7E61487F", "#B09C857F")) +
#     labs(x = "", y = "Concordance Index", color = "black", size = 14) +
#     geom_boxplot(outlier.shape = NA, width = 0.4, position = position_dodge(width = 0.1)) +
#     theme(axis.text.x = element_text(angle = 35, vjust = 0.5, hjust = 0.4, size = 09, color = "black"),
#           axis.text.y = element_text(size = 09, color = "black"),
#           axis.title.y = element_text(size = 10, color = "black", face = "bold")) +
#     # geom_hline(yintercept = 0.9, linetype = "dashed") +
#     theme(legend.position = "none") +
#     scale_y_continuous(limits = c(0, 0.4), oob = rescale_none, breaks = c(0, 0.2, 0.4)) +
#     # Add mean points
#     # stat_summary(fun = "mean", geom = "point", shape = 18, size = 2, color = "black") +
#     # Add mean labels above median line
#     stat_summary(fun.data = median_plus_mean_label, geom = "text", size = 3)

#   print(plt)
# }

# {
#   allValue <- pivot_longer(as.data.frame(ResTable),
#                            cols = everything(),
#                            names_to = "Method",
#                            values_to = "Value")
#   allValue <- hablar::retype(allValue)

#   # Ensure Value is numeric - force conversion and handle any issues
#   allValue$Value <- as.numeric(as.character(allValue$Value))

#   # Check for conversion issues
#   # print(paste("Number of NA values after conversion:", sum(is.na(allValue$Value))))

#   # Drop NA values if needed
#   # allValue <- allValue[!is.na(allValue$Value), ]

#   allValue$Method <- factor(allValue$Method, levels = methods)

#   # Calculate mean values for each method
#   mean_values <- aggregate(Value ~ Method, data = allValue, FUN = mean)

#   # Define a function to calculate median and add mean label
#   median_plus_mean_label <- function(x) {
#     m <- median(x)
#     mean_val <- mean(x)
#     return(data.frame(y = m + 0.02, label = sprintf("%.3f", mean_val)))
#   }

#   # Create plot
#   plt <- ggplot(allValue, aes(x = Method, y = Value, fill = Method)) +
#     theme_classic() +
#     scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF",
#                                  "#E64B357F", "#4DBBD57F", "#00A0877F", "#3C54887F", "#F39B7F7F", "#8491B47F", "#91D1C27F", "#DC00007F", "#7E61487F", "#B09C857F")) +
#     labs(x = "", y = "Concordance Index", color = "black", size = 14) +
#     geom_boxplot(outlier.shape = NA, width = 0.4, position = position_dodge(width = 0.1)) +
#     theme(axis.text.x = element_text(angle = 35, vjust = 0.5, hjust = 0.4, size = 09, color = "black"),
#           axis.text.y = element_text(size = 09, color = "black"),
#           axis.title.y = element_text(size = 10, color = "black", face = "bold")) +
#     # geom_hline(yintercept = 0.9, linetype = "dashed") +
#     theme(legend.position = "none") +
#     scale_y_continuous(limits = c(21, 30), oob = rescale_none, breaks = seq(21, 30, 4)) +
#     # Add mean points
#     # stat_summary(fun = "mean", geom = "point", shape = 18, size = 2, color = "black") +
#     # Add mean labels above median line
#     stat_summary(fun.data = median_plus_mean_label, geom = "text", size = 3)

#   print(plt)
# }

# {
#   cindex_cb_df <- as.data.frame(cindex_cb)
#   cindex_cb_df$dataset_datatype <- rownames(cindex_cb_df)
#   cindex_cb_df <- cindex_cb_df %>%
#     separate(dataset_datatype, into = c("Dataset", "DataType"), sep = "\\_")
#   # pva_cox_tab_df$DataType[pva_cox_tab_df$DataType == "Metabo"] <- "Metabolomics"
#   # Reshape the data to long format
#   allValue <- pivot_longer(cindex_cb_df,
#                            cols = -c(Dataset, DataType),
#                            names_to = "Method",
#                            values_to = "Value")
#   # Convert to appropriate types
#   allValue <- hablar::retype(allValue)
#   allValue$Value <- as.numeric(allValue$Value)
#   # If you have a specific order for methods, set it here
#   if (exists("methods")) {
#     allValue$Method <- factor(allValue$Method, levels = methods)
#   }

#   # Create the plot with Method on x-axis and DataType determining the fill color
#   # out_name <- "minuslog_pval_bydatatype"
#   # pdf(sprintf("/data/daotran/Cancer_RP/Subtyping/Metabolite_Analysis/Plots/%s.pdf", out_name), width = 5, height = 4)

#   # Calculate median values for each Method and DataType combination
#   median_values <- allValue %>%
#     group_by(Method, DataType) %>%
#     summarize(median_val = median(Value, na.rm = TRUE), .groups = "drop")

#   plt <- ggplot(allValue, aes(x = Method, y = Value, fill = DataType)) +
#     theme_classic() +
#     scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF",
#                                  "#E64B357F", "#4DBBD57F", "#00A0877F", "#3C54887F", "#F39B7F7F", "#8491B47F", "#91D1C27F", "#DC00007F", "#7E61487F", "#B09C857F")) +
#     labs(x = "", y = "Concordance Index", fill = NULL, size = 14) +
#     geom_boxplot(outlier.shape = NA, width = 0.6, position = position_dodge(width = 0.8)) +
#     # geom_boxplot(width = 0.6, position = position_dodge(width = 0.8)) +
#     theme(axis.text.x = element_text(angle = 35, vjust = 0.5, hjust = 0.4, size = 9, color = "black", face = "bold"),
#           axis.text.y = element_text(size = 9, color = "black"),
#           axis.title.y = element_text(size = 10, color = "black", face = "bold"),
#           legend.position = "top") +
#     geom_hline(yintercept = 1.3, linetype = "dashed", color = "gray50") +
#     scale_y_continuous(limits = c(0.4, 0.7), oob = rescale_none, breaks = seq(0.4,0.7, 0.1)) +
#     # Add median value labels above each boxplot
#     geom_text(data = median_values,
#               aes(x = Method, y = median_val, label = sprintf("%.2f", median_val), group = DataType),
#               position = position_dodge(width = 0.8),
#               vjust = -0.5,
#               size = 2.3)

#   print(plt)

# }

ggsave("/data/daotran/Cancer_RP/Subtyping/BiB_Subtyping/Plot/SP_compare.pdf", plot = combined_plot, width = 12, height = 9)

