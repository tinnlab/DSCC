library(tidyverse)
library(parallel)
library(limma)
library(ggplot2)
library(ggrepel)

mapmRNA_to_PwGene <- function(omicsdata) {
  Samples <- rownames(omicsdata)
  allGenes <- colnames(omicsdata)
  allGenes <- strsplit(allGenes, "___")
  allGenes <- lapply(c(1:length(allGenes)), function(i) {
    Gene <- allGenes[[i]][1]
    Gene
  }) %>% do.call(what = c)
  uniGenes <- as.list(unique(allGenes))
  omicsdata <- mclapply(uniGenes, mc.cores = 8, function(gene) {
    Mat <- as.matrix(omicsdata[, allGenes %in% gene])
    Gmean <- as.matrix(rowMeans(Mat, na.rm = T))
    Gmean
  }) %>% do.call(what = cbind)
  colnames(omicsdata) <- unique(allGenes)
  rownames(omicsdata) <- Samples
  return(as.data.frame(omicsdata))
}

datPath <- "/data/daotran/Cancer_RP/Subtyping/Data/TCGA_original/TCGA-ACC.rds"
rds <- readRDS(datPath)

subtypes <- readRDS("/data/daotran/Cancer_RP/Subtyping/BiB_Subtyping/Results/allbutClin/DSCCM-TCGA-ACC.rds")
subtypes <- subtypes$cluster

expdat <- rds$mRNATPM
expdat <- mapmRNA_to_PwGene(expdat)
expdat[is.na(expdat)] <- 0
expdat <- log2(expdat + 1)
length(unique(colnames(expdat)))

subtypes <- as.data.frame(as.matrix(subtypes))
# subtypes$IDs <- paste0("P", rownames(subtypes))
subtypes$IDs <- rownames(subtypes)
subtypes$V1[subtypes$V1 == 1] <- "c"
subtypes$V1[subtypes$V1 == 2] <- "d"
subtypes$V1[subtypes$V1 == 3] <- "d"
subtypes <- subtypes[, c("IDs", "V1")]

cSamples <- intersect(rownames(expdat), subtypes$IDs)
expdat <- t(expdat)
expdat <- expdat[, cSamples]
subtypes <- subtypes[cSamples, ]
counts <- as.matrix(expdat)
condition <- factor(subtypes$V1)

c_ids <- grep("c", condition)
d_ids <- grep("d", condition)
c_exp <- counts[, c_ids]
d_exp <- counts[, d_ids]

FC <- rowMeans(c_exp) - rowMeans(d_exp)
FC <- as.matrix(FC)
rownames(FC) <- toupper(rownames(FC))
FC <- as.data.frame(FC)
FC$Gene <- rownames(FC)
colnames(FC) <- c("Fold-Change", "Genes")
FC <- FC[, c("Genes", "Fold-Change")]
write.table(FC, file="/data/daotran/Cancer_RP/Subtyping/BiB_Subtyping/Data/TCGA-ACC_FoldChange.txt", sep="\t", quote=FALSE,
            col.names=TRUE, row.names=FALSE)

#### try limma
# Create a design matrix
design <- model.matrix(~condition)
colnames(design) <- c("Intercept", "conditionTreatment")
design[, 2] <- 1 - design[, 2]

fit <- lmFit(counts, design)

# Apply empirical Bayes smoothing
fit <- eBayes(fit)

# Extract results
results <- topTable(fit, coef = "conditionTreatment", n = Inf)

# Add a column indicating significance
results$significant <- ifelse(results$adj.P.Val < 0.05, "Yes", "No")

volcano_plot <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
  theme_bw() +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Volcano Plot of Differential Expression",
       x = "log2 Fold Change",
       y = "-log10 Adjusted P-value",
       color = "Significant")

print(volcano_plot)
significant_genes <- results[results$adj.P.Val < 0.05,]
sum(significant_genes$logFC >= 0)


### get clinical variables involved
# rds <- readRDS("/data/daotran/Cancer_RP/Subtyping/Data/TCGA-mapped-2/TCGA-ACC.rds")
# clin <- rds$clinicalImputedV3
# colnames(clin)[grep("stage", colnames(clin))]
# clin_stage <- clin[cSamples, grep("stage", colnames(clin))]
# stages <- c("pathologic_stage_Stage I", "pathologic_stage_Stage II",
#             "pathologic_stage_Stage III", "pathologic_stage_Stage IV")
# clin_stage <- clin_stage[, stages]
# stage_vec <- lapply(1:nrow(clin_stage), function(i){
#   vec <- clin_stage[i, ]
#   which(vec == 1)
# }) %>% do.call(what = c)
# stage_vec <- names(stage_vec)
# stage_vec <- gsub("pathologic_stage_", "", stage_vec)
#
# conti_race <- lapply(sort(unique(stage_vec)), function(stage){
#   table(subtypes[stage_vec == stage, "V1"])
# }) %>% do.call(what = rbind) %>% `rownames<-` (sort(unique(stage_vec)))
#
# race_chisq <- chisq.test(conti_race)
# std_resid_race <- race_chisq$stdres
#
# clin_gender <- clin[cSamples, grep("gender", colnames(clin))]
# gender_vec <- lapply(1:nrow(clin_gender), function(i){
#   vec <- clin_gender[i, ]
#   which(vec == 1)
# }) %>% do.call(what = c)
# gender_vec <- names(gender_vec)
# gender_vec <- gsub("gender_", "", gender_vec)
#
# conti_gender <- lapply(sort(unique(gender_vec)), function(gender){
#   table(subtypes[gender_vec == gender, "V1"])
# }) %>% do.call(what = rbind) %>% `rownames<-` (sort(unique(gender_vec)))
#
# gender_chisq <- chisq.test(conti_gender)
# std_resid_gender <- gender_chisq$stdres

### variant analysis
org_data <- readRDS("/data/share/dungp/gene-network/data/TCGA-ACC.rds")
data <- org_data[c("SNVMuseGeneFreq", "SNVVarscanGeneFreq", "SNVSniperGeneFreq", "SNVMutectGeneFreq")]

data <- lapply(data, function(table){
  table[rownames(table) %in% cSamples, ]
  table[is.na(table)] <- 0
  table
})

allGenes <- lapply(data, colnames) %>% unlist() %>% unique()

data.new <- lapply(data, function(table) {
  new.mat <- matrix(0, ncol = length(allGenes), nrow = nrow(table))
  colnames(new.mat) <- allGenes
  matching_cols <- intersect(colnames(table), allGenes)
  col_indices <- match(matching_cols, allGenes)
  new.mat[, col_indices] <- as.matrix(table[, matching_cols])
  new.mat
})
data.new <- Reduce("+", data.new)/length(data.new)

mutation <- lapply(1:ncol(data.new), function(i){
  ls_ids <- which(subtypes$V1 == "c")
  hg_ids <- which(subtypes$V1 == "d")
  c(sum(data.new[ls_ids, i] > 0), sum(data.new[hg_ids, i] > 0))
}) %>% do.call(what = rbind) %>% `rownames<-` (colnames(data.new))
mutation <- as.data.frame(mutation)
mutation$gene <- rownames(mutation)
highlight_genes <- rownames(mutation)[order(abs(mutation[, 2] - mutation[, 1]), decreasing = T)[1:6]]
colnames(mutation) <- c("lower_surviving_group", "higher_surviving_group", "gene")
mutation$highlight <- mutation$gene %in% highlight_genes

p <- ggplot(mutation, aes(x = lower_surviving_group, y = higher_surviving_group)) +
  # Add black points for non-highlighted genes
  geom_point(data = subset(mutation, !highlight),
             color = "black", size = 1) +
  # Add red points for highlighted genes
  geom_point(data = subset(mutation, highlight),
             color = "#E74C3C", size = 2) +
  # Add labels for highlighted genes with boxes
  geom_label_repel(data = subset(mutation, highlight),
                  aes(label = gene),
                  box.padding = 0.5,
                  point.padding = 0.5,
                  segment.color = 'transparent',
                  label.size = 0.5) +
  # Set axis labels
  labs(
    title = "Mutations in ACC subtypes",
    x = "short-term survival",
    y = "long-term survival"
  ) +
  # Set theme with white background
  theme_minimal() +
  # Remove grid lines
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray95")
  )

ggsave(paste0("/data/daotran/Cancer_RP/Subtyping/BiB_Subtyping/Plot/mutation_plot/mutation_plot.pdf"), p, width = 5, height = 4)

