# Differential expression analysis
library(parallel)
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(clusterProfiler)
library(dplyr)
library(ggrepel)

# import data:
res <- mclapply(dir(pattern="*.counts.txt", full.names=TRUE), function(fil){
  read.delim(fil, header=FALSE, stringsAsFactors=FALSE)
  }, mc.cores=8)

names(res) <- gsub("*.counts.txt", "" , dir(pattern="*.counts.txt"))
count_list <- lapply(res, function(x) setNames(x[, 2], x[, 1]))

# then we extract the additional info that HTSeq writes at the end of every file detailing 
addInfo <- c("__no_feature","__ambiguous",
             "__too_low_aQual","__not_aligned",
             "__alignment_not_unique")

# create dataframe
count_data <- do.call(cbind, count_list)
last_5_lines <- tail(count_data, 5)

# remove the addInfo part from the count_data
count_data <- count_data[!(rownames(count_data) %in% addInfo), ]

# set condition based on sample_id
sample_ids <- c("SRR1982462", "SRR1982463", "SRR1982464", "SRR1982465", "SRR1982466", "SRR1982467", "SRR1982468")
type <- c("TAF10KO", "TAF10HET", "WT", "TAF10HET", "TAF10HET", "WT", "TAF10KO")

condition_mapping <- data.frame(sample_id = sample_ids,
                                condition = type)

colData <- DataFrame(condition = factor(condition_mapping$condition))

# create DESeqDataSet from matrix
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = colData,
                              design = ~ condition)

dds$condition <- factor(dds$condition)

# filtering of lowly expressed genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# DESeq analysis
dds <- DESeq(dds, test = "Wald")


## PCA ##
# variance stabilization transformation
dds_vst <- vst(dds, blind = TRUE)

# Plot PCA
PCA <- plotPCA(dds_vst, intgroup = "condition")
PCA +  geom_point(aes(shape = dds_vst$type, color = dds_vst$type), size = 4.3) +
  scale_shape_manual(values = c("TAF10KO" = 15, "TAF10HET" = 16, "WT" = 17)) +
  scale_color_manual(values = c("TAF10KO" = "#46ACC8", "TAF10HET" = "#E58601", "WT"="#B40F20")) +
  theme_pubr() +
  guides(shape = guide_legend(title = "Sample shape"),
         color = guide_legend(title = "Sample colour")) +
  theme(legend.position = "right") +
  ggtitle("PCA of Fetal Liver Samples in Mice") +
  theme(plot.title = element_text(hjust = 0.5))




## Top 5 differentially expressed genes ##

# DEA results
res_KO_HET <- results(dds, contrast = c('condition', 'TAF10KO', 'TAF10HET'))
res_WT_KO <- results(dds, contrast = c('condition', 'WT', 'TAF10KO'))
res_WT_HET <- results(dds, contrast = c('condition', 'WT', 'TAF10HET'))
summary(res_KO_HET)
summary(res_WT_KO)
summary(res_WT_HET)



# define new dds for the DEA between WT/TAF10HET vs TAF10KO
# set condition based on sample_id
sample_ids <- c("SRR1982462", "SRR1982463", "SRR1982464", "SRR1982465", "SRR1982466", "SRR1982467", "SRR1982468")
type <- c("TAF10KO", "TAF10HET", "WT", "TAF10HET", "TAF10HET", "WT", "TAF10KO")
condition <- c("TAF10KO", "CONTROL", "CONTROL", "CONTROL", "CONTROL", "CONTROL", "TAF10KO")

condition_mapping2 <- data.frame(sample_id = sample_ids,
                                type = type,
                                condition = condition)

colData2 <- DataFrame(condition = factor(condition_mapping2$condition), type=type)

# create DESeqDataSet from matrix
dds2 <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = colData2,
                              design = ~ condition)

dds2$condition <- factor(dds2$condition)

# filtering of lowly expressed genes
keep2 <- rowSums(counts(dds2)) >= 10
dds2 <- dds2[keep2,]

# DESeq analysis
dds2 <- DESeq(dds2, test = "Wald")

# store results
results <- results(dds2)
summary(results)


# create function that will return a mapping df having all the gene names that correspond to the ensembl ids
ens_to_gene <- function(res) {
  z=c()
  for (i in row.names(res)){
    a=strsplit(i,"[[:punct:]]")[[1]][1]
    z=append(z, a)
  }
  
  eg = as.data.frame(clusterProfiler::bitr(z,
                                           fromType = "ENSEMBL",
                                           toType = "SYMBOL",
                                           OrgDb = "org.Mm.eg.db", drop = TRUE))
  
  return(eg)
}

map_df <- ens_to_gene(res_KO_HET)


# create a function that converts the res into a df and adds a new column of ensembl ids
transform_res <- function(res, map_df){
  # convert to DF
  res1 <- as.data.frame(res)
  
  # make column genes from rownames, by removing part after (.)
  res1$ensembl_id <- gsub("\\..*", "", rownames(res1))
  
  # map based on ensembl id
  res2 <- left_join(res1, map_df, by = c("ensembl_id" = "ENSEMBL"))
  
  # replace NA in SYMBOL with ensembl_id 
  res2 <- res2 %>%
    mutate(genes = ifelse(is.na(SYMBOL), ensembl_id, SYMBOL))
  
  return(res2)
}

# apply the function above
res_KO_HET <- transform_res(res_KO_HET, map_df)
res_WT_HET <- transform_res(res_WT_HET, map_df)
res_WT_KO <- transform_res(res_WT_KO, map_df)
results <- transform_res(results, map_df)

# now return the top 5 DEGs
top_KO_HET <- head(res_KO_HET$genes[order(res_KO_HET$padj, na.last = TRUE)], n = 5)
top_WT_HET <- head(res_WT_HET$genes[order(res_WT_HET$padj, na.last = TRUE)], n = 5)
top_WT_KO <- head(res_WT_KO$genes[order(res_WT_KO$padj, na.last = TRUE)], n = 5)
top_CTR_KO <- head(results$genes[order(results$padj, na.last = TRUE)], n = 5)

# Print the top 5 DEGs
cat("Top 5 Differentially Expressed Genes in TAF10KO vs TAF10HET:", top_KO_HET, sep = "\n")
cat("Top 5 Differentially Expressed Genes in WT vs TAF10HET:", top_WT_HET, sep = "\n")
cat("Top 5 Differentially Expressed Genes in WT vs TAF10KO:", top_WT_KO, sep = "\n")

cat("Top 5 Differentially Expressed Genes in CONTROL vs TAF10KO:", top_CTR_KO, sep = "\n")

top_KO_HET[1] <- "Gm6560"
top_WT_HET[2] <- "Gm24187"
top_WT_HET[3] <- "Gvin3"
top_WT_HET[5] <- "Gm29650"
top_WT_KO[1] <- "Gm6560"
top_CTR_KO[1] <- "Gm6560"

res_KO_HET$genes[res_KO_HET$genes=="ENSMUSG00000104913"] <- "Gm6560"
res_WT_HET$genes[res_WT_HET$genes=="ENSMUSG00000088609"] <- "Gm24187"
res_WT_HET$genes[res_WT_HET$genes=="ENSMUSG00000073902"] <- "Gvin3"
res_WT_HET$genes[res_WT_HET$genes=="ENSMUSG00000099876"] <- "Gm29650"
res_WT_KO$genes[res_WT_KO$genes=="ENSMUSG00000104913"] <- "Gm6560"
results$genes[results$genes=="ENSMUSG00000104913"] <- "Gm6560"

## Volcano-Plot ##
# construct the volcano plot function:
# create a dataframe out of the results of the DEA that will be used for the volcano plot
# add the gene symbol to the data frame in order to label the top 5 differentially expressed genes

volcano_plot <- function(res, top_genes, title){
  res_volcano <- as.data.frame(res)
  
  # for the data frame: only pvalue & log2FCm columns will be utilized
  res_volcano <- res_volcano[,c("log2FoldChange", "pvalue")]
  
  # add column for differentially expressed genes (up or down-regulated)
  res_volcano$diff_expressed <- 'NO'
  res_volcano$diff_expressed[res_volcano$log2FoldChange > 0.5 & res_volcano$pvalue < 0.05] <- 'UP'
  res_volcano$diff_expressed[res_volcano$log2FoldChange < -0.5 & res_volcano$pvalue < 0.05] <- 'DOWN'
  
  # add column for top 5 differentially expressed genes
  res_volcano$top5labels <- ifelse(res$genes %in% top_genes,
                                   res$genes,
                                   NA)
  # create the plot
  ggplot(data = res_volcano, aes(x=log2FoldChange, y=-log10(pvalue),
                                 color = diff_expressed,
                                 label = top5labels)) +
    geom_vline(xintercept = c(-0.5,0.5), col = '#a1a1a1', linetype = 'dashed') +
    geom_hline(yintercept = c(-log10(0.05)), col = '#a1a1a1', linetype = 'dashed') +
    geom_point(size=1.5) +
    theme_pubr() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = 'right') +
    ggtitle(title) +
    scale_color_manual(values = c('#698fc9', '#cfcfcf', '#c96971'),
                       labels = c('Downregulated', 'Not significant', 'Upregulated')) +
    guides(color = guide_legend(title = 'Diff. expressed')) +
    coord_cartesian (xlim = c(-10, 10)) +
    labs(x=expression ("log" [2] *"FC"), y = expression ("-log" [10] *"p-value")) +
    # geom_text(vjust = -0.5, hjust = 1) +
    geom_text_repel(alpha = 0.8, box.padding = 1, show.legend = FALSE, size = 4)
}

volcano_plot(res = res_KO_HET, top_genes = top_KO_HET, title = "Volcano Plot of TAF10KO vs TAF10HET")
volcano_plot(res = res_WT_KO, top_genes = top_WT_KO, title = "Volcano Plot of WT vs TAF10KO")
volcano_plot(res = res_WT_HET, top_genes = top_WT_HET, title = "Volcano Plot of WT vs TAF10HET")
volcano_plot(res = results, top_genes = top_CTR_KO, title = "Volcano Plot of CONTROL vs TAF10HET")




## GOterm gene enrichment analysis ##
library(org.Mm.eg.db)

GOterm_plots <- function(res, condition, file.name){
  # export the list of the significant genes
  sign_df <- as.data.frame(res)
  genes_sign <- sign_df$ensembl_id[sign_df$pvalue<0.05]
  genes_sign <- genes_sign[genes_sign!='NA']
  
  # run enrichGO function to perform enrichment analysis for BP
  GO_res_BP <- enrichGO(gene = genes_sign,
                        OrgDb = 'org.Mm.eg.db',
                        keyType = 'ENSEMBL',
                        ont = 'BP')
  
  # run enrichGO function to perform enrichment analysis for CC
  GO_res_CC <- enrichGO(gene = genes_sign,
                        OrgDb = 'org.Mm.eg.db',
                        keyType = 'ENSEMBL',
                        ont = 'CC')
  
  # run enrichGO function to perform enrichment analysis for MF
  GO_res_MF <- enrichGO(gene = genes_sign,
                        OrgDb = 'org.Mm.eg.db',
                        keyType = 'ENSEMBL',
                        ont = 'MF')
  
  # create barplot for GOterm enrichment results
  plot_BP <- plot(barplot(GO_res_BP, showCategory = 20) +
                    theme_light() +
                    ggtitle(paste('BP:', condition)))
  plot_CC <- plot(barplot(GO_res_CC, showCategory = 20) +
                    theme_light() +
                    ggtitle(paste('CC:', condition)))
  plot_MF <- plot(barplot(GO_res_MF, showCategory = 20) +
                    theme_light() +
                    ggtitle(paste('MF:', condition)))
  
  pdf(paste(file.name, "BP_plot.pdf"), width = 5.71, height = 5.45)
  print(plot_BP, newpage = F)
  dev.off()
  
  pdf(paste(file.name, "CC_plot.pdf"), width = 5.71, height = 5.45)
  print(plot_CC, newpage = F)
  dev.off()
  
  pdf(paste(file.name, "MF_plot.pdf"), width = 5.71, height = 5.45)
  print(plot_MF, newpage = F)
  dev.off()
}

# apply the function above to save the enrichment analysis plots
GOterm_plots(results, "CONTROL vs TAF10KO", "CTR-KO")
