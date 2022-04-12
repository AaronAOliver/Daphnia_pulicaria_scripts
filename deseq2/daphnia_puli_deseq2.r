# INSTALL REQUIRED PACKAGES
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("tximport")
#BiocManager::install("DESeq2")

# Import Dependencies
library("DESeq2")
library("ggplot2")

# Read input files
samples = read.delim("samples.info", header = TRUE, sep = "\t")
counts = as.matrix(read.delim("featureCounts.matrix", sep="\t", row.names = "Gene"))

# Factorize Data
samples$Lake = factor(samples$Lake)
samples$Temperature = factor(samples$Temperature)
samples$Fish = factor(samples$Fish)
samples$Treatment = factor(samples$Treatment)
samples$Condition = factor(samples$Condition)

# Run DESeq Analysis
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = samples,
                              design = ~ Treatment)

dds = DESeq(dds)

vsd1 <- vst(dds, blind=FALSE)

# Filter out low (normalized) readcounts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

vsd <- vst(dds, blind=FALSE)

# calculate the variance for each gene
rv <- rowVars(assay(vsd))

# select the ntop genes by variance
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]



#rld <- rlog(dds, blind=FALSE)

pdf("PCA_topall.pdf", width = 10, height = 10)
pcaData <- plotPCA(vsd, intgroup=c("Lake","Condition"), returnData=TRUE, ntop=200000)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Lake, shape=Condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + theme_bw() + scale_color_manual(values = c("blue", "orange"))
dev.off()

# Loop through all possible pairs of treatments
treatment_list = c("B15Y", "G15Y", "B15N", "G15N", "B25Y", "G25Y", "B25N", "G25N")
for (i in 1:(length(treatment_list) - 1)) {
  default_treatment = treatment_list[[i]]
  for (j in (i + 1):length(treatment_list)) {
    compare_treatment = treatment_list[[j]]
    
    # Output analysis for all genes from a pairwise comparison
    res.pairwise = results(dds, contrast=c("Treatment", default_treatment, compare_treatment))
    write.table(res.pairwise, file = paste("", default_treatment,"_", compare_treatment,"_all.tab", sep = ""), sep = "\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
    
    # Output analysis for just the DEGs from a pairwise comparison
    deg = res.pairwise[!is.na(res.pairwise$padj) & res.pairwise$padj < 0.05 & abs(res.pairwise$log2FoldChange) > 1,]
    write.table(deg, file = paste("", default_treatment,"_", compare_treatment,"_DEG.tab", sep = ""), sep = "\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
    write.table(rownames(deg), file = paste("", default_treatment,"_", compare_treatment,"_DEG.info", sep = ""), quote=FALSE, row.names=F, col.names=T)
    
  }
}

