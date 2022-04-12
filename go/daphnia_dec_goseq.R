if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("goseq")

library(goseq)

split_gos <- function(i) {
  if (i != "") {
    strsplit(i, split=";")
  }
  else {
    ""
  }
}

gene_gos = read.csv("pulicaria_go_len.info", sep = "\t")
gene_lengths = as.integer(gene_gos$length)
gene_names = gene_gos$gene
names(gene_lengths) = gene_names
gene_gos$GO = sapply(gene_gos$GO, split_gos)

gene_go_mapping = read.csv("pulicaria_go_many2many.info", sep = "\t")
  
##########3

deg_blue_temp = read.csv("B15Y_B25Y_DEG.info", sep = "\t")
deg_blue_temp = as.integer(gene_names %in% deg_blue_temp$x)
names(deg_blue_temp) = gene_names

pwf_blue_temp = nullp(deg_blue_temp, bias.data = gene_lengths)
goseq_blue_temp = goseq(pwf_blue_temp, gene2cat = gene_go_mapping, method = "Wallenius")
goseq_blue_temp = goseq_blue_temp[goseq_blue_temp$numDEInCat > 0,]
goseq_blue_temp$over_represented_padjust = p.adjust(goseq_blue_temp$over_represented_pvalue, method = "fdr")
enriched_blue_temp = goseq_blue_temp[goseq_blue_temp$over_represented_padjust < 0.05, ]

write.table(enriched_blue_temp, file = "blue_go_15Y25Y.tab", sep = "\t", quote = FALSE, row.names = FALSE)

#########

deg_blue_temp = read.csv("B15N_B25N_DEG.info", sep = "\t")
deg_blue_temp = as.integer(gene_names %in% deg_blue_temp$x)
names(deg_blue_temp) = gene_names

pwf_blue_temp = nullp(deg_blue_temp, bias.data = gene_lengths)
goseq_blue_temp = goseq(pwf_blue_temp, gene2cat = gene_go_mapping, method = "Wallenius")
goseq_blue_temp = goseq_blue_temp[goseq_blue_temp$numDEInCat > 0,]
goseq_blue_temp$over_represented_padjust = p.adjust(goseq_blue_temp$over_represented_pvalue, method = "fdr")
enriched_blue_temp = goseq_blue_temp[goseq_blue_temp$over_represented_padjust < 0.05, ]

write.table(enriched_blue_temp, file = "blue_go_15N25N.tab", sep = "\t", quote = FALSE, row.names = FALSE)

############

deg_gard_temp = read.csv("G15Y_G25Y_DEG.info", sep = "\t")
deg_gard_temp = as.integer(gene_names %in% deg_gard_temp$x)
names(deg_gard_temp) = gene_names
#head(deg_gard_temp)

pwf_gard_temp = nullp(deg_gard_temp, bias.data = gene_lengths)
goseq_gard_temp = goseq(pwf_gard_temp, gene2cat = gene_go_mapping, method = "Wallenius")
goseq_gard_temp = goseq_gard_temp[goseq_gard_temp$numDEInCat > 0,]
goseq_gard_temp$over_represented_padjust = p.adjust(goseq_gard_temp$over_represented_pvalue, method = "fdr")
enriched_gard_temp = goseq_gard_temp[goseq_gard_temp$over_represented_padjust < 0.05, ]

write.table(enriched_gard_temp, file = "gard_go_15Y25Y.tab", sep = "\t", quote = FALSE, row.names = FALSE)

########

deg_gard_temp = read.csv("G15N_G25N_DEG.info", sep = "\t")
deg_gard_temp = as.integer(gene_names %in% deg_gard_temp$x)
names(deg_gard_temp) = gene_names
#head(deg_gard_temp)

pwf_gard_temp = nullp(deg_gard_temp, bias.data = gene_lengths)
goseq_gard_temp = goseq(pwf_gard_temp, gene2cat = gene_go_mapping, method = "Wallenius")
goseq_gard_temp = goseq_gard_temp[goseq_gard_temp$numDEInCat > 0,]
goseq_gard_temp$over_represented_padjust = p.adjust(goseq_gard_temp$over_represented_pvalue, method = "fdr")
enriched_gard_temp = goseq_gard_temp[goseq_gard_temp$over_represented_padjust < 0.05, ]

write.table(enriched_gard_temp, file = "gard_go_15N25N.tab", sep = "\t", quote = FALSE, row.names = FALSE)

