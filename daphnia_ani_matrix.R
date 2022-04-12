if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

library(devtools)
library(ggtree)
library(ggtreeExtra)
library(ggfittext)
library(ggplot2)
library(tidyverse)
library(ggstance)
library(gplots)
library(mudata2)


blastani = data.frame(read.csv(file = "ani_matrix.tab", sep = "\t", row.names = 1))
busco = data.frame(read.csv(file = "busco_matrix.tab", sep = "\t", row.names = 1))

blastani_orig = blastani
reverse_ani = blastani
reverse_ani[] = lapply(blastani, function(x) 100 - x) 
reverse_ani = as.data.frame(scale(reverse_ani))

distance_matrix = dist(reverse_ani, method = 'euclidean')
hclusters = hclust(distance_matrix, method = "average")
dendro = as.dendrogram(hclusters)
dendro

tax_order = rev(c("Blue",	"Gardisky",	"D.pulicaria",	"D.pulex",	"D.galeata",	"D.magna",	"E.texana"))
              blastani = blastani %>% as.data.frame() %>% rownames_to_column("hit_tax") %>%
  pivot_longer(-c("hit_tax"), names_to = "source_tax", values_to = "id")
blastani$source_tax <- factor(blastani$source_tax, levels = rev(tax_order))
blastani$hit_tax <- factor(blastani$hit_tax, levels = tax_order)
blastani$id = as.numeric(blastani$id)

busco = busco %>% as.data.frame() %>% rownames_to_column("sample") %>%
   pivot_longer(-c("sample"), names_to = "hit_type", values_to = "pct")
busco$sample = factor(busco$sample, levels = (tax_order))
busco$hit_type[busco$hit_type == "Single.copy"] = "Single copy"
busco$hit_type = factor(busco$hit_type, levels = rev(c("Single copy", "Duplicated", "Fragmented", "Missing")))


pdf("dendro.pdf", width = 12, height = 11)
ggtree(dendro) + geom_tiplab(align=TRUE,size=12, color="black", as_ylab = TRUE)
dev.off()

pdf("busco.pdf", width = 12, height = 11)
blastani$hit_tax
ggplot(blastani, aes(x=source_tax, y=hit_tax, fill = id, label = id)) + geom_tile() +
  geom_fit_text(contrast = TRUE) + scale_fill_gradient(low = "lightgrey", high = "darkgreen") + 
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(size = 14, color = "black", face="italic"),
        axis.text.y = element_text(size = 14, color = "black", face="italic"),
        panel.background = element_blank(), legend.position = "none") + scale_x_discrete(position = "top")  
dev.off()


gheatmap(p, blastani, offset=5, width=2, font.size=3, 
         colnames_angle=-45, hjust=0, high = "green", low = "white", color = "black")



pdf("busco.pdf", width = 12, height = 11) #, colormodel = "cmyk")
ggplot(busco, aes(x=sample, y=pct, fill=hit_type)) +  geom_bar(stat="identity") + coord_flip() + 
  scale_fill_manual(name = "BUSCO Match Type", values = rev(c("blue", "lightblue", "yellow", "red")), guide = guide_legend(reverse = TRUE)) + 
  ylab("Percentage of Arthropod Marker Genes") + 
  scale_y_continuous(breaks = round(seq(0, 100, by = 10),1)) +
  theme(
        #axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size = 14),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black", face="italic"),
        legend.text=element_text(size=14), legend.title = element_text(size=14),
        panel.background = element_blank()) + scale_x_discrete(position = "top")
dev.off()
