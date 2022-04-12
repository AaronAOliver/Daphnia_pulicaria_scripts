library("VennDiagram")

G25N_G15N = as.vector(read.csv("G15N_G25N_DEG.info", header = T)$x)
G15Y_G15N = as.vector(read.csv("G15Y_G15N_DEG.info", header = T)$x)
G25Y_G25N = as.vector(read.csv("G25Y_G25N_DEG.info", header = T)$x)
G25Y_G15Y = as.vector(read.csv("G15Y_G25Y_DEG.info", header = T)$x)

B25N_B15N = as.vector(read.csv("B15N_B25N_DEG.info", header = T)$x)
B15Y_B15N = as.vector(read.csv("B15Y_B15N_DEG.info", header = T)$x)
B25Y_B25N = as.vector(read.csv("B25Y_B25N_DEG.info", header = T)$x)
B25Y_B15Y = as.vector(read.csv("B15Y_B25Y_DEG.info", header = T)$x)

B15N_G15N = as.vector(read.csv("B15N_G15N_DEG.info", header = T)$x)
B15Y_G15Y = as.vector(read.csv("B15Y_G15Y_DEG.info", header = T)$x)
B25N_G25N = as.vector(read.csv("B25N_G25N_DEG.info", header = T)$x)
B25Y_G25Y = as.vector(read.csv("B25Y_G25Y_DEG.info", header = T)$x)

VennDiagram::get.venn.partitions(list(G25N_G15N=G25N_G15N, G25Y_G15Y=G25Y_G15Y, G15Y_G15N=G15Y_G15N, G25Y_G25N=G25Y_G25N))
VennDiagram::get.venn.partitions(list(B25N_B15N=B25N_B15N, B25Y_B15Y=B25Y_B15Y, B15Y_B15N=B15Y_B15N, B25Y_B25N=B25Y_B25N))
VennDiagram::get.venn.partitions(list(B15N_G15N = B15N_G15N, B15Y_G15Y = B15Y_G15Y, B25N_G25N = B25N_G25N, B25Y_G25Y = B25Y_G25Y))

pdf("Gardisky_venn_nos7.pdf", width = 11, height = 11)
grid.newpage() #141 DEGs padj <0.001, 284 DEGs padj <0.005 vs. 700 original
grid::grid.draw(VennDiagram::venn.diagram(
  list(G25N_G15N=G25N_G15N, G25Y_G15Y=G25Y_G15Y, G15Y_G15N=G15Y_G15N, G25Y_G25N=G25Y_G25N), 
  main = c("Gardisky"), main.cex = 2.25, main.fontface = "bold",
  category.names = c("25N vs. 15N", "25Y vs. 15Y", "15Y vs. 15N", "25Y vs. 25N"),
  cat.cex = 2, cat.fontface = "bold", NULL,cex = 4, lwd = 2, lty = "blank", margin.table(c(1,1,1,1)),
  fill = c("Chartreuse", "Coral", "DarkOrchid1", "Cyan")))
dev.off()

pdf("Blue_venn_nos7.pdf", height = 11, width = 11)
grid.newpage() #130 DEGs padj < 0.001, 195 DEGs padj < 0.005 vs. 151 original
grid::grid.draw(VennDiagram::venn.diagram(
  list(B25N_B15N=B25N_B15N, B25Y_B15Y=B25Y_B15Y, B15Y_B15N=B15Y_B15N, B25Y_B25N=B25Y_B25N), 
  main = c("Blue"), main.cex = 2.25, main.fontface = "bold",
  category.names = c("25N vs. 15N", "25Y vs. 15Y", "15Y vs. 15N", "25Y vs. 25N"),
  cat.cex = 2, cat.fontface = "bold", NULL,cex = 4, lwd = 2, lty = "blank",
  fill = c("Chartreuse", "Coral", "DarkOrchid1", "Cyan")))
dev.off()

pdf("Blue-x-Gardisky_venn_nos7.pdf", height = 11, width = 11)
grid.newpage() #130 DEGs padj < 0.001, 195 DEGs padj < 0.005 vs. 151 original
grid::grid.draw(VennDiagram::venn.diagram(
  list(B15N_G15N = B15N_G15N, B25Y_G25Y = B25Y_G25Y, B15Y_G15Y = B15Y_G15Y, B25N_G25N = B25N_G25N), 
  main = c("Between Clones"), main.cex = 2.25, main.fontface = "bold",
  category.names = c("15N", "25Y","15Y", "25N"),
  cat.cex = 2, cat.fontface = "bold", NULL,cex = 4, lwd = 2, lty = "blank",
  fill = c("Chartreuse", "Coral", "DarkOrchid1", "Cyan")))
dev.off()
