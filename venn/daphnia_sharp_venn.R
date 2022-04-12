library("VennDiagram")

par(xpd = NA)
grid.newpage()
grid::grid.draw(venn.diagram(x = NULL, filename = NULL, direct.area = T, area.vector = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),
             category = c("test1", "test2", "test3", "test4"), cat.cex = 2, cex = 3,
             main = c("Title"), main.cex = 2, main.fontfamily = "serif"))

print_venn <- function(count_vec = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 
                       cat_vec = c("test1", "test2", "test3", "test4"),
                       fill_vec = c("Chartreuse", "Coral", "DarkOrchid1", "Cyan"),
                       title = "Title", filename = "test") {
  pdf(paste0(filename, "_venn.pdf"), width = 12, height = 11) #, colormodel = "cmyk")
  grid.newpage()
  grid::grid.draw(venn.diagram(x = NULL, filename = NULL, direct.area = T, area.vector = count_vec,
                               category = cat_vec, cat.cex = 1.7, cex = 4,
                               main = c(title), main.cex = 4,
                               main.fontfamily = "serif", main.fontface = "bold", main.pos = c(0.5, 0.8),
                               lty = "blank",cat.fontface = "bold", fill = fill_vec,
                               resolution = 600, margin = 0.75,
                               disable.logging = TRUE))
  dev.off()
}

print_venn()

# DEG counts - Blue
print_venn(count_vec = c(0,0,0,0,0,0,0,0,221,4,0,0,0,28,67), 
           cat_vec = c("15N vs. 25N", "15Y vs. 25Y", "15N vs. 15Y", "25N vs. 25Y"),
           title = "Blue", filename = "Blue")


# DEG counts - Gardisky
print_venn(count_vec = c(0,0,0,0,0,0,0,0,396,0,0,0,0,52,127), 
           cat_vec = c("15N vs. 25N", "15Y vs. 25Y", "15N vs. 15Y", "25N vs. 25Y"),
           title = "Gardisky", filename = "Gardisky")


# Gene orthology - Puli DEGs vs. Pulex DEGs
print_venn(count_vec = c(21,23,44,61,38,1,0,1,231,157,2,0,0,0,0), 
           cat_vec = c("D. pulex (Temperature)", "D. pulicaria (Kairomones)", "D. pulex (Kairomones)", "D. pulicaria (Temperature)"),
           title = "Orthology: DEGs vs. DEGs", filename = "DEG_Ortho", fill = c("greenyellow", "firebrick3", "darkgoldenrod4", "deeppink"))


# Gene orthology - Puli genome vs. Pulex DEGs
print_venn(count_vec = c(2922,1619,41,227,61,2,1), 
                      cat_vec = c("D. pulicaria (Genome)", paste(italic("D. pulex"), "(Temperature)"), "D. pulex (Kairomones)"),
                      title = "Orthology: Genome vs. DEGs", filename = "Genome_Ortho", fill = c("deeppink", "greenyellow", "darkgoldenrod4"))

