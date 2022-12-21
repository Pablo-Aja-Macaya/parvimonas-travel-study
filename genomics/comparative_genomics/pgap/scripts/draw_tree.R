suppressPackageStartupMessages(library("ape"))
suppressPackageStartupMessages(library("Biostrings"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggtree"))
suppressPackageStartupMessages(library("phytools"))

# Input path to txt tree and output path to png
args = commandArgs(trailingOnly=TRUE)
input_tree = args[1]
output_file = args[2]
input_metadata = NULL # TODO: Make this take the metadata input
# input_tree = "/home/usuario/Proyectos/Results/ARGA/ARGA_ALL_v2/kpneumoniae_comparative_analysis/03_pirate/raxml/raxml/RAxML_bestTree.tree"
# output_file = "/home/usuario/Proyectos/Results/ARGA/ARGA_ALL_v2/kpneumoniae_comparative_analysis/03_pirate/raxml/RAxML_bestTree.png"
# input_metadata = "/home/usuario/Proyectos/Results/ARGA/ARGA_ALL_v2/reports/03_multi_tabular_reports/multi_report.tsv" # TODO: Make this take the metadata input

# Read tree and save to output_file
tree <- read.tree(input_tree)
if ('Reference' %in% tree$tip.label){
  tree <- midpoint.root(tree)
} else {
  tree <- tree
}

# Draw base tree
p <- ggtree(tree) + 
  geom_treescale() +
  geom_tiplab(size=3) + hexpand(.2) +
  coord_cartesian(clip="off")

if (!is.null(input_metadata)){
  # Read metadata
  metadata <- read.csv2(input_metadata, sep='\t')
  metadata <- metadata[metadata$extraction %in% tree$tip.label,]
  metadata <- metadata[,c("extraction","biosample","ST","seq_run")]
  
  # Add STs to plot
  p <- p %<+% metadata + 
    geom_tiplab(aes(fill = factor(ST)),
                color = "black", # color for label font
                geom = "label",  # labels not text
                label.padding = unit(0.1, "lines"), # amount of padding around the labels
                label.size = 0.02) + # size of label border
    theme(legend.position = "bottom") + 
    guides(
      fill = guide_legend(
        title = "ST",
        override.aes = aes(label = "")
      )
    )
}

p
ggsave(output_file, width=2500, height=3000, units="px")






# 
# metadata$bc <- metadata$extraction
# metadata$ST <- as.factor(metadata[["ST"]])
# 
# x <- ggtree(tree) %<+% metadata + 
#   geom_tippoint(aes(color=ST)) +
#   geom_tiplab(aes(label=bc), align=T, size=3, offset=0.0001) +
#   geom_tiplab(aes(label=biosample), align=T, offset=0.001*1.5, linetype=NA, size=3) +
#   geom_tiplab(aes(label=ST), align=T, offset=0.001*2.7, linetype=NA, size=3)
# x

# 
# library(randomcoloR)
# color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)] # Creando vector de colores contrastantes
# heatmap.colours <- distinctColorPalette(length(levels(as.factor(metadata[["ST"]]))))
# names(heatmap.colours) <- levels(as.factor(metadata[["ST"]]))
# 
# metadata$st.color <- heatmap.colours[match(metadata$ST, names(heatmap.colours))]
# 
# gheatmap(x, metadata[,"ST",drop=F], offset = 0.001*3, width=0.05, color=NULL, 
#          colnames_position="top", 
#          colnames_angle=90, colnames_offset_y = 1, 
#          hjust=1, font.size=2) +
#   scale_fill_manual(values=metadata$st.color)
# 
# 

