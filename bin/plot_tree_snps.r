#!/usr/bin/env Rscript


library(argparse)

usage <- function() {
  usage.text <- '\nUsage: script.R --treefile <path to the iqtree file --snps <path to the snps csv file --prefix <prefix to the output filename>\n\n'
  return(usage.text)
}

parser <- ArgumentParser()
parser$add_argument("--treefile", default=NULL, help="path to the iqtree file (.treefile)")
parser$add_argument("--snps", default=NULL, help="path to the snps csv file")
parser$add_argument("--metadata", default=NULL, help="path to the metadata file having the column: lineage (with values letters A-O)")
parser$add_argument("--prefix", default=NULL, help="prefix to the output filename")
args <- parser$parse_args()

###########################################
#######        load packages        #######
###########################################


library(treeio)
library(ggtree)
library(ggplot2)
library(tidytree)
library(RColorBrewer)

###########################################
#######           checks           #######
###########################################

if (is.null(args$treefile)) {
  parser$print_help()
  stop("Please provide the iqtree generated tree file", call.=FALSE)
}
if (!file.exists(args$treefile)) {
  parser$print_help()
  stop(paste("The following tree file don't exist:", 
             paste(args$treefile, sep='', collapse=' '), sep=' '), call.=FALSE)
} else {
  tree_file <- args$treefile
}

if (is.null(args$snps)) {
  parser$print_help()
  stop("Please provide the SNPs csv file with headers as 'sequence identifiers' (similar to tree tiplabels)", call.=FALSE)
}
if (!file.exists(args$snps)) {
  parser$print_help()
  stop(paste("The following SNPs file don't exist:", 
             paste(args$snps, sep='', collapse=' '), sep=' '), call.=FALSE)
} else {
  snps_file <- args$snps
}

if (is.null(args$metadata)) {
  parser$print_help()
  stop("Please provide the metadata csv file with headers 'lineage' (values = letters A-O)", call.=FALSE)
}
if (!file.exists(args$metadata)) {
  parser$print_help()
  stop(paste("The following metadata file don't exist:", 
             paste(args$metadata, sep='', collapse=' '), sep=' '), call.=FALSE)
} else {
  metadata_file <- args$metadata
}

prefix <- args$prefix

# read the phylogenetic tree
tre <- read.iqtree(file = tree_file)

x <- as_tibble(tre)
x$label <- noquote(x$label)

xlin <- x[is.na(x$SH_aLRT) & is.na(x$UFboot),]
xnlin <- x[!is.na(x$SH_aLRT) & !is.na(x$UFboot),]


# split on label column
xlin$lineage <- sub(".*\\|", "", xlin$label)
xnlin$lineage <- "-"
xx <- rbind(xlin, xnlin)



tre <- groupOTU(tre, list(A = as.vector(xx[xx$lineage == "A",]$label),
                            B = as.vector(xx[xx$lineage == "B",]$label),
                            C = as.vector(xx[xx$lineage == "C",]$label),
                            D = as.vector(xx[xx$lineage == "D",]$label),
                            E = as.vector(xx[xx$lineage == "E",]$label),
                            F = as.vector(xx[xx$lineage == "F",]$label),
                            G = as.vector(xx[xx$lineage == "G",]$label),
                            H = as.vector(xx[xx$lineage == "H",]$label),
                            I = as.vector(xx[xx$lineage == "I",]$label),
                            J = as.vector(xx[xx$lineage == "J",]$label),
                            K = as.vector(xx[xx$lineage == "K",]$label),
                            L = as.vector(xx[xx$lineage == "L",]$label),
                            M = as.vector(xx[xx$lineage == "M",]$label),
                            N = as.vector(xx[xx$lineage == "N",]$label),
                            O = as.vector(xx[xx$lineage == "O",]$label)
))



# visualize the tree
nb.cols <- length(letters[1:15])
cols <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

colors <- c("#FFC125","#87CEFA","#7B68EE","#808080","#800080",
            "#9ACD32","#D15FEE","#FFC0CB","#EE6A50","#8DEEEE",
            "#006400","#800000","#B0171F","#191970", "#90BBD0")

#read the sampling info data
metadata <- read.csv(metadata_file, header = TRUE, sep = ",")
# remove duplicates and empty accession columns
metadata <- metadata[!duplicated(metadata$accession), ]
metadata <- metadata[!is.na(metadata$accession),]
metadata$id <- paste(metadata$accession, sep = "|", metadata$lineage)
metadata <- metadata[which(metadata$id %in% tre@phylo$tip.label),]

# read and process snps data
SNPs <- read.csv(snps_file, header = T, stringsAsFactor = F)
gapChar <- "?"
SNP <- t(SNPs)
lSNP <- apply(SNP, 1, function(x) {
  x != SNP[1,] & x != gapChar & SNP[1,] != gapChar
})
lSNP <- as.data.frame(lSNP)
lSNP$pos <- as.numeric(rownames(lSNP))
lSNP <- tidyr::gather(lSNP, name, value, -pos)
SNP_data <- lSNP[lSNP$value, c("name", "pos")]
SNP_data$name <- chartr(".", "|", SNP_data$name)


# visualize the tree
nb.cols <- length(unique(metadata$lineage))
cols <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

p <- ggtree(tre, ladderize = T, open.angle=15, aes(color=group)) %>% rotate(rootnode(tre)) +
  geom_tiplab(aes(label = paste0("italic('", label, "')")), parse = TRUE, size = 3.0, align=TRUE) +
  geom_tippoint(mapping=aes(color=group), size=1.5, show.legend=FALSE) +
  geom_treescale(x = 0.05, y = 1, width = 0.005, fontsize = 3) + 
  scale_color_manual(values = c(colors, "black"), breaks = c("A", "B", "C", "D", "E", "F", "G", "H",
                                                           "I", "J", "K", "L", "M", "N", "O")) +
  ggplot2::guides(color = guide_legend(title="Lineage", override.aes = list(size = 5, shape = 15))) +
  theme_tree2()

# add sampling information data set and add symbols colored by country
p <- p %<+% metadata + geom_tippoint(mapping = aes(color=group)) 

# visualize SNP data using dot and bar charts,
# and align them based on tree structure
outfile <- paste(prefix, ".tree.snps.pdf", sep = "")
pdf(file=outfile, w=12, h=8, bg='white')
p2 <- p + geom_facet(panel = "SNP", data = SNP_data, geom = geom_point,
                mapping=aes(x = pos, color = lineage), shape = '|') +
  theme_tree2(legend.position=c(.04, .30))
print(p2)
dev.off()