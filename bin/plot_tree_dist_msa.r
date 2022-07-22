#!/usr/bin/env Rscript


library(argparse)

usage <- function() {
  usage.text <- '\nUsage: script.R --treefile <path to the iqtree file --snps <path to the snps csv file --prefix <prefix to the output filename>\n\n'
  return(usage.text)
}

parser <- ArgumentParser()
parser$add_argument("--treefile", default=NULL, help="path to the iqtree file (.treefile)")
parser$add_argument("--alignment", default=NULL, help="path to the alignment FASTA file")
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
library(Biostrings)

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

if (is.null(args$alignment)) {
  parser$print_help()
  stop("Please provide the alignment FASTA file with sequence identifiers similar to tree tiplabels", call.=FALSE)
}
if (!file.exists(args$alignment)) {
  parser$print_help()
  stop(paste("The following alignment file don't exist:", 
             paste(args$alignment, sep='', collapse=' '), sep=' '), call.=FALSE)
} else {
  alignment_file <- args$alignment
}

prefix <- args$prefix

# read the phylogenetic tree
tree <- read.iqtree(file = tree_file)
x <- as_tibble(tree)
x$label <- noquote(x$label)

xlin <- x[is.na(x$SH_aLRT) & is.na(x$UFboot),]
xnlin <- x[!is.na(x$SH_aLRT) & !is.na(x$UFboot),]


# split on label column
xlin$lineage <- sub(".*\\|", "", xlin$label)
xnlin$lineage <- "-"
xx <- rbind(xlin, xnlin)



tree <- groupOTU(tree, list(A = as.vector(xx[xx$lineage == "A",]$label),
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


outfile <- paste(prefix, ".tree.pdf", sep = "")
pdf(file=outfile, w=12, h=8, bg='white')
pt <- ggtree(tree, ladderize = T, open.angle=15, aes(color=group)) %>% rotate(rootnode(tree)) +
  geom_tiplab(aes(label = paste0("italic('", label, "')")), parse = TRUE, size = 3.0, align=TRUE) +
  geom_tippoint(mapping=aes(color=group), size=1.5, show.legend=FALSE) +
  geom_treescale(x = 0.05, y = 1, width = 0.005, fontsize = 3) + 
  scale_color_manual(values = c(colors, "black"), breaks = c("A", "B", "C", "D", "E", "F", "G", "H",
                                                           "I", "J", "K", "L", "M", "N", "O")) +
  ggplot2::guides(color = guide_legend(title="Lineage", override.aes = list(size = 5, shape = 15))) +
  theme_tree2()
print(pt)
dev.off()




p <- ggtree(tree, aes(color = group), ladderize = TRUE) %>% rotate(rootnode(tree)) +
  geom_tiplab(aes(label = paste0("italic('", label, "')")), parse = TRUE, size = 2.5) +
  geom_tippoint(mapping=aes(color=group), size=1.5, show.legend=FALSE) +
  geom_treescale(x = 0.07, y = 1, width = 0.003, fontsize = 3) + 
  scale_color_manual(values = c(colors, "black"), na.value = "black", name = "Lineage",
                     breaks = c("A", "B", "C", "D", "E", "F", "G", "H",
                                "I", "J", "K", "L", "M", "N", "O")) +
  ggplot2::guides(color = guide_legend(override.aes = list(size = 5, shape = 15))) +
  theme_tree2()



# get tip labels
# tip_labels <- tree$tip.label
tip_labels <- tree@phylo$tip.label
accessions <- c(sub("\\|.$", "", tip_labels))
names(tip_labels) <- accessions

aln <- readDNAStringSet(filepath = alignment_file)

# tipseq_aln <- NULL
# for (i in aln@ranges@NAMES) {
#   if (width(aln[i]) > 815 | width(aln[i]) > 1306 ) {
#     tipseq_aln = subseq(aln, start = 816, end = 1305)
#   } else {
#     tipseq_aln = aln
#   }
# }
tipseq_aln <- DNAStringSet(aln)


# calculate pairwise hamming distances among sequences
tipseq_dist <- stringDist(tipseq_aln, method = "hamming")

# calculate percentage of differences
tipseq_d <- as.matrix(tipseq_dist) / width(tipseq_aln[1]) * 100

## convert the matrix to tidy data frame for facet_plot
dd <- as_tibble(tipseq_d)
dd$seq1 <- rownames(tipseq_d)

library("tidyr")
td <- gather(dd, seq2, dist, -seq1)
pdata <- p$data[p$data$group != 0,]
g <- pdata$group
names(g) <- pdata$label
td$clade <- g[td$seq2]

# visualize the sequence differences using dot plot and line plot
# and align the sequence difference plot to the tree using facet_plot
outfile <- paste(prefix, ".tree.dist.pdf", sep = "")
pdf(file=outfile, w=12, h=8, bg='white')
p2 <- facet_plot(p, panel = "Sequence Distance", data = td, geom_point, 
                 mapping = aes(x = dist, color = clade, shape = clade), 
                 alpha = .6) %>%
  facet_plot(panel = "Sequence Distance", data = td, geom = geom_path, 
             mapping=aes(x = dist, group = seq2, color = clade), alpha = .6) + 
  scale_shape_manual(values = 1:15, guide = "none")
print(p2)
dev.off()

################################################################################
p_aln <- ggtree(tree, aes(color = group), ladderize = TRUE) %>% rotate(rootnode(tree)) +
  geom_tiplab(aes(label = paste0("italic('", label, "')")), parse = TRUE, size = 2.5) +
  geom_treescale(x = 0.07, y = 1, width = 0.003, fontsize = 3) +
  scale_color_manual(values = c(colors, "black"), na.value = "black", name = "Lineage",
                     breaks = c("A", "B", "C", "D", "E", "F", "G", "H",
                                "I", "J", "K", "L", "M", "N", "O")) +
  ggplot2::guides(color = guide_legend(override.aes = list(size = 5, shape = 15))) +
  theme_tree2()


for (i in aln@ranges@NAMES) {
  if (width(aln[i]) > 815 | width(aln[i]) > 1306 ) {
    p <- msaplot(p_aln, alignment_file, offset=0.1, width=30, window=c(816, 1306))
  } else {
    p <- msaplot(p_aln, alignment_file, offset=0.1, width=30)
  }
}

outfile <- paste(prefix, ".tree.msa.pdf", sep = "")
pdf(file=outfile, w=12, h=8, bg='white')
print(p)
dev.off()



################################################################################