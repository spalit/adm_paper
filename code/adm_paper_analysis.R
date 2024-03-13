## Bulk RNA-Seq from whole lung tissue of mice
## treated with LPS [Lipopolysaccharide] or control PBS [Phosphate-buffered saline] 
## and are wild type or Adm (Adrenomedullin) haplo insufficient (only 1 copy left)
## There are 4 conditions (3 replicates each)
## Standard bulk analysis, differential expression, PCA, gene set enrichment

# RESOURCES
## Adding interaction term https://support.bioconductor.org/p/46153/
## https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
## https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html

# Load R libs ####
library(dplyr)
library(edgeR)
library(sva)
library(limma)
library(reshape2)
library(ggplot2)
library(viridis)
library(latex2exp)
library(stringr)
library(ggrepel)
library(pheatmap)
library(gridExtra)
library(RColorBrewer)

# Load lung cell marker lists from Tabula Muris ####
file_path <- "../data/Tabula_Muris.txt"
data_list <- list()

con <- file(file_path, "r")

while (length(line <- readLines(con, n = 1)) > 0) {
  fields <- strsplit(line, "\t")
  df <- as.data.frame(matrix(unlist(fields)[-c(1,2)]), nrow = 1, byrow = TRUE)
  
  df_list <- list(df)
  names(df_list) <- unlist(fields)[1]
  data_list <- c(data_list, df_list)
}

close(con)

lung_data_list <- data_list[grep('Lung',names(data_list),value = T)]

# Load bulk RNA seq data (Binoy) ####
metadata <- readxl::read_xlsx('../data/annotation.xlsx', sheet = 1)[,c(1,2)]
metadata$Group %>% unique #4

data <- read.csv('../data/combined-rnaseq.featureCounts-genes_Binoy.csv')
data$GeneBiotype %>% unique

# Fig 1 - PCA ####
data_pc <- data[which(data$GeneBiotype == 'protein_coding'),]
data_pc[,c(1,2)] %>% unique %>% dim

rownames(data_pc) <- data_pc$GeneID
data_pc <- data_pc[,c(4:dim(data_pc)[2])]

metadata$Sample <- gsub('-','',metadata$Sample)
metadata$Group <- gsub('-','_',metadata$Group)

# Both true continue with analysis
all(metadata$Sample %in% colnames(data_pc))
all(metadata$Sample == colnames(data_pc))

# which(rowSums(data_pc) == 0) %>% length
# 2807 genes (symbols) have no expression in any of the 12 samples
# 2817 ensembl IDs have no expression in any of the 12 samples

# Create DEGList object
# Use the group as vector to indicate each sample's group type
DEGL <- DGEList(counts=data_pc, group=metadata$Group)
dim(DEGL)
table(rowSums(DEGL$counts==0)==12)

# Remove genes without enough counts for statistical analysis
keep.exprs <- filterByExpr(DEGL, group=metadata$Group)
DEGL <- DEGL[keep.exprs,, keep.lib.sizes=FALSE]
dim(DEGL)

DEGL <- calcNormFactors(DEGL, method = "TMM")
DEGL$samples$norm.factors

lcpm <- cpm(DEGL, log=TRUE, normalized.lib.sizes = TRUE)

plot_pca <- function(input,title){
  input_pca <- prcomp(t(input))
  input_pca_out <- as.data.frame(input_pca$x)
  
  if(all(rownames(input_pca_out) == metadata$Sample) & all(rownames(input_pca_out) %in% metadata$Sample))
  { 
    input_pca_out$Group <- metadata$Group}
  
  percentage <- round(input_pca$sdev / sum(input_pca$sdev) * 100, 2)
  percentage <- paste(colnames(input_pca_out), "(", paste( as.character(percentage), "%", ")", sep=""))
  
  p <- ggplot(input_pca_out, aes(x=PC1,y=PC2, color = Group)) +
    geom_point() + 
    xlab(percentage[1]) + 
    ylab(percentage[2])  + 
    labs(title=title) +
    theme_classic()
  return(p)
}
plot_pca(lcpm,'PCA of gene expression')

ggsave(
  '/Users/lukas/OneDrive/Miko/THINC/projects/Shivanna/adm_paper/Fig1_pca.pdf', 
  width = 6, height = 5)

# Model with each coefficient corresponding to a group mean - NO INTERACTION
group <- as.factor(metadata$Group)
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design

contr.matrix <- makeContrasts(
  AMKOLPSvsPBS = AMKO_LPS - AMKO_PBS,
  WTLPSvsPBS = WT_LPS - WT_PBS,
  Diff = (AMKO_LPS - AMKO_PBS)-(WT_LPS - WT_PBS),
  levels = colnames(design))
contr.matrix

par(mfrow=c(1,2))
v <- voom(DEGL, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

resultz <- summary(decideTests(efit))
vennDiagram(resultz)

tfit <- treat(vfit, lfc=0.5)
dt <- decideTests(tfit)
summary(dt)

AMKOLPSvsPBS <- topTreat(tfit, coef=1, n=Inf)
head(AMKOLPSvsPBS)

WTLPSvsPBS <- topTreat(tfit, coef=2, n=Inf)
head(WTLPSvsPBS)

diff <- topTreat(tfit, coef=3, n=Inf)
head(diff)

# Fig 1 - top 50 heatmap ####
Top50 = topTable(efit,number=50,adjust="BH")
str(Top50)

Xtop = lcpm[rownames(Top50),] %>% as.data.frame

all(rownames(Xtop) == rownames(Top50))
all(rownames(Xtop) %in% rownames(Top50))

metadata01 <- as.data.frame(metadata)
rownames(metadata01) <- metadata01$Sample
metadata01 <- metadata01[colnames(Xtop),]
metadata01 <- metadata01[,-1,drop = F]

Xtop$gene <- data$GeneSymbol[match(rownames(Xtop),data$GeneID)]
rownames(Xtop) <- Xtop$gene
Xtop <- Xtop[,-which(colnames(Xtop)=='gene')]

pheatmap(Xtop, scale = "row", border_color = NA, fontsize_row = 6, 
         fontsize_col = 6, annotation_col = metadata01, 
         cutree_rows = 2, cutree_cols = 2,
         show_colnames = F, show_rownames = T,
         breaks = seq(-1.5, 1.5, length = 100))


# Fig 1 - Volcano plot Adm KO ####
AMKOLPSvsPBS$diffexpressed <- "NO"
AMKOLPSvsPBS$diffexpressed[AMKOLPSvsPBS$logFC > 0.5 & AMKOLPSvsPBS$adj.P.Val < 0.05] <- "UP"
AMKOLPSvsPBS$diffexpressed[AMKOLPSvsPBS$logFC < -0.5 & AMKOLPSvsPBS$adj.P.Val < 0.05] <- "DOWN"

AMKOLPSvsPBS$delabel <- data$GeneSymbol[match(rownames(AMKOLPSvsPBS),data$GeneID)]

temp <- AMKOLPSvsPBS[which(AMKOLPSvsPBS$logFC > 0),]
up_delabels <- temp[order(temp$adj.P.Val),] %>% rownames %>% head(10)
temp <- AMKOLPSvsPBS[which(AMKOLPSvsPBS$logFC < 0),]
down_delabels <- temp[order(temp$adj.P.Val),] %>% rownames %>% head(10)

AMKOLPSvsPBS$delabel[which(!rownames(AMKOLPSvsPBS) %in% c(up_delabels,down_delabels))] <- NA


ggplot(data=AMKOLPSvsPBS, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, 
                              label=delabel)) +
  geom_point(size = 0.95) + 
  geom_text_repel() +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept=c(-0.5, 0.5), col="red", lty = 23) +
  geom_hline(yintercept=-log10(0.05), col="red", lty = 23) + theme_classic() + 
  coord_equal(ratio = 4)


# Fig 1 - Volcano plot WT ####
WTLPSvsPBS$diffexpressed <- "NO"
WTLPSvsPBS$diffexpressed[WTLPSvsPBS$logFC > 0.5 & WTLPSvsPBS$adj.P.Val < 0.05] <- "UP"
WTLPSvsPBS$diffexpressed[WTLPSvsPBS$logFC < -0.5 & WTLPSvsPBS$adj.P.Val < 0.05] <- "DOWN"

WTLPSvsPBS$delabel <- data$GeneSymbol[match(rownames(WTLPSvsPBS),data$GeneID)]

temp <- WTLPSvsPBS[which(WTLPSvsPBS$logFC > 0),]
up_delabels <- temp[order(temp$adj.P.Val),] %>% rownames %>% head(10)
temp <- WTLPSvsPBS[which(WTLPSvsPBS$logFC < 0),]
down_delabels <- temp[order(temp$adj.P.Val),] %>% rownames %>% head(10)

WTLPSvsPBS$delabel[which(!rownames(WTLPSvsPBS) %in% c(up_delabels,down_delabels))] <- NA

ggplot(data=WTLPSvsPBS, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, label=delabel)) +
  geom_point(size = 0.95) + 
  geom_text_repel() +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept=c(-0.5, 0.5), col="red", lty = 23) +
  geom_hline(yintercept=-log10(0.05), col="red", lty = 23) + theme_classic() + 
  coord_equal(ratio = 4)


# Fig 1 - boxplots ####
plot_gene <- function(gene){
  geneid_from_publ <- data$GeneID[match(gene, data$GeneSymbol)]
  
  foo <- data.frame(
    group = metadata$Group,
    expr = lcpm[which(rownames(lcpm) == geneid_from_publ),])
  
  foo %>%
    ggplot(aes(x=group, y=expr, fill=group)) +
    geom_boxplot() + geom_point() +
    ggtitle(gene) + 
    theme_classic()
}
p1 <- plot_gene("Slit2")
p2 <- plot_gene("Mmp9")

p <- grid.arrange(p1, p2, ncol=2)

ggsave(
  p,
  filename = "/Users/lukas/OneDrive/Miko/THINC/projects/Shivanna/adm_paper/fig1_boxplots.pdf",
  width = 8, height = 4)


# Fig 2 - qPCR/protein data validation ####
df <- readxl::read_excel(
  '../data/Gene and protein expression results from previous publication.xlsx',col_types = "numeric", sheet = 1)
df <- subset(df, rowSums(is.na(df)) != ncol(df))

foo <- c(1:5)
my_list <- list()

genes <- c('CCL2','CCL3','Cxcl1','ICAM1','IL-1B','TNFa','stat3','stat1')
merged_old <- do.call(rbind, lapply(genes, function(each){
  my_list[[each]] <- df[foo,] %>% as.data.frame
  colnames(my_list[[each]]) <- c("Adm+/+mice_PBS","Adm+/-mice_PBS","Adm+/+mice_L6","Adm+/-mice_L6")
  foo <<- foo + 5
  
  tmp_df <- my_list[[each]]
  tmp_df <- tmp_df[,c("Adm+/-mice_L6","Adm+/-mice_PBS","Adm+/+mice_L6","Adm+/+mice_PBS")]
  
  tmp_df <- reshape2::melt(tmp_df)
  tmp_df$gene <- each
  
  tmp_df
}))
merged_old$gene <- factor(
  merged_old$gene, levels = genes)

genes_from_publ <- c('Ccl2','Ccl3','Cxcl1','Icam1','Il1b','Tnf','Stat3','Stat1')
geneid_from_publ <- data$GeneID[match(genes_from_publ,data$GeneSymbol)] #data$GeneID[which(data$GeneSymbol %in% genes_from_publ)]

merged_new <- do.call(rbind, lapply(geneid_from_publ, function(each){
  foo <- cbind.data.frame(group = metadata$Group, 
                          expr = lcpm[which(rownames(lcpm) == each),])
  foo$gene <- each
  foo
}))
merged_new$gene <- factor(
  merged_new$gene, levels = geneid_from_publ)

labels <- c(
  "WT_LPS","AMKO_LPS","WT_PBS","AMKO_PBS")
names(labels) <- c(
  "Adm+/+mice_L6","Adm+/-mice_L6","Adm+/+mice_PBS","Adm+/-mice_PBS")
merged_old$group <- labels[as.character(merged_old$variable)]

p2 <- ggplot(merged_new, aes(x=group, y=expr, fill=group)) +
  geom_boxplot() +
  facet_wrap(~gene, scales = "free_y", ncol = 1) +
  theme_classic() +
  theme(axis.text.x = NULL)

p1 <- ggplot(merged_old, aes(x=group, y=value, fill=group)) +
  geom_boxplot() +
  facet_wrap(~gene, scales = "free_y", ncol = 1) +
  theme_classic() +
  theme(axis.text.x = NULL)

p <- gridExtra::grid.arrange(p1, p2, ncol = 2)

ggsave(
  p,
  filename = "/Users/lukas/OneDrive/Miko/THINC/projects/Shivanna/adm_paper/figure_2_matched_boxplots.pdf",
  height = 14, width = 7)


# Fig 3 - double contrast scatter plot ####
ok <- intersect(rownames(AMKOLPSvsPBS),rownames(WTLPSvsPBS))
XY <- data.frame(
  AMKOLPSvsPBS[ok,],
  WTLPSvsPBS[ok,])

XY$group <- paste(
  XY$diffexpressed, XY$diffexpressed.1, sep = '_')
XY$group <- factor(
  XY$group, levels = c("NO_NO", "DOWN_DOWN", "UP_UP",
                       "DOWN_NO", "NO_DOWN", "NO_UP", "UP_NO"))

genes <- c("Pdlim3", "Mat1a",
           "Neb", "Syt6")
XY$label <- data$GeneSymbol[match(rownames(XY), data$GeneID)]
XY$label[-which(XY$label %in% genes)] <- NA

ggplot(XY[order(XY$group), ], aes(
  logFC.1, logFC, 
  col = group, 
  size = (group %in% c("DOWN_NO", "NO_DOWN", "NO_UP", "UP_NO")),
  label = label)) + 
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_point(shape = 19) +
  labs(
    title = 'Double-contrast scatter-plot',
    x = 'WT',
    y = 'Adm KO',
    size = 'Specific') +
  ggrepel::geom_label_repel() +
  scale_size_manual(values = c(1, 2)) +
  ggpubr::stat_cor(aes(group = NA)) +
  scale_color_manual(
    values = c("lightgrey", "grey50", "grey40",
               "#AF8DC3", "#1B7837",
               "#E9A3C9", "#A1D76A")) +
  theme_classic()

ggsave(
  filename = "/Users/lukas/OneDrive/Miko/THINC/projects/Shivanna/adm_paper/double_contrast.pdf",
  height = 9, width = 10)


# Fig 3 - boxplot of selected genes ####
p1 <- plot_gene('Pdlim3')
p2 <- plot_gene('Mat1a')
p3 <- plot_gene('Neb')
p4 <- plot_gene('Syt6')

p <- grid.arrange(p4, p2, p1, p3, ncol = 2)

ggsave(
  p,
  filename = "/Users/lukas/OneDrive/Miko/THINC/projects/Shivanna/adm_paper/double_contrast_boxplots.pdf",
  height = 10, width = 10)

# Fig 4 - boxplot of selected genes ####
p1 <- plot_gene('Ppm1j')
p2 <- plot_gene('S100a16')
p3 <- plot_gene('Lgi3')
p4 <- plot_gene('Tnfrsf9')
p5 <- plot_gene('Ppp1r16b')
p6 <- plot_gene('Myo5c')

p <- grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3)

ggsave(
  p,
  filename = "/Users/lukas/OneDrive/Miko/THINC/projects/Shivanna/adm_paper/cell_type_marker_boxplots.pdf",
  height = 4, width = 9)

