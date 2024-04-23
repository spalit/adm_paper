## Bulk RNA-Seq from whole lung tissue of mice
## treated with LPS [Lipopolysaccharide] or control PBS [Phosphate-buffered saline] 
## and are wild type or Adm (Adrenomedullin) haplo insufficient (only 1 copy left)
## There are 4 conditions (3 replicates each)
## Standard bulk analysis, differential expression, PCA, gene set enrichment

# RESOURCES
## Adding interaction term https://support.bioconductor.org/p/46153/
## https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
## https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html
## https://ycl6.github.io/GO-Enrichment-Analysis-Demo/4_enrichR.html#Session_information
## https://f1000research.com/articles/5-1438
## Stat17Stat3 downstream targets - https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/STAT1_01.html

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

##########----------------------##########
# Make lung cell marker lists from Tabula Muris
file_path <- "C:/Users/spalit/Downloads/Tabula_Muris.txt"
data_list <- list()

con <- file(file_path, "r")

while (length(line <- readLines(con, n = 1)) > 0) {
  fields <- strsplit(line, "\t")
  df <- as.data.frame(matrix(unlist(fields)[-c(1,2)]), nrow = 1, byrow = TRUE)
  # if you want list of lists, run the next line
  # df_list <- list(df$V1)
  df_list <- list(df)
  names(df_list) <- unlist(fields)[1]
  data_list <- c(data_list, df_list)
}

close(con)

lung_data_list <- data_list[grep('Lung',names(data_list),value = T)]

##########----------------------##########
# Load bulk RNA seq data (Binoy)
metadata <- readxl::read_xlsx('annotation.xlsx', sheet = 1)[,c(1,2)]
metadata$Group %>% unique #4

data <- read.csv('combined-rnaseq.featureCounts-genes_Binoy.csv')
data$GeneBiotype %>% unique

# Perform analysis with only protein-coding genes
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
  
  xy.limits <- range(c(-80,80)) # fix aspect ratio
  p <- ggplot(input_pca_out, aes(x=PC1,y=PC2, color = Group)) + geom_point(size = 0.8) + 
    # geom_label_repel(aes(label = rownames(input_pca_out),
    #                      fill = Group),
    #                  color = 'white', size = 2,
    #                  box.padding = unit(0.3, "lines"),
    #                  point.padding = unit(0.25, "lines")) +
    
    xlab(percentage[1]) + 
    ylab(percentage[2])  + 
    theme(legend.position = "top") + 
    labs(title=title) + theme_classic() + coord_fixed(ratio = 1) + 
    scale_x_continuous(limits = xy.limits) + scale_y_continuous(limits = xy.limits)
  return(p)
}
# FIGURE 1A
plot_pca(lcpm,'PCA of TMM-normalised log CPM counts')


# Model with each coefficient corresponding to a group mean - NO INTERACTION
group <- as.factor(metadata$Group)
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design

contr.matrix <- makeContrasts(
  AMKOLPSvsPBS = AMKO_LPS - AMKO_PBS,
  WTLPSvsPBS = WT_LPS - WT_PBS,
  AMKOvsWTLPS = AMKO_LPS - WT_LPS,
  # Diff = (AMKO_LPS - AMKO_PBS)-(WT_LPS - WT_PBS),
  levels = colnames(design))
contr.matrix

par(mfrow=c(1,1))
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

AMKOvsWTLPS <- topTreat(tfit, coef=3, n=Inf)
head(AMKOvsWTLPS)

AMKOLPSvsPBS$gene <- data$GeneSymbol[match(rownames(AMKOLPSvsPBS),data$GeneID)]
WTLPSvsPBS$gene <- data$GeneSymbol[match(rownames(WTLPSvsPBS),data$GeneID)]

write.csv(AMKOLPSvsPBS, 'C:/Users/spalit/Documents/differential_genes_LPSvsPBS_AMKO.csv')
write.csv(WTLPSvsPBS, 'C:/Users/spalit/Documents/differential_genes_LPSvsPBS_WT.csv')
# diff <- topTreat(tfit, coef=3, n=Inf)
# head(diff)

# FIGURE 1B
Top50 = topTable(efit,number=50,adjust="BH")
str(Top50)

## visualize
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
         breaks = seq(-1.5, 1.5, length = 100))

# FIGURE 1C & D
## Volcano plot of DE genes in AMKO-LPS vs PBS

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
  # theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept=c(-0.5, 0.5), col="red", lty = 23) +
  geom_hline(yintercept=-log10(0.05), col="red", lty = 23) + theme_classic() + 
  coord_equal(ratio = 4)


## Volcano plot of DE genes in WT-LPS vs PBS
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
  # theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept=c(-0.5, 0.5), col="red", lty = 23) +
  geom_hline(yintercept=-log10(0.05), col="red", lty = 23) + theme_classic() + 
  coord_equal(ratio = 4)

################################################################################

# FIGURE 2
# 01-03-2024
# Matching gene/protein expression results from previous studies

# qPCR/protein data in box-plots
df <- readxl::read_excel('Gene and protein expression results from previous publication.xlsx',col_types = "numeric", sheet = 1)
df <- subset(df, rowSums(is.na(df)) != ncol(df))

foo <- c(1:5)
my_list <- list()
# plot parameters
plot_list = list()
i = 1

for(each in c('CCL2','CCL3','Cxcl1','ICAM1','IL-1B','TNFa','stat3','stat1')){
  my_list[[each]] <- df[foo,] %>% as.data.frame
  colnames(my_list[[each]]) <- c("Adm+/+mice_PBS","Adm+/-mice_PBS","Adm+/+mice_L6","Adm+/-mice_L6")
  foo <- foo + 5
  
  tmp_df <- my_list[[each]]
  tmp_df <- tmp_df[,c("Adm+/-mice_L6","Adm+/-mice_PBS","Adm+/+mice_L6","Adm+/+mice_PBS")]
  
  p <- stack(tmp_df) %>%
    ggplot( aes(x=ind, y=values, fill=ind)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    ggtitle(each) +  theme_classic() + 
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1))
  # print(p)
  plot_list[[i]] = p
  i = i+1
  # ggsave(filename = paste0('Figure 2_qPCR_protein_val_',each,'_boxplot.pdf'),
  #        plot = p, width = 8, height = 5)
}

genes_from_publ <- c('Ccl2','Ccl3','Cxcl1','Icam1','Il1b','Tnf','Stat3','Stat1')
geneid_from_publ <- data$GeneID[match(genes_from_publ,data$GeneSymbol)] #data$GeneID[which(data$GeneSymbol %in% genes_from_publ)]

plot_list_RNA = list()
i = 1

for (each in geneid_from_publ) {
  foo <- cbind.data.frame(group = metadata$Group, 
                          expr = lcpm[which(rownames(lcpm) == each),])
  p <- foo %>%
    ggplot( aes(x=group, y=expr, fill=group)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ggtitle(paste(data$GeneSymbol[which(data$GeneID == each)])) + 
    theme_classic()
  # print(p)
  
  plot_list_RNA[[i]] = p
  i = i+1
  # ggsave(filename = paste0('Figure 2_',data$GeneSymbol[which(data$GeneID == each)],'_boxplot.pdf'),
  #        plot = p, width = 10, height = 7)
}

# makes 1 plot for each pair - columns will be genes and rows will be the qPCR or protein data matched w the RNA-seq
foo <- c(1:2)
for (i in 1:4) {
   ggsave(filename = paste0('Figure 2_matched_boxplots',i,'.pdf'),
         plot = marrangeGrob(append(plot_list[foo],plot_list_RNA[foo]), nrow = 2, ncol = 2),
         width = 14, height = 7)
  foo <- foo+2
}
# to make separate plots for each gene
my_vec <- c('CCL2','CCL3','Cxcl1','ICAM1','IL-1B','TNFa','stat3','stat1')
for (i in 1:length(plot_list)) {
  ggsave(filename = paste0('Figure 2_',my_vec[i],'_boxplot.pdf'),
         plot = marrangeGrob(list(plot_list[[i]],plot_list_RNA[[i]]), nrow = 1, ncol = 2),
         width = 14, height = 7)
}

dev.off()
###############################################################################
# FIGURE 3
# Double-contrast scatter-plot:
X <- AMKOLPSvsPBS[,c(1,5,6), drop = F]
Y <- WTLPSvsPBS[,c(1,5,6), drop = F]

Y <- Y[rownames(X),,drop = F]

XY <- data.frame(AMKO_lfc = X$logFC, AMKO_adjpval = X$adj.P.Val, 
                 WT_lfc = Y$logFC, WT_adjpval = Y$adj.P.Val)
rownames(XY) <- rownames(X)


XY_sig <- XY[which(abs(XY$AMKO_lfc) > 0.5 & XY$AMKO_adjpval < 0.05),] # get all genes significant in AMKO
XY_sig$group <- rep('not significant in WT',nrow(XY_sig))
XY_sig$group[which(XY_sig$WT_adjpval < 0.05)] <- 'significant in WT'


# ALL GENES SPECIFICALLY DIFFERENTIALLY EXPRESSED IN AMKO (not in WT)
list_of_s_AMKO_genes <- XY_sig[which(XY_sig$group == 'not significant in WT'),] # 752 genes
XY$group <- 'others'
XY$group[which(rownames(XY) %in% rownames(list_of_s_AMKO_genes))] <- 'significant'

ggplot(XY, aes(x=WT_lfc, y=AMKO_lfc, col=group)) + 
  geom_point(size = 1) + geom_abline(intercept = 0, slope = 1, col="red", size = 0.25, lwd = 0.5, lty = 4) + 
  ggtitle('Genes DE in AMKO but not in WT') + theme_classic() +
  scale_color_manual(values=c('blue','red')) + 
  xlim(-7,12) + ylim(-7,12) #, aes(shape=group)

list_of_s_AMKO_genes$GeneSymbol <- data$GeneSymbol[match(rownames(list_of_s_AMKO_genes),data$GeneID)]
write.csv(list_of_s_AMKO_genes, 'C:/Users/spalit/Documents/130324_genes_sign_in_AMKO_not_in_WT.csv')


# ALL GENES NOT DIFFERENTIALLY EXPRESSED IN AMKO
# Adm KO probably fails to recruit some specific celltypes/pathways - enrichment analysis with the genes below

XY_sig <- XY[which(abs(XY$WT_lfc) > 0.5 & XY$WT_adjpval < 0.05),] # genes significantly enriched in WT
XY_sig$group <- rep('not significant in KO',nrow(XY_sig))
XY_sig$group[which(XY_sig$AMKO_adjpval < 0.05)] <- 'significant in KO'

list_of_ns_AMKO_genes <- XY_sig[which(XY_sig$group == 'not significant in KO'),] # 155 genes
XY$group <- 'others'
XY$group[which(rownames(XY) %in% rownames(list_of_ns_AMKO_genes))] <- 'not significant'

ggplot(XY, aes(x=WT_lfc, y=AMKO_lfc, col=group)) + 
  geom_point(size = 1) + geom_abline(intercept = 0, slope = 1, col="red", size = 0.25, lwd = 0.5, lty = 4) + 
  ggtitle('Genes DE in WT but not in AMKO') + theme_classic() + scale_color_manual(values=c('blue','red')) + 
  xlim(-7,12) + ylim(-7,12) #, aes(shape=group)

list_of_ns_AMKO_genes$GeneSymbol <- data$GeneSymbol[match(rownames(list_of_ns_AMKO_genes),data$GeneID)]
write.csv(list_of_ns_AMKO_genes, 'C:/Users/spalit/Documents/130324_genes_sign_in_WT_not_in_AMKO.csv')

# BOXPLOTS OF SELECTED GENES

gene_ <- 'Pdlim3'#'Ncf1'
geneid_from_publ <- data$GeneID[match(gene_,data$GeneSymbol)]

foo <- cbind.data.frame(group = metadata$Group, 
                          expr = lcpm[which(rownames(lcpm) == geneid_from_publ),])
  p <- foo %>%
    ggplot( aes(x=group, y=expr, fill=group)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ggtitle(gene_) + 
    theme_classic()
  print(p)
  
  # Gene DE in AMKO not in WT
  # Mat1a
  foo <- cbind.data.frame(group = metadata$Group, 
                          expr = lcpm[which(rownames(lcpm) == 'ENSMUSG00000037798'),])
  p <- foo %>%
    ggplot( aes(x=group, y=expr, fill=group)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ggtitle(paste(data$GeneSymbol[which(data$GeneID == 'ENSMUSG00000037798')])) + 
    theme_classic()
  print(p)
  
  # Gene DE in WT not in AMKO
  # Neb
  foo <- cbind.data.frame(group = metadata$Group, 
                          expr = lcpm[which(rownames(lcpm) == 'ENSMUSG00000026950'),])
  p <- foo %>%
    ggplot( aes(x=group, y=expr, fill=group)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ggtitle(paste(data$GeneSymbol[which(data$GeneID == 'ENSMUSG00000026950')])) + 
    theme_classic()
  print(p)
  
  ######################################################################################
  # Enrichment analysis
  # install.packages("enrichR")
  library("enrichR")
  dbs <- listEnrichrDbs()
  up.genes.symbol <- list_of_s_AMKO_genes$GeneSymbol[which(list_of_s_AMKO_genes$AMKO_lfc > 0)]
  upEnriched_go <- enrichr(genes = up.genes.symbol, databases = 'Tabula_Muris')
  upEnriched_go$Tabula_Muris %>% View
  tmp <- upEnriched_go[[1]][grep('Lung',upEnriched_go[[1]]$Term),]
  plotEnrich(tmp, showTerms = 8, numChar = 50, 
             y = "Ratio", orderBy = "P.value") + theme_classic() + ggtitle('Positively DE in AMKO (not in WT)')
  
  dn.genes.symbol <- list_of_s_AMKO_genes$GeneSymbol[which(list_of_s_AMKO_genes$AMKO_lfc < 0)]
  dnEnriched_go <- enrichr(genes = dn.genes.symbol, databases = 'Tabula_Muris')
  dnEnriched_go$Tabula_Muris %>% View
  tmp <- dnEnriched_go[[1]][grep('Lung',dnEnriched_go[[1]]$Term),]
  plotEnrich(tmp, showTerms = 8, numChar = 50, 
             y = "Ratio", orderBy = "P.value") + theme_classic() + ggtitle('Negatively DE in AMKO (not in WT)')
  
  up.genes.symbol <- list_of_ns_AMKO_genes$GeneSymbol[which(list_of_ns_AMKO_genes$WT_lfc > 0)]
  upEnriched_go <- enrichr(genes = up.genes.symbol, databases = 'Tabula_Muris')
  upEnriched_go$Tabula_Muris %>% View
  tmp <- upEnriched_go[[1]][grep('Lung',upEnriched_go[[1]]$Term),]
  plotEnrich(tmp, showTerms = 8, numChar = 40, 
             y = "Ratio", orderBy = "P.value") + theme_classic() + ggtitle('Positively DE in WT (not in AMKO)')
  
  ##############################################################################
  #Figure 4 - lung atlas plots
  # MODIFIED
  # USE BOTH UP AS WELL AS DOWN REGULATED AMKO LUNG MARKER GENES
  # NK genes
  intersect(list_of_s_AMKO_genes$GeneSymbol,str_to_title(lung_data_list[['Natural Killer Cell Lung CL:0000623']]$V1))
  # endothelial genes
  intersect(list_of_s_AMKO_genes$GeneSymbol,str_to_title(lung_data_list$`Endothelial Cell Lung CL:0000115`$V1))
  list_of_s_AMKO_genes[which(list_of_s_AMKO_genes$GeneSymbol %in% str_to_title(lung_data_list$`Endothelial Cell Lung CL:0000115`$V1)),] %>% View
  # type II pneumocyte genes
  intersect(list_of_s_AMKO_genes$GeneSymbol,str_to_title(lung_data_list$`Type II Pneumocyte Lung CL:0002063`$V1))
  
  # Figure 4 - boxplots
  
  gene_ <- c('Lgi3','Myo5c','S100a16','Ppp1r16b','Ppm1j','Tnfrsf9')


  for (each in gene_) {
    geneid_from_publ <- data$GeneID[match(each,data$GeneSymbol)]
    foo <- cbind.data.frame(group = metadata$Group, 
                            expr = lcpm[which(rownames(lcpm) == geneid_from_publ),])
    
    p <- foo %>%
      ggplot( aes(x=group, y=expr, fill=group)) +
      geom_boxplot() +
      scale_fill_viridis(discrete = TRUE, alpha=0.6) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
      ggtitle(each) + 
      theme_classic()
    ggsave(
      filename = paste0("Figure 4_boxplots ",each,".pdf"), 
      plot = p, 
      width = 15, height = 9
    )
  }
  
#  FIGURE 4 Heatmaps
  ### Refit model as a two-factor model (rather than using the group variable)
  # INTERACTION
  vars <- colsplit(metadata$Group, "_", c("Condition", "Treatment"))
  metadata <- cbind.data.frame(metadata,vars)
  
  
  mm <- model.matrix(~metadata$Condition*metadata$Treatment)
  #mm <- model.matrix(~metadata$Treatment*metadata$Condition)
  colnames(mm)
  
  y <- voom(DEGL, mm, plot = TRUE)
  fit <- lmFit(y, mm)
  head(coef(fit))
  
  tmp <- contrasts.fit(fit, coef = 4)
  tmp <- eBayes(tmp)
  plotSA(tmp, main="Final model: Mean-variance trend")
  
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  head(top.table, 20)
  
  length(which(top.table$adj.P.Val < 0.5))
  # top.table$GeneSymbol <- data$GeneSymbol[match(rownames(top.table),data$GeneID)]
  # write.csv(top.table, 'C:/Users/spalit/Documents/130324_all_genes_interaction.csv')

  colnames(top.table) <- paste(colnames(top.table),'interact',sep = '_')
  colnames(WTLPSvsPBS) <- paste(colnames(WTLPSvsPBS),'WT_LPS_vs_PBS',sep = '_')
  colnames(AMKOLPSvsPBS) <- paste(colnames(AMKOLPSvsPBS),'AMKO_LPS_vs_PBS',sep = '_')
  colnames(AMKOvsWTLPS) <- paste(colnames(AMKOvsWTLPS),'AMKO_vs_WT_LPS',sep = '_')
  
  AMKOvsWTLPS <- AMKOvsWTLPS[rownames(WTLPSvsPBS),]
  AMKOLPSvsPBS <- AMKOLPSvsPBS[rownames(WTLPSvsPBS),]
  top.table <- top.table[rownames(WTLPSvsPBS),]
  
  final_df <- cbind.data.frame(top.table,AMKOLPSvsPBS,WTLPSvsPBS,AMKOvsWTLPS)
  final_df$GeneSymbol <- data$GeneSymbol[match(rownames(final_df),data$GeneID)]
  final_df <- final_df[,c("GeneSymbol",colnames(final_df)[1:length(colnames(final_df))-1])]# %>% View
  # write.csv(final_df,'all_genes_all_contrast.csv')
  final_df <- final_df[,-c(4,5,7,9,10,11,14,15,16,19,20,21)]
  
  ## FOR GENES TO BE SELECTED AS CANDIDATES 3 RULES MUST BE APPLIED:
  ## i. Expression above the 15th percentile of the data
  
  # compare the quantitative value at the 15th percentile to the mean expression
  # genes with expression above the 15 percentile
  rMean_lcpm <- lcpm %>% rowMeans() # calculate mean gene expression across all samples
  # WHICH GENES HAVE MEAN EXPRESSION GREATER THAN THE 15th PERCENTILE
  expressed_genes <- which(rMean_lcpm > quantile(rMean_lcpm, probs = 0.15)) %>% names
  final_df <- final_df[expressed_genes,]
  
  ## ii. is a lung marker (ref. Tabula Muris lung cell atlas)
  NK_markers_all <- final_df[which(final_df$GeneSymbol %in% str_to_title(lung_data_list[['Natural Killer Cell Lung CL:0000623']]$V1)),]
  
  ## iii. has significant interaction p-value
  ## finding genes with differential effects to LPS treatment in Adm KO or WT)
  Xx <- NK_markers_all[which(NK_markers_all$logFC_AMKO_LPS_vs_PBS > 0 & NK_markers_all$logFC_AMKO_vs_WT_LPS > 0),] #15 genes in total
  NK_markers_UP_AMKO_head10 <- Xx[order(Xx$adj.P.Val_interact),] %>% head(10)
  
  endo_markers_all <- final_df[which(final_df$GeneSymbol %in% str_to_title(lung_data_list[['Endothelial Cell Lung CL:0000115']]$V1)),]
  Xx.1 <- endo_markers_all[which(endo_markers_all$logFC_AMKO_LPS_vs_PBS < 0 & endo_markers_all$logFC_AMKO_vs_WT_LPS < 0),] #37
  endo_markers_DOWN_AMKO_head10 <- Xx.1[order(Xx.1$adj.P.Val_interact),] %>% head(10)
  
  
  type2_markers_all <- final_df[which(final_df$GeneSymbol %in% str_to_title(lung_data_list[['Type II Pneumocyte Lung CL:0002063']]$V1)),]
  Xx.2 <- type2_markers_all[which(type2_markers_all$logFC_AMKO_LPS_vs_PBS < 0 & type2_markers_all$logFC_AMKO_vs_WT_LPS < 0),] #27
  type2_markers_DOWN_AMKO_head10 <- Xx.2[order(Xx.2$adj.P.Val_interact),] %>% head(10)
  
  lcpm01 <- lcpm[which(rownames(lcpm) %in% c(rownames(NK_markers_UP_AMKO_head10),
                                             rownames(endo_markers_DOWN_AMKO_head10),
                                             rownames(type2_markers_DOWN_AMKO_head10))
  ),] %>% as.data.frame
  lcpm01 <- lcpm01[c(rownames(NK_markers_UP_AMKO_head10),rownames(endo_markers_DOWN_AMKO_head10),rownames(type2_markers_DOWN_AMKO_head10)),]
  
  lcpm01$GeneSymbol <- data$GeneSymbol[match(rownames(lcpm01),data$GeneID)]
  rownames(lcpm01) <- lcpm01$GeneSymbol
  lcpm01 <- as.matrix(lcpm01[,-dim(lcpm01)[2]])
  
  # WT LPS normalized to WT PBS
  mean_lcpm_WT_PBS <- lcpm01[,metadata$Sample[which(metadata$Treatment == 'PBS' & metadata$Condition == 'WT')]] %>% rowMeans()
  mean_lcpm_WT_LPS <- lcpm01[,metadata$Sample[which(metadata$Treatment == 'LPS' & metadata$Condition == 'WT')]] - mean_lcpm_WT_PBS
  
  
  # AMKO LPS normalized to AMKO PBS
  mean_lcpm_AMKO_PBS <- lcpm01[,metadata$Sample[which(metadata$Treatment == 'PBS' & metadata$Condition == 'AMKO')]] %>% rowMeans()
  mean_lcpm_AMKO_LPS <- lcpm01[,metadata$Sample[which(metadata$Treatment == 'LPS' & metadata$Condition == 'AMKO')]] - mean_lcpm_AMKO_PBS
  
  if(all(rownames(mean_lcpm_AMKO_LPS) == rownames(mean_lcpm_WT_LPS)) & all(rownames(mean_lcpm_AMKO_LPS) %in% rownames(mean_lcpm_WT_LPS)))
    mean_lcpm_LPS <- cbind.data.frame(mean_lcpm_AMKO_LPS,mean_lcpm_WT_LPS)
  
  metadata01 <- metadata
  rownames(metadata01) <- metadata01$Sample
  metadata01 <- metadata01[,-c(1,2)]
  # metadata01 <- metadata01[,-c(3,4)]
  
  X <- data.frame(gene = rownames(lcpm01),celltype = 'Type II pneumocytes')
  X$celltype[which(X$gene %in% NK_markers_UP_AMKO_head10$GeneSymbol)] <- 'NK cells'
  X$celltype[which(X$gene %in% endo_markers_DOWN_AMKO_head10$GeneSymbol)] <- 'Endothelial cells'
  
  rownames(X) <- X$gene
  X <- X[,2,drop = FALSE]
  
  cal_z_score <- function(x){
    (x - mean(x)) / sd(x)
  }
  
  mean_lcpm_LPS_norm <- t(apply(mean_lcpm_LPS, 1, cal_z_score))
  # https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
  pheatmap(mean_lcpm_LPS_norm, annotation_col = metadata01, annotation_row = X, cutree_cols = 2, cluster_rows = FALSE,# cutree_rows = 3,
           border_color = NA, cellheight = 10, cellwidth = 10, gaps_row = c(10,20))
  
#  Figure 4
  lcpm01_norm <- t(apply(lcpm01, 1, cal_z_score))
  pheatmap(lcpm01_norm[,metadata$Sample[order(metadata$Group)]], annotation_col = metadata01, cutree_cols = 2, annotation_row = X, cluster_rows = F, cluster_cols = F,
           border_color = NA, cellheight = 10, cellwidth = 10, gaps_row = c(10,20), gaps_col = 6)
  
  
################################################################################
