## https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html
## https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
## chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://bioconductor.org/packages/devel/bioc/manuals/fgsea/man/fgsea.pdf
## https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html#gsea-analysis


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
  df_list <- list(df$V1)
  #df_list <- list(df)
  names(df_list) <- unlist(fields)[1]
  data_list <- c(data_list, df_list)
}

close(con)

lung_data_list <- data_list[grep('Lung',names(data_list),value = T)]

## 	Tracheal aspirate transcriptomic and miRNA signatures of extreme premature 
## birth with bronchopulmonary dysplasia
## Homo sapiens
## Series GSE156028
## code availability - https://psilveyra.github.io/silveyralab/

data <- read.csv('C:/Users/spalit/Downloads/GSE156028_NGSCSVBPDwCfirst.csv/GSE156028_NGSCSVBPDwCfirst.csv')
data <- data[-which(rowSums(data[,c(2:ncol(data))]) == 0),] # removing genes with zero counts across all samples

table(data$gene) %>% table # 50539 genes only occur once (not duplicated)
# work with these for now!
data <- data[which(data$gene %in% names(which(table(data$gene) == 1))),]
# Replacing row names
rownames(data) <- data$gene
data <- data[,-1]
data_sum <- data
# REMOVE mitochondrial genes
# data_sum <- data_sum[-grep('^MT-',rownames(data_sum)),]

# Format metadata
metadata <- data.frame(Sample = colnames(data_sum),
                       Group = gsub('[[:digit:]]','',colnames(data_sum)))
metadata$Group[metadata$Group == 'B'] <- 'BPD'
metadata

all(metadata$Sample %in% colnames(data_sum))
all(metadata$Sample == colnames(data_sum))


# Create DEGList object
# Use the group as vector to indicate each sample's group type
DEGL <- DGEList(counts=data_sum, group=metadata$Group)
dim(DEGL)
table(metadata$Group) #21 BPD and 17 Control samples

DEGL <- calcNormFactors(DEGL, method = "TMM")
DEGL$samples$norm.factors

# Cleaning up
## Filtering to remove low counts
keep <- rowSums(cpm(DEGL) > 0.5) >= 10
table(keep)
DEGL <- DEGL[keep, , keep.lib.sizes=FALSE]
dim(DEGL)

lcpm <- cpm(DEGL, log=TRUE)

# Unsupervised PCA of BPD data
plot_pca <- function(input,title){
  input_pca <- prcomp(t(input))
  input_pca_out <- as.data.frame(input_pca$x)
  
  if(all(rownames(input_pca_out) == metadata$Sample) & all(rownames(input_pca_out) %in% metadata$Sample))
  { 
    input_pca_out$Group <- metadata$Group}
  
  percentage <- round(input_pca$sdev / sum(input_pca$sdev) * 100, 2)
  percentage <- paste(colnames(input_pca_out), "(", paste( as.character(percentage), "%", ")", sep=""))
  
  # xy.limits <- range(c(-150,150)) # fix aspect ratio
  p <- ggplot(input_pca_out, aes(x=PC1,y=PC2, color = Group)) + geom_point() + 
    geom_label_repel(aes(label = rownames(input_pca_out),
                         fill = Group),
                     color = 'white', size = 2,
                     box.padding = unit(0.3, "lines"),
                     point.padding = unit(0.25, "lines")) +
    
    xlab(percentage[1]) + 
    ylab(percentage[2])  + 
    theme(legend.position = "top") + 
    labs(title=title) + theme_classic() + coord_fixed(ratio = 1)# + 
  # scale_x_continuous(limits = xy.limits) + scale_y_continuous(limits = xy.limits)
  return(p)
}
plot_pca(lcpm,'PCA of TMM-normalised log CPM counts - BPD')

# Differential expression analysis
group <- as.factor(metadata$Group)
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design

contr.matrix <- makeContrasts(
  CvBPD=BPD-C,
  levels = colnames(design))
contr.matrix

par(mfrow=c(1,2))
v <- voom(DEGL, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

CvBPD=topTable(efit, coef = "CvBPD", adjust = "BH", sort.by = "p",n=Inf)
CvBPD
write.csv(CvBPD, 'C:/Users/spalit/Documents/differential_genes_BPD.csv')
# CvBPD$is_DE <- "NO"
# CvBPD$is_DE[CvBPD$logFC > 0.1 & CvBPD$adj.P.Val < 0.8] <- "UP"
# CvBPD$is_DE[CvBPD$logFC < -0.1 & CvBPD$adj.P.Val < 0.8] <- "DOWN"

CvBPD$is_NK_marker <- 'NO'
CvBPD$is_endo_marker <- 'NO'
CvBPD$is_type2_marker <- 'NO'

CvBPD$is_NK_marker[which(rownames(CvBPD) %in% lung_data_list$`Natural Killer Cell Lung CL:0000623`)] <- 'YES'
CvBPD$is_endo_marker[which(rownames(CvBPD) %in% lung_data_list$`Endothelial Cell Lung CL:0000115`)] <- 'YES'
CvBPD$is_type2_marker[which(rownames(CvBPD) %in% lung_data_list$`Type II Pneumocyte Lung CL:0002063`)] <- 'YES'

CvBPD$genes <- CvBPD %>% rownames
CvBPD$genes[setdiff(seq(1,nrow(CvBPD)),which(CvBPD[order(CvBPD$P.Value),]$is_NK_marker == 'YES') %>% head())] <- NA
# Volcano with NK markers highlighted
ggplot(data=CvBPD, aes(x=logFC, y=-log10(P.Value), col=is_NK_marker, label = genes)) +
  geom_point(data = filter(CvBPD, !is_NK_marker == "YES"), size = 0.5) + #, alpha = .05
  geom_point(data = filter(CvBPD, is_NK_marker == "YES"), size = 1.5) +
  geom_text_repel(max.overlaps = 30) +
  scale_color_manual(values=c("grey", "blue")) + theme_classic() + 
  xlab('log2FoldChange') + ylab('-log10(pvalue)')

foo <- cbind.data.frame(group = metadata$Group, 
                        expr = lcpm[which(rownames(lcpm) == 'SIDT1'),])
foo %>%
  ggplot( aes(x=group, y=expr, fill=group)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.6, alpha=0.9, width = 0.1) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme_classic() + ggtitle('SIDT1')

# Enrichment analysis using Lung markers (NK, type 2, endothelial gene lists) and all log fold changes from GEO BPD study

# gseaDat <- filter(CvBPD, is_NK_marker == 'YES')
# gseaDat <- CvBPD[which(CvBPD$is_NK_marker == 'YES' | CvBPD$is_endo_marker == 'YES' | CvBPD$is_type2_marker == 'YES'),]
# gseaDat$entrezid <- mapIds(org.Hs.eg.db, rownames(gseaDat), 'ENTREZID', 'SYMBOL') %>% as.vector
ranks <- CvBPD$logFC
names(ranks) <- CvBPD %>% rownames
# names(ranks) <- gseaDat$entrezid
head(ranks)
barplot(sort(ranks, decreasing = T))
ranks <- sort(ranks, decreasing=TRUE)

library(fgsea)
res <- fgseaMultilevel(pathways = lung_data_list[c(6,13,16)], 
                       stats = ranks)

plotEnrichment(lung_data_list$`Natural Killer Cell Lung CL:0000623`, ranks)

# ECDF plots
CvBPD_lfc <- CvBPD[,1,drop = F]
binoy$BPD_lfc <- CvBPD$logFC[match(binoy$gene_hsa,rownames(CvBPD_lfc))]

foo <- binoy[,c(12,11,13)]

foo[complete.cases(foo),]$term %>% table
# Endothelial          NK      Type 2 
# 28          25          29

# load library ggplot2
library(ggplot2)

# Basic ECDF plot using ggplot package
# col parameter is used to color plot 
# according to group
ggplot(foo[complete.cases(foo),], aes(x=BPD_lfc, col=term)) + 
  
  # stat_ecdf() function is used to plot ECDF plot
  stat_ecdf() + xlab('BPD - log2FoldChange') + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlim(-2,2) +
  theme_classic()
# ggplot(CvBPD, aes(x=logFC, color=is_NK_marker)) + 
#   
#   # color property for changing color of plot
#   # geom_density() function plots the density plot
#   geom_density()

# temp <- CvBPD[which(CvBPD$is_DE == 'UP'),]
# up_genes <- temp[order(temp$adj.P.Val),] %>% rownames %>% head(30)
# temp <- CvBPD[which(CvBPD$is_DE == 'DOWN'),]
# down_genes <- temp[order(temp$adj.P.Val),] %>% rownames %>% head(30)

# CvBPD$genes[which(!rownames(CvBPD) %in% lung_data_list$`Natural Killer Cell Lung CL:0000623`$V1)] <- NA
# CvBPD$genes[which(!rownames(CvBPD) %in% c(up_genes,down_genes))] <- NA

# ggplot(data=CvBPD, aes(x=logFC, y=-log10(adj.P.Val), col=is_DE, 
#                               label=genes)) +
#   geom_point(size = 0.85) + 
#   geom_text_repel(max.overlaps = 30) +
#   scale_color_manual(values=c("blue", "grey", "red")) +
#   geom_vline(xintercept=c(-0.1, 0.1), col="red", lty = 23) +
#   geom_hline(yintercept=-log10(0.8), col="red", lty = 23) + theme_classic()



###################################################################
# grep('B',colnames(data),value = T) %>% length #21

# # convert between human and mouse gene symbols (because Binoy's data using mice symbols)
# require("biomaRt")
# human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
# mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
# 
# genes_mmu = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = rownames(data_sum) , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
# 
# # meta is just the RNAseq data from Binoy - we are only going to use information
# # on which of the mice genes are protein coding and which are not
# # subset on only protein coding genes
# meta <- read.csv('C:/Users/spalit/Downloads/combined-rnaseq.featureCounts-genes_Binoy.csv')[,c('GeneSymbol','GeneBiotype')] %>% unique
# meta <- meta[which(meta$GeneBiotype == 'protein_coding'),]
# 
# genes_mmu <- genes_mmu[which(genes_mmu$MGI.symbol %in% meta$GeneSymbol),] # keep only protein coding genes
# 
# data_sum <- data_sum[which(rownames(data_sum) %in% genes_mmu$HGNC.symbol),]
# data_sum$gene_mmu <- genes_mmu$MGI.symbol[match(rownames(data_sum),genes_mmu$HGNC.symbol)]
# 
# data_sum <- data_sum %>% group_by(gene_mmu) %>% summarise(across(everything(), list(median))) %>% as.data.frame()
# rownames(data_sum) <- data_sum$gene_mmu
# data_sum <- data_sum[,-1]
# 
# colnames(data_sum) <- gsub('_.*$','',colnames(data_sum))
# # TRY
# # REMOVE mitochondrial genes
# # data_sum <- data_sum[-grep('^mt-',rownames(data_sum)),]
# # rm(mouse,human,meta,genes_mmu)

# samplenames <- colnames(DEGL)
# # Cutoffs
# L <- mean(DEGL$samples$lib.size) * 1e-6
# M <- median(DEGL$samples$lib.size) * 1e-6
# c(L, M)
# 
# lcpm.cutoff <- log2(10/M + 2/L)
# 
# nsamples <- ncol(DEGL)
# col <- brewer.pal(nsamples, "Paired")
# # par(mfrow=c(1,2))
# plot(density(lcpm[,1]), col=col[1], lwd=2, las=2, main="", xlab="")
# title(main="A. Raw data", xlab="Log-cpm")
# abline(v=lcpm.cutoff, lty=3)
# for (i in 2:nsamples){
#   den <- density(lcpm[,i])
#   lines(den$x, den$y, col=col[i], lwd=2)
# }
# legend("topright", samplenames, text.col=col, bty="n")
# 
# lcpm <- cpm(DEGL, log = TRUE)
# plot(density(lcpm[,1]), col=col[1], lwd=2, las=2, main="", xlab="")
# title(main="B. After filtering", xlab="Log-cpm")
# abline(v=lcpm.cutoff, lty=3)
# for (i in 2:nsamples){
#   den <- density(lcpm[,i])
#   lines(den$x, den$y, col=col[i], lwd=2)
# }
# legend("topright", samplenames, text.col=col, bty="n")

# resultz <- summary(decideTests(efit))
# vennDiagram(resultz)

# tfit <- treat(vfit, lfc=0.5)
# dt <- decideTests(tfit)
# summary(dt)

# DElist <- topTreat(tfit, coef=1, n=Inf)
# DElist %>% View
###################################################################
sessionInfo()
# R version 4.2.1 (2022-06-23 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United Kingdom.utf8 
# [2] LC_CTYPE=English_United Kingdom.utf8   
# [3] LC_MONETARY=English_United Kingdom.utf8
# [4] LC_NUMERIC=C                           
# [5] LC_TIME=English_United Kingdom.utf8    
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] enrichR_3.2         biomaRt_2.52.0      RColorBrewer_1.1-3 
# [4] gridExtra_2.3       pheatmap_1.0.12     ggrepel_0.9.3      
# [7] stringr_1.5.0       latex2exp_0.9.6     viridis_0.6.2      
# [10] viridisLite_0.4.1   ggplot2_3.4.2       reshape2_1.4.4     
# [13] sva_3.44.0          BiocParallel_1.30.4 genefilter_1.78.0  
# [16] mgcv_1.8-40         nlme_3.1-157        edgeR_3.38.4       
# [19] limma_3.52.4        dplyr_1.1.1        
# 
# loaded via a namespace (and not attached):
#   [1] bitops_1.0-7           matrixStats_0.63.0     bit64_4.0.5           
# [4] filelock_1.0.2         progress_1.2.2         httr_1.4.6            
# [7] GenomeInfoDb_1.32.4    tools_4.2.1            backports_1.4.1       
# [10] utf8_1.2.3             R6_2.5.1               DBI_1.1.3             
# [13] BiocGenerics_0.42.0    colorspace_2.1-0       withr_2.5.0           
# [16] tidyselect_1.2.0       prettyunits_1.1.1      curl_5.0.0            
# [19] bit_4.0.5              compiler_4.2.1         cli_3.6.1             
# [22] Biobase_2.56.0         xml2_1.3.4             labeling_0.4.2        
# [25] scales_1.2.1           rappdirs_0.3.3         digest_0.6.31         
# [28] XVector_0.36.0         pkgconfig_2.0.3        WriteXLS_6.5.0        
# [31] highr_0.10             dbplyr_2.3.2           fastmap_1.1.1         
# [34] rlang_1.1.0            readxl_1.4.2           rstudioapi_0.14       
# [37] RSQLite_2.3.0          generics_0.1.3         farver_2.1.1          
# [40] car_3.1-2              RCurl_1.98-1.12        magrittr_2.0.3        
# [43] GenomeInfoDbData_1.2.8 Matrix_1.5-4           Rcpp_1.0.10           
# [46] munsell_0.5.0          S4Vectors_0.34.0       fansi_1.0.4           
# [49] abind_1.4-5            lifecycle_1.0.3        stringi_1.7.12        
# [52] carData_3.0-5          zlibbioc_1.42.0        BiocFileCache_2.4.0   
# [55] plyr_1.8.8             grid_4.2.1             blob_1.2.4            
# [58] parallel_4.2.1         crayon_1.5.2           lattice_0.20-45       
# [61] Biostrings_2.64.1      splines_4.2.1          annotate_1.74.0       
# [64] hms_1.1.3              KEGGREST_1.36.3        locfit_1.5-9.7        
# [67] knitr_1.43             pillar_1.9.0           ggpubr_0.6.0          
# [70] rjson_0.2.21           ggsignif_0.6.4         codetools_0.2-18      
# [73] stats4_4.2.1           XML_3.99-0.14          glue_1.6.2            
# [76] evaluate_0.21          png_0.1-8              vctrs_0.6.1           
# [79] cellranger_1.1.0       gtable_0.3.3           purrr_1.0.1           
# [82] tidyr_1.3.0            cachem_1.0.7           xfun_0.39             
# [85] xtable_1.8-4           broom_1.0.4            rstatix_0.7.2         
# [88] survival_3.3-1         tibble_3.2.1           AnnotationDbi_1.58.0  
# [91] memoise_2.0.1          IRanges_2.30.1 