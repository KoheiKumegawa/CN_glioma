#-------------------------------------------------------------------------------
# 02_PCA_hg38.R
#-------------------------------------------------------------------------------
library(SummarizedExperiment)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggsignif)
library(ggrastr)
library(ggrepel)
"%ni%" <- Negate("%in%")
theme_cb <- function(color="black"){theme_classic() + theme(axis.text = element_text(color=color), axis.ticks = element_line(color=color))}

se <- readRDS("rds/01_se.rds")

#---------- Principal Component Analysis ----------#
row_vars1 <- rowVars(assays(se)$log2tpm)
upper25perc <- which(row_vars1 > quantile(row_vars1, 0.75)) # 6621 genes
pca1 <- prcomp(t(assays(se)$log2tpm[upper25perc,]))

summary(pca1)
# 0.3087  0.2624

df1 <- as.data.frame(pca1$x[,c(1:2)])
df1$treatment <- se$sampleType

treatment_colors <- c("PBS"="blue","CNO"="red")

p1 <- ggplot(df1, aes(x = PC1, y = PC2, color = treatment)) + geom_point(size=4, alpha=0.7) + 
  geom_vline(xintercept = 0, lty = "dotted", size = 0.5) + geom_hline(yintercept = 0, lty = "dotted", size = 0.5) +
  theme_cb() + scale_color_manual(values = treatment_colors) + 
  labs(x = "PC1 (30.9% variance explained)", y = "PC2 (26.2% variance explained)") +
  theme(legend.key.size = unit(0.4, 'cm'))
pdf("output/Plots/02_PCA_hg38.pdf", width = 4, height = 3)
p1
dev.off()

### PC loading ###
pca1_loading <- t(t(pca1$rotation)*pca1$sdev)[,c(1:2)] %>% as.data.frame()
pca1_loading$symbol <- rowData(se[rownames(pca1_loading),])$symbol
write.csv(pca1_loading, "output/Tables/02_PCA_loading_hg38.csv")
