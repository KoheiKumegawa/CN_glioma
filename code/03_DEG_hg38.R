#-------------------------------------------------------------------------------
# 03_DEG_hg38.R
#-------------------------------------------------------------------------------
source("code/edgeR_PairwiseFunction.R")
library(ggplot2)
library(ggrastr)
library(ggrepel)
theme_cb <- function(color="black"){theme_classic() + theme(axis.text = element_text(color=color), axis.ticks = element_line(color=color))}

se <- readRDS("rds/01_se.rds")

DiffTest <- edgeR_pairwise(se, compareCol = "sampleType", topGroup = "CNO", bottomGroup = "PBS")

#----- Volcano plot -----#
df1 <- assay(DiffTest) %>% as.data.frame()
df1$logFDR <- -log10(df1$FDR)
df1$Signif <- "NotSignif"
df1$Signif[which(df1$log2FC > 1 & df1$FDR < 0.01)] <- "CNO_UP"
df1$Signif[which(df1$log2FC < -1 & df1$FDR < 0.01)] <- "CNO_DN"
df1$Signif <- factor(df1$Signif, levels = c("NotSignif","CNO_UP","CNO_DN"))
df1 <- df1[order(df1$Signif),]
df1$label <- ""
df1$label[which(rownames(df1) %in% c("KLF5","CYR61","CTGF"))] <- rownames(df1)[which(rownames(df1) %in% c("KLF5","CYR61","CTGF"))]

p1 <- ggplot(df1, aes(x=log2FC,y=logFDR,color=Signif,label=label)) + geom_point_rast(size=2.5) + theme_cb() +
  geom_vline(xintercept = c(-1,1), lty = "dotted", size = 0.5) + 
  geom_hline(yintercept = 2, lty = "dotted", size = 0.5) +
  geom_vline(xintercept = 0, size = 0.5) +
  labs(x="log2FC",y="-log10(FDR)") +
  geom_label_repel(nudge_x = 1, nudge_y = 0.5) +
  scale_color_manual(values = c("CNO_UP"="red","CNO_DN"="blue","NotSignif"="darkgray"))
pdf("output/Plots/03_Volcano_hg38.pdf", width = 7, height = 6)
p1
dev.off()

write.csv(df1, "output/Tables/03_SourceData_Volcano_hg38.csv")
