# read in htseq-count result, which counts the reads per annotation element
S1 = read.csv("S1 .csv")
S2 = read.csv("S2 .csv")
S3 = read.csv("S3 .csv")
S4 = read.csv("S4 .csv")
S5 = read.csv("S5 .csv")
S6 = read.csv("S6 .csv")
S7 = read.csv("S7 .csv")

lut <- data.frame(S2$V2,S3$V2,S4$V2,S5$V2,S6$V2,S7$V2,S1$V2)
colnames(lut) <- c("Rc1","Rc2","Rc3","Rvi1","Rvi2","Rvi3","Nar")
rownames(lut) <- S1$V1

# metadata
genotype <- c("var","var","var","var","var","var","narWT")
condition <- c("control","control","control","rnai","rnai","rnai","outgroup")
lut_meta <- data.frame(genotype,condition)
rownames(lut_meta) <- c("Rc1","Rc2","Rc3","Rvi1","Rvi2","Rvi3","Nar")

# sample order 
all(rownames(lut_meta) == colnames(lut))

library(DESeq2)
# Create DESeq object
lut_object <- DESeqDataSetFromMatrix(countData = lut,
                                     colData = lut_meta,
                                     design = ~condition)

# Determine the size factors to use for normalization
lut_object <- estimateSizeFactors(lut_object)
sizeFactors(lut_object)

# Extract the normalized counts
lut_normalized_counts <- counts(lut_object, normalized=TRUE)
head(lut_normalized_counts)

# Vst transformation: 
vsd_lut <- vst(lut_object, blind = TRUE)

# Extract the matrix of transformed counts
vsd_mat_lut <- assay(vsd_lut)

# Compute the correlation values between samples
vsd_cor_lut <- cor(vsd_mat_lut) 
head(vsd_cor_lut)

# Load library for pheatmap,load library for tidyverse
library(pheatmap)
library(tidyverse)

# Plot the heatmap of sample correlation, well at least nar stands out
pheatmap(vsd_cor_lut,main ="sample correlation", annotation = select(lut_meta, condition))

# Plot the PCA of PC1 and PC2
plotPCA(vsd_lut, intgroup="condition")

# Run the DESeq2 analysis
lut_object <- DESeqDataSetFromMatrix(countData = lut,
                                     colData = lut_meta,
                                     design = ~condition)
deseq_lut <- DESeq(lut_object)

# Creating data frame with mean and variance for every gene
lut_matrix<- data.matrix(lut)
mean_count <- rowMeans(lut_matrix)
var_count <-rowVars(lut_matrix)
df <- data.frame(mean_count, var_count)
# load ggplot2
library(ggplot2)
  ggplot(df) +
  ggtitle("mean-variance") +
  geom_point(aes(x=mean_count, y=var_count)) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("Mean counts per gene") +
  ylab("Variance per gene")

# Dispersion 
# decrease with increasing mean and cluster around the maximum likelihood 
plotDispEsts(deseq_lut, main = "dispersion")

# Get the DESeq2 contrasts:
lut_result <- results(deseq_lut, alpha = 0.05)
lut_result
# MA plot 
plotMA(lut_result, ylim=c(-8,8),main = "MA plot")
# LFC shrinkage: not a lot of differentially expressed genes: 
lut_result <- lfcShrink(deseq_lut, 
                        contrast=c("condition", "rnai", "control"),
                        res=lut_result)
plotMA(lut_result, ylim=c(-6,6),main = "MA after shrinkage")

# DESeq2 results:
mcols(lut_result)
head(lut_result, n = 10)
# multiple test correction: (BH method), p-adjusted
summary(lut_result)
lut_all <- data.frame(S1$V1,lut_result)
View(lut_all)
# arrange by more significant DE
lut_sig <- subset(lut_all, padj <.05)
lut_sig <- lut_sig %>%  arrange(padj)
lut_sig_res <- lut_sig[,-1]
rownames(lut_sig_res) <- lut_sig[,1]


# heat map
sig_norm_counts_lut <- lut_normalized_counts[rownames(lut_sig_res), ]
pheatmap(sig_norm_counts_lut, 
         cluster_rows = F,
         cluster_cols = T, 
         show_rownames = F,
         annotation = select(lut_meta, condition), 
         scale = "row", main = "sample expression profile")

# Obtain logical vector regarding whether padj values are less than 0.05  
lut_logic_all <- data.frame(lut_all) %>% mutate(threshold = padj < 0.05)




# Volcano plot
# Create the volcano plot
v_basic <- ggplot(lut_logic_all,aes(x = log2FoldChange, y = -log10(padj), color = threshold)) + 
  geom_point(alpha = .5) + 
  ggtitle("vocano plot (threshold padj = .05)") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") + 
  theme(legend.position = "none", 
        plot.title = element_text(size = rel(1.5), hjust = 0.5), 
        axis.title = element_text(size = rel(1.25))) 
v_basic
v_line <- v_basic + geom_hline(yintercept=-log10(.05), color = "gray50") 
v_line 


# we have a list of myb bHLH that are differentially expressed, and we want to visualize them: 
library(ggrepel)
myb_bHLH <- read.csv("myb_bHLH.csv")
rownames(myb_bHLH) <- myb_bHLH[,1]
v_label <- ggplot(data = myb_bHLH, aes(x = -log2FoldChange, y = -log10(padj),color = annotation)) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") + 
  ggtitle("Mlu_12200 among DE mybs and bHLHs") +
  geom_point(size = 2) +
  scale_color_manual(values=c("#d8b365",  "#5ab4ac","#FC4E07")) +
  geom_text_repel(label = rownames(myb_bHLH),size=4) +
  geom_hline(yintercept=-log10(.05), color = "gray50") +
  geom_vline(xintercept= 0, color = "gray50") +  theme_linedraw()
v_label

myb_bHLH$annotation[7] = "x_myb"
v_label2 <- ggplot(data = myb_bHLH, aes(x = -log2FoldChange, y = -log10(padj),color = annotation)) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") + 
  ggtitle("DE mybs and bHLHs") +
  geom_point(size = 2) +
  scale_color_manual(values=c("#d8b365",  "#5ab4ac")) +
  geom_text_repel(label = rownames(myb_bHLH),size=4) +
  geom_hline(yintercept=-log10(.05), color = "gray50") +
  geom_vline(xintercept= 0, color = "gray50")  +  theme_linedraw()
v_label2


v_mybs <- ggplot(lut_logic_all,aes(x = log2FoldChange, y = -log10(padj), color = threshold) ) + 
  geom_point() +
  ggtitle("differetially expressed mybs and bHLHs") +
  geom_point(data = myb_bHLH,aes(x = -log2FoldChange, y = -log10(padj), color = annotation), shape =17, size =3) +
  scale_color_manual(values=c("#999999", "#E69F00","red","blue")) +
  geom_hline(yintercept=-log10(.05), color = "gray50") + theme_minimal()
  
v_mybs

# not log transformed read count plot:
plotCounts(deseq_lut, gene = "Mlu_12200", intgroup="condition",transform=F)
# log-transformed read count plot:
plotCounts(deseq_lut, gene = "Mlu_12200", intgroup="condition")


