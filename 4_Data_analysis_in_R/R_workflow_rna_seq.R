# ---------------- Libraries -----------------------------------
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(ggplot2)
library(forcats)
library(enrichplot)
library(ggstance)
library(ggupset)

# ---------------- Import --------------------------------------
# Import data from featureCounts
countdata <- read.table("feature_counts.txt", header=TRUE, row.names=1)

# ---------------- 5. Exploratory data analysis ----------------
# Check if the samples from the same experimental group show similar gene
# expression patterns Remove first five columns, by removing (chr, start, end,
# strand, length)
countdata <- countdata[ ,6:ncol(countdata)]

# Remove sorted.bam from filenames
colnames(countdata) <- gsub("\\.sorted.bam$", "", colnames(countdata))

# Convert to matrix
countdata <- as.matrix(countdata)

# specify the experimental group of each sample. 
group <- factor(c(rep("HER2+", 3), rep("NonTNBC", 3), rep("TNBC", 3), rep("Normal", 3)))
condition <- factor(c(rep("Sample", 9), rep("Normal", 3)))
sample <- factor(c("1", "2", "3", "1", "2", "3", "1", "2", "3", "1", "2", "3"))

# Create the DESeqDataSet object: Read in the counts from FeatureCounts 
(coldata <- data.frame(row.names=colnames(countdata), group, sample, condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~group+sample)

# Run the DESeq2::DESeq() function and save its output to a variable.
dds <- DESeq(dds)

# log transformation. Remove the dependence of the variance on the mean
rld <- rlog(dds, blind=TRUE)

# histogram of the log-transformation
hist(assay(rld))

# Assess how the samples cluster based on their gene expression profiles
plotPCA(rld, intgroup="group")


# ---------------- 6. Differential expression analysis ----------------
## Extract the results for at least one pairwise contrast
res1 <- results(dds, contrast=c("group", "TNBC", "NonTNBC"), alpha = 0.05)
res2 <- results(dds, contrast=c("group", "HER", "TNBC"), alpha = 0.05)

## Order by adjusted p-value
res1 <- res1[order(res1$padj), ]
res2 <- res2[order(res2$padj), ]

# print out a summary of the results 
summary(res1)
summary(res2)

## Merge with normalized count data
resdata1 <- merge(as.data.frame(res1), as.data.frame(counts(dds, normalized=TRUE)), by="row.names")
resdata2 <- merge(as.data.frame(res2), as.data.frame(counts(dds, normalized=TRUE)), by="row.names")
names(resdata1)[1] <- "Gene"
names(resdata2)[1] <- "Gene"

## write the csv output of all differential expressed genes 
write.csv(resdata1, file="diffexpr-results_1.csv")
write.csv(resdata2, file="diffexpr-results_2.csv")

## Plot of p-values
hist(res1$pvalue, breaks=50, col="grey")
hist(res2$pvalue, breaks=50, col="grey")

## Expression comparison of chosen genes 
par(mfrow=c(3,2))
HRAS <- plotCounts(dds,"ENSG00000174775", returnData = TRUE)
plot(count ~ group , data=HRAS, main = "Expression of HRAS")
plot(count ~ condition , data=HRAS, main = "Expression of HRAS")
as.data.frame(HRAS)

SPARC <- plotCounts(dds, "ENSG00000113140", returnData = TRUE)
boxplot(count ~ group , data=SPARC, main = "Expression of SPARC")
boxplot(count ~ condition , data=SPARC, main = "Expression of SPARC")
as.data.frame(SPARC)

OAZ1 <- plotCounts(dds,"ENSG00000104904", returnData = TRUE)
boxplot(count ~ group , data=OAZ1, main = "Expression of OAZ1")
boxplot(count ~ condition , data=OAZ1, main = "Expression of OAZ1")
as.data.frame(OAZ1)

## volcano plot of the log2fold change of pairwise comparison 
par(mfrow=c(2,1))
plot( res1$log2FoldChange, -log(res1$padj), ylab="-log", xlab="log2FoldChange", main = "TNBC vs. Non-TNBC")
plot( res2$log2FoldChange, -log(res2$padj), ylab="-log", xlab="log2FoldChange", main = "HER2+ vs. TNBC")

# ---------------- 7. Over-representation analysis ----------------
# ClusterProfiler to identify Gene Ontology terms that contain more differential
# expressed genes than expected by chance for the pairwise comparisons

# GO of pairwise comparison number 1 
data     <-  resdata1
geneList <- as.character(data$Gene)
data(geneList)
de       <- names(geneList)[abs(geneList) > 2]

gene.df1 <- bitr(de, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = "org.Hs.eg.db")

# GO of pairwise comparison number 2
data     <-  resdata2
geneList <- as.character(data$Gene)
data(geneList)
de       <- names(geneList)[abs(geneList) > 2]

gene.df2 <- bitr(de, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = "org.Hs.eg.db")

## GO over-representation analysis
# The clusterProfiler package for gene ontology over-representation test.

ego1  <- enrichGO(gene         = gene.df1$ENSEMBL,
                 OrgDb         = "org.Hs.eg.db",
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 readable      = TRUE)

ego2  <- enrichGO(gene         = gene.df2$ENSEMBL,
                 OrgDb         = "org.Hs.eg.db",
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 readable      = TRUE)

# ---------------- Visualization of functional enrichment result ---------------- 

## Mutate the results
x1 <- mutate(ego1, qscore = -log(p.adjust, base=10))
y1 <- mutate(ego1, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
x2 <- mutate(ego2, qscore = -log(p.adjust, base=10))
y2 <- mutate(ego2, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

## simplify the output
simplify(y1, cutoff = 0.7, by = "p.adjust")
simplify(y2, cutoff = 0.7, by = "p.adjust")

## dotPlot showing the over-representation 
dotplot(x1, showCategory=14) + ggtitle("Gene Ratio per GO term")
dotplot(x2, showCategory=14) + ggtitle("Gene Ratio per GO term")

## ggPlot for the enriched disease ontology 
ggplot(y1, showCategory = 14, 
       aes(richFactor, fct_reorder(Description, richFactor))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + 
  xlab("rich factor") +
  ylab(NULL) + 
  ggtitle("Enriched Ontology")

ggplot(y2, showCategory = 14, 
       aes(richFactor, fct_reorder(Description, richFactor))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + 
  xlab("rich factor") +
  ylab(NULL) + 
  ggtitle("Enriched Ontology")

