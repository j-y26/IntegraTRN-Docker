## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
# Load the package
library(IntegraTRN)
library(dplyr)

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
# Since devtools does not install the vignette dependencies, we need to install
# these Suggested packages manually
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
}
if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
  BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}
if (!requireNamespace("ggpubr", quietly = TRUE)) {
  install.packages("ggpubr")
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}

## ----message = FALSE, eval=FALSE----------------------------------------------
#  require("devtools")
#  devtools::install_github("j-y26/IntegraTRN", build_vignettes = TRUE)
#  library("IntegraTRN")

## ----message = FALSE, eval=FALSE----------------------------------------------
#  ls("package:IntegraTRN")

## ----message = FALSE, tidy=TRUE-----------------------------------------------
data("RNAseq_heart") # RNAseq count matrix
data("RNAseq_heart_samples") # RNAseq sample information
data("smallRNAseq_heart") # small RNAseq count matrix
data("smallRNAseq_heart_samples") # small RNAseq sample information

## ----message = FALSE----------------------------------------------------------
myMOList <- MOList(
  RNAseq = RNAseq_heart,
  RNAGroupBy = RNAseq_heart_samples$Age,
  smallRNAseq = smallRNAseq_heart,
  smallRNAGroupBy = smallRNAseq_heart_samples$Age
)

## ----message = FALSE, tidy=TRUE-----------------------------------------------
atacPeak1Path <- system.file("extdata", "peak1.bed", package = "IntegraTRN")
atacPeak2Path <- system.file("extdata", "peak2.bed", package = "IntegraTRN")

## ----message = FALSE----------------------------------------------------------
myMOList <- MOList(myMOList,
  ATACpeak1 = atacPeak1Path,
  ATACpeak2 = atacPeak2Path
)

## ----message = FALSE, warning=FALSE-------------------------------------------
myMOList <- diffOmics(myMOList,
  rnaseqBatch = RNAseq_heart_samples$Batch,
  program = "DESeq2"
)

## ----message = FALSE, eval=FALSE----------------------------------------------
#  # NOT run, just illustrate the usage of the batch argument
#  myMOList <- diffOmics(myMOList,
#    smallRnaBatch = smallRNAseq_heart_samples$Sex,
#    program = "edgeR"
#  )

## ----message = FALSE, tidy=TRUE, fig.align='center', fig.width=6, fig.height=4, fig.cap="Figure 1. Volcano plot of RNAseq data"----
plotVolcano(myMOList, omic = "RNAseq", log2FC = 0.1)

## ----message = FALSE, tidy=TRUE-----------------------------------------------
myRNADETag <- myMOList$DERNAseq

## ----message = FALSE, tidy=TRUE-----------------------------------------------
myRNADETag

## ----message = FALSE, tidy=TRUE-----------------------------------------------
exportDE(myRNADETag, original = TRUE) %>% head(5)

## ----message = FALSE, warning=FALSE-------------------------------------------
myRNATOPTag <- TOPTag(myRNADETag,
  pCutoff = 0.05,
  logFCCutoff = 0.1,
  topGenes = 200,
  direction = "both"
)

## ----message = FALSE, tidy=TRUE-----------------------------------------------
myRNATOPTag

## ----message = FALSE, tidy=TRUE-----------------------------------------------
# Export the desired results with rank information
myTopRNAResults <- exportDE(myRNATOPTag)
# Make numeric values easier to read
myTopRNAResults[, 1:4] <- format(myTopRNAResults[, 1:4],
  digits = 4,
  scientific = TRUE
)
head(myTopRNAResults)

## ----message = FALSE----------------------------------------------------------
genesToLabel <- exportDE(myRNADETag) %>%
  dplyr::filter(padj < 0.05, abs(logFC) > 0.1) %>%
  head(3) %>%
  rownames()

## ----message = FALSE, fig.align='center', fig.width=6, fig.height=4, fig.cap="Figure 2. Volcano plot of RNAseq data with some genes labelled"----
rnaVPlot <- plotVolcano(myMOList, omic = "RNAseq", log2FC = 0.1)
rnaVPlot <- rnaVPlot + ggrepel::geom_label_repel(
  data = rnaVPlot$data[rownames(rnaVPlot$data) %in% genesToLabel, ],
  mapping = ggplot2::aes(
    x = logFC,
    y = -log10(padj),
    label = rownames(rnaVPlot$data[rownames(rnaVPlot$data) %in% genesToLabel, ])
  ),
  box.padding = 0.5,
  size = 3
)
rnaVPlot

## ----message = FALSE----------------------------------------------------------
smallRNAAnnotation <- data.frame(
  transcript = c("hsa-miR-1-3p", "hsa-miR-1-5p", "hsa-miR-2-5p"),
  type = c("miRNA", "miRNA", "miRNA")
)

## ----message = FALSE, eval=FALSE----------------------------------------------
#  # NOT run, just illustrate the usage of the smallRNAAnnotation argument
#  myMOList <- annotateSmallRNA(myMOList, anno = smallRNAAnnotation)

## ----message = FALSE----------------------------------------------------------
myMOList <- annotateSmallRNA(myMOList, anno = "human")

## ----message = FALSE, tidy=TRUE, fig.align='center', fig.width=6, fig.height=4, fig.cap="Figure 3. Volcano plot of small RNAseq data"----
plotVolcano(myMOList, omic = "smallRNAseq", log2FC = 0)

## ----message = FALSE, tidy=TRUE, fig.align='center', fig.width=6, fig.height=4, fig.cap="Figure 4. Volcano plot of small RNAseq data with each type of small RNA coloured"----
# Use the RColorBrewer package to assign a color palette
library("RColorBrewer")
smallRNAColors <- brewer.pal(6, "BuPu")

# Plot Volcano plot of small RNAseq data
plotVolcanoSmallRNA(myMOList, log2FC = 0, colScheme = smallRNAColors)

## ----message = FALSE, tidy=TRUE, results='hide'-------------------------------
pcaPlotList <- plotSmallRNAPCAs(myMOList)

## ----message = FALSE, tidy=TRUE-----------------------------------------------
length(pcaPlotList)

## ----message = FALSE, tidy=TRUE, fig.align='center', fig.width=6, fig.height=4, fig.cap="Figure 5. PCA plot of miRNA"----
pcaPlotList$miRNA

## ----message = FALSE, tidy=TRUE, fig.align='center', fig.width=9, fig.height=6, fig.cap="Figure 6. Combined PCA plot of all small RNA types"----
ggpubr::ggarrange(pcaPlotList$miRNA,
  pcaPlotList$circRNA,
  pcaPlotList$piRNA,
  pcaPlotList$snoRNA,
  pcaPlotList$snRNA,
  pcaPlotList$tRNA,
  ncol = 3, nrow = 2
)

## ----message = FALSE, tidy=TRUE, warning=FALSE--------------------------------
# Annotate ATAC Peaks
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("org.Hs.eg.db")
library("BSgenome.Hsapiens.UCSC.hg38")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
annoDb <- "org.Hs.eg.db"
bsgenome <- BSgenome.Hsapiens.UCSC.hg38

## ----message = FALSE, tidy=TRUE-----------------------------------------------
# Load PWMs from JASPAR
data("jasparVertebratePWM")

## ----message = FALSE, warning=FALSE, results='hide'---------------------------
myMOList <- annotateATACPeaksMotif(myMOList,
  tssRegion = c(-3000, 3000),
  TxDb = txdb,
  annoDb = annoDb,
  bsgenome = bsgenome,
  pwmL = jasparVertebratePWM
)

## ----message = FALSE, tidy=TRUE, fig.align='center', fig.width=6, fig.height=4, fig.cap="Figure 7. Annotation of ATAC peaks"----
plotATACAnno(myMOList)

## ----message = FALSE, tidy=TRUE, fig.align='center', fig.width=6, fig.height=4, fig.cap="Figure 8. ATACseq coverage"----
plotATACCoverage(myMOList)

## ----message = FALSE, tidy=TRUE-----------------------------------------------
myMotifs <- as.data.frame(myMOList$DEATAC)
head(myMotifs[, 1:8]) # Many annotations, see the first 8 attributes

## ----message = FALSE, tidy=TRUE, fig.align='center', fig.width=8, fig.height=4, fig.cap="Figure 9. Motif enrichment on ATACseq peaks"----
plotATACMotifHeatmap(myMOList, pValue = 0.01)

## ----message = FALSE, tidy=FALSE, warning=FALSE-------------------------------
myMOList <- matchSamplesRNAsmallRNA(myMOList,
  sampleDFRNAseq = RNAseq_heart_samples,
  sampleDFSmallRNAseq = smallRNAseq_heart_samples
)

## ----message = FALSE, tidy=TRUE-----------------------------------------------
exportMatchResult(myMOList) %>% head(5) # For the first 5 matches

## ----message = FALSE, tidy=TRUE-----------------------------------------------
nrow(exportMatchResult(myMOList))

## ----message = FALSE, tidy=TRUE-----------------------------------------------
data("miR2Genes")
data("tf2Genes")
myMOList <- loadExtInteractions(myMOList,
  miR2Genes = miR2Genes,
  tf2Genes = tf2Genes
)

## ----message = FALSE----------------------------------------------------------
omiCutoffs <- setOmicCutoffs(
  rnaAdjPval = 0.05,
  rnaLogFC = 0.1,
  rnaTopGenes = 50,
  smallRNAAdjPval = 0.05,
  smallRNALogFC = 0.1,
  smallRNATopGenes = 50,
  atacMotifPval = 0.01
)

## ----message = FALSE, tidy=TRUE-----------------------------------------------
omiCutoffs

## ----message = FALSE, warning=FALSE-------------------------------------------
myTRNet <- constructTRN(myMOList,
  omiCutoffs = omiCutoffs,
  smallRNAtypes = "all",
  targetDirection = "both",
  predicted = TRUE,
  nthreads = 1,
  seed = 123
)

## ----message = FALSE, tidy=TRUE, warning=FALSE, fig.align='center', fig.width=6, fig.height=5, fig.cap="Figure 10. The Transcriptional Regulatory Network"----
plotNetwork(myTRNet, interactive = TRUE)

## ----message = FALSE----------------------------------------------------------
omiCutoffs2 <- omiCutoffs
omiCutoffs2$rnaAdjPval <- 0.01
omiCutoffs2$rnaLogFC <- 0.2

## ----message = FALSE, warning=FALSE-------------------------------------------
myTRNet2 <- constructTRN(myMOList,
  omiCutoffs = omiCutoffs2,
  smallRNAtypes = "miRNA",
  targetDirection = "up",
  predicted = TRUE
)

## ----message = FALSE, tidy=TRUE, warning=FALSE, fig.align='center', fig.width=6, fig.height=5, fig.cap="Figure 11. The Transcriptional Regulatory Network of Up-regulated Genes"----
plotNetwork(myTRNet2, vertex.size = 10)

## ----message = FALSE, tidy=TRUE, warning=FALSE--------------------------------
exportIgraph(myTRNet2)

## ----message = FALSE, tidy=TRUE-----------------------------------------------
sessionInfo("IntegraTRN")

