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
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  BiocManager::install("ComplexHeatmap")
}
if (!requireNamespace("circlize", quietly = TRUE)) {
  install.packages("circlize")
}
if (!requireNamespace("multiMiR", quietly = TRUE)) {
  BiocManager::install("multiMiR")
}

## ----message = FALSE, tidy=TRUE-----------------------------------------------
data("RNAseq_heart") # RNAseq count matrix
data("RNAseq_heart_samples") # RNAseq sample information
data("smallRNAseq_heart") # small RNAseq count matrix
data("smallRNAseq_heart_samples") # small RNAseq sample information
data("protein_heart") # Proteomics count matrix
data("protein_heart_samples") # Proteomics sample information
data("miR2Genes") # miRNA-target interactions, externally curated
data("tf2Genes") # TF-target interactions, externally curated
data("proteinGeneIDConvert") # Protein ID conversion table

## ----message = FALSE, tidy=TRUE-----------------------------------------------
atacPeak1Path <- system.file("extdata", "peak1.bed", package = "IntegraTRN")
atacPeak2Path <- system.file("extdata", "peak2.bed", package = "IntegraTRN")

## ----warning=FALSE------------------------------------------------------------
myMOList <- MOList(
  RNAseq = RNAseq_heart,
  RNAGroupBy = RNAseq_heart_samples$Age,
  smallRNAseq = smallRNAseq_heart,
  smallRNAGroupBy = smallRNAseq_heart_samples$Age,
  proteomics = protein_heart,
  proteomicsGroupBy = protein_heart_samples$Age,
  ATACpeak1 = atacPeak1Path,
  ATACpeak2 = atacPeak2Path
)

## ----message = FALSE, tidy=TRUE-----------------------------------------------
myMOList <- diffOmics(myMOList, program = "edgeR")

## ----message = FALSE, tidy=TRUE-----------------------------------------------
myRNADETag <- myMOList$DERNAseq

## ----message = FALSE, tidy=TRUE-----------------------------------------------
exportDE(myRNADETag) %>% head(5)

## ----message = FALSE, warning=FALSE-------------------------------------------
myRNATOPTag <- TOPTag(myRNADETag,
  pCutoff = 0.05,
  logFCCutoff = 0,
  topGenes = 20,
  direction = "both"
)

## ----message = FALSE, tidy=TRUE-----------------------------------------------
exportDE(myRNATOPTag) %>% head(5)

## ----message = FALSE, tidy=TRUE-----------------------------------------------
topRNANormCounts <- exportNormalizedCounts(myRNATOPTag)

## ----message = FALSE, fig.align='center', fig.width=6, fig.height=6, fig.cap="Figure 1. Heatmap of top 20 differentially expressed genes"----
heatmapMatrix <- topRNANormCounts %>%
  as.matrix() %>%
  t() %>%
  scale() %>%
  t()
exprRange <- range(heatmapMatrix)
col <- circlize::colorRamp2(
  c(exprRange[1], 0, exprRange[2]),
  c("blue", "white", "red")
)
ComplexHeatmap::Heatmap(heatmapMatrix,
  col = col,
  show_column_names = FALSE
)

## ----message = FALSE, warning=FALSE-------------------------------------------
rankedResult <- TOPTag(myRNADETag,
  pCutoff = 1,
  logFCCutoff = 0,
  topGenes = 1,
  direction = "both"
) %>% exportDE()

## ----message = FALSE, tidy=TRUE-----------------------------------------------
colnames(RNAseq_heart_samples)

## ----message = FALSE, warning=FALSE-------------------------------------------
# Create the edgeR DGEList object
dgeList <- edgeR::DGEList(
  counts = RNAseq_heart,
  group = RNAseq_heart_samples$Age
)
# Filter lowly expressed genes
keep <- edgeR::filterByExpr(dgeList)
dgeList <- dgeList[keep, ]
# Perform differential expression analysis
design <- model.matrix(~ RNAseq_heart_samples$Age +
  RNAseq_heart_samples$Sex +
  RNAseq_heart_samples$Batch)
dgeList <- edgeR::calcNormFactors(dgeList, method = "TMM") # or a compatible method
dgeList <- edgeR::estimateDisp(dgeList, design)
fit <- edgeR::glmQLFit(dgeList, design) # or other suitable models
qlf <- edgeR::glmQLFTest(fit, coef = ncol(RNAseq_heart_samples))

## ----message = FALSE, tidy=TRUE-----------------------------------------------
# Obtaining the DE results and normalized count
deResult <- edgeR::topTags(qlf, n = nrow(dgeList)) %>% as.data.frame()
normCounts <- edgeR::cpm(dgeList) %>% as.matrix()

# Construct the DETag object
myNewRNADETag <- DETag(
  DEResult = deResult,
  method = "edgeR",
  normalizedCounts = normCounts
)

## ----message = FALSE, tidy=TRUE-----------------------------------------------
myMOList$DERNAseq <- myNewRNADETag

## ----message = FALSE, warning=FALSE-------------------------------------------
set.seed(91711)
# Count-based data
RNAseq_heart <- RNAseq_heart[sample(nrow(RNAseq_heart), 500), ]
smallRNAseq_heart <- smallRNAseq_heart[sample(nrow(smallRNAseq_heart), 500), ]
protein_heart <- protein_heart[sample(nrow(protein_heart), 500), ]
# Chromatin accessibility data
peak1 <- read.table(atacPeak1Path, sep = "\t", header = FALSE)
peak2 <- read.table(atacPeak2Path, sep = "\t", header = FALSE)
peak1 <- peak1[sample(nrow(peak1), 500), ]
peak2 <- peak2[sample(nrow(peak2), 500), ]
set.seed(NULL)

## ----message = FALSE, warning=FALSE-------------------------------------------
myMOList <- MOList(
  RNAseq = RNAseq_heart,
  RNAGroupBy = RNAseq_heart_samples$Age,
  smallRNAseq = smallRNAseq_heart,
  smallRNAGroupBy = smallRNAseq_heart_samples$Age,
  proteomics = protein_heart,
  proteomicsGroupBy = protein_heart_samples$Age,
  ATACpeak1 = peak1,
  ATACpeak2 = peak2
)

## ----message = FALSE, warning=FALSE-------------------------------------------
myMOList <- diffOmics(myMOList)

## ----message = FALSE----------------------------------------------------------
setOmicCutoffs(
  rnaAdjPval = 0.05,
  rnaLogFC = 0.1,
  rnaTopGenes = 200
)

## -----------------------------------------------------------------------------
setOmicCutoffs(
  rnaAdjPval = 0.05,
  rnaLogFC = 0.1,
  rnaTopGenes = 200,
  proteomicsAdjPval = 0.05,
  proteomicsLogFC = 0
)

## -----------------------------------------------------------------------------
setOmicCutoffs(
  rnaAdjPval = 0.05,
  rnaLogFC = 0.1,
  rnaTopGenes = 200,
  proteomicsAdjPval = 1,
  proteomicsLogFC = 0
)

## ----message = FALSE----------------------------------------------------------
myMOList <- setGene2Protein(myMOList, proteinGeneIDConvert)

## ----message = FALSE----------------------------------------------------------
# Annotate the small RNAs
myMOList <- annotateSmallRNA(myMOList, anno = "human")

# Perform optimal matching between the RNAseq and small RNAseq data
myMOList <- matchSamplesRNAsmallRNA(myMOList,
  sampleDFRNAseq = RNAseq_heart_samples,
  sampleDFSmallRNAseq = smallRNAseq_heart_samples
)

# Define a cutoff for RNAseq and small RNAseq
setOmicCutoffs(
  rnaAdjPval = 0.05,
  rnaLogFC = 0,
  rnaTopGenes = 0.2,
  smallRNAAdjPval = 0.05,
  smallRNALogFC = 0,
  smallRNATopGenes = 0.2
)

## ----message = FALSE, tidy=TRUE-----------------------------------------------
geneList <- exportDE(myMOList$DERNAseq) %>% rownames()
geneList[1:12]

## ----message = FALSE----------------------------------------------------------
smallRNAs <- getDESmallRNA(myMOList,
  padj = 0.05,
  log2fc = 0,
  type = "miRNA"
)
head(smallRNAs, 5)

## ----message = FALSE, warning=FALSE, eval=FALSE-------------------------------
#  mir2geneMultiMiR <- multiMiR::get_multimir(
#    org = "hsa",
#    mirna = rownames(smallRNAs),
#    table = "validated",
#    summary = FALSE
#  )
#  mir2geneMultiMiR <- mir2geneMultiMiR@data

## ----message = FALSE, warning=FALSE-------------------------------------------
data("mir2geneMultiMiR")

## ----message = FALSE, warning=FALSE-------------------------------------------
mir2geneMultiMiR <- mir2geneMultiMiR %>%
  dplyr::select(mature_mirna_id, target_symbol) %>%
  dplyr::rename(
    regulator = mature_mirna_id,
    target = target_symbol
  )

## ----message = FALSE, warning=FALSE-------------------------------------------
# Combine the two sets of interactions
miR2Genes <- as.data.frame(miR2Genes)
miR2Genes <- rbind(miR2Genes, mir2geneMultiMiR)

# Remove any potential duplications across databases
miR2Genes <- miR2Genes %>%
  dplyr::distinct() %>%
  as.list()

## ----message = FALSE----------------------------------------------------------
myMOList <- loadExtInteractions(myMOList,
  miR2Genes = miR2Genes
)

## ----message = FALSE----------------------------------------------------------
omiCutoffs <- setOmicCutoffs(
  rnaAdjPval = 0.05,
  rnaLogFC = 0,
  rnaTopGenes = 0.5,
  smallRNAAdjPval = 0.05,
  smallRNALogFC = 0,
  smallRNATopGenes = 0.5
)

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
myTRN1 <- constructTRN(myMOList,
  omiCutoffs = omiCutoffs,
  useOmics = c("RNAseq", "smallRNAseq"),
  targetDirection = "both"
)

## ----message = FALSE, warning=FALSE, fig.align='center', fig.width=7, fig.height=5, fig.cap="Figure 2. TRN constructed with without considering protein expression"----
plotNetwork(myTRN1, interactive = TRUE)

## ----message = FALSE----------------------------------------------------------
omiCutoffs <- setOmicCutoffs(
  rnaAdjPval = 0.05,
  rnaLogFC = 0.1,
  rnaTopGenes = 0.5,
  smallRNAAdjPval = 0.05,
  smallRNALogFC = 0,
  smallRNATopGenes = 0.5,
  proteomicsAdjPval = 0.05,
  proteomicsLogFC = 0
)

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
myTRN2 <- constructTRN(myMOList,
  omiCutoffs = omiCutoffs,
  useOmics = c("RNAseq", "smallRNAseq", "proteomics"),
  targetDirection = "both"
)

## ----message = FALSE, warning=FALSE, fig.align='center', fig.width=7, fig.height=5, fig.cap="Figure 3. TRN constructed with coherent protein expression"----
plotNetwork(myTRN2, interactive = TRUE)

## ----message = FALSE, warning=FALSE, include=FALSE----------------------------
myTRN3 <- constructTRN(myMOList,
  omiCutoffs = omiCutoffs,
  useOmics = c("RNAseq", "smallRNAseq", "proteomics"),
  targetDirection = "both",
  smallRNAtypes = "miRNA"
)

## ----message = FALSE, warning=FALSE, fig.align='center', fig.width=7, fig.height=5, fig.cap="Figure 4. miRNA-target gene TRN constructed with coherent protein expression"----
plotNetwork(myTRN3, interactive = TRUE)

## ----message = FALSE----------------------------------------------------------
myMOList <- loadExtInteractions(myMOList,
  tf2Genes = tf2Genes
)

## ----message = FALSE, tidy=TRUE-----------------------------------------------
sessionInfo("IntegraTRN")

