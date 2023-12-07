FROM bioconductor/bioconductor_docker:RELEASE_3_18
LABEL maintainer="Jielin Yang"
LABEL email="jielin.yang@mail.utoronto.ca"
LABEL version="0.1.0"
LABEL description="Docker image for the IntegraTRN R package"

# Install R packages
RUN install2.r \
-d TRUE \
-r "https://cran.rstudio.com" \
RColorBrewer ggplot2 devtools \
rmarkdown knitr \
data.table \
tidyverse \
tidyr shiny \
shinyBS \
plotly \
kableExtra \
stringr stringi \
DT \
reshape2 here \
cowplot \
ggrepel \
readxl \
randomForest \
igraph \
networkD3 \
Rdpack \
MatchIt optmatch \
doParallel \
doRNG \
ggpubr \
testthat


# Install Bioconductor packages
RUN R --vanilla -e 'BiocManager::install(c("BiocGenerics", "Biostrings", \
"BSgenome", "BSgenome.Hsapiens.UCSC.hg38", \
"ChIPseeker", "ComplexHeatmap", "DESeq2", "edgeR", "GENIE3", "GenomicRanges", \
"monaLisa", "SummarizedExperiment", "TxDb.Hsapiens.UCSC.hg38.knownGene", \
"TxDb.Mmusculus.UCSC.mm10.knownGene", "VariantAnnotation", "vsn", "WGCNA", \
"JASPAR2022", "org.Hs.eg.db", "org.Mm.eg.db", "TFBSTools"))'

# Install IntegraTRN from GitHub
RUN R --vanilla -e 'devtools::install_github("j-y26/IntegraTRN", build_vignettes = TRUE)'
