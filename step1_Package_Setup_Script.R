
############################################################################################################################################
############################################### Load/Install Packages ######################################################################
############################################################################################################################################


# Load or install packages to your home directory
LIBLOC<-"/home/c.c2015805/R-Rstudio/x86_64-pc-linux-gnu-library/4.2"
.libPaths(LIBLOC)

if (length(old.packages(LIBLOC)[1])>0)
  update.packages(LIBLOC)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager",lib=LIBLOC)
BiocManager::install(version = "3.16")  ## update version

##install.packages("scCustomize",lib=LIBLOC)

BASIC <- c("devtools","BiocManager","remotes")

PACKAGES <- c("corrplot","data.table","Hmisc","plyr","dplyr","GGally","ggnetwork","ggplot2","ggpubr","ggrepel","ggthemes",
              "network","openxlsx","reshape2","Seurat","stringr","tidyr","rlist","parallel", "matrixStats", "lifecycle", "Rcpp", "cli", "rlang")

PACKAGES_TWO <- c("cli", "rlang", "tidygraph")

BIOCMANAGER_PACKAGES<-c("scater","scmap","edgeR",'SingleCellExperiment',"tximport","biomaRt","DESeq2","multtest",
                        "hypeR", "monocle","AUCell", "RcisTarget","GENIE3","zoo", "mixtools", "rbokeh","DT", "NMF","doMC", "doRNG",
                        "ComplexHeatmap","Nebulosa", "R2HTML", "Rtsne", "pcaMethods")

GITHUB_PACKAGES<-c("cole-trapnell-lab/leidenbase","cole-trapnell-lab/monocle3", "aertslab/SCopeLoomR", "velocyto-team/velocyto.R", "mojaveazure/seurat-disk")


install.packages("harmony", version = "3.8", lib=LIBLOC) # done

for (x in BASIC){
  if (!require(x,lib.loc=LIBLOC,character.only = TRUE)){
    install.packages(x,lib=LIBLOC)
    library(x,lib.loc=LIBLOC,character.only = TRUE,dependencies=TRUE)}}

### BASIC
library(BiocManager,lib.loc=LIBLOC)
library(devtools,lib.loc=LIBLOC)
library(remotes,lib.loc=LIBLOC)

for (x in PACKAGES){
  if (!require(x,lib.loc=LIBLOC,character.only = TRUE)){
    install.packages(x,lib=LIBLOC)
    library(x,lib.loc=LIBLOC,character.only = TRUE)}}

### PACKAGES
library(corrplot,lib.loc=LIBLOC)
library(data.table,lib.loc=LIBLOC)
library(dplyr,lib.loc=LIBLOC)
library(GGally,lib.loc=LIBLOC)
library(ggnetwork,lib.loc=LIBLOC)
library(ggplot2,lib.loc=LIBLOC)
library(ggpubr,lib.loc=LIBLOC)
library(ggrepel,lib.loc=LIBLOC)
library(ggthemes,lib.loc=LIBLOC)
library(Hmisc,lib.loc=LIBLOC)
library(lifecycle,lib.loc=LIBLOC)
library(matrixStats,lib.loc=LIBLOC)
library(network,lib.loc=LIBLOC)
library(openxlsx,lib.loc=LIBLOC)
library(parallel,lib.loc=LIBLOC)
library(plyr,lib.loc=LIBLOC)
library(Rcpp,lib.loc=LIBLOC)
library(reshape2,lib.loc=LIBLOC)
library(rlist,lib.loc=LIBLOC) 
library(Seurat,lib.loc=LIBLOC)
library(stringr,lib.loc=LIBLOC)
library(tidyr,lib.loc=LIBLOC)

for (x in PACKAGES_TWO){
  if (!require(x,lib.loc=LIBLOC,character.only = TRUE)){
    install.packages(x,lib=LIBLOC)
    library(x,lib.loc=LIBLOC,character.only = TRUE)}}

### PACKAGES_TWO
library(cli,lib.loc=LIBLOC)
library(rlang,lib.loc=LIBLOC)
library(tidygraph)

for (x in BIOCMANAGER_PACKAGES){
  if (!require(x,lib.loc=LIBLOC,character.only = TRUE)){
    BiocManager::install(x,lib=LIBLOC)
    library(x,lib.loc=LIBLOC,character.only = TRUE)}}

### BIOCMANAGER_PACKAGES
library(AUCell,lib.loc=LIBLOC)
library(biomaRt,lib.loc=LIBLOC)
library(ComplexHeatmap,lib.loc=LIBLOC)
library(DESeq2,lib.loc=LIBLOC)
library(doMC,lib.loc=LIBLOC)
library(doRNG,lib.loc=LIBLOC)
library(DT,lib.loc=LIBLOC)
library(edgeR,lib.loc=LIBLOC)
library(GENIE3,lib.loc=LIBLOC)
library(hypeR,lib.loc=LIBLOC)
library(mixtools,lib.loc=LIBLOC)
library(monocle,lib.loc=LIBLOC)
library(multtest,lib.loc=LIBLOC)
library(Nebulosa,lib.loc=LIBLOC)
library(NMF,lib.loc=LIBLOC)
library(pcaMethods,lib.loc=LIBLOC)
library(R2HTML,lib.loc=LIBLOC)
library(rbokeh,lib.loc=LIBLOC)
library(RcisTarget,lib.loc=LIBLOC)
library(Rtsne,lib.loc=LIBLOC)
library(scater,lib.loc=LIBLOC)
library(scmap,lib.loc=LIBLOC)
library(SingleCellExperiment,lib.loc=LIBLOC)
### library(tximeta,lib.loc=LIBLOC) could not install 14/07/2022
library(tximport,lib.loc=LIBLOC)
library(zoo,lib.loc=LIBLOC)

for (x in GITHUB_PACKAGES){
  if (!require(unlist(strsplit(x,split="/"))[2],character.only = TRUE)){
    ### BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats','limma', 'S4Vectors', 'SingleCellExperiment', 'SummarizedExperiment', 
    ###   'lme4', 'batchelor', 'Matrix.utils', 'HDF5Array', 'terra','Cairo', 'ggrastr'),lib=LIBLOC)
    devtools::install_github(x,lib=LIBLOC)
    library(x,lib.loc=LIBLOC,character.only = TRUE)}}


devtools::install_version('speedglm', '0.3-4', repos = 'https://packagemanager.rstudio.com/cran/2023-03-31')
devtools::install_github("cole-trapnell-lab/monocle3")

# monocle installation
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'), lib=LIBLOC, force=TRUE)

install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')

remotes::install_github("cvarrichio/Matrix.utils")
remotes::install_github("igraph/rigraph", ref = "v1.5.0", lib=LIBLOC)

library(monocle3)

remotes::install_github("mojaveazure/seurat-disk")
devtools::install_github('satijalab/seurat-data') # non-zero exit status

### GITHUB_PACKAGES
library(batchelor,lib.loc=LIBLOC)
library(BiocGenerics,lib.loc=LIBLOC)
library(Cairo,lib.loc=LIBLOC)
library(DelayedArray,lib.loc=LIBLOC)
library(DelayedMatrixStats,lib.loc=LIBLOC)
library(ggrastr,lib.loc=LIBLOC)
library(HDF5Array,lib.loc=LIBLOC)
library(limma,lib.loc=LIBLOC)
library(lme4,lib.loc=LIBLOC)
library(Matrix.utils,lib.loc=LIBLOC) #Tried to install manually, not available for this verison of R
library(S4Vectors,lib.loc=LIBLOC)
library(SummarizedExperiment,lib.loc=LIBLOC)
library(terra,lib.loc=LIBLOC)

devtools::install_github("aertslab/SCENIC", force = TRUE) 
packageVersion("SCENIC")
library(SCENIC)

library(R.utils)
remotes::install_github('satijalab/seurat-wrappers', lib=LIBLOC)

library(cowplot)
library(hdf5r)
library(Matrix)
library(patchwork)
library(RColorBrewer)
library(SeuratWrappers) 
library(tibble)



BiocManager::install("glmGamPoi")
#devtools::install_github("satijalab/sctransform", ref = "develop")
library(glmGamPoi,lib.loc=LIBLOC)
#library(sctransform,lib.loc=LIBLOC)

############################################################################################################################
######## from here check that everything has loaded and if not select what is missing in the Packages list manually ######## 
############################################################################################################################

.libPaths(LIBLOC)


library(AUCell,lib.loc=LIBLOC)
library(batchelor,lib.loc=LIBLOC)
library(BiocGenerics,lib.loc=LIBLOC)
library(BiocManager,lib.loc=LIBLOC)
library(biomaRt,lib.loc=LIBLOC)
library(Cairo,lib.loc=LIBLOC)
library(cli,lib.loc=LIBLOC)
library(ComplexHeatmap,lib.loc=LIBLOC)
library(corrplot,lib.loc=LIBLOC)
library(cowplot,lib.loc=LIBLOC)
library(data.table,lib.loc=LIBLOC)
library(DelayedArray,lib.loc=LIBLOC)
library(DelayedMatrixStats,lib.loc=LIBLOC)
library(DESeq2,lib.loc=LIBLOC)
library(devtools,lib.loc=LIBLOC)
library(doMC,lib.loc=LIBLOC)
library(doRNG,lib.loc=LIBLOC)
library(dplyr,lib.loc=LIBLOC)
library(DT,lib.loc=LIBLOC)
library(edgeR,lib.loc=LIBLOC)
library(GENIE3,lib.loc=LIBLOC)
library(GGally,lib.loc=LIBLOC)
library(ggnetwork,lib.loc=LIBLOC)
library(ggplot2,lib.loc=LIBLOC)
library(ggpubr,lib.loc=LIBLOC)
library(ggrastr,lib.loc=LIBLOC)
library(ggrepel,lib.loc=LIBLOC)
library(ggthemes,lib.loc=LIBLOC)
library(glmGamPoi,lib.loc=LIBLOC) 
library(harmony,lib.loc=LIBLOC)
library(HDF5Array,lib.loc=LIBLOC)
library(hdf5r,lib.loc=LIBLOC)
library(Hmisc,lib.loc=LIBLOC)
library(hypeR,lib.loc=LIBLOC)
library(leidenbase,lib.loc=LIBLOC)
library(lifecycle,lib.loc=LIBLOC)
library(limma,lib.loc=LIBLOC)
library(lme4,lib.loc=LIBLOC)
library(Matrix,lib.loc=LIBLOC)
library(Matrix.utils,lib.loc=LIBLOC) 
library(matrixStats,lib.loc=LIBLOC)
library(mixtools,lib.loc=LIBLOC)
library(monocle,lib.loc=LIBLOC)
library(monocle3,lib.loc=LIBLOC) 
library(multtest,lib.loc=LIBLOC)
library(Nebulosa,lib.loc=LIBLOC)
library(network,lib.loc=LIBLOC)
library(NMF,lib.loc=LIBLOC)
library(openxlsx,lib.loc=LIBLOC)
library(parallel,lib.loc=LIBLOC)
library(patchwork,lib.loc=LIBLOC)
library(pcaMethods,lib.loc=LIBLOC)
library(plyr,lib.loc=LIBLOC)
library(R.utils,lib.loc=LIBLOC)
library(R2HTML,lib.loc=LIBLOC)
library(rbokeh,lib.loc=LIBLOC)
library(RcisTarget,lib.loc=LIBLOC)
library(RColorBrewer,lib.loc=LIBLOC)
library(Rcpp,lib.loc=LIBLOC)
library(remotes,lib.loc=LIBLOC)
library(reshape2,lib.loc=LIBLOC)
library(rlang,lib.loc=LIBLOC)
library(rlist,lib.loc=LIBLOC)
library(Rtsne,lib.loc=LIBLOC)
library(S4Vectors,lib.loc=LIBLOC)
library(scater,lib.loc=LIBLOC)
library(SCENIC,lib.loc=LIBLOC)
library(scmap,lib.loc=LIBLOC)
library(SCopeLoomR,lib.loc=LIBLOC)
library(sctransform,lib.loc=LIBLOC)
library(Seurat,lib.loc=LIBLOC)
library(SeuratData,lib.loc=LIBLOC) 
library(SeuratDisk,lib.loc=LIBLOC)
library(SeuratWrappers,lib.loc=LIBLOC)
library(SingleCellExperiment,lib.loc=LIBLOC)
library(stringr,lib.loc=LIBLOC)
library(SummarizedExperiment,lib.loc=LIBLOC)
library(terra,lib.loc=LIBLOC) 
library(tibble,lib.loc=LIBLOC)
library(tidygraph,lib.loc=LIBLOC)
library(tidyr,lib.loc=LIBLOC)
### library(tximeta,lib.loc=LIBLOC) could not install 
library(tximport,lib.loc=LIBLOC)
library(velocyto.R,lib.loc=LIBLOC) 
library(zoo,lib.loc=LIBLOC)

# re-check
library(BiocManager)
library(cli)
library(devtools)
library(dplyr)
library(glmGamPoi,lib.loc=LIBLOC)
library(harmony)
library(lifecycle)
library(remotes)
library(rlang)
library(sctransform,lib.loc=LIBLOC)
