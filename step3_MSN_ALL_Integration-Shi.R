
############################################################################################################################################
############################################### Load/Install Packages ######################################################################
############################################################################################################################################

# Load or install packages
LIBLOC<-"/home/c.c2015805/R-Rstudio/x86_64-pc-linux-gnu-library/4.2/"
.libPaths(LIBLOC)


############################################################################################################################
######## from here check that everything has loaded and if not select what is missing in the Packages list manually ######## 
############################################################################################################################

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
library(Matrix.utils,lib.loc=LIBLOC) # no
library(matrixStats,lib.loc=LIBLOC)
library(mixtools,lib.loc=LIBLOC)
library(monocle,lib.loc=LIBLOC)
library(monocle3,lib.loc=LIBLOC) # no
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
library(SCopeLoomR,lib.loc=LIBLOC) # no
library(sctransform,lib.loc=LIBLOC)
library(Seurat,lib.loc=LIBLOC)
library(SeuratData,lib.loc=LIBLOC) # no
library(SeuratDisk,lib.loc=LIBLOC)
library(SeuratWrappers,lib.loc=LIBLOC) # no
library(SingleCellExperiment,lib.loc=LIBLOC)
library(stringr,lib.loc=LIBLOC)
library(SummarizedExperiment,lib.loc=LIBLOC)
library(terra,lib.loc=LIBLOC)
library(tibble,lib.loc=LIBLOC)
library(tidygraph,lib.loc=LIBLOC)
library(tidyr,lib.loc=LIBLOC)
library(tximport,lib.loc=LIBLOC)
library(velocyto.R,lib.loc=LIBLOC) # no
library(zoo,lib.loc=LIBLOC)

library(scCustomize,lib.loc=LIBLOC)
library(scales,lib.loc=LIBLOC)
library(ggalluvial,lib.loc=LIBLOC)
library(MAST,lib.loc=LIBLOC)
library(metap,lib.loc=LIBLOC)
library(future,lib=LIBLOC)
library(SCpubr,lib=LIBLOC)

setwd("/scratch/scw1229/Emily/Integration_Shi_FULL")

####################################################################################################################################
######################################### Import data & Create Seurat Object #######################################################
####################################################################################################################################

ALL_MSN_extended = LoadH5Seurat("/scratch/scw1229/Seurat_2023/Raw_iPSC_datasets_for_integration/ALL_MSN_extended_raw.h5Seurat")

Shi_all = LoadH5Seurat("/scratch/scw1229/Seurat_2023/Integration_Shi/Shi_all_raw.h5Seurat")


#######################################################################################################################
######################################### Assign Groups to data #######################################################
#######################################################################################################################

### Assign Sample Types
## iPSC - "KOLF2_D13", "KOLF2_D18", "KOLF2_D24_1", "KOLF2_D40_1", "KOLF2_D24_2", "KOLF2_D40_2", "i202_D24_1", "i202_D24_2", "i202_D40_1", "i202_D40_2"


ALL_MSN_extended[["SampleType"]] = case_when(ALL_MSN_extended@meta.data$orig.ident %in% c("KOLF2_D13", "KOLF2_D18", "KOLF2_D24_1", "KOLF2_D40_1", "KOLF2_D24_2", 
                                                                                          "KOLF2_D40_2", "i202_D24_1", "i202_D24_2", "i202_D40_1", "i202_D40_2") ~ 'iPSC')
Shi_all[["SampleType"]] = case_when(Shi_all@meta.data$orig.ident %in% c("GE_9W", "GE_12W_1", "GE_12W_2", "GE_13W", "GE_16W",
                                                                        "CGE_18W","GE_18W","LGE_18W","MGE_18W") ~ 'fetal_GE')

### Assign TimePoint
## D13 - "KOLF2_D13"
## D18 - "KOLF2_D18"
## D24 - "KOLF2_D24_1", "KOLF2_D24_2", "i202_D24_1", "i202_D24_2
## D40 - "KOLF2_D40_1", "KOLF2_D40_2", "i202_D40_1", "i202_D40_2"
ALL_MSN_extended[["TimePoint"]] = case_when(ALL_MSN_extended@meta.data$orig.ident %in% c("KOLF2_D13") ~ 'D13',
                                            ALL_MSN_extended@meta.data$orig.ident %in% c("KOLF2_D18") ~ 'D18',
                                            ALL_MSN_extended@meta.data$orig.ident %in% c("KOLF2_D24_1", "KOLF2_D24_2", "i202_D24_1", "i202_D24_2") ~ 'D24',
                                            ALL_MSN_extended@meta.data$orig.ident %in% c("KOLF2_D40_1", "KOLF2_D40_2", "i202_D40_1", "i202_D40_2") ~ 'D40')

Shi_all[["TimePoint"]] = case_when(Shi_all@meta.data$orig.ident %in% c("GE_9W") ~ '9PCW',
                                   Shi_all@meta.data$orig.ident %in% c("GE_12W_1") ~ '12PCW',
                                   Shi_all@meta.data$orig.ident %in% c("GE_12W_2") ~ '12PCW',
                                   Shi_all@meta.data$orig.ident %in% c("GE_13W") ~ '13PCW',
                                   Shi_all@meta.data$orig.ident %in% c("GE_16W") ~ '16PCW',
                                   Shi_all@meta.data$orig.ident %in% c("CGE_18W") ~ '18PCW',
                                   Shi_all@meta.data$orig.ident %in% c("GE_18W") ~ '18PCW',
                                   Shi_all@meta.data$orig.ident %in% c("LGE_18W") ~ '18PCW',
                                   Shi_all@meta.data$orig.ident %in% c("MGE_18W") ~ '18PCW')

###################################################################################################    
#################################### Filter Parameters ############################################
###################################################################################################

#### QC filters
ALL_MSN_extended = subset(ALL_MSN_extended, subset = nFeature_RNA >= 300 & nCount_RNA <= 150000 & percent.mt < 10)

Shi_all = subset(Shi_all, subset = nFeature_RNA >= 800 & nFeature_RNA <= 4000 & nCount_RNA <= 15000 & percent.mt < 5)

###################################################################################################    
#################################### Combine Datasets #############################################
###################################################################################################

## Create a list of Seurat objects 
# Control iPSCs + Shi

iPSC.Shi.list = list()
iPSC.Shi.list[["iPSC"]] = ALL_MSN_extended
iPSC.Shi.list[["Shi"]] = Shi_all

View(iPSC.Shi.list$iPSC@meta.data)
View(iPSC.Shi.list$Shi@meta.data)

########################################################################################################################################  
############################################   Run SCT and regress out 'percent.mt'   ##################################################
########################################################################################################################################

# Control iPSCs + Shi 

iPSC.Shi.list = lapply(X = iPSC.Shi.list, FUN = function(x) {
  x <- SCTransform(x, vst.flavor = "v2", vars.to.regress="percent.mt")})
DefaultAssay(iPSC.Shi.list$iPSC) <- "SCT"
DefaultAssay(iPSC.Shi.list$Shi) <- "SCT"
features = SelectIntegrationFeatures(object.list = iPSC.Shi.list, nfeatures = 3000)
iPSC.Shi.list = PrepSCTIntegration(object.list = iPSC.Shi.list, anchor.features = features)


######################################################################################################################
#################################### Full Integration of iPSC with Published Data ####################################
######################################################################################################################

# Control iPSCs + Shi 
iPSC.Shi.anchors = FindIntegrationAnchors(object.list = iPSC.Shi.list, normalization.method = "SCT",
                                          anchor.features = features)
iPSC.Shi.sct = IntegrateData(anchorset = iPSC.Shi.anchors, normalization.method = "SCT")

# Save fully integrated data
SaveH5Seurat(iPSC.Shi.sct,
             filename = paste0("MSNiPSC.Shi_FULLintegrated.h5Seurat"),
             overwrite = TRUE,
             verbose = TRUE)

# iPSC.Shi.sct = LoadH5Seurat("/scratch/scw1229/Emily/Integration_Shi_FULL/MSNiPSC.Shi_FULLintegrated.h5Seurat")


#####################################################################################

rm(list=ls())
setwd("/scratch/scw1229/Emily/Integration_Shi_SHORT")

####################################################################################################################################
######################################### Import data & Create Seurat Object #######################################################
####################################################################################################################################

ALL_MSN_extended = LoadH5Seurat("/scratch/scw1229/Seurat_2023/Raw_iPSC_datasets_for_integration/ALL_MSN_extended_raw.h5Seurat")

Shi_all = LoadH5Seurat("/scratch/scw1229/Seurat_2023/Integration_Shi/Shi_all_raw.h5Seurat")


#######################################################################################################################
######################################### Assign Groups to data #######################################################
#######################################################################################################################

### Assign Sample Types
## iPSC - "KOLF2_D13", "KOLF2_D18", "KOLF2_D24_1", "KOLF2_D40_1", "KOLF2_D24_2", "KOLF2_D40_2", "i202_D24_1", "i202_D24_2", "i202_D40_1", "i202_D40_2"


ALL_MSN_extended[["SampleType"]] = case_when(ALL_MSN_extended@meta.data$orig.ident %in% c("KOLF2_D13", "KOLF2_D18", "KOLF2_D24_1", "KOLF2_D40_1", "KOLF2_D24_2", 
                                                                                          "KOLF2_D40_2", "i202_D24_1", "i202_D24_2", "i202_D40_1", "i202_D40_2") ~ 'iPSC')
Shi_all[["SampleType"]] = case_when(Shi_all@meta.data$orig.ident %in% c("CGE_18W","LGE_18W","MGE_18W") ~ 'fetal_GE') ## this should create empty values in SampleType column for all other orig.ident values

## create smaller Shi seurat object with just LGE, CGE, MGE data

Shi_short <- subset(Shi_all, SampleType == "fetal_GE") ## this should only have orig.ident belonging to "CGE_18W","LGE_18W","MGE_18W"

### Assign TimePoint
## D13 - "KOLF2_D13"
## D18 - "KOLF2_D18"
## D24 - "KOLF2_D24_1", "KOLF2_D24_2", "i202_D24_1", "i202_D24_2
## D40 - "KOLF2_D40_1", "KOLF2_D40_2", "i202_D40_1", "i202_D40_2"
ALL_MSN_extended[["TimePoint"]] = case_when(ALL_MSN_extended@meta.data$orig.ident %in% c("KOLF2_D13") ~ 'D13',
                                            ALL_MSN_extended@meta.data$orig.ident %in% c("KOLF2_D18") ~ 'D18',
                                            ALL_MSN_extended@meta.data$orig.ident %in% c("KOLF2_D24_1", "KOLF2_D24_2", "i202_D24_1", "i202_D24_2") ~ 'D24',
                                            ALL_MSN_extended@meta.data$orig.ident %in% c("KOLF2_D40_1", "KOLF2_D40_2", "i202_D40_1", "i202_D40_2") ~ 'D40')

Shi_short[["TimePoint"]] = case_when(Shi_short@meta.data$orig.ident %in% c("CGE_18W") ~ '18PCW',
                                     Shi_short@meta.data$orig.ident %in% c("LGE_18W") ~ '18PCW',
                                     Shi_short@meta.data$orig.ident %in% c("MGE_18W") ~ '18PCW')

###################################################################################################    
#################################### Filter Parameters ############################################
###################################################################################################

#### QC filters
ALL_MSN_extended = subset(ALL_MSN_extended, subset = nFeature_RNA >= 300 & nCount_RNA <= 150000 & percent.mt < 10)

Shi_short = subset(Shi_short, subset = nFeature_RNA >= 800 & nFeature_RNA <= 4000 & nCount_RNA <= 15000 & percent.mt < 5)

###################################################################################################    
#################################### Combine Datasets #############################################
###################################################################################################

## Create a list of Seurat objects 
# Control iPSCs + Shi 

iPSC.Shi.list = list()
iPSC.Shi.list[["iPSC"]] = ALL_MSN_extended
iPSC.Shi.list[["Shi"]] = Shi_short

View(iPSC.Shi.list$iPSC@meta.data)
View(iPSC.Shi.list$Shi@meta.data)

########################################################################################################################################  
############################################   Run SCT and regress out 'percent.mt'   ##################################################
########################################################################################################################################

# Control iPSCs + Shi 

iPSC.Shi.list = lapply(X = iPSC.Shi.list, FUN = function(x) {
  x <- SCTransform(x, vst.flavor = "v2", vars.to.regress="percent.mt")})
DefaultAssay(iPSC.Shi.list$iPSC) <- "SCT"
DefaultAssay(iPSC.Shi.list$Shi) <- "SCT"
features = SelectIntegrationFeatures(object.list = iPSC.Shi.list, nfeatures = 3000)
iPSC.Shi.list = PrepSCTIntegration(object.list = iPSC.Shi.list, anchor.features = features)


######################################################################################################################
#################################### Full Integration of iPSC with Published Data ####################################
######################################################################################################################

# Control iPSCs + Shi 
iPSC.Shi.anchors = FindIntegrationAnchors(object.list = iPSC.Shi.list, normalization.method = "SCT",
                                          anchor.features = features)
iPSC.Shi.sct.short = IntegrateData(anchorset = iPSC.Shi.anchors, normalization.method = "SCT")

# Save fully integrated data
SaveH5Seurat(iPSC.Shi.sct.short,
             filename = paste0("MSNiPSC.Shi_SHORTintegrated.h5Seurat"),
             overwrite = TRUE,
             verbose = TRUE)

# iPSC.Shi.sct = LoadH5Seurat("/scratch/scw1229/Emily/Integration_Shi_SHORT/MSNiPSC.Shi_SHORTintegrated.h5Seurat")
