
# remove all the objects from the R session
rm(list=ls())

# Set working directory
setwd("/scratch/scw1229/Emily/Integration_Bocchi")

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
#library(BiocManager,lib.loc=LIBLOC) # issue
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
#library(monocle3,lib.loc=LIBLOC)
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
#library(SCopeLoomR,lib.loc=LIBLOC) # issue
library(sctransform,lib.loc=LIBLOC)
library(Seurat,lib.loc=LIBLOC)
#library(SeuratData,lib.loc=LIBLOC) # issue
library(SeuratDisk,lib.loc=LIBLOC)
#library(SeuratWrappers,lib.loc=LIBLOC) # issue
library(SingleCellExperiment,lib.loc=LIBLOC)
library(stringr,lib.loc=LIBLOC)
library(SummarizedExperiment,lib.loc=LIBLOC)
library(terra,lib.loc=LIBLOC)
library(tibble,lib.loc=LIBLOC)
library(tidygraph,lib.loc=LIBLOC)
library(tidyr,lib.loc=LIBLOC)
### library(tximeta,lib.loc=LIBLOC) could not install 14/07/2022
library(tximport,lib.loc=LIBLOC)
# library(velocyto.R,lib.loc=LIBLOC) # issue
library(zoo,lib.loc=LIBLOC)

library(scCustomize,lib.loc=LIBLOC)
library(scales,lib.loc=LIBLOC)
library(ggalluvial,lib.loc=LIBLOC)

# re-check
library(BiocManager)
library(cli)
library(devtools)
library(dplyr)
library(lifecycle)
library(remotes)
library(rlang)

#BiocManager::install("MAST",lib=LIBLOC)
library(MAST,lib.loc=LIBLOC)

#install.packages('metap',lib=LIBLOC)
library(metap,lib.loc=LIBLOC)

#install.packages('future',lib=LIBLOC)
library(future,lib=LIBLOC)

#install.packages('SCpubr',lib=LIBLOC)
library(SCpubr,lib=LIBLOC)

# attempting installation of velocyto
library(devtools)
devtools::install_github("velocyto-team/velocyto.R")



############################################
############################################

setwd("/scratch/scw1229/Emily/Integration_Bocchi")

res05_seurat = LoadH5Seurat("/scratch/scw1229/Emily/Integration_Bocchi/iPSCBocchi_res05_FULLintegrated_MF.h5Seurat")
view(res05_seurat@meta.data)

#res05_seurat[["clusters_IDa"]] = case_when(
#  res05_seurat@meta.data$seurat_clusters %in% c('9','17','18')~ 'Forebrain progenitors',
#  res05_seurat@meta.data$seurat_clusters %in% c('1','3','5','7','8',
#                                                '12','14','16',
#                                                '20','22','23')~ 'LGE progenitors',
#  res05_seurat@meta.data$seurat_clusters %in% c('2')~ 'Migrating interneurons',
#  res05_seurat@meta.data$seurat_clusters %in% c('19')~ 'MGE interneurons',
#  res05_seurat@meta.data$seurat_clusters %in% c('0')~ 'Nascent MSNs',
#  res05_seurat@meta.data$seurat_clusters %in% c('11') ~ 'dMSN and iMSN precursors',
#  res05_seurat@meta.data$seurat_clusters %in% c('4','10') ~ 'dMSNs',
#  res05_seurat@meta.data$seurat_clusters %in% c('6','15') ~ 'iMSNs',
#  res05_seurat@meta.data$seurat_clusters %in% c('13')~ 'Ventral neocortical neurons',
#  res05_seurat@meta.data$seurat_clusters %in% c('21', '24')~ 'Endothelial cells'
#)

levels(res05_seurat) 
Idents(res05_seurat) <- "clusters_IDa"
levels(res05_seurat)

# Reorder the legend
res05_seurat$clusters_IDa = factor(res05_seurat$clusters_IDa,
                                   level = c('Forebrain progenitors',
                                             'LGE progenitors',
                                             'Nascent MSNs',
                                             'dMSN and iMSN precursors',
                                             'dMSNs',
                                             'iMSNs',
                                             'Migrating interneurons',
                                             'MGE interneurons',
                                             'Ventral neocortical neurons',
                                             'Endothelial cells'))

levels(res05_seurat)

SaveH5Seurat(res05_seurat, 
             filename = paste0("/scratch/scw1229/Emily/Integration_Bocchi"),
            overwrite = TRUE,
            verbose = TRUE)

#################################################################################
################# UMAP & alluvial plots with  cluter IDs  #######################
#################################################################################

setwd("/scratch/scw1229/Emily/Integration_Bocchi")

col = show_col(hue_pal()(10))

cluster_cols <- c('Forebrain progenitors' = '#F8766D',
                  'LGE progenitors' = '#DB8E00',
                  'Nascent MSNs' = '#00C1A7',
                  'dMSN and iMSN precursors' = '#00A6FF',
                  'dMSNs' = '#EF67EB',
                  'iMSNs' = '#B385FF',
                  'Migrating interneurons' = '#00BD5C',
                  'MGE interneurons' = '#AEA200',  
                  'Ventral neocortical neurons' = '#64B200',
                  'Endothelial cells' = '#00BADE')

umap3 = DimPlot(res05_seurat, reduction="umap", group.by = "clusters_IDa", raster = FALSE,
                label = FALSE, repel = TRUE, cols = cluster_cols)+ ggtitle(paste0("iPSC-Fetal LGE cluster IDs")) 
umap3

ggsave(umap3, filename = "UMAP_iPSCBocchi_res05_by_clusterIDa.tiff",width=10,height=6, units="in",
       dpi=700,device="tiff",limitsize = TRUE)

umap4 = DimPlot(res05_seurat, reduction="umap", group.by = "clusters_IDa", 
                split.by = "clusters_IDa", raster = FALSE,
                label = FALSE, repel = TRUE, cols = cluster_cols)+ ggtitle(paste0("iPSC-Fetal LGE cluster IDs"))+
  facet_wrap(~clusters_IDa, ncol=4) 
umap4

ggsave(umap4, filename = "UMAP_iPSCBocchi_res05_splitby_clusterIDa.tiff",width=14,height=8, units="in",
       dpi=700,device="tiff",limitsize = TRUE)

umap5 = DimPlot(res05_seurat, reduction="umap", group.by = "clusters_IDa", 
                split.by = "SampleType", raster = FALSE,
                label = FALSE, repel = TRUE, cols = cluster_cols)+ ggtitle(paste0("iPSC-Fetal LGE cluster IDs")) 
umap5

ggsave(umap5, filename = "UMAP_iPSCBocchi_res05_SampleType_clusterIDa.tiff",width=14,height=6, units="in",
       dpi=700,device="tiff",limitsize = TRUE)

#########################################################################
# Produce curve plots 
View(res05_seurat@meta.data)

map_summary_by.SampleType = res05_seurat@meta.data %>% group_by(SampleType,clusters_IDa) %>% dplyr::summarise(n())

map_summary_by.SampleType$SampleType = factor(x = map_summary_by.SampleType$SampleType, 
                                              levels = c('fetal_LGE',
                                                         'iPSC'))

map_summary_by.SampleType$clusters_IDa = factor(map_summary_by.SampleType$clusters_IDa, 
                                                level=c('Forebrain progenitors',
                                                        'LGE progenitors',
                                                        'Nascent MSNs',
                                                        'dMSN and iMSN precursors',
                                                        'dMSNs',
                                                        'iMSNs',
                                                        'Migrating interneurons',
                                                        'MGE interneurons',  
                                                        'Ventral neocortical neurons',
                                                        'Endothelial cells'))

ggplot(data = map_summary_by.SampleType,
       aes(axis1 = SampleType, axis2 = clusters_IDa, y = `n()`)) +
  geom_alluvium(aes(fill = clusters_IDa),
                curve_type = "sigmoid") +
  geom_stratum(aes(fill = clusters_IDa)) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)), size = 5) +
  theme_void()+
  guides(fill = guide_legend(title = "SampleType by Cluster_IDa"))+
  scale_fill_discrete(limits = c('Forebrain progenitors',
                                 'LGE progenitors',
                                 'Nascent MSNs',
                                 'dMSN and iMSN precursors',
                                 'dMSNs',
                                 'iMSNs',
                                 'Migrating interneurons',
                                 'MGE interneurons',  
                                 'Ventral neocortical neurons',
                                 'Endothelial cells')) +
  scale_fill_manual(
    breaks = c('Forebrain progenitors',
               'LGE progenitors',
               'Nascent MSNs',
               'dMSN and iMSN precursors',
               'dMSNs',
               'iMSNs',
               'Migrating interneurons',
               'MGE interneurons',  
               'Ventral neocortical neurons',
               'Endothelial cells'),
    
    values = c('Forebrain progenitors' = '#F8766D',
               'LGE progenitors' = '#DB8E00',
               'Nascent MSNs' = '#00C1A7',
               'dMSN and iMSN precursors' = '#00A6FF',
               'dMSNs' = '#EF67EB',
               'iMSNs' = '#B385FF',
               'Migrating interneurons' = '#00BD5C',
               'MGE interneurons' = '#AEA200',  
               'Ventral neocortical neurons' = '#64B200',
               'Endothelial cells' = '#00BADE'))

filename <- "CurvePlot_iPSCBocchi_SmapleType_to_clusterIDa.tiff"
ggsave(filename,
       plot = last_plot(),
       device = "tiff",
       ### path = "HighRes/",
       scale = 1,
       width = 28,
       height = 18,
       units = c("cm"),
       dpi = 700,
       limitsize = TRUE)

#################################################################################
#############  FindConservedMarkers between datasets by cluster IDs  ############
#################################################################################

res05_seurat = LoadH5Seurat("/scratch/scw1229/Emily/Integration_Bocchi/iPSCBocchi_res05_FULLintegrated_MF.h5Seurat")

setwd("/scratch/scw1229/Seurat_2023/Integration_Bocchi/Full_integration/analysis_sct/marker_analysis_MF/conserved_markers")

conserved.markers <- data.frame()
all.clusters <- levels(res05_seurat$clusters_IDa)

# loop through each cluster
for (i in all.clusters) {
  
  # print the cluster you're on
  print(i)
  
  # find conserved marker in cluster vs all other clusters
  markers <- FindConservedMarkers(res05_seurat,
                                  ident.1 = i, # subset to single cluster
                                  grouping.var = "SampleType", # compare by dataset
                                  only.pos = TRUE, # default
                                  min.pct = 0.25, logfc.threshold = 0.25,
                                  test.use = "MAST"
  )
  
  # skip if none
  if(nrow(markers) == 0) {
    next
  }
  
  # make rownames a column
  markers <- rownames_to_column(markers, var = "gene")
  
  # make cluster number a column
  markers$cluster <- i
  
  # add to final table
  conserved.markers <- rbind(conserved.markers, markers)
}

# create delta_pct
conserved.markers$iPSC_delta_pct <- abs(conserved.markers$iPSC_pct.1 - 
                                          conserved.markers$iPSC_pct.2)
conserved.markers$fetal_LGE_delta_pct <- abs(conserved.markers$fetal_LGE_pct.1 - 
                                               conserved.markers$fetal_LGE_pct.2)

# more stringent filtering
markers.strict <- conserved.markers[
  conserved.markers$iPSC_delta_pct > summary(conserved.markers$iPSC_delta_pct)[5],]
markers.strict <- conserved.markers[
  conserved.markers$fetal_LGE_delta_pct > summary(conserved.markers$fetal_LGE_delta_pct)[5],]
markers.strict$gene_name <- markers.strict$gene
markers.strict$row.num <- 1:nrow(markers.strict)

# compare 
#table(as.numeric(conserved.markers$cluster))
#table(as.numeric(markers.strict$cluster))

#save 

saveRDS(conserved.markers, "/scratch/scw1229/Emily/Integration_Bocchi/iPSCBocchi_conserved_markers_byClusterIDs.rds")
saveRDS(markers.strict, "/scratch/scw1229/Emily/Integration_Bocchi/iPSCBocchi_conserved_markers_byClusterIDs_strict.rds")


write.xlsx(conserved.markers, "/scratch/scw1229/Emily/Integration_Bocchi/iPSCBocchi_conserved_markers_byClusterIDs.xlsx")
write.xlsx(markers.strict, "/scratch/scw1229/Emily/Integration_Bocchi/iPSCBocchi_conserved_markers_byClusterIDs_strict.xlsx")

########################################################
####### Dotplot for Conserved Markers ##################
########################################################

setwd("/scratch/scw1229/Emily/Integration_Bocchi")

Idents(res05_seurat) <- "clusters_IDa"
levels(res05_seurat)
view(res05_seurat@meta.data)

# NEEDS TO BE FIXED AND SET ORDER CORRECTLY
Idents(res05_seurat) = factor(Idents(res05_seurat), 
                              levels = c('Forebrain progenitors',
                                         'LGE progenitors',
                                         'Nascent MSNs',
                                         'dMSN and iMSN precursors',
                                         'dMSNs',
                                         'iMSNs',
                                         'Migrating interneurons',
                                         'MGE interneurons',
                                         'Ventral neocortical neurons',
                                         'Endothelial cells'))

# where does this info come from
markers.to.plot <- c("CENPF","TOP2A","PAX6","FOXG1","KIF22",   ###  forebrain progenitors
                     "HIRIP3","FABP7","DLX2","GSX2","ASCL1",  ###  LGE progenitors
                     "DLX6-AS1","DCX","MEIS2","BCL11B","FOXP1","FOXP2",
                     "ZFHX3","TAC1","PCSK1N","AKAP9","SIX3","SP9", "PRRT2","SEZ6L2","TLCD3B",  ###  dMSN/iMSN, "NGRN" not present in integrated data slot
                     "PLS3","GAD1","GAD2","SCG3",   ###  Migrating interneurons & MGE interneurons
                     "NEUROD6","TBR1",   ###  Ventral neocortical neurons
                     "S100A11")   ###  Endothelial cells
###  16p genes - "BOLA2","MAZ","PPP4C" not present in integrated data slot

dot = DotPlot(res05_seurat, features = markers.to.plot, cols = c("#00BFC4", "#F8766D"),     ###   assay = "SCT", 
              dot.scale = 8, split.by = "SampleType") +
  RotatedAxis()

dot

ggsave(dot, filename = "dotplot_iPSCBocchi_conserved_markersDEG_by_clusterID_SampleType.tiff",
       width=14,height=7, units="in",dpi=700,device="tiff")

###########################################################################
#############  FindAllMarkers between datasets by cluster IDs  ############
###########################################################################

plan("future::multisession")
options(future.globals.maxSize= 5943718400)

setwd("/scratch/scw1229/Emily/Integration_Bocchi/DEG_by_Cluster_ID")


Idents(res05_seurat) <- "celltype.sampletype"
levels(res05_seurat)

view(res05_seurat)

###############   Forebrain progenitors   ################
## iPSC vs Bocchi 
temp <- FindMarkers(res05_seurat, ident.1 = "Forebrain progenitors_iPSC", 
                    ident.2 = "Forebrain progenitors_fetal_LGE", 
                    min.pct = 0.25, logfc.threshold = 0.25,
                    only.pos = TRUE,
                    verbose = FALSE)
temp$pctratio = temp$pct.1/temp$pct.2
temp$gene = row.names(temp)
temp = inner_join(temp,gene_match_list[,c(2:3)],by="gene") 
write.xlsx(temp,file=paste0("DEG_forebrainprogenitor_iPSCvsBocchi.xlsx"))

## Bocchi vs iPSC
temp <- FindMarkers(res05_seurat, ident.1 = "Forebrain progenitors_fetal_LGE", 
                    ident.2 = "Forebrain progenitors_iPSC", 
                    min.pct = 0.25, logfc.threshold = 0.25,
                    only.pos = TRUE,
                    verbose = FALSE)
temp$pctratio = temp$pct.1/temp$pct.2
temp$gene = row.names(temp)
temp = inner_join(temp,gene_match_list[,c(2:3)],by="gene")
write.xlsx(temp,file=paste0("DEG_forebrainprogenitor_BocchivsiPSC.xlsx"))
# done 

###############   LGE progenitors   ################
## iPSC vs Bocchi 
temp <- FindMarkers(res05_seurat, ident.1 = "LGE progenitors_iPSC", 
                    ident.2 = "LGE progenitors_fetal_LGE", 
                    min.pct = 0.25, logfc.threshold = 0.25,
                    only.pos = TRUE,
                    verbose = FALSE)
temp$pctratio = temp$pct.1/temp$pct.2
temp$gene = row.names(temp)
temp = inner_join(temp,gene_match_list[,c(2:3)],by="gene")
write.xlsx(temp,file=paste0("DEG_LGEprogenitor_iPSCvsBocchi.xlsx"))


## Bocchi vs iPSC
temp <- FindMarkers(res05_seurat, ident.1 = "LGE progenitors_fetal_LGE", 
                    ident.2 = "LGE progenitors_iPSC", 
                    min.pct = 0.25, logfc.threshold = 0.25,
                    only.pos = TRUE,
                    verbose = FALSE)
temp$pctratio = temp$pct.1/temp$pct.2
temp$gene = row.names(temp)
temp = inner_join(temp,gene_match_list[,c(2:3)],by="gene")
write.xlsx(temp,file=paste0("DEG_LGEprogenitor_BocchivsiPSC.xlsx"))
#done

############### Nascent MSNs   ################
## iPSC vs Bocchi 
temp <- FindMarkers(res05_seurat, ident.1 = "Nascent MSNs_iPSC", 
                    ident.2 = "Nascent MSNs_fetal_LGE", 
                    min.pct = 0.25, logfc.threshold = 0.25,
                    only.pos = TRUE,
                    verbose = FALSE)
temp$pctratio = temp$pct.1/temp$pct.2
temp$gene = row.names(temp)
temp = inner_join(temp,gene_match_list[,c(2:3)],by="gene")
write.xlsx(temp,file=paste0("DEG_nascentMSNs_iPSCvsBocchi.xlsx"))


## Bocchi vs iPSC
temp <- FindMarkers(res05_seurat, ident.1 = "Nascent MSNs_fetal_LGE", 
                    ident.2 = "Nascent MSNs_iPSC", 
                    min.pct = 0.25, logfc.threshold = 0.25,
                    only.pos = TRUE,
                    verbose = FALSE)
temp$pctratio = temp$pct.1/temp$pct.2
temp$gene = row.names(temp)
temp = inner_join(temp,gene_match_list[,c(2:3)],by="gene")
write.xlsx(temp,file=paste0("DEG_nascentMSNs_BocchivsiPSC.xlsx"))
# done

##################   dMSN and iMSN precursors   ################
## iPSC vs Bocchi
temp <- FindMarkers(res05_seurat, ident.1 = "dMSN and iMSN precursors_iPSC", 
                    ident.2 = "dMSN and iMSN precursors_fetal_LGE", 
                    min.pct = 0.25, logfc.threshold = 0.25,
                    only.pos = TRUE,
                    verbose = FALSE)
temp$pctratio = temp$pct.1/temp$pct.2
temp$gene = row.names(temp)
temp = inner_join(temp,gene_match_list[,c(2:3)],by="gene")
write.xlsx(temp,file=paste0("DEG_dMSNandiMSNprecursors_iPSCvsBocchi.xlsx"))
# done

## Bocchi vs iPSC
temp <- FindMarkers(res05_seurat, ident.1 = "dMSN and iMSN precursors_fetal_LGE", 
                    ident.2 = "dMSN and iMSN precursors_iPSC", 
                    min.pct = 0.25, logfc.threshold = 0.25,
                    only.pos = TRUE,
                    verbose = FALSE)
temp$pctratio = temp$pct.1/temp$pct.2
temp$gene = row.names(temp)
temp = inner_join(temp,gene_match_list[,c(2:3)],by="gene")
write.xlsx(temp,file=paste0("DEG_dMSNandiMSNprecursors_BocchivsiPSC.xlsx"))



##################   dMSNs   ################
## iPSC vs Bocchi
temp <- FindMarkers(res05_seurat, ident.1 = "dMSNs_iPSC", 
                    ident.2 = "dMSNs_fetal_LGE", 
                    min.pct = 0.25, logfc.threshold = 0.25,
                    only.pos = TRUE,
                    verbose = FALSE)
temp$pctratio = temp$pct.1/temp$pct.2
temp$gene = row.names(temp)
temp = inner_join(temp,gene_match_list[,c(2:3)],by="gene")
write.xlsx(temp,file=paste0("DEG_dMSNs_iPSCvsBocchi.xlsx"))


## Bocchi vs iPSC
temp <- FindMarkers(res05_seurat, ident.1 = "dMSNs_fetal_LGE", 
                    ident.2 = "dMSNs_iPSC", 
                    min.pct = 0.25, logfc.threshold = 0.25,
                    only.pos = TRUE,
                    verbose = FALSE)
temp$pctratio = temp$pct.1/temp$pct.2
temp$gene = row.names(temp)
temp = inner_join(temp,gene_match_list[,c(2:3)],by="gene")
write.xlsx(temp,file=paste0("DEG_dMSNs_BocchivsiPSC.xlsx"))


##################   iMSNs   ################
## iPSC vs Bocchi
temp <- FindMarkers(res05_seurat, ident.1 = "iMSNs_iPSC", 
                    ident.2 = "iMSNs_fetal_LGE", 
                    min.pct = 0.25, logfc.threshold = 0.25,
                    only.pos = TRUE,
                    verbose = FALSE)
temp$pctratio = temp$pct.1/temp$pct.2
temp$gene = row.names(temp)
temp = inner_join(temp,gene_match_list[,c(2:3)],by="gene")
write.xlsx(temp,file=paste0("DEG_iMSNs_iPSCvsBocchi.xlsx"))



## Bocchi vs iPSC
temp <- FindMarkers(res05_seurat, ident.1 = "iMSNs_fetal_LGE", 
                    ident.2 = "iMSNs_iPSC", 
                    min.pct = 0.25, logfc.threshold = 0.25,
                    only.pos = TRUE,
                    verbose = FALSE)
temp$pctratio = temp$pct.1/temp$pct.2
temp$gene = row.names(temp)
temp = inner_join(temp,gene_match_list[,c(2:3)],by="gene")
write.xlsx(temp,file=paste0("DEG_iMSNs_BocchivsiPSC.xlsx"))


##################  Migrating interneurons   ################
## iPSC vs Bocchi
temp <- FindMarkers(res05_seurat, ident.1 = "Migrating interneurons_iPSC", 
                    ident.2 = "Migrating interneurons_fetal_LGE", 
                    min.pct = 0.25, logfc.threshold = 0.25,
                    only.pos = TRUE,
                    verbose = FALSE)
temp$pctratio = temp$pct.1/temp$pct.2
temp$gene = row.names(temp)
temp = inner_join(temp,gene_match_list[,c(2:3)],by="gene")
write.xlsx(temp,file=paste0("DEG_Migratinginterneurons_iPSCvsBocchi.xlsx"))



## Bocchi vs iPSC
temp <- FindMarkers(res05_seurat, ident.1 = "Migrating interneurons_fetal_LGE", 
                    ident.2 = "Migrating interneurons_iPSC", 
                    min.pct = 0.25, logfc.threshold = 0.25,
                    only.pos = TRUE,
                    verbose = FALSE)
temp$pctratio = temp$pct.1/temp$pct.2
temp$gene = row.names(temp)
temp = inner_join(temp,gene_match_list[,c(2:3)],by="gene")
write.xlsx(temp,file=paste0("DEG_Migratinginterneurons_BocchivsiPSC.xlsx"))


##################   MGE interneurons   ################
## iPSC vs Bocchi
temp <- FindMarkers(res05_seurat, ident.1 = "MGE interneurons_iPSC", 
                    ident.2 = "MGE interneurons_fetal_LGE", 
                    min.pct = 0.25, logfc.threshold = 0.25,
                    only.pos = TRUE,
                    verbose = FALSE)
temp$pctratio = temp$pct.1/temp$pct.2
temp$gene = row.names(temp)
temp = inner_join(temp,gene_match_list[,c(2:3)],by="gene")
write.xlsx(temp,file=paste0("DEG_MGEinterneurons_iPSCvsBocchi.xlsx"))


## Bocchi vs iPSC
temp <- FindMarkers(res05_seurat, ident.1 = "MGE interneurons_fetal_LGE", 
                    ident.2 = "MGE interneurons_iPSC", 
                    min.pct = 0.25, logfc.threshold = 0.25,
                    only.pos = TRUE,
                    verbose = FALSE)
temp$pctratio = temp$pct.1/temp$pct.2
temp$gene = row.names(temp)
temp = inner_join(temp,gene_match_list[,c(2:3)],by="gene")
write.xlsx(temp,file=paste0("DEG_MGEinterneurons_BocchivsiPSC.xlsx"))


################   Ventral neocortical neurons   ################
## iPSC vs Bocchi
temp <- FindMarkers(res05_seurat, ident.1 = "Ventral neocortical neurons_iPSC", 
                    ident.2 = "Ventral neocortical neurons_fetal_LGE", 
                    min.pct = 0.25, logfc.threshold = 0.25,
                    only.pos = TRUE,
                    verbose = FALSE)
temp$pctratio = temp$pct.1/temp$pct.2
temp$gene = row.names(temp)
temp = inner_join(temp,gene_match_list[,c(2:3)],by="gene")
write.xlsx(temp,file=paste0("DEG_Ventralneocorticalneurons_iPSCvsBocchi.xlsx"))


## Bocchi vs iPSC
temp <- FindMarkers(res05_seurat, ident.1 = "Ventral neocortical neurons_fetal_LGE", 
                    ident.2 = "Ventral neocortical neurons_iPSC", 
                    min.pct = 0.25, logfc.threshold = 0.25,
                    only.pos = TRUE,
                    verbose = FALSE)
temp$pctratio = temp$pct.1/temp$pct.2
temp$gene = row.names(temp)
temp = inner_join(temp,gene_match_list[,c(2:3)],by="gene")
write.xlsx(temp,file=paste0("DEG_Ventralneocorticalneurons_BocchivsiPSC.xlsx"))


################  Endothelial cells   ################
## iPSC vs Bocchi
temp <- FindMarkers(res05_seurat, ident.1 = "Endothelial cells_iPSC", 
                    ident.2 = "Endothelial cells_fetal_LGE", 
                    min.pct = 0.25, logfc.threshold = 0.25,
                    only.pos = TRUE,
                    verbose = FALSE)
temp$pctratio = temp$pct.1/temp$pct.2
temp$gene = row.names(temp)
temp = inner_join(temp,gene_match_list[,c(2:3)],by="gene")
write.xlsx(temp,file=paste0("DEG_Endothelialcells_iPSCvsBocchi.xlsx"))



## Bocchi vs iPSC
temp <- FindMarkers(res05_seurat, ident.1 = "Endothelial cells_fetal_LGE", 
                    ident.2 = "Endothelial cells_iPSC", 
                    min.pct = 0.25, logfc.threshold = 0.25,
                    only.pos = TRUE,
                    verbose = FALSE)
temp$pctratio = temp$pct.1/temp$pct.2
temp$gene = row.names(temp)
temp = inner_join(temp,gene_match_list[,c(2:3)],by="gene")
write.xlsx(temp,file=paste0("DEG_Endothelialcells_BocchivsiPSC.xlsx"))


###########################################################################
############  FindAllMarkers Upregulated genes by cluster IDa  ############
###########################################################################

setwd("/scratch/scw1229/Emily/Integration_Bocchi/DEG_by_Cluster_ID")

View(res05_seurat@meta.data)
Idents(res05_seurat)= "clusters_IDa"
levels(res05_seurat)

### run FindALLMarkers for clusters_IDa  ###
### and plot signature genes as heatmap  ###

gene_match_list <- read.table("/scratch/scw1229/Seurat_analysis/gene_match_list.txt", 
                              quote = "", header = T, check.names = F, sep = "\t") 
names(gene_match_list)[2]<-"gene"

temp<-FindAllMarkers(res05_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
temp$pctratio<-temp$pct.1/temp$pct.2
temp<-inner_join(temp,gene_match_list[,c(2:3)],by="gene")
write.xlsx(temp,file=paste0("DEG_iPSCBocchi_UP_by_clusterIDa.xlsx"))

# dont overwrite temp
### use celltype.sampletype instead

Idents(res05_seurat)= "celltype.sampletype"
levels(res05_seurat)

temp<-FindAllMarkers(res05_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
temp$pctratio<-temp$pct.1/temp$pct.2
temp<-inner_join(temp,gene_match_list[,c(2:3)],by="gene")
write.xlsx(temp,file=paste0("DEG_iPSCBocchi_UP_by_clusterIDa_sampletype.xlsx"))


#################################################################################
########################  Heatmap with  cluster IDa  ############################
#################################################################################

marker_list_clustersIDa_up = read.xlsx("/scratch/scw1229/Emily/Integration_Bocchi/DEG_by_Cluster_ID/DEG_iPSCBocchi_UP_by_clusterIDa.xlsx",sheet=1)

setwd("/scratch/scw1229/Emily/Integration_Bocchi/Heatmaps")

##### DoHeatmap() generates an expression heatmap for given cells and features. 
##### In this case, we are plotting the top 10 markers (or all markers if less than 10) for each cluster.

############ TOP 10 GENES ############
marker_list_clustersIDa_up <- as.data.frame(read.xlsx("/scratch/scw1229/Emily/Integration_Bocchi/DEG_by_Cluster_ID/DEG_iPSCBocchi_UP_by_clusterIDa.xlsx", sheet = 1))

# Ensure that the column names are correct and there is no name conflict
colnames(marker_list_clustersIDa_up)

marker_list_clustersIDa_up %>%     #### file with DEGs in marker_list
 group_by(cluster) %>%     
  top_n(n = 10, wt = avg_log2FC) -> top10 #### see how this table is arranged and you can plot this way any table
 
 # Other way of getting top 10 avg log:
 # Select the top 10 genes by avg_log2FC for each cluster using dplyr functions explicitly
 # top10 <- dplyr::top_n(dplyr::group_by(marker_list_clustersIDa_up, cluster), n = 10, wt = avg_log2FC)

### genes chosen by higher FC but then sorted by p_adj
pheatmap = DoHeatmap(res05_seurat, assay = "SCT", cells = sample(Cells(res05_seurat)),
                     group.by = "clusters_IDa", features = top10$gene, 
                     angle = 60, size = 4.5,
                     draw.lines = TRUE)+ NoLegend() + 
  theme(axis.text.y = element_text(size = 6),
        plot.margin = unit(c(0.5,3,0.5,0.5), "cm"))
pheatmap

ggsave(pheatmap, filename = "Heatmap_iPSCBocchi_top10genes_UP_by_clustersIDa.tiff",width=16,height=14, 
       units="in",dpi=500,device="tiff")
#DONE
############ TOP 100 GENES ############
#res05_seurat_2 <- subset(res05_seurat, cluster !="Endothelial cells")

marker_list_clustersIDa_up %>%       #### file with DEGs in marker_list
 group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top100 #### see how this table is arranged and you can plot this way any table

pheatmap = DoHeatmap(res05_seurat, assay = "SCT", cells = sample(Cells(res05_seurat)),
                     group.by = "clusters_IDa", features = top100$gene, 
                     angle = 60, size = 4.5,
                     draw.lines = TRUE)+ NoLegend() + 
  theme(axis.text.y = element_text(size = 6),
        plot.margin = unit(c(0.5,3,0.5,0.5), "cm"))
pheatmap

ggsave(pheatmap, filename = "Heatmap_iPSCBocchi_top100genes_UP_by_clustersIDa.tiff",width=16,height=14, 
       units="in",dpi=500,device="tiff")

##################################################################################
########################## Heatmap by  cluster IDs  ##############################
##################################################################################

Idents(res05_seurat) = "clusters_IDa" 
levels(res05_seurat)
View(res05_seurat)

DefaultAssay(res05_seurat) = 'SCT'

###############   Forebrain progenitors   ################
marker_list_clustersIDa_up <- as.data.frame(read.xlsx("/scratch/scw1229/Emily/Integration_Bocchi/DEG_by_Cluster_ID/DEG_iPSCBocchi_UP_by_clusterIDa_sampletype.xlsx", sheet = 1))

forebrainIPS <- subset(marker_list_clustersIDa_up, cluster == "Forebrain progenitors_iPSC")
forebrainLGE <- subset(marker_list_clustersIDa_up, cluster == "Forebrain progenitors_fetal_LGE")

# Ensure that the column names are correct and there is no name conflict
colnames(marker_list_clustersIDa_up)

forebrainIPS %>%     #### file with DEGs in marker_list
  top_n(n = 50, wt = avg_log2FC) -> top50iPSC #### see how this table is arranged and you can plot this way any table
## subset by cluster ID

# Ensure that the column names are correct and there is no name conflict
colnames(marker_list_clustersIDa_up)

forebrainLGE %>%     #### file with DEGs in marker_list
  top_n(n = 50, wt = avg_log2FC) -> top50LGE #### see how this table is arranged and you can plot this way any table
## subset by cluster ID

top100 <- rbind(top50iPSC, top50LGE)

forebrain_progenitors_seurat <- subset(res05_seurat, clusters_IDa == "Forebrain progenitors")
view(res05_seurat)

# Heatmap
pheatmap = DoHeatmap(forebrain_progenitors_seurat, assay = "SCT", cells = sample(Cells(forebrain_progenitors_seurat)),
                     group.by = "SampleType", features = top100$gene, 
                     angle = 60, size = 4.5,
                     draw.lines = TRUE)+ NoLegend() + 
  theme(axis.text.y = element_text(size = 6),
        plot.margin = unit(c(0.5,3,0.5,0.5), "cm"))
pheatmap

ggsave(pheatmap, filename = "Heatmap_iPSCBocchi_top50ForebrainProgenitors.tiff",width=16,height=14, 
       units="in",dpi=500,device="tiff")


###############   LGE progenitors   ################
marker_list_clustersIDa_up <- as.data.frame(read.xlsx("/scratch/scw1229/Emily/Integration_Bocchi/DEG_by_Cluster_ID/DEG_iPSCBocchi_UP_by_clusterIDa_sampletype.xlsx", sheet = 1))

LGEprogiPSC <- subset(marker_list_clustersIDa_up, cluster == "LGE progenitors_iPSC")
LGEprogLGE <- subset(marker_list_clustersIDa_up, cluster == "LGE progenitors_fetal_LGE")

# Ensure that the column names are correct and there is no name conflict
colnames(marker_list_clustersIDa_up)

LGEprogiPSC %>%     #### file with DEGs in marker_list
  top_n(n = 50, wt = avg_log2FC) -> top50iPSC #### see how this table is arranged and you can plot this way any table
## subset by cluster ID

# Ensure that the column names are correct and there is no name conflict
colnames(marker_list_clustersIDa_up)

LGEprogLGE %>%     #### file with DEGs in marker_list
  top_n(n = 50, wt = avg_log2FC) -> top50LGE #### see how this table is arranged and you can plot this way any table


top100 <- rbind(top50iPSC, top50LGE)

LGE_progenitors_seurat <- subset(res05_seurat, clusters_IDa == "LGE progenitors")

# Heatmap
pheatmap = DoHeatmap(LGE_progenitors_seurat, assay = "SCT", cells = sample(Cells(LGE_progenitors_seurat)),
                     group.by = "SampleType", features = top100$gene, 
                     angle = 60, size = 4.5,
                     draw.lines = TRUE)+ NoLegend() + 
  theme(axis.text.y = element_text(size = 6),
        plot.margin = unit(c(0.5,3,0.5,0.5), "cm"))
pheatmap

ggsave(pheatmap, filename = "Heatmap_iPSCBocchi_top50LGEProgenitors.tiff",width=16,height=14, 
       units="in",dpi=500,device="tiff")

# Plot the heatmap for "LGE progenitors"
print(pheatmap_fb)



############### Nascent MSNs   ################
marker_list_clustersIDa_up <- as.data.frame(read.xlsx("/scratch/scw1229/Emily/Integration_Bocchi/DEG_by_Cluster_ID/DEG_iPSCBocchi_UP_by_clusterIDa_sampletype.xlsx", sheet = 1))

NascMSNiPSC <- subset(marker_list_clustersIDa_up, cluster == "Nascent MSNs_iPSC")
NascMSNLGE <- subset(marker_list_clustersIDa_up, cluster == "Nascent MSNs_fetal_LGE")

# Ensure that the column names are correct and there is no name conflict
colnames(marker_list_clustersIDa_up)

NascMSNiPSC %>%     #### file with DEGs in marker_list
  top_n(n = 50, wt = avg_log2FC) -> top50iPSC #### see how this table is arranged and you can plot this way any table
## subset by cluster ID

# Ensure that the column names are correct and there is no name conflict
colnames(marker_list_clustersIDa_up)

NascMSNLGE %>%     #### file with DEGs in marker_list
  top_n(n = 50, wt = avg_log2FC) -> top50LGE #### see how this table is arranged and you can plot this way any table


top100 <- rbind(top50iPSC, top50LGE)

Nascent_MSNs_seurat <- subset(res05_seurat, clusters_IDa == "Nascent MSNs")

# Heatmap
pheatmap = DoHeatmap(Nascent_MSNs_seurat, assay = "SCT", cells = sample(Cells(Nascent_MSNs_seurat)),
                     group.by = "SampleType", features = top100$gene, 
                     angle = 60, size = 4.5,
                     draw.lines = TRUE)+ NoLegend() + 
  theme(axis.text.y = element_text(size = 6),
        plot.margin = unit(c(0.5,3,0.5,0.5), "cm"))
pheatmap


ggsave(pheatmap, filename = "Heatmap_iPSCBocchi_top50NascentMSNs.tiff",width=16,height=14, 
       units="in",dpi=500,device="tiff")

##################   dMSN and iMSN precursors   ################
marker_list_clustersIDa_up <- as.data.frame(read.xlsx("/scratch/scw1229/Emily/Integration_Bocchi/DEG_by_Cluster_ID/DEG_iPSCBocchi_UP_by_clusterIDa_sampletype.xlsx", sheet = 1))

dMSNiMSNiPSC <- subset(marker_list_clustersIDa_up, cluster == "dMSN and iMSN precursors_iPSC")
dMSNiMSNLGE <- subset(marker_list_clustersIDa_up, cluster == "dMSN and iMSN precursors_LGE")

# Ensure that the column names are correct and there is no name conflict
colnames(marker_list_clustersIDa_up)

dMSNiMSNiPSC %>%     #### file with DEGs in marker_list
  top_n(n = 50, wt = avg_log2FC) -> top50iPSC #### see how this table is arranged and you can plot this way any table
## subset by cluster ID

# Ensure that the column names are correct and there is no name conflict
colnames(marker_list_clustersIDa_up)

dMSNiMSNLGE %>%     #### file with DEGs in marker_list
  top_n(n = 50, wt = avg_log2FC) -> top50LGE #### see how this table is arranged and you can plot this way any table


top100 <- rbind(top50iPSC, top50LGE)

dMSNiMSN_seurat <- subset(res05_seurat, clusters_IDa == "dMSN and iMSN precursors")

# Heatmap
pheatmap = DoHeatmap(dMSNiMSN_seurat, assay = "SCT", cells = sample(Cells(dMSNiMSN_seurat)),
                     group.by = "SampleType", features = top100$gene, 
                     angle = 60, size = 4.5,
                     draw.lines = TRUE)+ NoLegend() + 
  theme(axis.text.y = element_text(size = 6),
        plot.margin = unit(c(0.5,3,0.5,0.5), "cm"))
pheatmap


ggsave(pheatmap, filename = "Heatmap_iPSCBocchi_top50dMSNiMSN.tiff",width=16,height=14, 
       units="in",dpi=500,device="tiff")


##################   dMSNs   ################
marker_list_clustersIDa_up <- as.data.frame(read.xlsx("/scratch/scw1229/Emily/Integration_Bocchi/DEG_by_Cluster_ID/DEG_iPSCBocchi_UP_by_clusterIDa_sampletype.xlsx", sheet = 1))

dMSNsiPSC <- subset(marker_list_clustersIDa_up, cluster == "dMSNs_iPSC")
dMSNsLGE <- subset(marker_list_clustersIDa_up, cluster == "dMSNs_LGE")

# Ensure that the column names are correct and there is no name conflict
colnames(marker_list_clustersIDa_up)

dMSNsiPSC %>%     #### file with DEGs in marker_list
  top_n(n = 50, wt = avg_log2FC) -> top50iPSC #### see how this table is arranged and you can plot this way any table
## subset by cluster ID

# Ensure that the column names are correct and there is no name conflict
colnames(marker_list_clustersIDa_up)

dMSNsLGE %>%     #### file with DEGs in marker_list
  top_n(n = 50, wt = avg_log2FC) -> top50LGE #### see how this table is arranged and you can plot this way any table


top100 <- rbind(top50iPSC, top50LGE)

dMSNs_seurat <- subset(res05_seurat, clusters_IDa == "dMSNs")

# Heatmap
pheatmap = DoHeatmap(dMSNs_seurat, assay = "SCT", cells = sample(Cells(dMSNs_seurat)),
                     group.by = "SampleType", features = top100$gene, 
                     angle = 60, size = 4.5,
                     draw.lines = TRUE)+ NoLegend() + 
  theme(axis.text.y = element_text(size = 6),
        plot.margin = unit(c(0.5,3,0.5,0.5), "cm"))
pheatmap


ggsave(pheatmap, filename = "Heatmap_iPSCBocchi_top50dMSNs.tiff",width=16,height=14, 
       units="in",dpi=500,device="tiff")

##################   iMSNs   ################
marker_list_clustersIDa_up <- as.data.frame(read.xlsx("/scratch/scw1229/Emily/Integration_Bocchi/DEG_by_Cluster_ID/DEG_iPSCBocchi_UP_by_clusterIDa_sampletype.xlsx", sheet = 1))

iMSNsiPSC <- subset(marker_list_clustersIDa_up, cluster == "iMSNs_iPSC")
iMSNsLGE <- subset(marker_list_clustersIDa_up, cluster == "iMSNs_LGE")

# Ensure that the column names are correct and there is no name conflict
colnames(marker_list_clustersIDa_up)

iMSNsiPSC %>%     #### file with DEGs in marker_list
  top_n(n = 50, wt = avg_log2FC) -> top50iPSC #### see how this table is arranged and you can plot this way any table
## subset by cluster ID

# Ensure that the column names are correct and there is no name conflict
colnames(marker_list_clustersIDa_up)

iMSNsLGE %>%     #### file with DEGs in marker_list
  top_n(n = 50, wt = avg_log2FC) -> top50LGE #### see how this table is arranged and you can plot this way any table


top100 <- rbind(top50iPSC, top50LGE)

iMSNs_seurat <- subset(res05_seurat, clusters_IDa == "iMSNs")

# Heatmap
pheatmap = DoHeatmap(iMSNs_seurat, assay = "SCT", cells = sample(Cells(iMSNs_seurat)),
                     group.by = "SampleType", features = top100$gene, 
                     angle = 60, size = 4.5,
                     draw.lines = TRUE)+ NoLegend() + 
  theme(axis.text.y = element_text(size = 6),
        plot.margin = unit(c(0.5,3,0.5,0.5), "cm"))
pheatmap


ggsave(pheatmap, filename = "Heatmap_iPSCBocchi_top50iMSNs.tiff",width=16,height=14, 
       units="in",dpi=500,device="tiff")

##################  Migrating interneurons   ################ NOT DONE YET
marker_list_clustersIDa_up <- as.data.frame(read.xlsx("/scratch/scw1229/Emily/Integration_Bocchi/DEG_by_Cluster_ID/DEG_iPSCBocchi_UP_by_clusterIDa_sampletype.xlsx", sheet = 1))

MigIntiPSC <- subset(marker_list_clustersIDa_up, cluster == "Migrating Interneurons_iPSC")
MigIntLGE <- subset(marker_list_clustersIDa_up, cluster == "Migrating Interneurons_LGE")

# Ensure that the column names are correct and there is no name conflict
colnames(marker_list_clustersIDa_up)

MigIntiPSC %>%     #### file with DEGs in marker_list
  top_n(n = 50, wt = avg_log2FC) -> top50iPSC #### see how this table is arranged and you can plot this way any table
## subset by cluster ID

# Ensure that the column names are correct and there is no name conflict
colnames(marker_list_clustersIDa_up)

MigIntLGE %>%     #### file with DEGs in marker_list
  top_n(n = 50, wt = avg_log2FC) -> top50LGE #### see how this table is arranged and you can plot this way any table


top100 <- rbind(top50iPSC, top50LGE)

Migrating_interneurons_seurat <- subset(res05_seurat, clusters_IDa == "Migrating interneurons")

# Heatmap
pheatmap = DoHeatmap(Migrating_interneurons_seurat, assay = "SCT", cells = sample(Cells(Migrating_interneurons_seurat)),
                     group.by = "SampleType", features = top100$gene, 
                     angle = 60, size = 4.5,
                     draw.lines = TRUE)+ NoLegend() + 
  theme(axis.text.y = element_text(size = 6),
        plot.margin = unit(c(0.5,3,0.5,0.5), "cm"))
pheatmap #Error here, try again later


ggsave(pheatmap, filename = "Heatmap_iPSCBocchi_top50MigInt.tiff",width=16,height=14, 
       units="in",dpi=500,device="tiff")


##################   MGE interneurons   ################
marker_list_clustersIDa_up <- as.data.frame(read.xlsx("/scratch/scw1229/Emily/Integration_Bocchi/DEG_by_Cluster_ID/DEG_iPSCBocchi_UP_by_clusterIDa_sampletype.xlsx", sheet = 1))

MGEIntiPSC <- subset(marker_list_clustersIDa_up, cluster == "MGE interneurons_iPSC")
MGEIntLGE <- subset(marker_list_clustersIDa_up, cluster == "MGE interneurons_fetal_LGE")

# Ensure that the column names are correct and there is no name conflict
colnames(marker_list_clustersIDa_up)

MGEIntiPSC %>%     #### file with DEGs in marker_list
  top_n(n = 50, wt = avg_log2FC) -> top50iPSC #### see how this table is arranged and you can plot this way any table
## subset by cluster ID

# Ensure that the column names are correct and there is no name conflict
colnames(marker_list_clustersIDa_up)

MGEIntLGE %>%     #### file with DEGs in marker_list
  top_n(n = 50, wt = avg_log2FC) -> top50LGE #### see how this table is arranged and you can plot this way any table


top100 <- rbind(top50iPSC, top50LGE)

MGE_interneurons_seurat <- subset(res05_seurat, clusters_IDa == "MGE interneurons")

# Heatmap
pheatmap = DoHeatmap(MGE_interneurons_seurat, assay = "SCT", cells = sample(Cells(MGE_interneurons_seurat)),
                     group.by = "SampleType", features = top100$gene, 
                     angle = 60, size = 4.5,
                     draw.lines = TRUE)+ NoLegend() + 
  theme(axis.text.y = element_text(size = 6),
        plot.margin = unit(c(0.5,3,0.5,0.5), "cm"))
pheatmap


ggsave(pheatmap, filename = "Heatmap_iPSCBocchi_top50MGEinterneurons.tiff",width=16,height=14, 
       units="in",dpi=500,device="tiff")

################   Ventral neocortical neurons   ################ also not done, no observations in VentNeo - check metadata
marker_list_clustersIDa_up <- as.data.frame(read.xlsx("/scratch/scw1229/Emily/Integration_Bocchi/DEG_by_Cluster_ID/DEG_iPSCBocchi_UP_by_clusterIDa_sampletype.xlsx", sheet = 1))

VentNeoiPSC <- subset(marker_list_clustersIDa_up, cluster == "Ventral neocortical neurons_iPSC")
VentNeoLGE <- subset(marker_list_clustersIDa_up, cluster == "Ventral neocortical neurons_LGE") #no observations?

# Ensure that the column names are correct and there is no name conflict
colnames(marker_list_clustersIDa_up)

VentNeoiPSC %>%     #### file with DEGs in marker_list
  top_n(n = 50, wt = avg_log2FC) -> top50iPSC #### see how this table is arranged and you can plot this way any table
## subset by cluster ID

# Ensure that the column names are correct and there is no name conflict
colnames(marker_list_clustersIDa_up)

VentNeoLGE %>%     #### file with DEGs in marker_list
  top_n(n = 50, wt = avg_log2FC) -> top50LGE #### see how this table is arranged and you can plot this way any table


top100 <- rbind(top50iPSC, top50LGE)

Ventral_neocortical_neurons_seurat <- subset(res05_seurat, clusters_IDa == "Ventral neocortical neurons")

# Heatmap
pheatmap = DoHeatmap(Ventral_neocortical_neurons_seurat, assay = "SCT", cells = sample(Cells(Ventral_neocortical_neurons_seurat)),
                     group.by = "SampleType", features = top100$gene, 
                     angle = 60, size = 4.5,
                     draw.lines = TRUE)+ NoLegend() + 
  theme(axis.text.y = element_text(size = 6),
        plot.margin = unit(c(0.5,3,0.5,0.5), "cm"))
pheatmap


ggsave(pheatmap, filename = "Heatmap_iPSCBocchi_top50Ventralneocortical.tiff",width=16,height=14, 
       units="in",dpi=500,device="tiff")

################  Endothelial cells   ################
marker_list_clustersIDa_up <- as.data.frame(read.xlsx("/scratch/scw1229/Emily/Integration_Bocchi/DEG_by_Cluster_ID/DEG_iPSCBocchi_UP_by_clusterIDa_sampletype.xlsx", sheet = 1))

EndoiPSC <- subset(marker_list_clustersIDa_up, cluster == "Endothelial cells_iPSC")
EndoLGE <- subset(marker_list_clustersIDa_up, cluster == "Endothelial cells_LGE")

# Ensure that the column names are correct and there is no name conflict
colnames(marker_list_clustersIDa_up)

EndoiPSC %>%     #### file with DEGs in marker_list
  top_n(n = 50, wt = avg_log2FC) -> top50iPSC #### see how this table is arranged and you can plot this way any table
## subset by cluster ID

# Ensure that the column names are correct and there is no name conflict
colnames(marker_list_clustersIDa_up)

EndoLGE %>%     #### file with DEGs in marker_list
  top_n(n = 50, wt = avg_log2FC) -> top50LGE #### see how this table is arranged and you can plot this way any table


top100 <- rbind(top50iPSC, top50LGE)

Endothelial_cells_seurat <- subset(res05_seurat, clusters_IDa == "Endothelial cells")

# Heatmap
pheatmap = DoHeatmap(Endothelial_cells_seurat, assay = "SCT", cells = sample(Cells(Endothelial_cells_seurat)),
                     group.by = "SampleType", features = top100$gene, 
                     angle = 60, size = 4.5,
                     draw.lines = TRUE)+ NoLegend() + 
  theme(axis.text.y = element_text(size = 6),
        plot.margin = unit(c(0.5,3,0.5,0.5), "cm"))
pheatmap


ggsave(pheatmap, filename = "Heatmap_iPSCBocchi_top50Endo.tiff",width=16,height=14, 
       units="in",dpi=500,device="tiff")