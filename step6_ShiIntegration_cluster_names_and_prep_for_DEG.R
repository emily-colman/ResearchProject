
### you have to consult xlsx output file in ClusterIDs folder *_IDenrichment_clusters.xlsx for FULL and then SHORT integration analysis
### in a new Excel table have all cluster numbers from 0 to 23 (FULL) or 21 (SHORT) and next to it the name of the cluster
### keep the same names as below for cluster groups
### all clusters identified as 'AP1-3' are Forebrain progenitors
### all clusters identified as 'BP1-3' are LGE progenitors
### all clusters identified as 'BP1-3' are LGE progenitors
### all clusters identified as 'Migr.int.' are Migrating interneurons
### all clusters identified as 'pre.MSNs_gene' from Bocchi and 'post-mitotic' from Shi are Nascent MSNs
### all clusters identified as 'D1.MSNs.mat._gene' or 'D1.MSNs.imm._gene' are dMSNs
### all clusters identified as 'D2.MSNs_gene' are iMSNs
### if any clusters are identified as either 'D1.MSNs' and 'D2.MSNs_gene' simultaneously then they are dMSN and iMSN precursors
### all clusters identified as 'D2.MSNs_gene' are iMSNs
### all clusters identified as 'v.Cx_gene' are Ventral neocortical neurons
### all clusters identified as 'Endo_gene' are Endothelial cells
### all clusters identified as 'v.CGEint._gene' are CGE interneurons


#################################################################################
################# Full integration #################################################
#####################################################################################
### creating new variable for clusters_IDa = to group individual clusters into groups

### open your seurat object
res05_seurat = LoadH5Seurat("/scratch/scw1229/Emily/Integration_Shi_FULL/MSNiPSC.Shi_FULLintegrated.h5Seurat")
view(res05_seurat@meta.data)

### update below with relevant cluster numbers next to cluster group names (info below is from Bocchi script) - do with shi full and short
### add new groups as required, e.g. CGE interneurons (they do not always show up)

res05_seurat[["clusters_IDa"]] = case_when(
  res05_seurat@meta.data$seurat_clusters %in% c('17','4','11')~ 'Forebrain progenitors', 
  res05_seurat@meta.data$seurat_clusters %in% c('2','5','7','15','16')~ 'LGE progenitors',
  res05_seurat@meta.data$seurat_clusters %in% c('3','8')~ 'Nascent MSNs',
  res05_seurat@meta.data$seurat_clusters %in% c('9')~ 'dMSN and iMSN precursors',
  res05_seurat@meta.data$seurat_clusters %in% c('14')~ 'dMSNs',
  res05_seurat@meta.data$seurat_clusters %in% c('12')~ 'iMSNs',
  res05_seurat@meta.data$seurat_clusters %in% c('0','1')~ 'MGE interneurons',
  res05_seurat@meta.data$seurat_clusters %in% c('10')~ 'CGE interneurons',
  res05_seurat@meta.data$seurat_clusters %in% c('13','6','20')~ 'Excitatory neurons',
  res05_seurat@meta.data$seurat_clusters %in% c('21','22','23')~ 'Endothelial cells',
  res05_seurat@meta.data$seurat_clusters %in% c('18')~ 'Microglia',
  res05_seurat@meta.data$seurat_clusters %in% c('19')~ 'OPC'
)


### Changes to clusters (generally Bocchi have more diversity in MSN mature clusters so use this paper first please): 
### I moved cluster 1 to MGE interneurons from CGE as overlap with Bocchi MGE is larger than Shi CGE 
### cluster 10 is nascent MSNs according to Shi but overlaps more closely with CGE interneurons (most significantly), dMSNs and cortical neurons from Bocchi
### cluster 12 is LGE identity from Shi but more significantly D2 MSN or iMSN from Bocchi.
### cluster 14 is LGE identity from Shi but more significantly D1 MSN or dMSN from Bocchi.
### clusters 17,4 AP1-3 from Bocchi so Forebrain progenitors
### clusters 2,5,7 BP1-3 from Bocchi so LGE progenitors
### for simplistic purposes Cx, Thalamic and Excitatory neuron clusters 13 & 6 can just be named Excitatory neurons
### cluster 11 based on gene expression is RGC which is the same as Forebrain progenitors
### clusters 15 & 16 do not appear in overlap tables but based on DEG table gene expression they are LGE progenitors 

### I have ordered the cluster groups below in the correct order to use for plots

levels(res05_seurat) 
Idents(res05_seurat) <- "clusters_IDa"
levels(res05_seurat)

### Reorder the legend
res05_seurat$clusters_IDa = factor(res05_seurat$clusters_IDa,
                                   level = c('Forebrain progenitors',
                                             'LGE progenitors',
                                             'Nascent MSNs',
                                             'dMSN and iMSN precursors',
                                             'dMSNs',
                                             'iMSNs',
                                             'MGE interneurons',
                                             'CGE interneurons',
                                             'Excitatory neurons',
                                             'Endothelial cells',
                                             'Microglia',
                                             'OPC'))

levels(res05_seurat)

### re-save the seurat object to update the old one
SaveH5Seurat(res05_seurat, 
             filename = paste0("/scratch/scw1229/Emily/Integration_Shi_FULL/MSNiPSC.Shi_FULLintegrated.h5Seurat"),
             overwrite = TRUE,
             verbose = TRUE)

### then repeat the steps for SHORT

### next you can run the script from step9_MSN_ALL_Fullintegration_Bocchi_Emily.R
### to perform DEG analysis

### you will have to update colour number and values for MAP & alluvial plots
### if you have a different number of cluster groups than in Bocchi analysis

col = show_col(hue_pal()(12)) ### put number of cluster groups here

### update below with new cluster group names and values

cluster_cols <- c('Forebrain progenitors' = '#F8766D',
                  'LGE progenitors' = '#DB8E00',
                  'Nascent MSNs' = '#B79f00',
                  'dMSN and iMSN precursors' = '#00BA38',
                  'dMSNs' = '#00c08B',
                  'iMSNs' = '#00BFC4',
                  'MGE interneurons' = '#00B4F0',
                  'CGE interneurons' = '#619CFF', 
                  'Excitatory neurons' = '#C77CFF',
                  'Endothelial cells' = '#F564E3',
                  'Microglia' = '#FF64B0',
                  'OPC' = '#00BADE')

### you can follow the Bocchi script next to plot graphs and perform DEG analysis and then plot Heatmaps


##########################################################################################
##################### Short integration ###############################################
##########################################################################################

### creating new variable for clusters_IDa = to group individual clusters into groups

# Short
res05_seurat = LoadH5Seurat("/scratch/scw1229/Emily/Integration_Shi_SHORT/MSNiPSC.Shi_SHORTintegrated.h5Seurat")
view(res05_seurat@meta.data)

### all updated
res05_seurat[["clusters_IDa"]] = case_when(
  res05_seurat@meta.data$seurat_clusters %in% c('10','2','7','9') ~ 'Forebrain progenitors',
  res05_seurat@meta.data$seurat_clusters %in% c('12','14','4') ~ 'LGE progenitors',
  res05_seurat@meta.data$seurat_clusters %in% c('1','3','13')~ 'Nascent MSNs',
  res05_seurat@meta.data$seurat_clusters %in% c('0')~ 'dMSN and iMSN precursors',
  res05_seurat@meta.data$seurat_clusters %in% c('6','8')~ 'dMSNs',
  res05_seurat@meta.data$seurat_clusters %in% c('5')~ 'iMSNs',
  res05_seurat@meta.data$seurat_clusters %in% c('11','17')~ 'Migrating interneurons',
    res05_seurat@meta.data$seurat_clusters %in% c('18') ~ 'Excitatory neurons',
  res05_seurat@meta.data$seurat_clusters %in% c('19','20','21') ~ 'Endothelial cells',
  res05_seurat@meta.data$seurat_clusters %in% c('15','22')~ 'Microglia',
  res05_seurat@meta.data$seurat_clusters %in% c('16')~ 'OPC'
)

levels(res05_seurat) 
Idents(res05_seurat) <- "clusters_IDa"
levels(res05_seurat)

### Reorder the legend
res05_seurat$clusters_IDa = factor(res05_seurat$clusters_IDa,
                                   level = c('Forebrain progenitors',
                                             'LGE progenitors',
                                             'Nascent MSNs',
                                             'dMSN and iMSN precursors',
                                             'dMSNs',
                                             'iMSNs',
                                             'Migrating interneurons',
                                             'Excitatory neurons',
                                             'Endothelial cells',
                                             'Microglia',
                                             'OPC'))

levels(res05_seurat)

### re-save the seurat object to update the old one
SaveH5Seurat(res05_seurat, 
             filename = paste0("/scratch/scw1229/Emily/Integration_Shi_SHORT/MSNiPSC.Shi_SHORTintegrated.h5Seurat"),
             overwrite = TRUE,
             verbose = TRUE)

### then repeat the steps for SHORT

### next you can run the script from step9_MSN_ALL_Fullintegration_Bocchi_Emily.R
### to perform DEG analysis

### you will have to update colour number and values for MAP & alluvial plots
### if you have a different number of cluster groups than in Bocchi analysis

col = show_col(hue_pal()(11)) ### put number of cluster groups here

### update below with new cluster group names and values

cluster_cols <- c('Forebrain progenitors' = '#F8766D',
                  'LGE progenitors' = '#DB8E00',
                  'Nascent MSNs' = '#B79f00',
                  'dMSN and iMSN precursors' = '#00BA38',
                  'dMSNs' = '#00c08B',
                  'iMSNs' = '#00BFC4',
                  'Migrating interneurons' = '#00BD5C',
                  'Excitatory neurons' = '#C77CFF', 
                  'Endothelial cells' = '#F564E3',
                  'Microglia' = '#FF64B0',
                  'OPC' = '#00BADE')

### you can follow the Bocchi script next to plot graphs and perform DEG analysis and then plot Heatmaps



