

## Shi integration script part-2 17/07/2024

#### integrated seurat objects:
##   iPSC.Shi.sct
##   iPSC.Shi.sct.short

#### next steps include 
##  (i)    to run PCA function on different settings of dimensions: 10, 20, 30, 40, 50.
##  (ii)   choosing the right number of PCA 'ndims' for further analysis
##  (iii)  running cell cycle analysis
##  (iv)   running cluster analysis and choosing the right resolution 

#### the script below is focused on Bocchi integration pipeline will have to change file paths
#### to run this twice for FULL and SHORT Shi objects


Seurat_filtered_MTsct = LoadH5Seurat("/scratch/scw1229/Emily/Integration_Shi_FULL/MSNiPSC.Shi_FULLintegrated.h5Seurat")
# Seurat_filtered_MTsct = LoadH5Seurat("/scratch/scw1229/Emily/Integration_Shi_SHORT/MSNiPSC.Shi_SHORTintegrated.h5Seurat")



########################################################################################################################################  
########################################## Run PCA & UMAPs from Fully Integrated Dataset ###############################################
########################################################################################################################################

setwd("/scratch/scw1229/Emily/Integration_Shi_SHORT")
Seurat_filtered_MTsct = LoadH5Seurat("/scratch/scw1229/Emily/Integration_Shi_SHORT/MSNiPSC.Shi_SHORTintegrated.h5Seurat")


Seurat_filtered_MTsct = RunPCA(Seurat_filtered_MTsct)

## Plot the standard deviations of the PCA for easy identification of an elbow in the graph
Elbowplot_iPSC.Bocchi.sct = ElbowPlot(Seurat_filtered_MTsct,ndims=50) + ylim (0, 25)
Elbowplot_iPSC.Bocchi.sct
ggsave(Elbowplot_iPSC.Bocchi.sct, filename = "Elbowplot_iPSC_Shi_integratedSHORT_sct.tiff",
       width=8,height=6, units="in",dpi=600,device="tiff",limitsize = TRUE)


##### Test run which PCA dims for UMAP ######

if(!dir.exists(paste0("PCA_analysis")))
  dir.create(paste0("PCA_analysis"))

setwd("/scratch/scw1229/Emily/Integration_Shi_SHORT/PCA_analysis")

##### 10 PCA
PCA10 = RunUMAP(Seurat_filtered_MTsct,dims = 1:10)
UMAP_byID_PCA10 = DimPlot(PCA10,reduction="umap",group.by = "orig.ident")+ggtitle(paste0("10 PCA"))
UMAP_byID_PCA10
ggsave(UMAP_byID_PCA10, filename = "UMAP_byID_PCA10_iPSC-ShiSHORT.tiff",width=8,height=6, units="in",
       dpi=600,device="tiff",limitsize = TRUE)

##### 20 PCA
PCA20 = RunUMAP(Seurat_filtered_MTsct,dims = 1:20)
UMAP_byID_PCA20 = DimPlot(PCA20,reduction="umap",group.by = "orig.ident")+ggtitle(paste0("20 PCA"))
UMAP_byID_PCA20
ggsave(UMAP_byID_PCA20, filename = "UMAP_byID_PCA20_iPSC-ShiSHORT.tiff",width=8,height=6, units="in",
       dpi=600,device="tiff",limitsize = TRUE)

##### 30 PCA
PCA30 = RunUMAP(Seurat_filtered_MTsct,dims = 1:30)
UMAP_byID_PCA30 = DimPlot(PCA30,reduction="umap",group.by = "orig.ident")+ggtitle(paste0("30 PCA"))
UMAP_byID_PCA30
ggsave(UMAP_byID_PCA30, filename = "UMAP_byID_PCA30_iPSC-ShiSHORT.tiff",width=8,height=6, units="in",
       dpi=600,device="tiff",limitsize = TRUE)

##### 40 PCA
PCA40 = RunUMAP(Seurat_filtered_MTsct,dims = 1:40)
UMAP_byID_PCA40 = DimPlot(PCA40,reduction="umap",group.by = "orig.ident")+ggtitle(paste0("40 PCA"))
UMAP_byID_PCA40
ggsave(UMAP_byID_PCA40, filename = "UMAP_byID_PCA40_iPSC-ShiSHORT.tiff",width=8,height=6, units="in",
       dpi=600,device="tiff",limitsize = TRUE)

##### 50 PCA
PCA50 = RunUMAP(Seurat_filtered_MTsct,dims = 1:50)
UMAP_byID_PCA50 = DimPlot(PCA50,reduction="umap",group.by = "orig.ident")+ggtitle(paste0("50 PCA"))
UMAP_byID_PCA50
ggsave(UMAP_byID_PCA50, filename = "UMAP_byID_PCA50_iPSC-ShiSHORT.tiff",width=8,height=6, units="in",
       dpi=600,device="tiff",limitsize = TRUE)

rm(PCA10,PCA20,PCA30,PCA40,PCA50,
   UMAP_byID_PCA10, UMAP_byID_PCA20, UMAP_byID_PCA30, UMAP_byID_PCA40, UMAP_byID_PCA50)

View(Seurat_filtered_MTsct@meta.data)

#############################

setwd("/scratch/scw1229/Emily/Integration_Shi_SHORT")

## Select dims 1:50

Seurat_filtered_MTsct = RunUMAP(Seurat_filtered_MTsct,dims = 1:50)
UMAP_byID = DimPlot(Seurat_filtered_MTsct,reduction="umap",group.by = "orig.ident",raster=FALSE)+
  ggtitle(paste0("50 PCA orig.ident"))
UMAP_byID
ggsave(UMAP_byID, filename = "UMAP_by_orig.ident_PCA50_iPSC-ShiSHORT.tiff",width=8,height=6, units="in",
       dpi=600,device="tiff",limitsize = TRUE)

#Save seurat with new data
SaveH5Seurat(Seurat_filtered_MTsct, 
             filename = paste0("/scratch/scw1229/Emily/Integration_Shi_SHORT/MSNiPSC.Shi_SHORTintegrated.h5Seurat"),
             overwrite = TRUE,
             verbose = TRUE)

Seurat_filtered_MTsct = LoadH5Seurat("/scratch/scw1229/Emily/Integration_Shi_SHORT/MSNiPSC.Shi_SHORTintegrated.h5Seurat")

### UMAP by SampleType
UMAP_bySampleType = DimPlot(Seurat_filtered_MTsct,reduction="umap",group.by = "SampleType",raster=FALSE)+
  ggtitle(paste0("50 PCA SampleType"))
UMAP_bySampleType
ggsave(UMAP_bySampleType, filename = "UMAP_bySampleType_PCA50_iPSC-ShiSHORT.tiff",width=8,height=6, units="in",
       dpi=600,device="tiff",limitsize = TRUE)

# Split up by sample type
UMAP_bySampleType = DimPlot(Seurat_filtered_MTsct,reduction="umap",group.by = "SampleType",
                            split.by = "SampleType", raster=FALSE)+
  ggtitle(paste0("50 PCA SampleType"))
UMAP_bySampleType
ggsave(UMAP_bySampleType, filename = "UMAP_bySampleType_PCA50_iPSC-Shi_splitSHORT.tiff",width=8,height=6, units="in",
       dpi=600,device="tiff",limitsize = TRUE)

### UMAP by TimePoint - check Shi Seurat objects to see what levels are within TimePoint variable for this dataset and update values below


levels(Seurat_filtered_MTsct)   # checks active variable of seurat object, usually it is orig.ident
Idents(Seurat_filtered_MTsct) <- "TimePoint"  ## change to TimePoint
levels(Seurat_filtered_MTsct)   # check levels within TimePoint to update data for graph below
Idents(Seurat_filtered_MTsct) <- "orig.ident"  ## change back to what the active variable was before

# Reorder TimePoint
Seurat_filtered_MTsct$TimePoint = factor(x = Seurat_filtered_MTsct$TimePoint, 
                                         levels = c("D13","D18","D24","D40","9PCW","12PCW","13PCW","16PCW","18PCW"))


UMAP_byTimePoint = DimPlot(Seurat_filtered_MTsct,reduction="umap",group.by = "TimePoint",
                           raster=FALSE)+
  ggtitle(paste0("50 PCA TimePoint"))
UMAP_byTimePoint
ggsave(UMAP_byTimePoint, filename = "UMAP_byTimePoint_PCA50_iPSC-ShiSHORT.tiff",width=8,height=6, units="in",
       dpi=600,device="tiff",limitsize = TRUE)


UMAP_byTimePoint = DimPlot(Seurat_filtered_MTsct,reduction="umap",group.by = "TimePoint",
                           split.by = "TimePoint", raster=FALSE, combine = TRUE, ncol = 4)+
  ggtitle(paste0("50 PCA TimePoint"))
UMAP_byTimePoint
ggsave(UMAP_byTimePoint, filename = "UMAP_byTimePoint_PCA50_iPSC-Shi_splitSHORT.tiff",width=8,height=6, units="in",
       dpi=600,device="tiff",limitsize = TRUE)


#######################################################################################################
########################################## analysis_sct ###############################################
#######################################################################################################

setwd("/scratch/scw1229/Emily/Integration_Shi_FULL")

if(!dir.exists(paste0("analysis_sct")))
  dir.create(paste0("analysis_sct"))

setwd("/scratch/scw1229/Emily/Integration_Shi_FULL/analysis_sct")
Seurat_filtered_MTsct <- LoadH5Seurat("/scratch/scw1229/Emily/Integration_Shi_SHORT/MSNiPSC.Shi_SHORTintegrated.h5Seurat")


##########################################################################
####################### Cell Cycle Scoring for SCT #######################
##########################################################################

if(!dir.exists(paste0("CellCycle")))
  dir.create(paste0("CellCycle"))

setwd("/scratch/scw1229/Emily/Integration_Shi_FULL/analysis_sct/CellCycle")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

Seurat_filtered_MTsct = CellCycleScoring(Seurat_filtered_MTsct, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = FALSE)
phase_summary = Seurat_filtered_MTsct@meta.data %>% group_by(orig.ident,Phase) %>% dplyr::summarise(n()) %>% 
  dcast(orig.ident~Phase)
Idents(Seurat_filtered_MTsct) = Seurat_filtered_MTsct$orig.ident

View(Seurat_filtered_MTsct@meta.data)

### Visualize the distribution of cell cycle markers across
RidgePlot(Seurat_filtered_MTsct, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
filename <- "Proliferation_iPSC-ShiSHORT.tiff"
ggsave(plot = last_plot(), filename = filename,width=30,height=30, units="cm",dpi=600,device="tiff",limitsize = TRUE)

### Visualise the cell cycle

cell_cycle_sct <- RunPCA(Seurat_filtered_MTsct, features = c(s.genes, g2m.genes))
View(cell_cycle_sct@meta.data)

CellCycle_UMAP = DimPlot(cell_cycle_sct, reduction = 'umap', group.by = 'Phase', shuffle = TRUE)
CellCycle_UMAP
ggsave(CellCycle_UMAP, filename = "CellCycle_UMAP_iPSC-ShiSHORT.tiff",width=8,height=6, 
       units="in",dpi=600,device="tiff",limitsize = TRUE)

# ID
CellCycle_SplitByID = DimPlot(cell_cycle_sct, reduction = 'umap', group.by = 'Phase', shuffle = TRUE, 
                              split.by = 'orig.ident')+facet_wrap(~orig.ident, ncol=4)
CellCycle_SplitByID
ggsave(CellCycle_SplitByID, filename = "CellCycle_SplitBy_orig.ident_iPSC-ShiSHORT.tiff",width=30,height=27, 
       units="in",dpi=600,device="tiff",limitsize = TRUE)

# Sample type
CellCycle_SplitBySampleType = DimPlot(cell_cycle_sct, reduction = 'umap', group.by = 'Phase', shuffle = TRUE, 
                                      split.by = 'SampleType')
CellCycle_SplitBySampleType
ggsave(CellCycle_SplitBySampleType, filename = "CellCycle_SplitBySampleType_iPSC-ShiSHORT.tiff",width=24,height=6, 
       units="in",dpi=600,device="tiff",limitsize = TRUE)


####################################################################################################################
cell_cycle_sct$TimePoint = factor(x = cell_cycle_sct$TimePoint, 
                                  levels = c("D13","D18","D24","D40","9PCW","12PCW","13PCW","16PCW","18PCW"))

CellCycle_SplitByTimePoint = DimPlot(cell_cycle_sct, reduction = 'umap', group.by = 'Phase', shuffle = TRUE, 
                                     split.by = 'TimePoint', combine = TRUE, ncol = 4)
CellCycle_SplitByTimePoint
ggsave(CellCycle_SplitByTimePoint, filename = "CellCycle_SplitBy_TimePoint_iPSC-ShiFULL.tiff",width=24,height=6, 
       units="in",dpi=600,device="tiff",limitsize = TRUE)


########################################################################
############ Plotting Cell Cycle Summary for each SampleID #############
########################################################################

## SampleID ##
x<-cell_cycle_sct@meta.data %>% group_by(orig.ident,Phase) %>% dplyr::summarize(N=n())
x2<-cell_cycle_sct@meta.data %>% group_by(orig.ident) %>% dplyr::summarize(Total=n())
phase.summary = inner_join(x, x2, by = "orig.ident")
write.table(phase.summary, file = paste("CellCycle_count_by_SampleID_iPSC-ShiFULL.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")
write.xlsx(phase.summary, file = paste("CellCycle_count_by_SampleID_iPSC-ShiFULL.xlsx"), rowNames = FALSE, colNames = TRUE, showNA = FALSE)

## SampleID - Percentage bar graph
p<-ggplot(data=phase.summary,
          aes(x=orig.ident,y=N/Total*100,fill=Phase))+
  geom_bar(stat="identity")+
  scale_color_manual(values = c("G1" = "#F8766D", "G2M"="00BA38", "S"="#619CFF"))+
  theme_bw()+
  labs(y="Percentage of cells",fill="Cell cycle phase")+
  theme(legend.position = "bottom",
        text=element_text(size=8))
p
ggsave(plot=p,filename = "CellCycle_count%_by_SampleID_iPSC-ShiFULL.tiff",
       width=40,height=8,units="cm",dpi=600,device="tiff",limitsize = TRUE)

## SampleID - Count bar graph
p<-ggplot(data=phase.summary,
          aes(x=orig.ident,y=N,fill=Phase))+
  geom_bar(stat="identity")+
  scale_color_manual(values = c("G1" = "#F8766D", "G2M"="00BA38", "S"="#619CFF"))+
  theme_bw()+
  labs(y="Number of cells",fill="Cell cycle phase")+
  theme(legend.position = "bottom",
        text=element_text(size=8))
p
ggsave(plot=p,filename = "CellCycle_count_by_SampleID_iPSC-ShiFULL.tiff",
       width=40,height=8,units="cm",dpi=600,device="tiff",limitsize = TRUE)

## TimePoint ##

x<-cell_cycle_sct@meta.data %>% group_by(TimePoint,Phase) %>% dplyr::summarize(N=n())
x2<-cell_cycle_sct@meta.data %>% group_by(TimePoint) %>% dplyr::summarize(Total=n())
phase.summary = inner_join(x, x2, by = "TimePoint")
write.table(phase.summary, file = paste("CellCycle_count_by_TimePoint_iPSC-ShiFULL.txt", sep = ""), quote = FALSE, row.names = FALSE, sep = "\t")
write.xlsx(phase.summary, file = paste("CellCycle_count_by_TimePoint_iPSC-ShiFULL.xlsx"), rowNames = FALSE, colNames = TRUE, showNA = FALSE)

## TimePoint - Percentage bar graph
p<-ggplot(data=phase.summary,
          aes(x=TimePoint,y=N/Total*100,fill=Phase))+
  geom_bar(stat="identity")+
  scale_color_manual(values = c("G1" = "#F8766D", "G2M"="00BA38", "S"="#619CFF"))+
  theme_bw()+
  labs(y="Percentage of cells",fill="Cell cycle phase")+
  theme(legend.position = "bottom",
        text=element_text(size=8))
p
ggsave(plot=p,filename = "CellCycle_count%_by_TimePoint_iPSC-ShiFULL.tiff",
       width=40,height=8,units="cm",dpi=600,device="tiff",limitsize = TRUE)


SaveH5Seurat(Seurat_filtered_MTsct, 
             filename = paste0("/scratch/scw1229/Emily/Integration_Shi_SHORT/MSNiPSC.Shi_FULLintegrated.h5Seurat"),
             overwrite = TRUE,
             verbose = TRUE)



########################################################################################################################################   
############################################### Clustering for Integration #############################################################
########################################################################################################################################


### need to test different setting for clustering to choose appropriate parameters

Seurat_filtered_MTsct <- LoadH5Seurat("/scratch/scw1229/Emily/Integration_Shi_FULL/MSNiPSC.Shi_SHORTintegrated.h5Seurat")

setwd("/scratch/scw1229/Emily/Integration_Shi_FULL/analysis_sct")

if(!dir.exists(paste0("Clusters")))
  dir.create(paste0("Clusters"))

setwd("/scratch/scw1229/Emily/Integration_Shi_FULL/analysis_sct/Clusters")


### Define the number of colors you want
nb.cols <- 36
mycolors <- colorRampPalette(brewer.pal(11, "Spectral"))(nb.cols)

FUN_seurat_clustering<-function(seurat,prefix){
  lapply(seq(1,6)/10,FUN=function(x){
    print (paste("Finding clusters res=",x,sep=" "))
    seurat<-FindNeighbors(seurat, dims = 1:50)
    seurat<-FindClusters(seurat, verbose = FALSE, res=x)
    p<-DimPlot(seurat, reduction = "umap", group.by = "seurat_clusters", 
               label = TRUE, repel = TRUE)+
      scale_fill_manual(values = mycolors)+ggtitle(paste0("res=",x))+
      theme(text=element_text(size=8))+guides(fill=guide_legend(ncol=2))
    ggsave(plot=p,filename=paste0(prefix,"/res",x,".tiff"),width=8,height=6,
           units="in",dpi=600,device="tiff",limitsize = TRUE)
    seurat
  })
}

Seurat_filtered_MTsct_cluster_list = FUN_seurat_clustering(Seurat_filtered_MTsct,"Clustering_UMAP_MSN-Shi")
names(Seurat_filtered_MTsct_cluster_list)<-seq(1,6)/10



### Get the stats from different resolutions
clustering_summary<-as.data.frame(t(rbindlist(lapply(1:6,FUN=function(i){
  temp<-as.data.frame(t(data.frame(NumCell=summary(Seurat_filtered_MTsct_cluster_list[[i]]$seurat_clusters))))
}),fill=TRUE)))
names(clustering_summary)<-names(Seurat_filtered_MTsct_cluster_list)
write.table(clustering_summary, file = paste("MSN_Shi_int_clustering_summary.txt", sep = ""), quote = FALSE, row.names = TRUE, sep = "\t")
write.xlsx(clustering_summary, file = paste("MSN_Shi_int_clustering_summary.xlsx"), rowNames = TRUE, colNames = TRUE, showNA = FALSE)

##################################################################################  
####################### Plotting Cell Cycle Summary in SCT #######################
##################################################################################

### Get bar graphs of Cell Phase for each cluster for each res
sapply(1:4,FUN=function(i){
  if (!dir.exists(paste0("Phase_summary_by_res")))
    dir.create(paste0("Phase_summary_by_res"))
  temp<-Seurat_filtered_MTsct_cluster_list[[i]]
  x<-temp@meta.data %>% group_by(seurat_clusters,Phase) %>% dplyr::summarize(N=n())
  write.xlsx(x, file = paste0("Phase_summary_by_res/",names(Seurat_filtered_MTsct_cluster_list)[i],"_Phase_count_MSN_Shi_MTsct.xlsx"), rowNames = FALSE, colNames = TRUE, showNA = FALSE)
  p<-ggplot(data=x,
            aes(x=seurat_clusters,y=N,fill=Phase))+
    geom_bar(stat="identity")+
    scale_fill_tableau(palette = "Tableau 10")+
    theme_bw()+
    labs(y="Number of cells",fill="Cell cycle phase")+
    scale_color_manual(values = c("G1" = "#F8766D", "G2M"="00BA38", "S"="#619CFF"))+
    theme(legend.position = "bottom",
          text=element_text(size=8))
  ggsave(plot=p,filename=paste0("Phase_summary_by_res/",names(Seurat_filtered_MTsct_cluster_list)[i],"_Phase_MSN_Shi_MTsct.tiff"),
         width=13,height=8,units="cm",dpi=600,device="tiff",limitsize = TRUE)
})

################################################################################  
####################### Plotting Sample Summary in SCT #######################
################################################################################

### Get bar graphs of SampleID for each cluster for each res in sublist
sapply(1:6,FUN=function(i){
  if (!dir.exists(paste0("SampleID_summary_by_res")))
    dir.create(paste0("SampleID_summary_by_res"))
  temp<-Seurat_filtered_MTsct_cluster_list[[i]]
  x<-temp@meta.data %>% group_by(seurat_clusters,orig.ident) %>% dplyr::summarize(N=n())
  write.xlsx(x, file = paste0("SampleID_summary_by_res/",names(Seurat_filtered_MTsct_cluster_list)[i],"_SampleID_count_MSN_Shi.xlsx"), rowNames = FALSE, colNames = TRUE, showNA = FALSE)
  p<-ggplot(data=x,
            aes(x=seurat_clusters,y=N,fill=orig.ident))+
    geom_bar(stat="identity")+
    scale_fill_manual(values = mycolors)+
    theme_bw()+
    labs(y="Number of cells",fill="SampleID")+
    theme(legend.position = "bottom",
          text=element_text(size=8))
  ggsave(plot=p,filename=paste0("SampleID_summary_by_res/",names(Seurat_filtered_MTsct_cluster_list)[i],"_SampleID_MSN_Shi.tiff"),
         width=13,height=8,units="cm",dpi=300,device="tiff",limitsize = TRUE)
})


### Get bar graphs of SampleType for each cluster for each res in sublist
sapply(1:6,FUN=function(i){
  if (!dir.exists(paste0("SampleType_summary_by_res")))
    dir.create(paste0("SampleType_summary_by_res"))
  temp<-Seurat_filtered_MTsct_cluster_list[[i]]
  x<-temp@meta.data %>% group_by(seurat_clusters,SampleType) %>% dplyr::summarize(N=n())
  write.xlsx(x, file = paste0("SampleID_summary_by_res/",names(Seurat_filtered_MTsct_cluster_list)[i],"_SampleType_count_MSN_Shi.xlsx"), 
             rowNames = FALSE, colNames = TRUE, showNA = FALSE)
  p<-ggplot(data=x,
            aes(x=seurat_clusters,y=N,fill=SampleType))+
    geom_bar(stat="identity")+
    theme_bw()+
    labs(y="Number of cells",fill="SampleType")+
    theme(legend.position = "bottom",
          text=element_text(size=8))
  ggsave(plot=p,filename=paste0("SampleType_summary_by_res/",names(Seurat_filtered_MTsct_cluster_list)[i],"_SampleType_MSN_Shi.tiff"),
         width=13,height=8,units="cm",dpi=300,device="tiff",limitsize = TRUE)
})
