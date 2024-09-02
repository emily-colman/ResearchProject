############################################################################
#### run similar to before for SHORT and FULL Shi integration
## all parts apart from trajectory analysis completed for both short and full integration
########################################################################################################################################
############################################### Clustering for Integration #############################################################
########################################################################################################################################

#Seurat_filtered_MTsct <- LoadH5Seurat("/scratch/scw1229/Emily/Integration_Shi_SHORT/MSNiPSC.Shi_SHORTintegrated.h5Seurat")
#View(Seurat_filtered_MTsct)
Seurat_filtered_MTsct <- LoadH5Seurat("/scratch/scw1229/Emily/Integration_Shi_FULL/MSNiPSC.Shi_FULLintegrated.h5Seurat")

### Define cluster resolution, we chose res = 0.3 for both SHORT and FULL Shi integration
### Also define dims
Seurat_filtered_MTsct<-FindNeighbors(Seurat_filtered_MTsct, dims = 1:50)
Seurat_filtered_MTsct<-FindClusters(Seurat_filtered_MTsct, verbose = FALSE, resolution = 0.3)

### save updated seurat object
SaveH5Seurat(Seurat_filtered_MTsct, 
             filename = paste0("/scratch/scw1229/Emily/Integration_Shi_FULL/MSNiPSC.Shi_FULLintegrated.h5Seurat"),
             overwrite = TRUE,
             verbose = TRUE)

setwd("/scratch/scw1229/Emily/Integration_Shi_FULL/analysis_sct")

### Define the number of colors
nb.cols <- 36
mycolors <- colorRampPalette(brewer.pal(11, "Spectral"))(nb.cols)

p<-DimPlot(Seurat_filtered_MTsct, reduction="umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)+scale_fill_manual(values = mycolors)+ggtitle(paste0("res=0.3"))+
  theme(text=element_text(size=8))+guides(fill=guide_legend(ncol=2))
p
ggsave(plot=p,filename = "UMAP_MSN_ALL_ShiFULL_by_SeuratClusters.tiff",
       width=8,height=6,units="in",dpi=700,device="tiff",limitsize = TRUE)


#########################################
############ save data.frames ###########
#########################################

setwd("/scratch/scw1229/Emily/Integration_Shi_FULL/analysis_sct")

if(!dir.exists(paste0("DataFrames")))
  dir.create(paste0("DataFrames"))

setwd("/scratch/scw1229/Emily/Integration_Shi_FULL/analysis_sct/DataFrames")

cluster03<-as.data.frame(Seurat_filtered_MTsct@meta.data)
write.xlsx(cluster03, file = paste("MSN_ShiFULL_cluster_res0.3.xlsx"), rowNames = TRUE, colNames = TRUE, showNA = FALSE)

rm(cluster03)

###########################################################################################################################################          
########################################## DEG - FindMarkers for seurat_cluster ###########################################################
###########################################################################################################################################     
setwd("/scratch/scw1229/Emily/Integration_Shi_FULL/analysis_sct/Clusters")

gene_match_list <- read.table("/scratch/scw1229/Seurat_analysis/gene_match_list.txt", quote = "", 
                              header = T, check.names = F, sep = "\t") 

names(gene_match_list)[2]<-"gene"

############# FindMarkers only.pos = TRUE  #############

## by clusters

levels(Seurat_filtered_MTsct) ## check that active levels are seurat_clusters
Idents(Seurat_filtered_MTsct) <- "seurat_clusters"  ## change to seurat_clusters
levels(Seurat_filtered_MTsct) ## check that active levels are seurat_clusters

temp<-FindAllMarkers(Seurat_filtered_MTsct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
temp$pctratio<-temp$pct.1/temp$pct.2
temp<-inner_join(temp,gene_match_list[,c(2:3)],by="gene")
write.xlsx(temp,file=paste0("DEG_UP_by_cluster_ShiFULL",".xlsx"))

################################################################################    
######################## Marker expressions in SCT #############################
################################################################################

if(!dir.exists(paste0("Marker_expression_by_cluster")))
  dir.create(paste0("Marker_expression_by_cluster"))

setwd("/scratch/scw1229/Emily/Integration_Shi_FULL/analysis_sct/Clusters/Marker_expression_by_cluster")

### Define the number of colors you want
nb.cols <- 36
mycolors <- colorRampPalette(brewer.pal(11, "Spectral"))(nb.cols)

marker_list_clusters<-read.xlsx("/scratch/scw1229/Seurat_analysis/KeyMarker_list_clusters.xlsx",sheet=1)

### Plot marker expression as violin plots
sapply(levels(factor(marker_list_clusters$CellClass)),FUN=function(x){
  temp<-marker_list_clusters$gene[marker_list_clusters$CellClass==x]
  p<-VlnPlot(Seurat_filtered_MTsct, features = unique(temp),pt.size = 0.2, ncol = 4) 
  #,assay="SCT",slot="data",log=TRUE)+
  #  ylab("Logarithmic Normalised Expression")+scale_fill_tableau()
  ggsave(plot=p,filename=paste0(x,".tiff"),
         width=36,height=ceiling(length(temp)/3)*7,units="cm",dpi=300,device="tiff",limitsize = FALSE) ### dpi=600 for publication
})

################################################################################
######################## Trajectories in Integration ###########################
################################################################################

setwd("/scratch/scw1229/Emily/Integration_Shi_SHORT/analysis_sct")

if(!dir.exists(paste0("trajectory")))
  dir.create(paste0("trajectory"))


setwd("/scratch/scw1229/Emily/Integration_Shi_SHORT/analysis_sct/trajectory")  
library(monocle3)
### Building trajectories with Monocle3 
Seurat_MTsct_batch.cds <- as.cell_data_set(Seurat_filtered_MTsct)

#### need to assign all cells to 1 partition
recreate.partition <- c(rep(1,length(Seurat_MTsct_batch.cds@colData@rownames)))
names(recreate.partition) <- Seurat_MTsct_batch.cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)
Seurat_MTsct_batch.cds@clusters$UMAP$partitions <- recreate.partition

Seurat_MTsct_batch.cds <- cluster_cells(cds = Seurat_MTsct_batch.cds, reduction_method = "UMAP")
Seurat_MTsct_batch.cds <- learn_graph(Seurat_MTsct_batch.cds, use_partition = FALSE)  #if forcing 1 partition, then FALSE

Seurat_MTsct_batch.cds <- order_cells(Seurat_MTsct_batch.cds, reduction_method = "UMAP")

### Plotting trajectories colored by pseudotime
monocle_trajectories_Seurat = plot_cells(cds = Seurat_MTsct_batch.cds,color_cells_by = "pseudotime", 
                                         show_trajectory_graph = TRUE)
monocle_trajectories_Seurat
ggsave(monocle_trajectories_Seurat, filename = "UMAP_ShiSHORT_trajectories.tiff",width=8,height=6, 
       units="in",dpi=700,device="tiff")

### Add Pseudotime info to Seurat metadata of selected resolution 
Seurat_filtered_MTsct <- AddMetaData(
  object = Seurat_filtered_MTsct,
  metadata = Seurat_MTsct_batch.cds@principal_graph_aux@listData$UMAP$pseudotime, 
  col.name = "Pseudotime")

View(Seurat_filtered_MTsct@meta.data)

UMAP_Seurat_pseudotime = FeaturePlot(Seurat_filtered_MTsct, c("Pseudotime"), pt.size = 0.4) & scale_color_viridis_c()
UMAP_Seurat_pseudotime
ggsave(UMAP_Seurat_pseudotime, filename = "UMAP_ShiSHORT_pseudotime.tiff",width=8,height=6, 
       units="in",dpi=700,device="tiff")

### save updated seurat object
SaveH5Seurat(Seurat_filtered_MTsct, 
             filename = paste0("/scratch/scw1229/Emily/Integration_Shi_SHORT/MSNiPSC.Shi_SHORTintegrated.h5Seurat"),
             overwrite = TRUE,
             verbose = TRUE)


###############################################
###############################################
#################################################################################
############################  Identify cluster IDs  #############################
#################################################################################

setwd("/scratch/scw1229/Emily/Integration_Shi_FULL/analysis_sct/ClusterIDs") 

Seurat_filtered_MTsct_markers <- read.xlsx("../Clusters/DEG_UP_by_cluster_ShiFULL.xlsx",sheet=1)

Seurat_filtered_MTsct_markers_list<-list(res0.3=Seurat_filtered_MTsct_markers)

# if needed to re-open the seurat object
Seurat_filtered_MTsct <- LoadH5Seurat("/scratch/scw1229/Emily/Integration_Shi_FULL/MSNiPSC.Shi_FULLintegrated.h5Seurat")
background=nrow(Seurat_filtered_MTsct) ### should return length of genenames vector to insert on line 261 

Bocchi_clusters <- read.xlsx("Bocchi_clusters.xlsx", sheet=1)
Bocchi<-list(AP1_gene=Bocchi_clusters[,1],
             AP2_gene=Bocchi_clusters[,2],
             AP3_gene=Bocchi_clusters[,3],
             BP1_gene=Bocchi_clusters[,4],
             BP2_gene=Bocchi_clusters[,5],
             BP3_gene=Bocchi_clusters[,6],
             pre.MSNs_gene=Bocchi_clusters[,7],
             D2.MSNs_gene=Bocchi_clusters[,8],
             D1.MSNs.imm._gene=Bocchi_clusters[,9],
             D1.MSNs.mat._gene=Bocchi_clusters[,10],
             MGE.int_gene=Bocchi_clusters[,11],
             Migr.int._gene=Bocchi_clusters[,12],
             v.CGEint._gene=Bocchi_clusters[,13],
             v.Cx_gene=Bocchi_clusters[,14],
             Endo_gene=Bocchi_clusters[,15])

Shi_clusters <- read.xlsx("Shi_scRNAseq_clusters.xlsx", sheet=1)
Shi<-list(MGE=Shi_clusters[,1],
          CGE=Shi_clusters[,2],
          LGE=Shi_clusters[,3],
          progenitor=Shi_clusters[,4],
          OPC=Shi_clusters[,5],
          Excitatory_neurons=Shi_clusters[,6],
          Excitatory_IPC=Shi_clusters[,7],
          Microglia=Shi_clusters[,8],
          Endothelial=Shi_clusters[,9],
          Thalamic_neurons=Shi_clusters[,10],
          post_mitotic_cells=Shi_clusters[,11],
          progenitor=Shi_clusters[,12],
          RGC=Shi_clusters[,13],
          IPC=Shi_clusters[,14],
          pC1=Shi_clusters[,15],
          pC2=Shi_clusters[,16],
          pC3=Shi_clusters[,17],
          pL1=Shi_clusters[,18],
          pL2=Shi_clusters[,19],
          pL3=Shi_clusters[,20],
          pM1=Shi_clusters[,21],
          pM2=Shi_clusters[,22],
          pM3=Shi_clusters[,23],
          pM4=Shi_clusters[,24],
          MGE_b1=Shi_clusters[,25],
          CLGE_b1=Shi_clusters[,26],
          LGE_b2=Shi_clusters[,27],
          CGE_b2=Shi_clusters[,28],
          LGE_cells_with_OB_potential=Shi_clusters[,29],
          LGE_cells_with_STR_potenital=Shi_clusters[,30])

genesets<-list(Bocchi=Bocchi,
               Shi=Shi)

FUN_Geneset_Overlap<-function(query,reference,background,prefix){
  output_dir<-paste0(prefix)
  if(!dir.exists(output_dir))
    dir.create(output_dir)
  lapply(1:length(query),FUN=function(i){
    temp<-hypeR(query[[i]],reference,test="hypergeometric",background=background, fdr=0.05)
    hyp_to_excel(temp, file_path=paste0(output_dir,"/",names(query)[i],".xlsx"))
    temp<-read.xlsx(xlsxFile=paste0(output_dir,"/",names(query)[i],".xlsx"),sheet=1)
    if (nrow(temp)>0)
      temp$cluster<-names(query)[i]
    write.xlsx(temp,file=paste0(output_dir,"/",names(query)[i],".xlsx"))
  })
  sig_overlap_summary<-as.data.frame(rbindlist(lapply(list.files(output_dir,full.names=TRUE,pattern=".xlsx"),FUN=function(x){
    temp<-read.xlsx(x)
    temp$pct_overlap<-temp$overlap/temp$signature*100
    temp}),fill=TRUE))
  sig_overlap_summary
}

#########################
seurat_sct_cluster_sublist_marker_IDenrichment<-lapply(1:length(Seurat_filtered_MTsct_markers_list),FUN=function(index){
  if(!dir.exists(paste0("DEG_after_clustering")))
    dir.create(paste0("DEG_after_clustering"))
  if(!dir.exists(paste0("DEG_after_clustering/ID_RawOutput")))
    dir.create(paste0("DEG_after_clustering/ID_RawOutput"))
  
  temp<-Seurat_filtered_MTsct_markers_list[[index]]
  DEG_sig<-subset(temp,p_val_adj<0.05,select=c(gene,cluster))
  DEG_sig_list<-lapply(levels(factor(DEG_sig$cluster)),FUN=function(x){DEG_sig$gene[DEG_sig$cluster==x]})
  names(DEG_sig_list)<-levels(factor(DEG_sig$cluster))
  DEG_sig_overlap_summary<-lapply(1:length(genesets),FUN=function(i){
    #print(genesets[i])
    DEG_sig_overlap_summary<-FUN_Geneset_Overlap(DEG_sig_list,genesets[[i]],3000, #### update this number for each dataset from nrow function
                                                 paste0("DEG_after_clustering/ID_RawOutput/",names(Seurat_filtered_MTsct_markers_list)[[index]]))
    write.xlsx(DEG_sig_overlap_summary,file=paste0("DEG_after_clustering/",names(Seurat_filtered_MTsct_markers_list)[[index]],"_Summary_Enrichment_",names(genesets)[i],".xlsx"))
    DEG_sig_overlap_summary$genesets<-rep(names(genesets)[i],nrow(DEG_sig_overlap_summary))
    DEG_sig_overlap_summary
  })
  DEG_sig_overlap_summary<-as.data.frame(rbindlist(DEG_sig_overlap_summary,use.names=TRUE,idcol=TRUE))
  temp1<-subset(DEG_sig_overlap_summary,fdr<0.05)
  temp1
})

names(seurat_sct_cluster_sublist_marker_IDenrichment)<-names(Seurat_filtered_MTsct_markers_list)

#### ClusterID enrichment summary ####

seurat_sct_cluster_sublist_marker_IDenrichment_summary<-as.data.frame(rbindlist(lapply(
  list.files(paste0("DEG_after_clustering"),pattern="_Summary_Enrichment_",full.names = TRUE),
  FUN=function(x){
    Geneset<-unlist(strsplit(last(unlist(strsplit(x,split="_",fixed=TRUE))),split=".",fixed=TRUE))[1]
    res<-last(unlist(strsplit(head(unlist(strsplit(x,split="_",fixed=TRUE))),split="/",fixed=TRUE))[4])
    temp<-read.xlsx(xlsxFile=x)
    temp$genesets<-Geneset
    temp$res<-res
    temp_list<-as.data.frame(rbindlist(lapply(levels(factor(temp$cluster)),FUN=function(y){
      temp_list_temp<-subset(temp,fdr<0.05 & cluster==y)
      temp_list_temp<-top_n(temp_list_temp,10,pct_overlap)
    }),use.names=TRUE))
    temp_list
  }),use.names=TRUE))

seurat_sct_cluster_sublist_marker_IDenrichment_summary$order_pct<-sapply(seurat_sct_cluster_sublist_marker_IDenrichment_summary$label,FUN=function(y){sum(seurat_sct_cluster_sublist_marker_IDenrichment_summary$pct_overlap[seurat_sct_cluster_sublist_marker_IDenrichment_summary$label==y])})

sapply(levels(factor(seurat_sct_cluster_sublist_marker_IDenrichment_summary$res)),FUN=function(x){
  if(!dir.exists(paste0("DEG_after_clustering")))
    dir.create(paste0("DEG_after_clustering"))
  sapply(levels(factor(seurat_sct_cluster_sublist_marker_IDenrichment_summary$genesets)),FUN=function(y){
    temp<-subset(seurat_sct_cluster_sublist_marker_IDenrichment_summary,res==x & genesets==y)
    p<-ggplot(data=temp,
              aes(x=cluster,y=reorder(label,order_pct),colour=pct_overlap,size=-log10(fdr)))+
      geom_point()+
      scale_colour_gradient2_tableau()+
      theme_bw()+
      theme(axis.text.y = element_text(size=12))+
      ylab("ID term")
    ggsave(plot=p,filename=paste0("DEG_after_clustering/Enrichment_",x,"_",y,".tiff"),
           width=(length(levels(factor(temp$cluster)))*0.5+30),height=(length(levels(factor(temp$label)))*0.5+3),
           units="cm",dpi=600,device="tiff",limitsize = TRUE)
  })
})

write.xlsx(seurat_sct_cluster_sublist_marker_IDenrichment_summary, file = paste("ShiFULL_IDenrichment_summary.xlsx"), rowNames = TRUE, colNames = TRUE, showNA = FALSE)
