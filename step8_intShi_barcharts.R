x<- res05_seurat@meta.data %>% group_by(clusters_IDa,SampleType) %>% dplyr::summarize(N=n())

x2<- res05_seurat@meta.data %>% group_by(clusters_IDa) %>% dplyr::summarize(Total=n())

SampleType.summary = inner_join(x, x2, by = "clusters_IDa")

write.xlsx(SampleType.summary, file = paste("SampleType_count_by_clusters_IDaSHORT.xlsx"), rowNames = FALSE, colNames = TRUE, showNA = FALSE)



## SampleType - Percentage bar graph

p<-ggplot(data=SampleType.summary,
          
          aes(x=clusters_IDa,y=N/Total*100,fill=SampleType))+
  
  geom_bar(stat="identity")+
  
  theme_bw()+
  
  labs(y="Percentage of cells",fill="SampleType")+
  
  theme(legend.position = "bottom",
        
        text=element_text(size=8))

p

ggsave(plot=p,filename = "SampleType_count%_by_clusters_IDa_iPSC-ShiSHORT.tiff",
       
       width=40,height=8,units="cm",dpi=600,device="tiff",limitsize = TRUE)



## SampleID - Count bar graph

p<-ggplot(data=SampleType.summary,
          
          aes(x= clusters_IDa,y=N,fill=SampleType))+
  
  geom_bar(stat="identity")+
  
  theme_bw()+
  
  labs(y="Number of cells",fill="SampleType")+
  
  theme(legend.position = "bottom",
        
        text=element_text(size=8))

p

ggsave(plot=p,filename = "SampleType_count_by_clusters_IDa_iPSC-ShiSHORT.tiff",
       
       width=40,height=8,units="cm",dpi=600,device="tiff",limitsize = TRUE)