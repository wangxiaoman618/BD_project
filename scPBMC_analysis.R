
# Part1: integrate 10X data

library(Seurat)
library(tidyverse)
library(ggpubr)

setwd("/home/hp/Documents/projects/BD/data/scRNA_seq_data/")
BD_PBMC.list<- list()
filename<- list.dirs("./",recursive = F)
for (i in filename) {
  tmp<-Read10X(i)
  print("finish read10x")
  tmp.obj<- CreateSeuratObject(counts = tmp, project =basename(i),
                               min.cells = 5,min.features = 200 )
  #remove cells with low quality
  tmp.obj[["percent.mt"]]<- PercentageFeatureSet(tmp.obj,pattern="^MT")
  
  tmp.obj<- subset(tmp.obj, 
                   subset= nFeature_RNA>200 & nFeature_RNA<4000 & 
                     percent.mt<10 )
  #normalization
  tmp.obj <- NormalizeData(tmp.obj)
  tmp.obj <- FindVariableFeatures(tmp.obj,selection.method = "vst")
  BD_PBMC.list[[basename(i)]]<- tmp.obj
}

BD_PBMC.anchors<- FindIntegrationAnchors(BD_PBMC.list,dims = 1:30)
BD_PBMC.comb <- IntegrateData(anchorset = BD_PBMC.anchors, dims = 1:30)
DefaultAssay(BD_PBMC.comb)<- "integrated"
saveRDS(BD_PBMC.bj, file= "./BD_PBMC.bj.rds")
