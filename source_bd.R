library(Seurat)
library(data.table)
library(tidyverse)
library(ggsci)
library(ggthemes)
library(clusterProfiler)
library(ggpubr)
library(ggplot2)
library(pheatmap)
library(cowplot)
library(ggsci)
library(RColorBrewer)
library(readxl)

enrichDEgene<- function(genedata = DEseq2_output,
                        species = c("Homo sapiens","Mus musculus"),
                        category=c("H","C1","C7"),subcategory = "CP:KEGG",
                        method=c("GO","KEGG","GSEA","enrichr"),
                        genelist_input = NULL, # when it is not null, enrich this gene list
                        direction =c("Up","Down")){
  if (species=="Homo sapiens") {
    OrgDb_value="org.Hs.eg.db"
    organism_value ="hsa"
  }else{
    if (species=="Mus musculus") {
      OrgDb_value="org.Mm.eg.db"
      organism_value ="mmu"
      require(org.Mm.eg.db)
    }else{
      stop("Wrong species! Either Homo sapiens or Mus musculus")
    }
    
  }
  
  #defined <- ls()
  passed.arg <- names(as.list(match.call())[-1])
  if (method %in% c("GSEA","enrichr")) {
    if (!("category" %in% passed.arg)) {
      stop("missing values for category")
    }
  }
  
  library(msigdbr)
  library(clusterProfiler)
  if (is.null(genelist_input)) {
    if (direction=="Up") {
      genelist = rownames(genedata[genedata$log2FoldChange>0,])
    } else{
      genelist = rownames(genedata[genedata$log2FoldChange<0,])
    }
  }else{
    genelist = genelist_input
  }
  
  
  if (method =="KEGG") {
    
    eg = bitr(genelist,fromType="SYMBOL", toType="ENTREZID", OrgDb=OrgDb_value)
    mkk <- enrichKEGG(gene = eg$ENTREZID,organism = organism_value)
  } 
  
  if (method=="GO") {
    mkk<-enrichGO(gene     = genelist,
                  OrgDb         = OrgDb_value,
                  keyType       = 'SYMBOL',
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
  } 
  
  if (method =="enrichr") {
    m_t2g <- msigdbr(species = species, category = category,subcategory = subcategory) %>% 
      dplyr::select(gs_name, gene_symbol)
    mkk<- enricher(genelist, TERM2GENE=m_t2g)
  }
  
  if (method =="GSEA") {
    
    geneList = setNames(genedata$log2FoldChange,rownames(genedata))
    geneList <- sort(geneList, decreasing = TRUE)
    m_t2g <- msigdbr(species = species, category = category,subcategory = subcategory) %>% 
      dplyr::select(gs_name, gene_symbol)
    mkk<- GSEA(geneList, TERM2GENE = m_t2g)
  }
  
  return(mkk)
}



#volcanplot
plotVolcano <- function(data = data,
                        x= "log2FoldChange",y="Padj",
                        fill = "Threshold", color = "Threshold",
                        color_palette = c("#8DA0CB","grey","#FC8D62"),
                        label="gene",label_data = subset_data,
                        alpha=0.5,plotLab = T,
                        xlab = "log2 Fold Change",ylab="-log10 Padj",
                        title="Volcano Plot of genes",xline=c(-1,1)){
  require(ggrepel)
  require(ggplot2)
  p=ggplot(data,aes_string(x=x,y=y,fill=fill,color=color))+
    geom_point(alpha=alpha)+
    scale_color_manual(values=color_palette)+
    labs(x=xlab,y=ylab,title=title)+theme_pubr()+
    geom_hline(yintercept=1.3, linetype="dashed", color = "grey")+
    geom_vline(xintercept=xline, linetype="dashed", color = "grey")
  if (plotLab==T) {
    p<- p+geom_text_repel(data=label_data, aes_string(label=label))
  }
  plot(p)
}

makeMarkers.list <- function(degene){
  require(tidyverse)
  require(clusterProfiler)
  Degenelist=list()
  
  for (cells in unique(degene$cluster)) {
    
    tmp.gene<-  degene %>% 
      dplyr::filter(cluster==cells, p_val_adj<.001) %>% .$gene
    
    eg = bitr(tmp.gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    Degenelist[[cells]] <- eg$ENTREZID
  }
  return(Degenelist)
}


#'Description: This is for gene signature heatmap figure
#' @author https://github.com/yujijun/Code_library
#' @param obj Seurat object
#' @param feature a vector of marker gene list
#' @return a stackVlnplot 
#' @export
modify_vlnplot<- function(obj,
                          feature,
                          col_cluster=col_cluster,
                          pt.size = 0,
                          plot.margin = unit(c(0, 0, 0, 0), "cm"),
                          ...) {
  require(Seurat)
  require(ggplot2)
  p<- VlnPlot(obj, features = feature, pt.size = pt.size,cols = col_cluster, ... )  +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_line(),
          axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
          plot.margin = plot.margin)+
    rremove("xlab")
  return(p)
}

## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0,
                          col_cluster=col_cluster,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  require(Seurat)
  require(ggplot2)
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x,
                                                              col_cluster=col_cluster, ...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 90,hjust=1), 
          axis.ticks.x = element_line()
    )
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}





