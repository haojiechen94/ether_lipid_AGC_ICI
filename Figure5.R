#-------------------------------------------------------------------------------------------------------------
# Figure 5 Oral PE-P supplementation enhances chemo-ICI efficacy by remodeling the TIME
#-------------------------------------------------------------------------------------------------------------

#Dissecting PE-P associated TIME changesbased on scRNA-seq data
library(Seurat)
library(dplyr)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)


obj_list<-c()
for(path in Sys.glob('./*/outs/filtered_feature_bc_matrix/')){
  name<-strsplit(path,'/')[[1]][10]
  print(name)
  data <- Read10X(data.dir = path)
  obj <- CreateSeuratObject(counts=data,project=name,min.cells=3,min.features=200)
  obj<-RenameCells(obj, new.names =unlist(lapply(colnames(obj),function(x){return(paste0(name,'_',x))})))
  obj_list<-c(obj_list,obj)
}

merge.obj<-merge(obj_list[[1]],obj_list[2:length(obj_list)])
merge.obj[["percent.mt"]] <- PercentageFeatureSet(merge.obj, pattern = "^mt-")
VlnPlot(merge.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.01,raster=FALSE)

sub.obj <- subset(merge.obj,  cells = rownames(merge.obj@meta.data)[
  c(merge.obj@meta.data$nFeature_RNA> 1000)&
    c(merge.obj@meta.data$nFeature_RNA< 6000)&
    c(merge.obj@meta.data$nCount_RNA< 40000)&
    c(merge.obj@meta.data$nCount_RNA> 2000)&
    c(merge.obj@meta.data$percent.mt<5)])

VlnPlot(sub.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.01)

dim(merge.obj)
dim(sub.obj)


unique(sub.obj@meta.data$orig.ident)

condation_map<-list( "CIT_PBMC1"='CIT_PBMC',"CIT_PBMC2"='CIT_PBMC',"CIT_PBMC3"='CIT_PBMC',
                     "CIT_PEP_PBMC1"='CIT_PEP_PBMC',"CIT_PEP_PBMC2"='CIT_PEP_PBMC',"CIT_PEP_PBMC3"='CIT_PEP_PBMC',
                     "CIT_PEP_T1"='CIT_PEP_T',"CIT_PEP_T2"='CIT_PEP_T',"CIT_PEP_T3"='CIT_PEP_T',
                     "CIT_T1"='CIT_T',"CIT_T2"='CIT_T',"CIT_T3"='CIT_T')

sub.obj@meta.data$condition<-as.character(condation_map[sub.obj@meta.data$orig.ident])
sub.obj<-JoinLayers(sub.obj)

sub.obj@meta.data$orig.ident

sub.obj[["RNA"]] <- split(sub.obj[["RNA"]], f = sub.obj@meta.data$orig.ident)

sub.obj <- NormalizeData(sub.obj, normalization.method = "LogNormalize", scale.factor = 10000)
sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 2000)
sub.obj <- ScaleData(sub.obj)
sub.obj <- RunPCA(sub.obj, features = VariableFeatures(object = sub.obj))



sub.obj <- IntegrateLayers(
  object = sub.obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)



sub.obj <- FindNeighbors(sub.obj, reduction = "harmony", dims = 1:30)
sub.obj <- FindClusters(sub.obj, resolution = 1, cluster.name = "harmony_cluster")

sub.obj <- RunUMAP(sub.obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

DimPlot(sub.obj,reduction = "umap.harmony",label=T,group.by = 'harmony_cluster')

DimPlot(sub.obj,reduction = "umap.harmony",label=T,group.by = 'orig.ident')

DimPlot(sub.obj,reduction = "umap.harmony",label=T,group.by = 'condition')


annotations<-list('0'='Monocyte_derived_macrophages','3'='Tissue_resident_macrophages','14'='Mmp12_macrophages',
                  '13'='Monocytes',
                  '31'='Dendritic_cells',
                  '25'='Proliferating_macrophages',
                  '10'='Neutrophils',
                  
                  '2'='B_cells','16'='B_cells',
                  '38'='Proliferating_B_cells',
                  
                  '1'='Naive_CD4_T_cells','17'='Naive_CD4_T_cells',
                  '23'='Tregs','37'='ISGs_CD4_T_cells',
                  '33'='Platelet',
                  
                  '7'='Naive_CD8_T_cells','26'='Naive_CD8_T_cells',
                  '9'='Effector_CD8_T_cells',
                  '21'='Proliferating_CD8_T_cells',
                  '20'='NK_cells',
                  
                  '4'='Malignant_cells','6'='Malignant_cells','8'='Malignant_cells',
                  '15'='Malignant_cells','30'='Malignant_cells',
                  '5'='Proliferating_malignant_cells',
                  
                  '18'='Fibroblasts','32'='Endothelial_cells',
                  
                  '12'='Doublets',
                  '27'='Doublets','11'='Doublets',
                  
                  '24'='Unknown','29'='Unknown','34'='Unknown','35'='Unknown','39'='Unknown',
                  '19'='Unknown','22'='Unknown','36'='Unknown','28'='Unknown')

sub.obj@meta.data[['celltype']]<-unlist(lapply(sub.obj@meta.data$harmony_cluster, 
                                               function(x){return(annotations[gettextf('%s',x)])}))


DimPlot(sub.obj, reduction = "umap.harmony",raster=FALSE,label = F, group.by = 'celltype')

DimPlot(sub.obj, reduction = "umap.harmony",raster=FALSE,label = F, group.by = 'condition')


unique(sub.obj@meta.data$celltype)

Seurat::Idents(object=sub.obj)<-sub.obj@meta.data$celltype

levels(sub.obj) <- c(
  "Naive_CD4_T_cells","ISGs_CD4_T_cells","Tregs",
  "Naive_CD8_T_cells","Effector_CD8_T_cells","Proliferating_CD8_T_cells",
  "NK_cells","B_cells","Proliferating_B_cells","Platelet",
  
  "Monocytes","Monocyte_derived_macrophages","Tissue_resident_macrophages","Mmp12_macrophages","Proliferating_macrophages",
  "Dendritic_cells",
  "Neutrophils",
  
  "Endothelial_cells","Fibroblasts",
  
  
  "Malignant_cells","Proliferating_malignant_cells",
  "Unknown","Doublets"
)




DotPlot(sub.obj, 
        features = c('Cd3d','Cd3e','Cd3g','Cd4',
                     'Ccr7','Lef1','Tcf7',
                     'Ifit1','Ifit3b','Ifit3',
                     'Foxp3','Ctla4','Il2ra',
                     'Cd8a','Cd8b1',
                     'Nkg7','Gzmb','Ccl5',
                     'Havcr2','Lag3','Pdcd1',
                     'Ncr1','Prf1','Klrd1',
                     'Cd79a','Cd79b','Ighd',
                     'Gp9','Ppbp','F5',
                     'Ccr2','Ly6c2','Chil3',
                     'C1qa','C1qb','C1qc',
                     'Il1b','Tnf','Fcgr1',
                     'Mertk','Mrc1','Folr2',
                     'Arg1','Clec4d','Nos2',
                     'Flt3','Clec9a','Clec10a',
                     'Retnlg','S100a9','S100a8',
                     'Cdh5','Cdh13','Adgrf5',
                     'Mfap5','Postn','Dcn',
                     'Sox9','Krt10','Wnt7b',
                     'Mki67','Top2a','Stmn1'
        )) + RotatedAxis()+
  scale_color_gradient2(low='darkblue',mid='white',high='darkred')


