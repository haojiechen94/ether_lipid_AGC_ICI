#-------------------------------------------------------------------------------------------------------------
# Figure 4 Serum ether-linked phospholipid enrichment is associated with an immune-activated and
# cytotoxic tumor immune microenvironment
#-------------------------------------------------------------------------------------------------------------

# 1) Dissecting TIME based on bulk RNA-seq data
library(readxl)
library(DESeq2)
library(xCell)
library(GSVA)


selected_samples<-res_m[c('EG160','EG225','EG232','EG25','EG6','EG49','EG71','EG79','EG44', 'EG51','EG178','EG19','EG223',
                                          'EG61','EG133','EG145','EG230','EG253', 'EG127','EG144','EG31','EG220','EG170','EG152', 'EG106',
                                          'EG221','EG209','EG136','EG11','EG20'),]
selected_samples$RNA_seq_id<-c('G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8', 'G9', 'G10', 'G11',
                               'G12', 'G13', 'G14', 'G15', 'G16', 'G17', 'G18', 'G19', 'G20', 'G21',
                               'G22', 'G23', 'G24', 'G25', 'G26', 'G27', 'G28', 'G29', 'G30')
selected_samples

TPMs<-read_excel('./TPM.annot.xlsx')
rownames(TPMs)<-TPMs$gene_name
d<-as.data.frame(log2(TPMs[,rownames(temp_df)])+1)
rownames(d)<-rownames(TPMs)
sig<-xCellAnalysis(d,cell.types.use=c('aDC','B-cells','Basophils','CD4+ memory T-cells',
                                      'CD4+ naive T-cells','CD4+ T-cells','CD4+ Tcm',
                                      'CD4+ Tem','CD8+ naive T-cells','CD8+ T-cells',
                                      'CD8+ Tcm','CD8+ Tem','cDC','Class-switched memory B-cells',
                                      'DC','Eosinophils','iDC','Macrophages',
                                      'Macrophages M1','Macrophages M2','Mast cells','Memory B-cells',
                                      'Monocytes','naive B-cells',
                                      'Neutrophils','NK cells','NKT','pDC','Plasma cells',
                                      'pro B-cells','Tgd cells','Th1 cells','Th2 cells','Tregs'))

e<-selected_samples

result1<-list()
for(i in rownames(sig)){
  res<-cor.test(sig[i,e$RNA_seq_id],e[,'ether'],method = 'spearman')
  result1[[i]]<-c(res$estimate,
                  res$p.value)
}

result_df1<-as.data.frame(result1)
result_df1<-as.data.frame(t(result_df1))
colnames(result_df1)<-c('scc','pval')
result_df1_ether<-result_df1
rownames(result_df1_ether)<-rownames(sig)
result_df1_ether[result_df1_ether$pval<0.2,]

result_df1_ether$celltype<-rownames(result_df1_ether)



ggplot(data=result_df1_ether,aes(x = scc, y = -log10(pval),size=-log10(pval))
)+geom_point(aes(fill=scc),shape = 21)+
  geom_hline(yintercept = -log10(0.05), color = "grey", linewidth = 1,linetype = "dashed")+
  geom_vline(xintercept = 0, color = "grey", linewidth = 1,linetype = "dashed")+
  scale_size_continuous(range = c(1,6))+
  scale_fill_gradient2(low = "#0099ff",
                       mid = "white",
                       high = "#ff004c") +
  geom_text_repel(data = result_df1_ether[result_df1_ether$pval<0.5,],
                  aes(label = celltype),size = 5,col = 'black',max.overlaps=30)+
  labs(x='SCC',y="-Log10(p-value)",title="Ether linked PLs")+
  #xlim(-0.5,0.5)+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5,size = 20),
                        axis.text = element_text(color = "black", size = 20),
                        axis.title = element_text(size = 20),
                        legend.text = element_text(size = 15, color = "black"),
                        legend.title = element_text(size = 15, face = "bold"))


library(ComplexHeatmap)
library(circlize)



col_fun1 <- colorRamp2(seq(-1,1, length = 3),  
                       c("#0099ff",'white',"#ff004c"))

col_fun2 <- colorRamp2(seq(-1,1, length = 2),  
                       c('lightgrey','#90e0ef'))

col_fun4 <- colorRamp2(seq(-1,1, length = 2),  
                       c('lightgrey','#9095d2'))


col_fun3 <- colorRamp2(seq(-1,1, length = 3),  
                       c('#B39DDB','#c2c2c2','#fee440'))


heatmap_annotation<-HeatmapAnnotation(ether=e[order(e$ether),'ether'],
                                      ceramide=e[order(e$ether),'Cer'],
                                      col = list(ether = col_fun2,
                                                 ceramide = col_fun4))

row_annotation<-rowAnnotation(SCC=c(result_df1_ether[rownames(result_df1_ether[result_df1_ether$pval<0.05,]),'scc']),
                              col=list(SCC = col_fun3))



Heatmap(t(scale(t(sig[c(rownames(result_df1_ether[result_df1_ether$pval<0.05,])),
                      e[order(e$ether),'RNA_seq_id']]))),
        cluster_columns=F,name = "Relative abundance1",top_annotation=heatmap_annotation,
        col = col_fun1,
        left_annotation = row_annotation)



# 2) Dissecting TIME based on scRNA-seq data
library(Seurat)
library(dplyr)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)



obj_list<-c()
for(path in Sys.glob('./cellranger_ouput/*')){
  name<-strsplit(strsplit(path,'/')[[1]][9],'_')[[1]][5]
  if(name %in% rownames(sample_id_map)){
    print(sample_id_map[name,'Metabolite_ID'])
    data <- Read10X(data.dir = path)
    obj <- CreateSeuratObject(counts=data,project=sample_id,min.cells=3,min.features=200)
    obj<-RenameCells(obj, new.names =unlist(lapply(colnames(obj),function(x){return(paste0(sample_id,'_',x))})))
    obj_list<-c(obj_list,obj)
  }
}

merge.obj<-merge(obj_list[[1]],obj_list[2:length(obj_list)])
merge.obj[["percent.mt"]] <- PercentageFeatureSet(merge.obj, pattern = "^MT-")
VlnPlot(merge.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.01)


#subseting
sub.obj <- subset(merge.obj,  cells = rownames(merge.obj@meta.data)[
  c(merge.obj@meta.data$nFeature_RNA> 2000)&
    c(merge.obj@meta.data$nFeature_RNA< 7000)&
    c(merge.obj@meta.data$nCount_RNA< 40000)&
    c(merge.obj@meta.data$nCount_RNA> 4000)&
    c(merge.obj@meta.data$percent.mt<20)])

VlnPlot(sub.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.01)


#integrating using Harmony
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



#Cell type annotation
sub.obj<-JoinLayers(sub.obj)

cluster_markers<-FindAllMarkers(sub.obj,only.pos=TRUE)




annotations<-list('0'='T_NK_cells','1'='T_NK_cells','4'='T_NK_cells','8'='T_NK_cells','10'='T_NK_cells',
                  
                  '13'='B_cells','5'='Plasma_cells','24'='Plasma_cells',
                  
                  '6'='Myeloid_cells','11'='Myeloid_cells','15'='Myeloid_cells','22'='Myeloid_cells',
                  '31'='Myeloid_cells',
                  
                  '17'='Mast_cells',
                  
                  '3'='Fibroblasts','9'='Fibroblasts',
                  
                  '2'='Endothelial_cells','29'='Endothelial_cells',
                  
                  '7'='Epithelial_cells','14'='Epithelial_cells','18'='Epithelial_cells','19'='Epithelial_cells',
                  '20'='Epithelial_cells','21'='Epithelial_cells','23'='Epithelial_cells','25'='Epithelial_cells',
                  '26'='Epithelial_cells','27'='Epithelial_cells','32'='Epithelial_cells',
                  '28'='Epithelial_cells',
                  
                  '12'='Basal_epithelial_cells','16'='Basal_epithelial_cells',
                  
                  '30'='Unknown','33'='Unknown')

sub.obj@meta.data[['celltype']]<-unlist(lapply(sub.obj@meta.data$harmony_cluster, 
                                               function(x){return(annotations[gettextf('%s',x)])}))


DimPlot(sub.obj,raster=FALSE,label = F, group.by = 'celltype',reduction = "umap.harmony",
        cols=c('#E6A6AC','#4fc3f7','#3366ff','#98DBC6','#9DCBDE','#D48982',
               '#A7E8E2','#F4CD8C','#E19764','grey'))






#T/NK cell subpopulation analysis

T_NK_cells<-subset(sub.obj,
                   cells=rownames(sub.obj@meta.data)[
                     c(sub.obj@meta.data$celltype=='T_NK_cells')])

T_NK_cells <- FindNeighbors(T_NK_cells, reduction = "harmony", dims = 1:30)
T_NK_cells <- FindClusters(T_NK_cells, resolution = 1, cluster.name = "harmony_cluster2")
T_NK_cells <- RunUMAP(T_NK_cells, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

DimPlot(T_NK_cells,reduction = "umap.harmony",label=T,group.by = 'harmony_cluster2')
DimPlot(T_NK_cells,reduction = "umap.harmony",label=T,group.by = 'group')

T_NK_cells<-JoinLayers(T_NK_cells)

cluster_markers<-FindAllMarkers(T_NK_cells,only.pos=TRUE)




annotations<-list('0'='Naive_CD4_cells','4'='Treg17',
                  '1'='Tregs','0'='Naive_CD4_cells',
                  
                  '2'='GZMK+effector_CD8_T_cells',
                  '3'='GZMH+effector_memory_CD8_T_cells',
                  '6'='GZMK+GZMH+effector_memory_CD8_T_cells',
                  
                  '5'='MAIT_cells',
                  
                  '7'='NK_cells','8'='Gamma_delta_T_cells',
                  
                  '9'='Proliferating_CD8_T_cells','10'='Proliferating_CD8_T_cells',
                  '11'='Proliferating_CD4_T_cells',
                  
                  '12'='Unknown')

T_NK_cells@meta.data[['celltype2']]<-unlist(lapply(T_NK_cells@meta.data$harmony_cluster2, 
                                                   function(x){return(annotations[gettextf('%s',x)])}))


DimPlot(T_NK_cells,raster=FALSE,label = F, group.by = 'celltype2',reduction = "umap.harmony",
        cols=c('#E6A6AC','#4fc3f7','#3366ff','#98DBC6','#9DCBDE','#D48982',
               '#A7E8E2','#F4CD8C','#E19764','#DAA5FF','#FFA5AB','grey'))





Seurat::Idents(object=T_NK_cells)<-T_NK_cells@meta.data$celltype2

levels(T_NK_cells) <- c("Naive_CD4_cells","Treg17","Tregs","Proliferating_CD4_T_cells" ,
                        "MAIT_cells",
                        "GZMK+effector_CD8_T_cells",
                        "GZMK+GZMH+effector_memory_CD8_T_cells",
                        "GZMH+effector_memory_CD8_T_cells",
                        "Gamma_delta_T_cells",
                        "Proliferating_CD8_T_cells",
                        "NCAM1+_NK_cells",
                        "Unknown")


DotPlot(T_NK_cells, 
        features = c('CD3D','CD3E','CD3G',
                     'CD4',
                     'LEF1','CCR7','TCF7','SELL',
                     'IL17A','IL17F','IL23R','BATF','STAT3',
                     'FOXP3','IL2RA','TNFRSF4',
                     'CD8A','CD8B',
                     'SLC4A10','NCR3','KLRB1',
                     'GZMK','CXCR4','CCL5',
                     'GZMH','CX3CR1',
                     'GZMB','PRF1',
                     'TRDV1','TRGV8','KIR2DL4',
                     'NCR1','NCAM1','NCR2',
                     'MKI67','TOP2A','STMN1')) + RotatedAxis()+
  scale_color_gradient2(low="#9095d2",mid='white',high="#E6A6AC")




T_NK_cells<-AddModuleScore(T_NK_cells,list(c('GZMK','GZMA','GZMB','PRF1','NKG7')),name='cytotoxic')

T_NK_cells<-AddModuleScore(T_NK_cells,list(c('TIGIT','CTLA4','LAG3','PDCD1','FOXP3')),name='inhibitory')

DotPlot(T_NK_cells,c('cytotoxic1','inhibitory1'))+
  scale_color_gradient2(low="#4fc3f7",mid='white',high="#ff6f61", midpoint = 0,limits = c(-1, 1),oob = scales::squish)


T_NK_cells@meta.data[['group']]<-ULM_scores[T_NK_cells@meta.data$orig.ident,'group']


cell_numbers<-table(T_NK_cells$celltype2, T_NK_cells$group)
cell_numbers

oddratios<-c()
pvals<-c()
celltypes<-c()
for(i in rownames(cell_numbers)){
  manual_matrix <- matrix(c(cell_numbers[i,'High'], cell_numbers[i,'Low'], 
                            sum(cell_numbers[,'High'])-cell_numbers[i,'High'],
                            sum(cell_numbers[,'Low'])-cell_numbers[i,'Low']), nrow = 2, byrow = TRUE)
  
  print(i)
  res<-fisher.test(manual_matrix)
  oddratios<-c(oddratios,res$estimate)
  pvals<-c(pvals,res$p.value)
  celltypes<-c(celltypes,i)
}

fisher_res_high<-data.frame('OR'=oddratios,'pval'=pvals,'padj'=p.adjust(pvals,method='BH'),
                            'celltype'=celltypes)
fisher_res_high



oddratios<-c()
pvals<-c()
celltypes<-c()
for(i in rownames(cell_numbers)){
  manual_matrix <- matrix(c(cell_numbers[i,'Low'], cell_numbers[i,'High'], 
                            sum(cell_numbers[,'Low'])-cell_numbers[i,'Low'],
                            sum(cell_numbers[,'High'])-cell_numbers[i,'High']), nrow = 2, byrow = TRUE)
  
  print(i)
  res<-fisher.test(manual_matrix)
  oddratios<-c(oddratios,res$estimate)
  pvals<-c(pvals,res$p.value)
  celltypes<-c(celltypes,i)
}

fisher_res_low<-data.frame('OR'=oddratios,'pval'=pvals,'padj'=p.adjust(pvals,method='BH'),
                           'celltype'=celltypes)
fisher_res_low


or_matrix<-cbind(fisher_res_high[,'OR',drop=F],fisher_res_low[,'OR',drop=F])
rownames(or_matrix)<-fisher_res_high$celltype
colnames(or_matrix)<-c('High','Low')

padj_matrix<-cbind(fisher_res_high[,'padj',drop=F],fisher_res_low[,'padj',drop=F])
rownames(padj_matrix)<-fisher_res_high$celltype
colnames(padj_matrix)<-c('High','Low')
padj_matrix

log2_or_matrix <- log2(or_matrix)
log2_or_matrix


sig_matrix <- matrix(rep('ns',24), nrow = 12, byrow = TRUE)

sig_matrix


rownames(sig_matrix) <- rownames(or_matrix)
colnames(sig_matrix) <- colnames(or_matrix)


for(i in rownames(sig_matrix)){
  for(j in colnames(sig_matrix)){
    if(padj_matrix[i,j]<=0.05 & padj_matrix[i,j]>0.01){
      sig_matrix[i,j]<-'*'
    }else if(padj_matrix[i,j]<=0.01 & padj_matrix[i,j]>0.001){
      sig_matrix[i,j]<-'**'
    }else if(padj_matrix[i,j]<0.001){
      sig_matrix[i,j]<-'***'
    }
  }
}





library(ComplexHeatmap)
library(circlize)


col_fun <- colorRamp2(c(-2, 0, 2), c("#9095d2", "white", "#E6A6AC"))

Heatmap(
  matrix = log2_or_matrix, 
  col = col_fun,
  name = "log2(OR)",
  cluster_columns = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sig_matrix[i, j], x, y, gp = gpar(fontsize = 10))
  },
  rect_gp = gpar(col = "white", lwd = 1))










