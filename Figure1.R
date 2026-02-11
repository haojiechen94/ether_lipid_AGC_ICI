#-------------------------------------------------------------------------------------------------------------
# Figure 1 Distinct serum metabolic signatures associated with therapeutic response and survival 
# outcomes in AGC patients treated with chemo-ICI
#-------------------------------------------------------------------------------------------------------------

library(ggrepel)
library(ggplot2)
library(factoextra)
library(survival)
library(survminer)
library(GSVA)
library(GSEABase)
library(ropls)
library(sva)



#gastric cancer patient serum metabolic profiles
metabolites<-read.table('../metabolites.txt',sep='\t',
                        stringsAsFactors = F,check.names = F,header=T)
rownames(metabolites)<-metabolites$Index
metabolites


#patient clinical information
patient_info<-read.table('../patient_infos.txt',sep='\t',
                         stringsAsFactors = F,check.names = F,header=T,row.names = 1)
patient_info


#log transform
metabolites<-log2(metabolites+5)

data<-cbind(metabolites,patient_info)
data

#response
data$response<-'NA'
data[data$treatment_response =='PR','response']<-'DCB'
data[data$treatment_response =='PD','response']<-'NCB'
data[data$treatment_response =='SD','response']<-'DCB'
data[data$treatment_response =='CR','response']<-'DCB'

#survival
median_surv <- 18*30
data$surv_group<-'NA'
data[c(data$OS_time<=median_surv & data$OS_event==1),'surv_group']<-'short_survival'
data[c(data$OS_time>median_surv),'surv_group']<-'long_survival'

#-------------------------------------------------------------------------------------------
# 1) Responder vs Non-responder

#differential abundance analysis using t-test
metabolites<-c()
lfcs<-c()
pvals<-c()
tvals<-c()
for(i in colnames(data)){
  r<-t.test(data[data$response=='R',i],
            data[data$response=='NR',i])
  metabolites<-c(metabolites,i)
  pvals<-c(pvals,r$p.value)
  tvals<-c(tvals,r$statistic)
  lfcs<-c(lfcs,mean(data[data$response=='R',i])-mean(data[data$response=='NR',i]))
}


res_df<-data.frame('metabolite'=metabolites,
                   'LFC'=lfcs,
                   'tval'=tvals,
                   'pval'=pvals)
res_df$padj<-p.adjust(res_df$pval,method='BH')


#metabolites ranked by t-statistics
ranks <- setNames(res_df$tval, res_df$metabolite)
ranks <- sort(ranks, decreasing = TRUE) 




#metabolite set enrichment analysis
metabolite_classification_df<-read.table('./metabolite_classification.txt',sep='\t',
                       stringsAsFactors = F,check.names = F,header=T,row.names = 1)

metabolite_classification_sets <- split(metabolite_classification_df$metabolite,
                                        metabolite_classification_df$classification)
metabolite_classification_sets



library(fgsea)


fgsea_res <- fgseaSimple(
  pathways = metabolite_classification_sets,
  stats = ranks,
  minSize = 10,
  nperm = 10000
)


#visualization
plotGseaTable(
  pathways = pathway_list[c('GP','GL','Bile acids','FA')],
  stats = ranks,
  fgseaRes = fgsea_res,
  gseaParam = 0.5)


#KEGG metabolism pathway enrichment analysis

kegg_pathway<-getGmt('./KEGG_pathways.gmt')

kegg_pathway_list <- list()
for(i in c(1:length(kegg_pathway@.Data))){
  kegg_pathway_list[[kegg_pathway@.Data[[i]]@setName]]<-kegg_pathway@.Data[[i]]@geneIds
  
}

kegg_fgsea_res <- fgsea(
  pathways = kegg_pathway_list,       
  stats = ranks,            
  minSize = 10,             
  maxSize = 1000
)




kegg_fgsea_res <- kegg_fgsea_res[order(pval), ] 

#visualization
library(dplyr)

bar_df <- kegg_fgsea_res %>%
  filter(!is.na(padj)) %>%
  arrange(desc(NES)) %>%
  slice_head(n = 5) %>%
  bind_rows(
    kegg_fgsea_res %>%
      arrange(NES) %>%
      slice_head(n = 5)
  ) %>%
  arrange(NES) %>%
  mutate(
    pathway = factor(pathway, levels = pathway),
    Direction = ifelse(NES > 0, "Upregulated", "Downregulated")
  )

bar_df



library(ggplot2)

ggplot(bar_df,
       aes(x = pathway, y = NES, fill = Direction)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Upregulated"   = "#ff004c",
      "Downregulated" = "#0099ff"
    )
  ) +
  theme_classic()+ theme(plot.title = element_text(hjust = 0.5,size = 20),
                         axis.text = element_text(color = "black", size = 20),
                         axis.title = element_text(size = 20),
                         legend.text = element_text(size = 15, color = "black"),
                         legend.title = element_text(size = 15, face = "bold"),
                         axis.text.x = element_text(angle=0, vjust=1, hjust=0.5))+
  labs(
    x = NULL,
    y = "Normalized Enrichment Score (NES)"
  )



#-------------------------------------------------------------------------------------------
# 2) Long survival group vs Short survival group

#differential abundance analysis using t-test
metabolites<-c()
lfcs<-c()
pvals<-c()
tvals<-c()

for(i in colnames(data.df2)){
  r<-t.test(data.df2[rownames(pheno_data2)[pheno_data2$surv_group=='long_survival'],i],
            data.df2[rownames(pheno_data2)[pheno_data2$surv_group=='short_survival'],i])
  
  metabolites<-c(metabolites,i)
  pvals<-c(pvals,r$p.value)
  tvals<-c(tvals,r$statistic)
  lfcs<-c(lfcs,mean(data.df2[rownames(pheno_data2)[pheno_data2$surv_group=='long_survival'],i])-mean(data.df2[rownames(pheno_data2)[pheno_data2$surv_group=='short_survival'],i]))
  
}


res_df<-data.frame('metabolite'=metabolites,
                   'LFC'=lfcs,
                   'tval'=tvals,
                   'pval'=pvals)


res_df$padj<-p.adjust(res_df$pval,method='BH')



#metabolites ranked by t-statistics
ranks <- setNames(res_df$tval, res_df$metabolite)
ranks <- sort(ranks, decreasing = TRUE) 




#metabolite set enrichment analysis
metabolite_classification_df<-read.table('./metabolite_classification.txt',sep='\t',
                                         stringsAsFactors = F,check.names = F,header=T,row.names = 1)

metabolite_classification_sets <- split(metabolite_classification_df$metabolite,
                                        metabolite_classification_df$classification)
metabolite_classification_sets



library(fgsea)


fgsea_res <- fgseaSimple(
  pathways = metabolite_classification_sets,
  stats = ranks,
  minSize = 10,
  nperm = 10000
)


#visualization
plotGseaTable(
  pathways = pathway_list[c('GP','GL','Bile acids','FA')],
  stats = ranks,
  fgseaRes = fgsea_res,
  gseaParam = 0.5)


#KEGG metabolism pathway enrichment analysis

kegg_pathway<-getGmt('./KEGG_pathways.gmt')

kegg_pathway_list <- list()
for(i in c(1:length(kegg_pathway@.Data))){
  kegg_pathway_list[[kegg_pathway@.Data[[i]]@setName]]<-kegg_pathway@.Data[[i]]@geneIds
  
}

kegg_fgsea_res <- fgsea(
  pathways = kegg_pathway_list,       
  stats = ranks,            
  minSize = 10,             
  maxSize = 1000
)




kegg_fgsea_res <- kegg_fgsea_res[order(pval), ] 

#visualization
library(dplyr)

bar_df <- kegg_fgsea_res %>%
  filter(!is.na(padj)) %>%
  arrange(desc(NES)) %>%
  slice_head(n = 5) %>%
  bind_rows(
    kegg_fgsea_res %>%
      arrange(NES) %>%
      slice_head(n = 5)
  ) %>%
  arrange(NES) %>%
  mutate(
    pathway = factor(pathway, levels = pathway),
    Direction = ifelse(NES > 0, "Upregulated", "Downregulated")
  )

bar_df



library(ggplot2)

ggplot(bar_df,
       aes(x = pathway, y = NES, fill = Direction)) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Upregulated"   = "#ff004c",
      "Downregulated" = "#0099ff"
    )
  ) +
  theme_classic()+ theme(plot.title = element_text(hjust = 0.5,size = 20),
                         axis.text = element_text(color = "black", size = 20),
                         axis.title = element_text(size = 20),
                         legend.text = element_text(size = 15, color = "black"),
                         legend.title = element_text(size = 15, face = "bold"),
                         axis.text.x = element_text(angle=0, vjust=1, hjust=0.5))+
  labs(
    x = NULL,
    y = "Normalized Enrichment Score (NES)"
  )

