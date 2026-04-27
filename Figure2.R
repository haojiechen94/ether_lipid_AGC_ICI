#-------------------------------------------------------------------------------------------------------------
# Figure 2 Metabolite set–based profiling identifies prognostic lipid metabolic serotypes associated 
# with chemo-ICI response in AGC
#-------------------------------------------------------------------------------------------------------------


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


#Using univariate linear regression method to estimate the expression activity of serum metabolite module signatures
#import metabolite module
metabolite_sets<-getGmt('./metabolons_hmdb.gmt')


ulm_scores<-list()
for(i in rownames(data)){
  a_list_of_metabolites<-unlist(data[i,])
  temp<-c()
  for(j in names(metabolite_sets)){
    temp_df<-data.frame('weight'=0,
                        'metabolite'=as.vector(a_list_of_metabolites))
    rownames(temp_df)<-names(a_list_of_metabolites)
    temp_df[rownames(temp_df) %in% metabolite_sets[[j]]@geneIds,'weight']<- 1
    r<-summary(lm(metabolite ~ weight,temp_df))
    temp<-c(temp,r$coefficients['weight','t value'])
  }
  ulm_scores[[i]]<-temp
}


ulm_res<-scale(t(as.data.frame(ulm_scores)))

colnames(ulm_res)<-names(metabolite_sets)
rownames(ulm_res)<-rownames(ulm_scores)



pvals_df<-as.data.frame(cbind(ulm_res1[,c('HR','p')],ulm_res2[,c('HR','p')],ulm_res3[,c('HR','p')]))

colnames(pvals_df)<-c('SYSUCC_ICI_1_HR','SYSUCC_ICI_1_p',
                      'GDGH_ICI_2_HR','GDGH_ICI_2_p',
                      'SYSUCC_ICI_2_HR','SYSUCC_ICI_2_p')
pvals_df


#Stoufer's method
stoufer_method<-function(HR,p){
  if(HR>1){
    z<-qnorm(1-p/2, mean = 0, sd = 1, lower.tail = TRUE)
    return(z)
  }else if(HR<=1){
    z<-qnorm(p/2, mean = 0, sd = 1, lower.tail = TRUE)
    return(z)        
  }
}

combined_p<-function(x,weights){
  weights<-weights/sum(weights)
  z<-sum(c(stoufer_method(x[1],x[2])*weights[1],
           stoufer_method(x[3],x[4])*weights[2],
           stoufer_method(x[5],x[6])*weights[3]))/sqrt(sum(weights**2))
  if(z>0){
    p<- 1-pnorm(z, mean = 0, sd = 1, lower.tail = TRUE)
    return(p)
  }else{
    p<- pnorm(z, mean = 0, sd = 1, lower.tail = TRUE)
    return(p)
  }
}


cps_s<-c()
for(i in rownames(pvals_df)){
  cps_s<-c(cps_s,combined_p(unlist(pvals_df[i,]),c(241,49)))
}

pvals_df$stouffer<-cps_s
pvals_df[order(pvals_df$stouffer),]



pvals_df$stouffer_rank<-rank(pvals_df$stouffer)
pvals_df$stouffer_sig<-'ns'
pvals_df[pvals_df$stouffer<0.05,'stouffer_sig']<-'sig'


pvals_df$metabolite<-rownames(pvals_df)

pvals_df


pvals_df$HR<-rowMeans(pvals_df[,c('SYSUCC_ICI_1_HR','SYSUCC_ICI_3_HR')])
pvals_df

#lipid metabolite modules vs non-lipid metabolite modules
#rank plot
pvals_df$type<-'lipids'
pvals_df[rownames(pvals_df) %in% c('Nucleotide_and_Its_metabolites','Sugars',
                                   'Amino_acid_derivatives','Phenolic_acids',
                                   'Organic_acid_and__Its_derivatives','Amines',
                                   'Organic_acid_and_Its_derivatives',
                                   'Heterocyclic_compounds','Small_Peptide','Polyamines','CoEnzyme_and_vitamins',
                                   'Alcohols','CAR','Indole_and_Its_derivatives','Phenolics'),
         'type']<-'non-lipids'

ggplot(data=pvals_df,aes(x = stouffer_rank, y = -log10(stouffer)))+geom_point(aes(fill=type),size=5,shape = 21)+
  scale_fill_manual(values = c("#ff004c", "#0099ff"))+
  scale_color_manual(values = c("#ff004c","#0099ff"))+
  geom_text_repel(data = pvals_df[pvals_df$stouffer<0.05,],
                  aes(label = metabolite),size = 5,col = 'black',max.overlaps=30)+
  labs(x='Rank',y="-Log10(combined p-value)",title="Stouffer's method")+
  theme_classic()+ theme(plot.title = element_text(hjust = 0.5,size = 20),
                         axis.text = element_text(color = "black", size = 20),
                         axis.title = element_text(size = 20),
                         legend.text = element_text(size = 15, color = "black"),
                         legend.title = element_text(size = 15, face = "bold"),
                         axis.text.x = element_text(angle=0, vjust=1, hjust=0.5))





#focus on lipid metabolite modules

#volcano plot
#Stouffer's method
ggplot(data=pvals_df[pvals_df$type=='lipids',],aes(x = log(HR), y = -log10(stouffer),size=-log10(stouffer))
)+geom_point(aes(fill=log(HR)),shape = 21)+
  geom_hline(yintercept = -log10(0.05), color = "grey", linewidth = 1,linetype = "dashed")+
  geom_vline(xintercept = 0, color = "grey", linewidth = 1,linetype = "dashed")+
  scale_size_continuous(range = c(1,6))+
  scale_fill_gradient2(low = "#0099ff",
                       mid = "white",
                       high = "#ff004c") +
  geom_text_repel(data = pvals_df[pvals_df$stouffer<0.05 & pvals_df$type=='lipids',],
                  aes(label = metabolite),size = 5,col = 'black',max.overlaps=30)+
  labs(x='Log(HR)',y="-Log10(combined p-value)",title="Stouffer's method")+
  xlim(-0.5,0.5)+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5,size = 20),
                        axis.text = element_text(color = "black", size = 20),
                        axis.title = element_text(size = 20),
                        legend.text = element_text(size = 15, color = "black"),
                        legend.title = element_text(size = 15, face = "bold"))




#---------------------------------------------------------------------------------------------------
#Clustering GC patients based on serum metabolon enrichment scores, serum metabolon-defined subtypes, serotypes
#applying consensus clustering
library(ConsensusClusterPlus)


setwd('./serum_types/')


data<-t(all_data[,c('PE-P, 'PE-O', 'LPE','PC','DG', 'Cer-NDS','Cer-NP', 'Cer-NS', 'Cer-AS',''Cer-AP')])
head(data)




results <- ConsensusClusterPlus(SYSUCC_ICI_1_and_2_data,maxK=5,reps=50,pItem=0.8,pFeature=1,
                                title='Discovery_cohort',clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="pdf")


#choosing k=2
all_data$group<-as.vector(results[[2]][["consensusClass"]])
all_data$serotype<-unlist(lapply(all_data$group, function(x){return(gettextf('Serotype%s',x))}))


#PCA visualization
library(pcaMethods)
pcs<-pca(all_data[,c('PE_P','PE_O','LPC','LPC_O','DG','Cer_NDS','Cer_AP','Cer_AS')],nPcs=3)


library(factoextra)
library(ggplot2)

pca_df<-as.data.frame(list(PC1=scores(pcs)[,'PC1'],
                           PC2=scores(pcs)[,'PC2'],
                           Serotype=unlist(lapply(all_data$group, function(x){return(gettextf('Serotype%s',x))}))))


ggplot(pca_df,
       aes(x = PC1,
           y = PC2)) +
  geom_point(aes(fill=Serotype),size=3,shape=21) +
  scale_fill_manual(values = c("#ff004c", "#0099ff"))+
  scale_color_manual(values = c("#ff004c","#0099ff"))+
  ggforce::geom_mark_ellipse(aes(fill = Serotype,
                                 color = Serotype))+
  labs(title='PCA on serum metabolite profiles')+xlab(gettextf('PC 1 (%.1f%%)',pcs@R2[1]*100))+
  ylab(gettextf('PC 2 (%.1f%%)',pcs@R2[2]*100))+
  coord_cartesian(xlim = c(-7.5, 7.5), ylim = c(-5,5))+ theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 18))




#survival differences
library(patchwork)

custom_theme <- function(base_size = 14) {
  theme_survminer() %+replace%
    theme(
      # 主图设置
      plot.title = element_text(size = base_size + 2, face = "bold"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 2),
      axis.text.y = element_text(size = base_size - 2, color = "black"),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 2),
      legend.position = "top",  
      
     
      table = list(
        theme(
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = base_size - 2, color = "black"),
          plot.title = element_text(size = base_size, hjust = 0)
        )
      )
    )
}



fit <- survfit(formula('Surv(OS_time , OS_event)~ serotype'), data = all_data)

all_data$serotype<-factor(all_data$serotype,levels = c('Serotype1','Serotype2'))
cox_model <- coxph(Surv(OS_time, OS_event) ~ serotype, data = all_data)
summary_cox <- summary(cox_model)

HR <- round(summary_cox$coef[2], 2)

CI_lower <- round(summary_cox$conf.int[,"lower .95"], 2)
CI_upper <- round(summary_cox$conf.int[,"upper .95"], 2)

ggsurv<-ggsurvplot(fit, data = all_data, title = "Discovery cohort: Serotype",conf.int = TRUE,
                   palette = c("#ff004c", "#0099ff"),
                   pval = T,
                   font.title = c(16, "bold", "darkblue"),
                   risk.table = TRUE,pval.size = 5,
                   risk.table.height = 0.25,risk.table.y.text = T,
                   surv.plot.height = 0.75,
                   fontsize =6,
                   ggtheme = custom_theme(base_size = 18),
                   tables.theme = custom_theme(base_size = 18)$table)

ggsurv$plot<-ggsurv$plot +
  annotate("text", 
           x = 25, 
           y = 0.4,
           label = paste0("HR = ", HR, " (95% CI: ", CI_lower, "-", CI_upper, ")"),
           size = 5)

(ggsurv$plot + theme(axis.text.x=element_blank())) / 
  ggsurv$table + 
  plot_layout(heights=c(3,1))




#heatmap visulization
library(ComplexHeatmap)
library(circlize)


col_fun <- colorRamp2(seq(-2,2, length = 3),  
                      c("#0099ff",'white',"#ff004c"))

col_fun1 <- colorRamp2(seq(15,35, length = 2),  
                       c('white','orchid'))

col_fun2 <- colorRamp2(seq(20,80, length = 2),  
                       c('white','blue'))


col_fun3 <- colorRamp2(seq(0,100, length = 2),  
                       c('white','#1A5B5B'))

col_fun4 <- colorRamp2(seq(0,150, length = 2),  
                       c('white','#009E73'))


samples<-c(rownames(all_data)[all_data$serotype=='Serotype1'],
           rownames(all_data)[all_data$serotype=='Serotype2'])


all_data$sex<-patient_info[samples,'sex']

heatmap_annotation<-HeatmapAnnotation(Group=c(rep('Serotype1',sum(all_data$serotype=='Serotype1')),
                                              rep('Serotype2',sum(all_data$serotype=='Serotype2'))),
                                      smocking_history=all_data[samples,'smocking'],
                                      drinking_history=all_data[samples,'drinking'],
                                      BMI=all_data[samples,'BMI'],
                                      Age=all_data[samples,'age'],
                                      Sex=all_data[samples,'sex'],
                                      
                                      HER2=all_data[samples,'HER2'],
                                      
                                      CPS=all_data[samples,'CPS'],
                                      
                                      TMD=all_data[samples,'TMD(mm)'],
                                      
                                      differentiation=all_data[samples,'differentiation(L=1;M=2;H=3)'],
                                      col=list(Group=c('Serotype1'="#ff004c",'Serotype2'="#0099ff"),
                                               smocking_history=c('0'="white",'1'="black"),
                                               drinking_history=c('0'="white",'1'="black"),
                                               Sex=c('0'="white",'1'="black"),
                                               BMI=col_fun1,
                                               Age=col_fun2,
                                               HER2=c('0'="white",'1'="black"),
                                               CPS=col_fun3,
                                               TMD=col_fun4,
                                               differentiation=c('1'="#8ecae6",'2'="#ffb703",'3'="#fb8500")),
                                      na_col = "grey")


Heatmap(t(all_data[samples,
                   c('PE_P','PE_O','LPC','LPC_O','DG','Cer_NDS','Cer_AP','Cer_AS')]),
        cluster_rows=F,
        cluster_columns=F,
        heatmap_legend_param=list(title='Rel.Ab.'),
        top_annotation=heatmap_annotation,
        show_row_names=T,show_column_names=F,col = col_fun)




#bulid a simple LASSO logistic regression model
library(glmnet)


set.seed(12345)

index <- createDataPartition(all_data$group, p=.8, list=FALSE, times=1)


train_df <- all_data[index,c('PE-P, 'PE-O', 'LPE','PC','DG', 'Cer-NDS','Cer-NP', 'Cer-NS', 'Cer-AS',''Cer-AP','group')]
test_df <- all_data[-index,c('PE-P, 'PE-O', 'LPE','PC','DG', 'Cer-NDS','Cer-NP', 'Cer-NS', 'Cer-AS',''Cer-AP','group')]

train_df$group[train_df$group==1] <- "Serotype1"
train_df$group[train_df$group==2] <- "Serotype2"

test_df$group[test_df$group==1] <- "Serotype1"
test_df$group[test_df$group==2] <- "Serotype2"


train_df$group <- as.factor(train_df$group)
test_df$group <- as.factor(test_df$group)


x <- model.matrix(group ~ ., train_df)[, -1]
y <- train_df$group

library(glmnet)

lasso <- cv.glmnet(x, y, family="binomial", alpha=1,nfolds=5,type.measure="deviance")



coefficients<-as.data.frame(as.matrix(coef(lasso, s = "lambda.1se")))

coefficients<-coefficients[2:11,,drop=F]

coefficients<-coefficients[order(coefficients$s1),,drop=F]



coefficients$m<-rownames(coefficients)

coefficients$m<-factor(coefficients$m,levels = rownames(coefficients))


ggplot(data=coefficients)+
  geom_point(aes(x=s1,y=m), size=10,shape=16,color='#E76F51')+
  geom_vline(xintercept = 0, linetype="dashed")+labs(x="Coefficient", y="")+ggtitle('LASSO logistic regression')+
  theme_classic()+  theme(plot.title = element_text(hjust = 0.5,size = 20),
                          axis.text = element_text(color = "black", size = 20),
                          axis.title = element_text(size = 20),
                          legend.text = element_text(size = 15, color = "black"),
                          legend.title = element_text(size = 15, face = "bold"))


predictions <- predict(lasso,
                       model.matrix(group ~ ., test_df)[, -1], s = "lambda.1se", type = "response")

pred_class <- ifelse(predictions > 0.5, "Serotype2", "Serotype1")


confusionMatrix(data=predictions, test_df$group)

library(pROC)


res_train<-as.data.frame(predict(model1, newdata=train_df,type = "prob"))
res_train$true<-train_df$group
res_train

roc_obj_train <- roc(res_train$true, res_train$Serotype2)
auc_value_train<- auc(roc_obj_train)
auc_value_train
plot(roc_obj_train, col = "red", lwd = 2, cex.lab=2, cex.axis=2,cex.main=2)


res_test<-as.data.frame(predict(model1, newdata=test_df,type = "prob"))
res_test$true<-test_df$group

roc_obj_test <- roc(res_test$true, res_test$Serotype2)
auc_value_test<- auc(roc_obj_test)
auc_value_test
plot(roc_obj_test, col = "red", lwd = 2, cex.lab=2, cex.axis=2,cex.main=2)










SYSUCC_ICI_ESCA_res_m<-read.table('./SYSUCC_ICI_ESCA_res_m.txt',
                                  sep='\t',header=T,stringsAsFactors = F,check.names = F,row.names = 1)


SYSUCC_ICI_ESCA_res_m




predictions <- predict(model1, newdata=SYSUCC_ICI_ESCA_res_m[,c('PE-P, 'PE-O', 'LPE','PC','DG', 'Cer-NDS','Cer-NP', 'Cer-NS', 'Cer-AS',''Cer-AP')])


SYSUCC_ICI_ESCA_res_m$group<-as.vector(predictions)


fit <- survfit(formula('Surv(OS_time , OS_event)~ group'), data = SYSUCC_ICI_ESCA_res_m)



SYSUCC_ICI_ESCA_res_m$group<-factor(SYSUCC_ICI_ESCA_res_m$group,levels = c('Serotype1','Serotype2'))
cox_model <- coxph(Surv(OS_time, OS_event) ~ group, data = SYSUCC_ICI_ESCA_res_m)
summary_cox <- summary(cox_model)

HR <- round(summary_cox$coef[2], 2)
CI_lower <- round(summary_cox$conf.int[,"lower .95"], 2)
CI_upper <- round(summary_cox$conf.int[,"upper .95"], 2)

ggsurv<-ggsurvplot(fit, data = SYSUCC_ICI_ESCA_res_m, title = "External cohort: Serotype",conf.int = TRUE,
                   palette = c("#ff004c", "#0099ff"),
                   pval = T,
                   font.title = c(16, "bold", "darkblue"),
                   risk.table = TRUE,pval.size = 5,
                   risk.table.height = 0.25,risk.table.y.text = T,
                   surv.plot.height = 0.75,
                   fontsize =6,
                   ggtheme = custom_theme(base_size = 18),
                   tables.theme = custom_theme(base_size = 18)$table)

ggsurv$plot<-ggsurv$plot +
  annotate("text", 
           x = 15, 
           y = 0.4,
           label = paste0("HR = ", HR, " (95% CI: ", CI_lower, "-", CI_upper, ")"),
           size = 5)


(ggsurv$plot + theme(axis.text.x=element_blank())) / 
  ggsurv$table + 
  plot_layout(heights=c(3,1))




#---------------------------------------------------------------------------------------------------


#import metabolite module
kegg_metabolism_pathways<-getGmt('./KEGG_pathways.gmt')

ulm_scores<-list()
for(i in rownames(data)){
  a_list_of_metabolites<-unlist(data[i,])
  temp<-c()
  for(j in names(kegg_metabolism_pathways)){
    temp_df<-data.frame('weight'=0,
                        'metabolite'=as.vector(a_list_of_metabolites))
    rownames(temp_df)<-names(a_list_of_metabolites)
    temp_df[rownames(temp_df) %in% kegg_metabolism_pathways[[j]]@geneIds,'weight']<- 1
    r<-summary(lm(metabolite ~ weight,temp_df))
    temp<-c(temp,r$coefficients['weight','t value'])
  }
  ulm_scores[[i]]<-temp
}


ulm_res<-scale(t(as.data.frame(ulm_scores)))

colnames(ulm_res)<-names(kegg_metabolism_pathways)
rownames(ulm_res)<-rownames(data)

ulm_res<-as.data.frame(ulm_res)
ulm_res

ulm_res$serotype<-all_data[rownames(ulm_res),'serotype']
head(ulm_res)

ms<-c()
diffs<-c()
ps<-c()

for(i in colnames(ulm_res)[1:(length(colnames(ulm_res))-1)]){
  print(i)
  r<-t.test(ulm_res[i][ulm_res['serotype']=='Serotype1'],
            ulm_res[i][ulm_res['serotype']=='Serotype2'])
  ms<-c(ms,i)
  diffs<-c(diffs,mean(ulm_res[i][ulm_res['serotype']=='Serotype1'])-mean(ulm_res[i][ulm_res['serotype']=='Serotype2']))
  ps<-c(ps,r$p.value)
}

DPA_df<-data.frame('pathway'=ms,'pval'=ps,'diff'=diffs)
DPA_df$padj<-p.adjust(DPA_df$pval,method='BH')


ggplot(data=DPA_df,aes(x = diff, y = -log10(padj),size=-log10(padj))
)+geom_point(aes(fill=diff),shape = 21)+
  geom_hline(yintercept = -log10(0.05), color = "grey", linewidth = 1,linetype = "dashed")+
  geom_vline(xintercept = 0, color = "grey", linewidth = 1,linetype = "dashed")+
  scale_size_continuous(range = c(1,6))+
  scale_fill_gradientn(colours = c( "#0099ff",'white',"#ff004c")) +
  geom_text_repel(data = DPA_df[DPA_df$padj<1e-10,],
                  aes(label = pathway),size = 5,col = 'black',max.overlaps=30)+
  labs(x='Mean differences',y="-Log10(adjusted p-value)",title="KEGG metabolism pathway")+
  #xlim(-0.5,0.5)+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5,size = 20),
                        axis.text = element_text(color = "black", size = 20),
                        axis.title = element_text(size = 20),
                        legend.text = element_text(size = 15, color = "black"),
                        legend.title = element_text(size = 15, face = "bold"))


ggboxplot(ulm_res,x='serotype',y="Ether_lipid_metabolism",color='serotype',size=2)+ stat_compare_means()+
  scale_fill_manual(values = rev(c('#0099ff',"#ff004c")))+
  scale_color_manual(values = rev(c('#0099ff',"#ff004c")))+
  theme_classic()+  theme(plot.title = element_text(hjust = 0.5,size = 20),
                          axis.text = element_text(color = "black", size = 18),
                          axis.title = element_text(size = 20),
                          legend.text = element_text(size = 15, color = "black"),
                          legend.title = element_text(size = 15, face = "bold"),
                          axis.text.x = element_text(angle=25, vjust=0.5, hjust=0.5))


ggboxplot(ulm_res,x='serotype',y="Sphingolipid_metabolism",color='serotype',size=2)+ stat_compare_means()+
  scale_fill_manual(values = rev(c('#0099ff',"#ff004c")))+
  scale_color_manual(values = rev(c('#0099ff',"#ff004c")))+
  theme_classic()+  theme(plot.title = element_text(hjust = 0.5,size = 20),
                          axis.text = element_text(color = "black", size = 18),
                          axis.title = element_text(size = 20),
                          legend.text = element_text(size = 15, color = "black"),
                          legend.title = element_text(size = 15, face = "bold"),
                          axis.text.x = element_text(angle=25, vjust=0.5, hjust=0.5))

