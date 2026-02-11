#-------------------------------------------------------------------------------------------------------------
# Figure 3 Serum ether-linked phospholipids are associated with stratify prognosis and clinical benefit
# to chemo-ICI in AGC
#-------------------------------------------------------------------------------------------------------------


library(ggrepel)
library(ggplot2)
library(factoextra)
library(survival)
library(survminer)
library(GSVA)
library(GSEABase)
library(pheatmap)
library(scales)


#Optimizing cutoff to group patients into high- and low- expression groups
c<-res_m
best_cutoff<-0.1
best_cutoff_p<-1
for(i in seq(0.1,0.9,0.001)){
  c[,'group']<-'High'
  c[c$ether < quantile(c$ether,i),'group']<-'Low'
  r<-survdiff(formula('Surv(OS_time , OS_event)~ group'), data = c)
  if(r$pvalue<best_cutoff_p){
    best_cutoff_p<-r$pvalue
    best_cutoff<-i
  }
}
best_cutoff

c[,'group']<-'High'
c[c$ether< quantile(c$ether,best_cutoff),'group']<-'Low'


fit <- survfit(formula('Surv(OS_time , OS_event)~ group'), data = c)


c$group<-factor(c$group,levels = c('Low','High'))
cox_model <- coxph(Surv(OS_time, OS_event) ~ group, data = c)
summary_cox <- summary(cox_model)

HR <- round(summary_cox$coef[2], 2)
CI_lower <- round(summary_cox$conf.int[,"lower .95"], 2)
CI_upper <- round(summary_cox$conf.int[,"upper .95"], 2)

ggsurv<-ggsurvplot(fit, data = c, title = "Ether-linked PLs",conf.int = TRUE,
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





c<-res_m
best_cutoff<-0.1
best_cutoff_p<-1
for(i in seq(0.1,0.9,0.001)){
  c[,'group']<-'High'
  c[c$Ceramide < quantile(c$Ceramide,i),'group']<-'Low'
  r<-survdiff(formula('Surv(OS_time , OS_event)~ group'), data = c)
  if(r$pvalue<best_cutoff_p){
    best_cutoff_p<-r$pvalue
    best_cutoff<-i
  }
}
best_cutoff

c[,'group']<-'High'
c[c$Ceramide< quantile(c$Ceramide,best_cutoff),'group']<-'Low'



c$group<-factor(c$group,levels = c('Low','High'))
cox_model <- coxph(Surv(OS_time, OS_event) ~ group, data = c)
summary_cox <- summary(cox_model)

HR <- round(summary_cox$coef[2], 2)
CI_lower <- round(summary_cox$conf.int[,"lower .95"], 2)
CI_upper <- round(summary_cox$conf.int[,"upper .95"], 2)

fit <- survfit(formula('Surv(OS_time , OS_event)~ group'), data = c)

ggsurv<-ggsurvplot(fit, data = c, title = "Ceramide",conf.int = TRUE,
                   palette = c("#0099ff","#ff004c"),
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








#-------------------------------------------------------------------------------------------------------------------------------
library(pROC)


temp_dat<-res_m[res_m$treatment_response %in% c('PR','PD','SD','CR'),]

temp_dat<-temp_dat[,c('CPS','ether','Ceramide','treatment_response')]
temp_dat <- na.omit(temp_dat)

temp_dat$treatment_response2<-'NCB'
temp_dat[temp_dat$treatment_response %in% c('SD','PR','CR'),'treatment_response2']<-'DCB'
temp_dat$treatment_response2<-factor(temp_dat$treatment_response2,levels=c('NCB','DCB'))

roc_obj <- roc(temp_dat[['treatment_response2']], temp_dat[['ether']],smooth = TRUE, smooth.method = "binormal")
auc_value_ether<- auc(roc_obj)

specificities<-c()
sensitivities<-c()
methods<-c()


specificities<-c(specificities,rev(1-roc_obj$specificities))
sensitivities<-c(sensitivities,rev(roc_obj$sensitivities))
methods<-c(methods,rep(gettextf('ether=%.3f',auc_value_ether),length(roc_obj$sensitivities)))


roc_obj <- roc(temp_dat[['treatment_response2']], temp_dat[['Ceramide']],smooth = TRUE, smooth.method = "binormal")
auc_value_ceramide<- auc(roc_obj)


specificities<-c(specificities,rev(1-roc_obj$specificities))
sensitivities<-c(sensitivities,rev(roc_obj$sensitivities))
methods<-c(methods,rep(gettextf('ceramide=%.3f',auc_value_ceramide),length(roc_obj$sensitivities)))




roc_obj <- roc(temp_dat[['treatment_response2']], temp_dat[['CPS']],smooth = TRUE, smooth.method = "binormal")
auc_value_cps<- auc(roc_obj)

specificities<-c(specificities,rev(1-roc_obj$specificities))
sensitivities<-c(sensitivities,rev(roc_obj$sensitivities))
methods<-c(methods,rep(gettextf('CPS=%.3f',auc_value_cps),length(roc_obj$sensitivities)))




auc_values_df<-data.frame(
  'Method'=methods,
  'FPR'=specificities,
  'TPR'=sensitivities
)

library(ggplot2)

ggplot(auc_values_df, aes(x = FPR, y = TPR, group = Method, color = Method)) +
  geom_line(size=1.5) +ggtitle('AUC')+
  scale_fill_manual(values = c('#0099ff','#eda1b9','#ff004c'))+
  scale_color_manual(values = c('#0099ff','#eda1b9','#ff004c'))+
  labs(color = "AUC")+
  theme_classic()+  theme(plot.title = element_text(hjust = 0.5,size = 20),
                          axis.text = element_text(color = "black", size = 20),
                          axis.title = element_text(size = 20),
                          legend.text = element_text(size = 15, color = "black"),
                          legend.title = element_text(size = 15, face = "bold"))+
  geom_abline()













