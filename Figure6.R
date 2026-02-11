#-------------------------------------------------------------------------------------------------------------
# Figure 6 Oral PE-P supplementation synergizes with chemo-ICI in huCD34+ humanized mouse
# models
#-------------------------------------------------------------------------------------------------------------


library(ggplot2)
library(ggprism)
library(dplyr)
library(ggpubr)

#--------------------------------------------------------------------------------------------------------------------------------
df_G11<-read.table('./G1_volumns.txt',
                   sep='\t',header=T,stringsAsFactors = F,check.names = F)
df_G11

df_G11$Condition<-factor(df_G11$Condition,levels = c('Vehicle','PE_P','antiPD1_Oxa','antiPD1_Oxa_PE_P'))


ggplot(df_G11,
       aes(x = Days, y = Volumn,color=Sample_id)) +
  geom_line(lwd=1)+facet_wrap(~ Condition, ncol = 4)+
  scale_fill_manual(values = rev(c(rep('#706F6F',4),rep('#793C92',4),rep('#89A7CB',4),rep('#EA5412',4))))+
  scale_color_manual(values = rev(c(rep('#706F6F',4),rep('#793C92',4),rep('#89A7CB',4),rep('#EA5412',4))))+
  xlim(0,15)+ylim(0,1000)+
  theme_prism() 


summary_data <- df_G11 %>% group_by(Condition,Days) %>% summarise(
  Mean = mean(Volumn, na.rm = TRUE),     
  SD = sd(Volumn, na.rm = TRUE),         
  N = n(),                              
  SE = SD / sqrt(N),                    
  .groups = 'drop'                       
)

summary_data


ggplot(summary_data, aes(x = Days, y = Mean, color = Condition, fill = Condition,shape=Condition)) +
  geom_errorbar(aes(ymin = Mean , ymax = Mean + SE), width = 0.5)+
  geom_line(size = 1)+  xlim(0,15)+ylim(0,800)+
  theme_prism() +
  geom_point(size = 4)+
  scale_fill_manual(values = c('#706F6F','#793C92','#89A7CB','#EA5412'))+
  scale_color_manual(values = c('#706F6F','#793C92','#89A7CB','#EA5412'))+
  scale_shape_manual(values = c(21,24,22,25))




df_G11<-read.table('./G1_weights.txt',
                   sep='\t',header=T,stringsAsFactors = F,check.names = F)
df_G11

df_G11$Condition<-factor(df_G11$Condition,levels = c('Vehicle','PE_P','antiPD1_Oxa','antiPD1_Oxa_PE_P'))

df_G11$Weight<-df_G11$Weight*1000

ggplot(df_G11,
       aes(x = Condition, y = Weight, fill  = Condition, shape = Condition,color=Condition)) +
  geom_point(
    position = position_jitter(width = 0.2, height = 0),
    size = 2.5,
    alpha = 1
  )+
  stat_summary(
    fun = mean,
    geom = "errorbar",
    aes(ymin = ..y.., ymax = ..y..),
    size = 1,
    alpha = 1,
    width=0.5
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    size = 1
  ) +scale_fill_manual(values = c('#706F6F','#793C92','#89A7CB','#EA5412'))+
  scale_color_manual(values = c('#706F6F','#793C92','#89A7CB','#EA5412'))+
  scale_shape_manual(values = c(21,24,22,25))+
  theme_prism(axis_text_angle = 45)+ylim(0,800)



t.test(df_G11$Weight[df_G11$Condition=='antiPD1_Oxa'],
       df_G11$Weight[df_G11$Condition=='antiPD1_Oxa_PE_P'])


t.test(df_G11$Weight[df_G11$Condition=='PE_P'],
       df_G11$Weight[df_G11$Condition=='antiPD1_Oxa_PE_P'])

t.test(df_G11$Weight[df_G11$Condition=='Vehicle'],
       df_G11$Weight[df_G11$Condition=='antiPD1_Oxa_PE_P'])


#--------------------------------------------------------------------------------------------------------------------------------
df_G11<-read.table('./G2_volumns.txt',
                   sep='\t',header=T,stringsAsFactors = F,check.names = F)
df_G11

df_G11$Condition<-factor(df_G11$Condition,levels = c('Vehicle','PE_P','antiPD1_Oxa','antiPD1_Oxa_PE_P'))


ggplot(df_G11,
       aes(x = Days, y = Volumn,color=Sample_id)) +
  geom_line(lwd=1)+facet_wrap(~ Condition, ncol = 4)+
  scale_fill_manual(values = rev(c(rep('#706F6F',4),rep('#793C92',4),rep('#89A7CB',4),rep('#EA5412',4))))+
  scale_color_manual(values = rev(c(rep('#706F6F',4),rep('#793C92',4),rep('#89A7CB',4),rep('#EA5412',4))))+
  xlim(0,20)+ylim(0,1500)+
  theme_prism() 


summary_data <- df_G11 %>% group_by(Condition,Days) %>% summarise(
  Mean = mean(Volumn, na.rm = TRUE),      
  SD = sd(Volumn, na.rm = TRUE),          
  N = n(),                               
  SE = SD / sqrt(N),                     
  .groups = 'drop'                      
)

summary_data


ggplot(summary_data, aes(x = Days, y = Mean, color = Condition, fill = Condition,shape=Condition)) +
  geom_errorbar(aes(ymin = Mean , ymax = Mean + SE), width = 0.5)+
  geom_line(size = 1)+  xlim(0,20)+ylim(0,1500)+
  theme_prism() +
  geom_point(size = 4)+
  scale_fill_manual(values = c('#706F6F','#793C92','#89A7CB','#EA5412'))+
  scale_color_manual(values = c('#706F6F','#793C92','#89A7CB','#EA5412'))+
  scale_shape_manual(values = c(21,24,22,25))


t.test(df_G11$Volumn[df_G11$Condition=='antiPD1_Oxa' & df_G11$Days==17],
       df_G11$Volumn[df_G11$Condition=='antiPD1_Oxa_PE_P' & df_G11$Days==17])

t.test(df_G11$Volumn[df_G11$Condition=='Vehicle' & df_G11$Days==17],
       df_G11$Volumn[df_G11$Condition=='antiPD1_Oxa_PE_P' & df_G11$Days==17])

t.test(df_G11$Volumn[df_G11$Condition=='PE_P' & df_G11$Days==17],
       df_G11$Volumn[df_G11$Condition=='antiPD1_Oxa_PE_P' & df_G11$Days==17])


df_G11<-read.table('./G2_weights.txt',
                   sep='\t',header=T,stringsAsFactors = F,check.names = F)
df_G11

df_G11$Condition<-factor(df_G11$Condition,levels = c('Vehicle','PE_P','antiPD1_Oxa','antiPD1_Oxa_PE_P'))

df_G11$Weight<-df_G11$Weight*1000

ggplot(df_G11,
       aes(x = Condition, y = Weight, fill  = Condition, shape = Condition,color=Condition)) +
  geom_point(
    position = position_jitter(width = 0.2, height = 0),
    size = 2.5,
    alpha = 1
  )+
  stat_summary(
    fun = mean,
    geom = "errorbar",
    aes(ymin = ..y.., ymax = ..y..),
    size = 1,
    alpha = 1,
    width=0.5
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    size = 1
  ) +scale_fill_manual(values = c('#706F6F','#793C92','#89A7CB','#EA5412'))+
  scale_color_manual(values = c('#706F6F','#793C92','#89A7CB','#EA5412'))+
  scale_shape_manual(values = c(21,24,22,25))+
  theme_prism(axis_text_angle = 45)+ylim(0,1000)



