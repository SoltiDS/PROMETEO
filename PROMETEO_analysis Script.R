library(ggplot2)
library(dplyr)
library(readr)
library(tidyverse)
library(impute)
library(readxl)

theme_basic<-function(){
  ggplot2::theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill=NA),
    legend.position = "bottom",
    text = element_text(family = "Arial", size = 10),
    axis.title = element_text(family = "Arial", size = 10, face = "bold"),
    axis.text = element_text(family = "Arial", size = 10, face="bold"),
    legend.text = element_text(family = "Arial", size = 10),
    legend.title = element_text(family = "Arial", size = 10),
    plot.title = element_text(family = "Arial", size = 10, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(family = "Arial", size = 10),
    plot.caption = element_text(family = "Arial", size = 10))
}


############## Figure 4a ###########################

Plotdata <- read_excel("SourceData.xlsx", 
                       sheet = "Figure 4a")

Plotdata<-Plotdata%>%mutate(TILs=as.numeric(TILs))

summarySE <- function(data = NULL, measurevar, groupvars = NULL, na.rm = FALSE, conf.interval = .95, .drop = TRUE) {
  require(plyr)
  length2 <- function (x, na.rm = FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  datac <- ddply(data, groupvars, .drop = .drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm = na.rm),
                     mean = mean(xx[[col]], na.rm = na.rm),
                     sd   = sd(xx[[col]], na.rm = na.rm)
                   )
                 },
                 measurevar
  )
  
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)
  ciMult <- qt(conf.interval / 2 + .5, datac$N - 1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# Crear resumen estadístico de los datos
tgc <- summarySE(Plotdata, measurevar="TILs", groupvars=c( "Visit"), na.rm=TRUE)

# Renombrar columnas según lo que se espera en ggplot
colnames(tgc) <- c("Visit", "N", "TILs", "sd", "se", "ci")


pdf("TILs_5timepoints.pdf", width=12, height=7)
ggplot(Plotdata, aes(x = Visit, y = TILs)) +
  geom_boxplot(fill = "grey", outlier.shape = NA, color = "black") +  # box without default outliers
  geom_jitter(shape = 3, width = 0.2, size = 2, alpha = 0.7) +        # points as '+' symbols
  theme_minimal(base_size = 12) +
  labs(x = "Group", y = "Value")+
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 10))+
  theme_classic()
dev.off()

write.table(tgc,"fig4a.stats.txt",sep = "\t",row.names = F)



############## Figure 4c ###########################
Plotdata <- read_excel("SourceData.xlsx", 
                       sheet = "Figure 4c")

Plotdata<-Plotdata%>%mutate(`PDL1%`=as.numeric(`PDL1%`))%>%
  mutate(Visit = recode(Visit, "PRENAC"="Baseline"))%>%
  filter(Visit%in%c("Baseline","SCR","SURGERY" ))%>%
  mutate(Visit = factor(Visit,levels = c("Baseline","SCR","SURGERY" )))

# Crear resumen estadístico de los datos
tgc <- summarySE(Plotdata, measurevar="PDL1%", groupvars=c( "Visit"), na.rm=TRUE)

# Renombrar columnas según lo que se espera en ggplot
colnames(tgc) <- c("Visit", "N", "TILs", "sd", "se", "ci")


pdf("pdl1_5timepoints.pdf", width=12, height=7)
ggplot(Plotdata, aes(x = Visit, y = `PDL1%`)) +
  geom_boxplot(fill = "grey", outlier.shape = NA, color = "black") +  # box without default outliers
  geom_jitter(shape = 3, width = 0.2, size = 2, alpha = 0.7) +        # points as '+' symbols
  theme_minimal(base_size = 12) +
  labs(x = "Group", y = "Value")+
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, by = 5))+
  theme_classic()
dev.off()

write.table(tgc,"fig4a.stats.txt",sep = "\t",row.names = F)

############## Figure 4f ###########################
Plotdata <- openxlsx::read.xlsx("SourceData.xlsx", 
                                sheet = "Figure 4e", startRow = 2)[1:201,c(1,2)]

Plotdata2<-read_excel("SourceData.xlsx", 
                      sheet = "Figure 4e", skip = 205)

Plotdata$Signature_ID_paper<-stringr::str_trim(Plotdata$Signature_ID_paper,"both")
Plotdata2$`Immune Gene Signature`<-trimws(Plotdata2$`Immune Gene Signature`,"both")

Plotdata$Signature_ID_paper<-gsub("^\\s+|\\s+$", "", Plotdata$Signature_ID_paper)
colnames(Plotdata2)[1]="Signature_ID_paper"

Plotdata1<-full_join(Plotdata,Plotdata2
                    )


matrix<-Plotdata1%>%select(Baseline,SCR,C2D1,C3D1,SUR)

matrix<-as.matrix(matrix)

matrix<-apply(matrix,2,as.numeric)

heatmap<-pheatmap(matrix, cluster_rows = TRUE, cluster_cols = TRUE, scale = "row") 

############## Figure 4i ###########################
library(ggplot2)
Plotdata <- read_excel("SourceData.xlsx", 
                         sheet = "Figure 4i")

colors_subtype<-c("HR_POS"= "#4169E1",
                  "TNBC"="#d01352")

Plotdata<-Plotdata%>%
  mutate(Visit = recode(Visit,"PRENAC"="Baseline"))%>%
  mutate(IGG_Signature_Mean = as.numeric(IGG_Signature_Mean))%>%
  filter(!is.na(IGG_Signature_Mean))

Plotdata$IGG_Signature_Mean<-scales::rescale(Plotdata$IGG_Signature_Mean, to=c(0,5))

plot<- ggplot2::ggplot(data=Plotdata,
                       
                       aes(x=Visit,y=IGG_Signature_Mean))+
  
  ggplot2::geom_boxplot(aes(fill=Subtype),outlier.shape = NA)+
  
  geom_jitter(aes(color="black"), position = position_jitterdodge(jitter.width = 0, dodge.width = 0.9), size = 1.5, alpha = 0.85, shape=3)+
  stat_summary(fun.y = median,
               
               geom = 'line',
               
               aes(group = Subtype, colour = Subtype),
               
               position = position_dodge(width = 0.9), linewidth=0.8)+
  
  theme_basic()+
  
  scale_fill_manual(values=colors_subtype)+
  scale_color_manual(values=colors_subtype)

ggsave("figure4i.pdf",plot, device=cairo_pdf)

###stats ###
box_core <- Plotdata %>%
  group_by(Subtype,Visit) %>%
  summarise(
    N = n(),
    lower = quantile(IGG_Signature_Mean, 0.25) - 1.5 * IQR(IGG_Signature_Mean),
    upper = quantile(IGG_Signature_Mean, 0.75) + 1.5 * IQR(IGG_Signature_Mean),
    lower_whisker = max(min(IGG_Signature_Mean), quantile(IGG_Signature_Mean, 0.25) - 1.5 * IQR(IGG_Signature_Mean)),
    upper_whisker = min(max(IGG_Signature_Mean), quantile(IGG_Signature_Mean, 0.75) + 1.5 * IQR(IGG_Signature_Mean)),
    centre = median(IGG_Signature_Mean)
  ) %>%
  select(lower_whisker, centre, upper_whisker, N, Visit)%>%
  arrange(Visit)

box_core
write.table(box_core,"fig4i.stats.txt",sep = "\t",row.names = F)


############## Figure 5a ######################
Plotdata <- read_excel("SourceData.xlsx", 
                       sheet = "Figure 5a")

colors_pCR<-c("No"= "#f9a22c",
                  "Yes"="#a42c33")

Plotdata<-Plotdata%>%
  mutate(Visit = recode(Visit,"PRENAC"="Baseline"))%>%
  mutate(TILs= as.numeric(TILs))%>%
  mutate(Visit=factor(Visit,levels=c("Baseline","SCR","C2D1","C3D1","Surgery")))%>%
  filter(!is.na(Visit))
  

#Plotdata$IGG_Signature_Mean<-scales::rescale(Plotdata$IGG_Signature_Mean, to=c(0,5))

plot1<- ggplot2::ggplot(data=Plotdata %>%filter(pCR=="No"),
                       
                       aes(x=Visit,y=TILs))+
  
  ggplot2::geom_boxplot(aes(fill=pCR),outlier.shape = NA)+
  
  geom_jitter(aes(color="black"), size = 1.5, alpha = 1, shape=3,width = 0)+
  
  theme_basic()+
  
  scale_fill_manual(values=colors_pCR)+
  scale_color_manual(values=colors_pCR)+
  
  ggnewscale::new_scale_colour()+
  geom_line(data = Plotdata %>% filter(pCR=="No"), aes(group=UNIQUEID,color=UNIQUEID))+
  theme(legend.position = "none")+
  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + theme(legend.position = "none")


plot2<- ggplot2::ggplot(data=Plotdata %>%filter(pCR=="Yes"),
                        
                        aes(x=Visit,y=TILs))+
  
  ggplot2::geom_boxplot(aes(fill=pCR),outlier.shape = NA)+
  
  geom_jitter(aes(color="black"), size = 1.5, alpha = 1, shape=3,width = 0)+
  
  theme_basic()+
  
  scale_fill_manual(values=colors_pCR)+
  scale_color_manual(values=colors_pCR)+
  
  ggnewscale::new_scale_colour()+
  geom_line(data = Plotdata %>% filter(pCR=="Yes"), aes(group=UNIQUEID,color=UNIQUEID))+
  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  theme(legend.position = "none")

plot<-ggpubr::ggarrange(plotlist = list(plot2,plot1),nrow=1)

ggsave("figure5a.pdf",plot, device=cairo_pdf)

###stats ###
box_core <- Plotdata %>%
  group_by(Subtype,Visit) %>%
  summarise(
    N = n(),
    lower = quantile(IGG_Signature_Mean, 0.25) - 1.5 * IQR(IGG_Signature_Mean),
    upper = quantile(IGG_Signature_Mean, 0.75) + 1.5 * IQR(IGG_Signature_Mean),
    lower_whisker = max(min(IGG_Signature_Mean), quantile(IGG_Signature_Mean, 0.25) - 1.5 * IQR(IGG_Signature_Mean)),
    upper_whisker = min(max(IGG_Signature_Mean), quantile(IGG_Signature_Mean, 0.75) + 1.5 * IQR(IGG_Signature_Mean)),
    centre = median(IGG_Signature_Mean)
  ) %>%
  select(lower_whisker, centre, upper_whisker, N, Visit)%>%
  arrange(Visit)

box_core
write.table(box_core,"fig4i.stats.txt",sep = "\t",row.names = F)



####################### Figure 5b #########################################
Plotdata <- read_excel("SourceData.xlsx",
                       sheet = "Figure 5b")

colors_pCR<-c("No"= "#f9a22c",
              "Yes"="#a42c33")

Plotdata1<-Plotdata%>%
  mutate(Visit = recode(Visit,"PRENAC"="Baseline"))%>%
  mutate(`PDL1%`= as.numeric(`PDL1%`))%>%
  mutate(Visit=factor(Visit,levels=c("Baseline","SCR","C2D1","C3D1","SURGERY")))%>%
  filter(!is.na(Visit))%>%
  filter(`PDL1%`>=1)


#Plotdata$IGG_Signature_Mean<-scales::rescale(Plotdata$IGG_Signature_Mean, to=c(0,5))

plot1<- ggplot2::ggplot(data=Plotdata,

                        aes(x=Visit,y=`PDL1%`))+

  ggplot2::geom_boxplot(aes(fill=pCR),outlier.shape = NA)+

  geom_jitter(aes(color="black"), size = 1.5, alpha = 1, shape=3,width = 0)+

  theme_basic()+

  scale_fill_manual(values=colors_pCR)+
  scale_color_manual(values=colors_pCR)+

  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + theme(legend.position = "none")

################### Figure 5c ##################################

colors_pCR<-c("No"= "#f9a22c",
              "Yes"="#a42c33")
classes<-c("B cells ", "T cells ", "B cells/T cells ", "Ig  ", "NK cells ")


Plotdata <- read_excel("SourceData.xlsx",
                       sheet = "Figure 5c",skip=1)[1:111,]

Plotdata1 <- read_excel("SourceData.xlsx",
                       sheet = "Figure 5c",skip=117)[,c(1,3)]

Plotdata1$`Signature_ID_paper `<-stringi::stri_trim(Plotdata1$`Signature_ID_paper `,"both")

colnames(Plotdata1)[1]<-"variable"

data<-as_tibble(reshape2::melt(Plotdata,id.var="TIME_PCR"))%>%full_join(Plotdata1)%>%mutate(value=as.numeric(value))

data$`Immune_Class `<-stringi::stri_trim(data$`Immune_Class `,"both")


ggplot2::ggplot(data=data%>%filter(`Immune_Class `=="B cells"),

                aes(x=TIME_PCR,y=value))+

  ggplot2::geom_boxplot(aes(fill=TIME_PCR),outlier.shape = NA)+

  geom_jitter(aes(color="black"), size = 1.5, alpha = 1, shape=3,width = 0)+

  theme_basic()+
  facet_wrap(vars(`Immune_Class `))

  scale_fill_manual(values=colors_pCR)+
  scale_color_manual(values=colors_pCR)+

  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + theme(legend.position = "none")


#######################EXTENDED FIGURE 3######################################################
library(networkD3)
library(dplyr)
library(tidyr)

#Data
data <- data.frame(
  Timepoint1 = c("Basal1", "Basal1", "Basal1", "Basal1", "Basal1", "Basal1", "Basal1", "Basal1", "Her21", "LumA1", "LumA1", "LumA1", "LumA1", "LumA1", "LumA1", "LumB1", "LumB1", "LumB1", "LumB1", "LumB1", "LumB1", "LumB1", "LumB1", "LumB1", "Normal1", "Normal1", "Normal1"),
  Timepoint2 = c("Basal2", "Basal2", "Basal2", "Basal2", "Basal2", "Normal2", "Normal2", "Normal2", "Normal2", "LumA2", "LumA2", "LumA2", "Normal2", "Normal2", "Normal2", "LumA2", "LumA2", "LumA2", "LumA2", "LumA2", "LumA2", "LumB2", "Normal2", "Normal2", "Basal2", "Normal2", "Normal2"),
  Timepoint3 = c("Basal3", "Basal3", "Basal3", "Basal3", "Normal3", "Normal3", "Normal3", "Normal3", "Normal3", "LumA3", "Normal3", "Normal3", "Normal3", "Normal3", "Normal3", "LumA3", "LumA3", "Normal3", "Normal3", "Normal3", "LumB3", "NA3", "Normal3", "Normal3", "Basal3", "Normal3", "Normal3"),
  Timepoint4 = c("Basal4", "Basal4", "Basal4", "Her24", "Normal4", "Normal4", "Normal4", "Normal4", "NA4", "Normal4", "LumA4", "Normal4", "Normal4", "Normal4", "Normal4", "NA4", "Normal4", "Normal4", "Basal4", "LumA4", "Normal4", "NA4", "Normal4", "Normal4", "Basal4", "Normal4", "Normal4"),
  Timepoint5 = c("Basal5", "Basal5", "Basal5", "Basal5", "Normal5", "Basal5", "Normal5", "Normal5", "Normal5", "Normal5", "LumB5", "Normal5", "Normal5", "Normal5", "Normal5", "Normal5", "Normal5", "Normal5", "Basal5", "Normal5", "Normal5", "NA5", "Normal5", "Normal5", "Basal5", "Basal5", "LumA5")
)

links <- data.frame()
for (i in 1:(ncol(data) - 1)) {
  temp_links <- data %>%
    select(starts_with("Timepoint")) %>%
    transmute(source = .[[i]], target = .[[i + 1]]) %>%
    group_by(source, target) %>%
    summarise(value = n(), .groups = 'drop')
  links <- bind_rows(links, temp_links)
}

nodes <- data.frame(name = unique(c(as.character(links$source), as.character(links$target))))

nodes$id <- 0:(nrow(nodes) - 1)

links <- links %>%
  left_join(nodes, by = c("source" = "name")) %>%
  rename(source_id = id) %>%
  left_join(nodes, by = c("target" = "name")) %>%
  rename(target_id = id)

sankeyNetwork(Links = links, Nodes = nodes,
              Source = 'source_id', Target = 'target_id',
              Value = 'value', NodeID = 'name',
              units = 'Observations')


########################### Extended figure 4 #############################
Plotdata <- read_excel("SourceData.xlsx", 
                        sheet = "Extended Fig4")
data_plot<-reshape2::melt(Plotdata,id.var=c("Visit","SUBTYPE"))

colors_subtype<-c("HR_POS"= "#4169E1",
                  "TNBC"="#d01352")

data_plot<-data_plot%>%
  mutate(Visit = recode(Visit,"PRENAC"="Baseline"))%>%
  mutate(value = as.numeric(value))


plot<- ggplot2::ggplot(data=data_plot,
                       
                       aes(x=Visit,y=value,fill=SUBTYPE))+
  
  ggplot2::geom_boxplot(aes(fill=SUBTYPE),outlier.shape = NA)+
  
  geom_jitter(aes(color="black"), position = position_jitterdodge(jitter.width = 0, dodge.width = 0.9), size = 1.5, alpha = 0.85, shape=3)+
  stat_summary(fun.y = median,
               
               geom = 'line',
               
               aes(group = SUBTYPE, colour = SUBTYPE),
               
               position = position_dodge(width = 0.9), linewidth=0.8)+
  facet_wrap(vars(variable),scales="free")+
  
  theme_basic()+
  
  scale_fill_manual(values=colors_subtype)+
  scale_color_manual(values=colors_subtype)

ggsave("Ext.fig4.pdf",plot, device=cairo_pdf)

###stats ###
data_plot<-as_tibble(data_plot)

data_plot<-data_plot%>%mutate(Visit=factor(Visit),
                            SUBTYPE= factor(SUBTYPE),
                            variable = factor(variable))

box_core <- data_plot%>%
  dplyr::group_by(Visit,variable,SUBTYPE) %>%
  summarise(
    lower = quantile(value, 0.25,na.rm=T) - 1.5 * IQR(value,na.rm = T),
    upper = quantile(value, 0.75,na.rm=T) + 1.5 * IQR(value,na.rm = T),
    lower_whisker = max(min(value,na.rm=T), quantile(value, 0.25,na.rm=T) - 1.5 * IQR(value,na.rm=T),na.rm=T),
    upper_whisker = min(max(value,na.rm = T), quantile(value, 0.75,na.rm=T) + 1.5 * IQR(value,na.rm=T),na.rm = T),
    centre = median(value,na.rm = T)
  )

box_core
write.table(box_core,"Exfig4i.stats.txt",sep = "\t",row.names = F)


################################ Extended Figure 5a ##################
Plotdata <-  read_excel("SourceData.xlsx", 
                       sheet = "Extended Fig5a")

Plotdata<-Plotdata%>%
  mutate(Visit = recode(Visit,"PRENAC"="Baseline"))%>%
  mutate(TILs= as.numeric(TILs))%>%
  filter(!is.na(Visit))


#Plotdata$IGG_Signature_Mean<-scales::rescale(Plotdata$IGG_Signature_Mean, to=c(0,5))

plot1<- ggplot2::ggplot(data=Plotdata%>%filter(pCR!="ND"),
                        
                        aes(x=pCR,y=TILs))+
  
  ggplot2::geom_boxplot(aes(fill="lightgrey"),outlier.shape = NA)+
  
  geom_jitter(aes(color="black"), size = 1.5, alpha = 1, shape=3,width = 0.2)+
  
  theme(legend.position = "none")+
  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + theme(legend.position = "none")

ggsave("Extfigure5a.pdf",plot, device=cairo_pdf)

###stats ###
box_core <- Plotdata %>%
  group_by(pCR) %>%
  summarise(
    N = n(),
    lower = quantile(TILs, 0.25) - 1.5 * IQR(TILs),
    upper = quantile(TILs, 0.75) + 1.5 * IQR(TILs),
    lower_whisker = max(min(TILs), quantile(TILs, 0.25) - 1.5 * IQR(TILs)),
    upper_whisker = min(max(TILs), quantile(TILs, 0.75) + 1.5 * IQR(TILs)),
    centre = median(TILs)
  )

box_core
write.table(box_core,"Extfig5a.stats.txt",sep = "\t",row.names = F)

################################ Extended Figure 5b ##################
Plotdata <-  read_excel("SourceData.xlsx", 
                        sheet = "Extended Fig5b")

Plotdata<-Plotdata%>%
  mutate(TILs= as.numeric(TILs))%>%
  filter(!is.na(Visit))


#Plotdata$IGG_Signature_Mean<-scales::rescale(Plotdata$IGG_Signature_Mean, to=c(0,5))

plot1<- ggplot2::ggplot(data=Plotdata%>%filter(pCR!="ND"),
                        
                        aes(x=pCR,y=TILs))+
  
  ggplot2::geom_boxplot(aes(fill="lightgrey"),outlier.shape = NA)+
  
  geom_jitter(aes(color="black"), size = 1.5, alpha = 1, shape=3,width = 0.2)+
  
  theme(legend.position = "none")+
  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + theme(legend.position = "none")

ggsave("Extfigure5b.pdf",plot, device=cairo_pdf)

###stats ###
box_core <- Plotdata %>%
  group_by(pCR) %>%
  summarise(
    N = n(),
    lower = quantile(TILs, 0.25) - 1.5 * IQR(TILs),
    upper = quantile(TILs, 0.75) + 1.5 * IQR(TILs),
    lower_whisker = max(min(TILs), quantile(TILs, 0.25) - 1.5 * IQR(TILs)),
    upper_whisker = min(max(TILs), quantile(TILs, 0.75) + 1.5 * IQR(TILs)),
    centre = median(TILs)
  )

box_core
write.table(box_core,"Extfig5b.stats.txt",sep = "\t",row.names = F)

################################ Extended Figure 5c ##################
Plotdata <-  read_excel("SourceData.xlsx", 
                        sheet = "Extended Fig5c")

Plotdata<-Plotdata%>%
  mutate(`PDL1%`= as.numeric(`PDL1%`))%>%
  filter(!is.na(Visit))


#Plotdata$IGG_Signature_Mean<-scales::rescale(Plotdata$IGG_Signature_Mean, to=c(0,5))

plot1<- ggplot2::ggplot(data=Plotdata%>%filter(pCR!="ND"),
                        
                        aes(x=pCR,y=`PDL1%`))+
  
  ggplot2::geom_boxplot(aes(fill="lightgrey"),outlier.shape = NA)+
  
  geom_jitter(aes(color="black"), size = 1.5, alpha = 1, shape=3,width = 0.2)+
  
  theme(legend.position = "none")+
  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + theme(legend.position = "none")

ggsave("Extfigure5c.pdf",plot, device=cairo_pdf)

###stats ###
box_core <- Plotdata %>%
  group_by(pCR) %>%
  na.omit()%>%
  summarise(
    N = n(),
    lower = quantile(`PDL1%`, 0.25) - 1.5 * IQR(`PDL1%`),
    upper = quantile(`PDL1%`, 0.75) + 1.5 * IQR(`PDL1%`),
    lower_whisker = max(min(`PDL1%`), quantile(`PDL1%`, 0.25) - 1.5 * IQR(`PDL1%`)),
    upper_whisker = min(max(`PDL1%`), quantile(`PDL1%`, 0.75) + 1.5 * IQR(`PDL1%`)),
    centre = median(`PDL1%`)
  )

box_core
write.table(box_core,"Extfig5c.stats.txt",sep = "\t",row.names = F)


################################ Extended Figure 5d ##################
Plotdata <-  read_excel("SourceData.xlsx", 
                        sheet = "Extended Fig5d")

Plotdata<-Plotdata%>%
  mutate(`PDL1%`= as.numeric(`PDL1%`))%>%
  filter(!is.na(Visit))


#Plotdata$IGG_Signature_Mean<-scales::rescale(Plotdata$IGG_Signature_Mean, to=c(0,5))

plot1<- ggplot2::ggplot(data=Plotdata%>%filter(pCR!="ND"),
                        
                        aes(x=pCR,y=`PDL1%`))+
  
  ggplot2::geom_boxplot(aes(fill="lightgrey"),outlier.shape = NA)+
  
  geom_jitter(aes(color="black"), size = 1.5, alpha = 1, shape=3,width = 0.2)+
  
  theme(legend.position = "none")+
  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + theme(legend.position = "none")

ggsave("Extfigure5c.pdf",plot, device=cairo_pdf)

###stats ###
box_core <- Plotdata %>%
  group_by(pCR) %>%
  na.omit()%>%
  summarise(
    N = n(),
    lower = quantile(`PDL1%`, 0.25) - 1.5 * IQR(`PDL1%`),
    upper = quantile(`PDL1%`, 0.75) + 1.5 * IQR(`PDL1%`),
    lower_whisker = max(min(`PDL1%`), quantile(`PDL1%`, 0.25) - 1.5 * IQR(`PDL1%`)),
    upper_whisker = min(max(`PDL1%`), quantile(`PDL1%`, 0.75) + 1.5 * IQR(`PDL1%`)),
    centre = median(`PDL1%`)
  )

box_core
write.table(box_core,"Extfig5d.stats.txt",sep = "\t",row.names = F)

######################################### Extended Figure 6 ################
data<-read.table("SCR_IMMUNEMODULESvsPCR.txt", as.is=T, head=T, quote= "\"'", sep="\t")
dim(data)

standardize<-function(x){
  annAll<-dimnames(x)
  x<-scale(x)
  dimnames(x)<-annAll
  return(x)
}

C <- data
B<-standardize(C[,c(5:205)])
dim(B)


B <- data.frame(PCR=C$PCR, B)
dim(B)


prob<-glm(PCR~Bcells.Cluster_Iglesia_CCR.2014_PMID.24916698,data=B,family="binomial")
summary(prob)

colnames(B)
resB  <-  data.frame(matrix(nrow=202,ncol=9))
colnames(resB)=c("gen", "AIC", "Estimate", "Std.Error", "tvalue", "OR", "Low","High", "pvalue")
line=1
for (i in c(2:202)){ 
  resB[line,1] <-  names(B)[i]
  resB[line,2] <-  glm(PCR~B[,i],data=B,family="binomial")$aic
  resB[line,3] <-  glm(PCR~B[,i],data=B,family="binomial")$coef[[2]]
  resB[line,4] <-  summary(glm(PCR~B[,i],data=B,family="binomial"))$coef[[4]]
  resB[line,5] <-  summary(glm(PCR~B[,i],data=B,family="binomial"))$coef[[6]]
  resB[line,6] <-  exp(cbind(OR = coef(glm(PCR~B[,i],data=B,family=binomial)), confint(glm(PCR~B[,i],data=B,family=binomial))))[[2]]
  resB[line,7] <- exp(cbind(OR = coef(glm(PCR~B[,i],data=B,family=binomial)), confint(glm(PCR~B[,i],data=B,family=binomial))))[[4]]
  resB[line,8] <- exp(cbind(OR = coef(glm(PCR~B[,i],data=B,family=binomial)), confint(glm(PCR~B[,i],data=B,family=binomial))))[[6]]
  resB[line,9] <-  summary(glm(PCR~B[,i],data=B,family="binomial"))$coef[[8]]
  line=line+1
}
resB 

write.table(resB, file="output_pCR_univariat.txt", sep="\t", col.names = NA)

install.packages("forestplot")
library(forestplot)

data <- read.table("Forrestplotoutput_pCR_univariat_immunemodules_PRENAC.txt", header = TRUE, sep = "\t")

print(data)

tabletext <- cbind(
  c("Variable", as.character(data$gen)),  # Nombres de las variables
  c("OR (95% CI)", paste0(round(data$OR, 2), " (", round(data$Low, 2), "-", round(data$High, 2), ")")),  # Intervalos de confianza
  c("p-value", round(data$pvalue, 4))  # Valores p
)

if (nrow(tabletext) != nrow(data) + 1) {
  stop("Número de filas en tabletext y en los datos no coincide.")
}

png("Forestplot_PRENAC_Results_Immune_V1.png", width = 5000, height = 2500, res = 300)

par(mar = c(5, 25, 4, 2) + 0.1)

forestplot(
  labeltext = tabletext,  # Incluye la fila de encabezados
  mean = c(NA, data$OR),  # Odds ratio, con NA al inicio para el encabezado
  lower = c(NA, data$Low),  # Límite inferior del intervalo de confianza, con NA al inicio para el encabezado
  upper = c(NA, data$High),  # Límite superior del intervalo de confianza, con NA al inicio para el encabezado
  new_page = TRUE,  # Crear una nueva página
  is.summary = c(TRUE, rep(FALSE, nrow(data))),  # Indica si la primera fila es un resumen
  xlab = "Odds Ratio (Log Scale)",  # Etiqueta del eje X
  zero = 1,  # Línea central en OR = 1
  boxsize = 0.1,  # Tamaño de los cuadros
  col = fpColors(box = "royalblue", line = "darkblue", summary = "royalblue"),  # Colores
  txt_gp = fpTxtGp(label = gpar(fontsize = 12), xlab = gpar(fontsize = 14), ticks = gpar(fontsize = 12)),  # Tamaño del texto
  line.margin = 0.2,  # Margen entre líneas
  xlog = TRUE  # Escala logarítmica en el eje X
)

dev.off()


