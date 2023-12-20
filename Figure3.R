# Bueno data
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [3A] 
rm(list=ls())

myinf1 = #path to gene expression
myinf2 = #path to clinical infomation
myinf3 = #path to PRACME

data = read.table(myinf1, sep="\t", header=T, row.names=1, quote="")
xx = apply(data, 2, sum)
mymed = median(xx)
for(k in 1:ncol(data))
{
  data[,k] = data[,k]*(mymed/xx[k])
}
data = log2(data+1)
mygen = read.table(myinf3, sep="\t", header=T, row.names=1, stringsAsFactors=F, quote="")
comxx = intersect(row.names(data), row.names(mygen))
data = data[comxx,]
mygen = mygen[comxx,]

#Calculating signature score
score = t(data) %*% mygen /length(mygen)

#----------------------
#Processing clinical information
info = read.table(myinf2, sep="\t", header=T, stringsAsFactors=F)
se = which(!is.na(info$RNA.Tumor.ID))
info = info[se,]
row.names(info) = paste("X", info$RNA.Tumor.ID, sep="")

xx = info$Stage
xx = gsub("Stage ", "", xx)
xx[grep("T4", xx)] ="IV"
xx[grep("T3", xx)] ="III"
xx[grep("T2", xx)] ="II"
xx[grep("T1", xx)] ="I"
stage.sim = xx
se = which(stage.sim%in%c("I", "II", "III", "IV", ""))

t.surv = info$Survival.from.surgery
e.surv = ifelse(info$Status=="d", 1, 0)
info = cbind(t.surv, e.surv, stage.sim, info)
row.names(info) = paste("X", info$RNA.Tumor.ID, sep="")

comxx = intersect(row.names(score), row.names(info))
score = score[comxx,]
info = info[comxx, ]

# Proceed with no scale
xx=cbind(score,info)


#Figure 3A (xx includes score here)
#-------------------------------
myfac = xx$Histology.reduced
tmp = table(myfac)
tmp = tmp[tmp>=5]
tmp
# Biphasic Epithelioid Sarcomatoid 
# 62         141           7 

myList = list(NULL)

for(k in 1:length(tmp))
{
  se = which(myfac==names(tmp)[k])
  myList[[k]] = xx$score[se]
}
names(myList) = names(tmp)

wilcox.test(myList[[1]], myList[[2]])
wilcox.test(myList[[1]], myList[[3]])
wilcox.test(myList[[2]], myList[[3]])

## Epithelioid < Biphasic < Sarcomatoid
# B:0.6110972
# E: 0.5942314
#S: 0.6396106
se= which(xx$Histology.reduced=="Desmoplastic")
xx= xx[-se,]
model_fit <- survfit(Surv(t.surv, e.surv)~Histology.reduced, data = xx)
ggsurvplot(
  model_fit,
  data = xx,
  size = 1,                 # change line size
  pval = TRUE,  
  pval.coord = c(0, 0.03),
  palette  =c("#6CB0D6","#F8766D","#EEB479"),
  legend.labs =
    c("Biphasic","Epithelioid","Sarcomatoid"),    # Change legend labels
  ggtheme = theme_classic(),    
  xlab = "\n Survival Time (Years)",
  ylab ="Overall Survival",
)
library(tidyverse)
data <- data.frame(
  name=c( rep("Sarcomatioid",7),rep("Biphasic",62), rep("Epithelioid",141) ),
  sig_score=c( myList$Sarcomatoid, myList$Biphasic, myList$Epithelioid )
)
my_xlab <- paste(levels(as.factor(data$name)),"\n(N=",table(data$name),")",sep="")

#reorder by sig score from high to low
blab=my_xlab[1]
elab=my_xlab[2]
slab=my_xlab[3]

# Plot
data %>%
  ggplot( aes(x=reorder(name,-sig_score), y=sig_score, fill=name)) +
  geom_boxplot(width=0.5/length(unique(data$name))) +
  scale_fill_manual(values=c("#6CB0D6","#F8766D","#EEB479"))+
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_classic() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("")+
  scale_x_discrete(labels=my_xlab)

#141 patients
se= which(xx$Histology.reduced=="Epithelioid")
tmp = xx[se,]
se = which(xx$score>= median(xx$score))
xx$range =xx$stage.sim
xx[se,]$range ="high-value"
xx[-se,]$range ="low-value"
model_fit <- survfit(Surv(t.surv, e.surv)~range, data = xx)

ggsurvplot(
  model_fit,
  data = xx,
  size = 1,                 # change line size
  pval = TRUE,    
  palette  = c("#F8766D","#4d4d4d"),
  legend.labs =
    c("High-value", "Low-value"),    # Change legend labels
  xlab = "\n Survival Time (Years)",
  ylab ="Overall Survival",
)
# [3C]
xx$asbesto = xx$Asbestos.exposure
table(xx$Asbestos.exposure)
levels(as.factor((xx$Asbestos.exposure)))
myList = list(NULL)
se = which(xx$Asbestos.exposure=="yes")
xx$asbesto[se] = "yes"
myList[[1]] = xx$score[se]
se = which(xx$Asbestos.exposure=="none" | xx$Asbestos.exposure=="none known")
xx$asbesto[se] = "no"
myList[[2]] = xx$score[se]
names(myList) = c("yes", "no")

data <- data.frame(
  name=c( rep("Has asbestos exposure",143),rep("No asbestos exposure",53) ),
  sig_score=c( myList$yes, myList$no )
)
my_xlab <- paste(levels(as.factor(data$name)),"\n(N=",table(data$name),")",sep="")


# Plot boxplot
data %>%
  ggplot( aes(x=reorder(name,-sig_score), y=sig_score, fill=name)) +
  geom_boxplot(width=0.5/3) +
  # scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  # palette = c("#6CB0D6","#6CBA7D","#ED85B0")+
  scale_fill_manual(values=c("#F8766D","#4d4d4d"))+
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_classic() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("")+
  scale_x_discrete(labels=my_xlab)

# [3D]
se = which(xx$asbesto=="yes"|xx$asbesto=="no")
yy= xx[se,]
model_fit <- survfit(Surv(t.surv, e.surv)~asbesto, data = yy)
library(ggplot2)
library(ggfortify)
library(survminer)
ggsurvplot(
  model_fit,
  data = yy,
  size = 1,                 # change line size
  pval = TRUE,    
  palette  = c("#F8766D","#4d4d4d"),
  legend.labs =
    c("No Exposure", "Has Exposure"),    # Change legend labels
  xlab = "\n Survival Time (Years)",
  ylab ="Overall Survival",
)
levels(yy$asbesto)
wilcox.test(myList[[1]], myList[[2]])

####################################################################################
#Immune infiltration 
rm(list=ls())
myinf1 = #path to bueno score
myinf2 = #path to immune infiltration

score = read.table(myinf1, sep=",", header=T)
row.names(score)=score$X

data = read.table(myinf2, sep="\t", header=T)
cnum = ncol(data)/4
data = data[, 1:cnum] - data[, (1+cnum):(cnum*2)]
tmp = colnames(data)
tmp = gsub("_up\\.ES", "", tmp)
colnames(data) = tmp

comxx = intersect(row.names(data), row.names(score))
data = data[comxx,]
score = score[comxx,]

cor(data, score, method="s")
df=cbind(data,score)
#scatterplot
ggplot(df, aes(x=df$x, y=df$Monocyte), pval.coord = c(0, 0.03)) + 
  geom_point(size=3) +
  geom_smooth(method=lm , color="red", se=FALSE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                         panel.background = element_blank(), axis.line = element_line(colour = "black"), 
                                                         axis.text.x = element_text(angle=90, hjust=1))+stat_cor(r.digits  = 2)


# tcga----------------------------

rm(list=ls())
myinf1 = #path to tcga score
myinf2 = # path to tcga infiltration

score = read.table(myinf1, sep="\t", header=T)


data = read.table(myinf2, sep="\t", header=T)
cnum = ncol(data)/4
data = data[, 1:cnum] - data[, (1+cnum):(cnum*2)]
tmp = colnames(data)
tmp = gsub("_up\\.ES", "", tmp)
colnames(data) = tmp

comxx = intersect(row.names(data), row.names(score))
data = data[comxx,]
score = score[comxx,]

df=cbind(data,score)
#scatterplot
ggplot(df, aes(x=df$PRACME, y=df$Monocyte), pval.coord = c(0, 0.03)) + 
  geom_point(size=3) +
  geom_smooth(method=lm , color="red", se=FALSE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                         panel.background = element_blank(), axis.line = element_line(colour = "black"), 
                                                         axis.text.x = element_text(angle=90, hjust=1))+stat_cor(r.digits  = 2)
################################################################################################
# [3.4] correltion  PRACME/ORACLE score with Immune infiltration  -- Thronson 
rm(list=ls())
myinf1 = #path to PRACME/ORACLE SCORE
myinf2 = #path to Throson table

score = read.table(myinf1, sep="\t", header=T)

imm= read.table(myinf2, sep=",", header=T, row.names=1, stringsAsFactors=F)
se = which(imm$TCGA.Study=="MESO")
data = imm[se,]

comxx = intersect(row.names(data), row.names(score))
data = data[comxx,]
score = score[comxx,]
dim(data)		## n=87 samples

data = data[, c(4:31, 36:54)]

xx=cor(data, score, method="s", use="pair")

xx1 = xx[order(apply(xx,1, max), decreasing=T), ]
xx2 = xx[order(apply(xx,1, max), decreasing=F), ]
xx2


#[3G] correlation with Immune genes
rm(list=ls())

myinf1 = #path to gene expresison
myinf2 = #path to score


imm.inh = c("CTLA4", "PDCD1", "LAG3", "BTLA", "CD160", "IDO1", "IL10", "IL10RB", "TGFB1", "TGFBR1", "VTCN1", "CD244", "LGALS9", "HAVCR2", "ADORA2A", "TIGIT", "CSF1R", "KIR2DL1", "KIR2DL2", "KIR2DL3", "KDR", "CD96", "PVRL2", "C10orf54")
imm.sti = c("MICA", "MICB", "CD27", "CD274", "CD28", "CD40", "CD40LG", "CD70", "CD80", "CD86", "ICOS", "ICOSLG", "IL6", "IL6R", "PDCD1LG2", "TMEM173", "TNFRSF13B", "TNFRSF13C", "TNFRSF14", "TNFRSF17", "TNFRSF18", "TNFRSF4", "TNFRSF9", "TNFSF13", "TNFSF13B", "TNFSF18", "TNFSF4", "TNFSF9", "TNFSF15", "TNFRSF25", "HHLA2", "TMIGD2", "BTNL2", "CD276", "CD48", "TNFSF14", "TNFRSF8", "PVR", "LTA",  "IL2RA", "ENTPD1", "NT5E", "CXCR4", "CXCL12", "KLRK1", "NKG2A", "RAET1E", "ULBP1")
imm.oth = c("GZMA", "PRF1")
mygen = c(imm.inh, imm.sti, imm.oth)
# mygen ="C10orf54"
#--------------
load(myinf1)
data = mydata
dim(data)
xx = apply(data>0,1, sum)
se = which(xx>ncol(data)*0.1)
length(se)
data = data[se,]
data = log2(data+1)

se = which(row.names(data)%in%mygen)
data = data[se,]
data = t(data)
#--------------
score = read.table(myinf2, sep="\t", header=T)

comxx = intersect(row.names(data), row.names(score))
data = data[comxx,]
score = score[comxx,]
dim(data)		## n=87 samples


xx = cor(data, score, method="s", use="pair")
xx1 = xx[order(apply(xx,1, max), decreasing=T), ]
xx2 = xx[order(apply(xx,1, max), decreasing=F), ]
xx2[1:5,]
data = as.data.frame(data)
library(ggplot2)
library(dplyr)
library(ggpubr)
df = cbind(data$C10orf54,score$PRACME)
# df.vista= as.data.frame(df$C10orf54)
colnames(df) = c("gexp", "pracme")
df = as.data.frame(df)
ggplot(df, aes(x=pracme, y=gexp), pval.coord = c(0, 0.03)) + 
  geom_point(size=3) +
  geom_smooth(method=lm , color="red", se=FALSE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                         panel.background = element_blank(), axis.line = element_line(colour = "black"), 
                                                         axis.text.x = element_text(angle=90, hjust=1))+stat_cor(r.digits  = 2)

#               ORACLE       PRACME
# CD276      0.566432164  0.552635416			***
#   TGFB1      0.560180798  0.490632062			***
#   TNFSF4     0.440001458  0.508055697
# TGFBR1     0.307538091  0.432164467			***
#   PVR        0.361777357  0.428464679
# KDR       -0.444011081 -0.326802508
# IL6R      -0.371400452 -0.379565503
# TNFSF13   -0.439035503 -0.565375082
# C10orf54  -0.551396078 -0.577932493			***  VISTA (associated with good prognosis)
# RAET1E    -0.591783918 -0.576492673

colnames(data)

# N1  N2         diff     tscore         t.P         w.P
# BAP1  46 151  0.006681115  0.9400833 0.350008674 0.349840732
# NF2   37 160  0.022881411  2.7992970 0.007175163 0.002458266
# SETD2 17 180 -0.012943562 -1.3993111 0.176266226 0.181111724
# TP53  16 181  0.002274205  0.2371993 0.814994567 0.943471599


#VSIR patient prognosis
rm(list=ls())

myinf1 = #path to gene expression
myinf2 = #path to score
# mygen ="C10orf54"
mygen ="CD276"
#--------------
load(myinf1)
data = mydata
dim(data)
xx = apply(data>0,1, sum)
se = which(xx>ncol(data)*0.1)
length(se)
data = data[se,]
data = log2(data+1)


se = which(row.names(data)%in%mygen)
data = data[se,]
data = as.data.frame(t(data))



library(ggplot2)
library(dplyr)
library(ggpubr)

#----------------------
myinf3 = #path to clinical info
info = read.table(myinf3, sep="\t", header=T, row.names=1, stringsAsFactors=F, quote="")
se = c("days_to_birth", "gender",  "stage_event.pathologic_stage", "vital_status","days_to_death", "days_to_last_followup", "history_asbestos_exposure")
info = info[,se]
colnames(info)[1:3] =c("age", "gender", "stage")
xx = ifelse(!is.na(info$days_to_death), info$days_to_death, info$days_to_last_followup)
t.surv = as.numeric(xx)
e.surv = ifelse(info[, "vital_status"]=="dead", 1, 0)
info = cbind(t.surv, e.surv, info)
info = info[!is.na(info$t.surv), ]
xx = info$stage
xx = gsub("stage ", "", xx)
xx[grep("iv", xx)] ="IV"
xx[grep("iii", xx)] ="III"
xx[grep("ii", xx)] ="II"
xx[grep("i", xx)] ="I"
info$stage = xx
info$age = -info$age/365
se = which(info$t.surv>=30)
info = info[se,]

comxx = intersect(row.names(data), row.names(info))
info = info[comxx, ]
data = data[comxx,]
xx = cbind(data, info)


library(survival)
mycox = coxph(Surv(t.surv/365, e.surv)~data, data = xx) 
mycox
# coef exp(coef) se(coef)      z        p
# data -0.4295    0.6508   0.1163 -3.694 0.000221
# 
# Likelihood ratio test=12.91  on 1 df, p=0.0003276
# n= 84, number of events= 57 
xx$range =xx$data
q= quantile(data, probs =c(1/3,2/3))
se = which(xx$data<=median(xx$data))
xx[se,]$range ="low-value"
xx[-se,]$range ="high-value"


library(ggplot2)
library(ggfortify)
library(survminer)
model_fit <- survfit(Surv(t.surv/365, e.surv)~range, data = xx)
ggsurvplot(
  model_fit,
  data = xx,
  size = 1,                 # change line size
  pval = TRUE,  
  pval.coord = c(0, 0.03),
  # Add p-value
  palette  = c("#F8766D","#4d4d4d"),
  legend.labs =
    c("High-expression","Low-expression"),    # Change legend labels
  xlab = "\n Survival Time (Years)",
  ylab ="Overall Survival",
)


#[3.6] correlation with GSVA scores (pathways)
rm(list=ls())
myinf1 = #path to gene expression
myinf2 = #../c2.all.v4.0.symbols.gmt
myinf3 = #path to score
myinf4 = #path to gene expression

#--------------
conIn = file(myinf2, "r")
data = readLines(conIn)
close(conIn)
se = grep("KEGG_|REACTOME_|PID_|BIOCARTA_", data)
data = data[se]

myList = list(NULL)
for(k in 1:length(data))
{
  xx = unlist(strsplit(data[k], "\t"))
  myList[[k]] = xx[3:length(xx)]
  names(myList)[k] = xx[1]
}


#--------------
load(myinf1)
data = mydata
dim(data)
xx = apply(data>0,1, sum)
se = which(xx>ncol(data)*0.1)
length(se)
data = data[se,]
data = log2(data+1)


data = as.matrix(data)
xx = apply(data, 1, mean)
data = data-xx



library(GSVA)
library(GSEABase)
library(GSVAdata)


gs.list = myList
gsva.es <- gsva(data, gs.list, verbose=FALSE, min.sz=20, max.sz=500)
data = t(gsva.es)

gsva.df= as.data.frame(gsva.es)
write.csv(as.data.frame(gsva.es),"pathwayq2.csv")
se=which(gsva.es>0.5)

#---------------
#each quadrant
mpm = read.table("/home/ylin/PRACME/pathway/MPM_with_quadrant.txt")
se=which(mpm$quadrant=="Q2")
se=intersect(row.names(data),row.names(mpm[se,]))
data=data[se,]
data = as.matrix(data)
xx = apply(data, 1, mean)
data = data-xx


library(GSVA)
library(GSEABase)
library(GSVAdata)


gs.list = myList
gsva.es <- gsva(data, gs.list, verbose=FALSE, min.sz=20, max.sz=500)
data = t(gsva.es)

gsva.df= as.data.frame(gsva.es)
write.csv(as.data.frame(gsva.es),"pathwayq2.csv")
se=which(gsva.es>0.5)


#--------------
score = read.table(myinf3, sep="\t", header=T)

comxx = intersect(row.names(data), row.names(score))
data = data[comxx,]
score = score[comxx,]
dim(data)		## n=87 samples

xx = cor(data, score, method="s", use="pair")
xx1 = xx[order(apply(xx,1, max), decreasing=T), ]
xx2 = xx[order(apply(xx,1, max), decreasing=F), ]

