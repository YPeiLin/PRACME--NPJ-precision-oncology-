# Bueno data 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2A] 
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
# xx=cbind(score,info)

# Scale method 1; used for univariate
scale.score = (score -mean(score))/sd(score)
xx = cbind(scale.score, info)

# Scale method 2 (0-100); used for multivariate analysis
# library(scales)
# normalized.score =rescale(score, to = c(0, 100))

# Scale method 3 (0-1)
# normalized.score = (x-min(x))/(max(x)-min(x))
# xx = cbind(normalized.score,info)

library(survival)
library(ggplot2)
library(ggfortify)
library(survminer)

mycox = coxph(Surv(t.surv, e.surv)~scale.score, data = xx) 
summary(mycox)


#Divide risk score into 3 groups
scale.score = (score -mean(score))/sd(score)
xx = cbind(scale.score, info)

xx$range =xx$scale.score
q= quantile(scale.score, probs =c(1/3,2/3))
se = which(xx$scale.score<=q[1])
xx[se,]$range ="low-value"
se = which(xx$scale.score>=q[2])
xx[se,]$range ="high-value"
se = which(xx$scale.score<q[2]&xx$scale.score>q[1])
xx[se,]$range ="intermediate-value"

model_fit <- survfit(Surv(t.surv, e.surv)~range, data = xx)
ggsurvplot(
  model_fit,
  data = xx,
  size = 1,                 # change line size
  pval = TRUE,  
  pval.coord = c(0, 0.03),
  palette  = c("#F8766D","#7CAE00","#4d4d4d"),
  legend.labs =
    c("High-value","Intermediate-value","Low-value"),    # Change legend labels
  ggtheme = theme_classic(),    
  xlab = "\n Survival Time (Years)",
  ylab ="Overall Survival",
)


###############################################################################
# 2C: Bueno Multivariate
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(survival)
library(forestmodel)
library(scales)
library(survminer)

normalized.score =rescale(score, to = c(0, 100))
xx=cbind(info,scale.score)

#183 patients
se = which(xx$stage.sim%in%c("I","II","III", "IV"))
yy = xx[se, ]
#remove desmoplastic
yy= yy[-which(yy$Histology.reduced=="Desmoplastic"),]
yy$stage.sim = as.factor(as.character(yy$stage.sim))
se= which(yy$Asbestos.exposure=="39 years known")
yy[se,]$Asbestos.exposure="yes"
se= which(yy$Asbestos.exposure=="none known" |yy$Asbestos.exposure=="none mentioned" )
yy[se,]$Asbestos.exposure="none"
se = which(yy$Asbestos.exposure=="yes" | yy$Asbestos.exposure=="none")
yy=yy[se,]

yy <- within(yy, {
  "Age" = Age.at.surgery
  "Stage"=stage.sim
  "Histology" = Histology.reduced
})

yy=within(yy, Histology <- relevel(as.factor(Histology), ref = "Epithelioid"))
mycox = coxph(Surv(t.surv, e.surv)~score+Age+Stage+Sex+Histology+Asbestos.exposure, data = yy) 
summary(mycox)

print(ggforest(mycox,data=yy))
summary(mycox)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [2B] Bott_GSE29354 data
rm(list=ls())

myinf1 = #path to gene expression
myinf2 = #path to clinical infomation
myinf3 = #path to PRACME
  
data = read.table(myinf1, sep="\t", header=T, row.names=1, quote="")
mygen = read.table(myinf3, sep="\t", header=T, row.names=1, stringsAsFactors=F, quote="")

comxx = intersect(row.names(data), row.names(mygen))
data = data[comxx,]
mygen = mygen[comxx,]
dim(data)

score = t(data) %*% mygen /length(mygen)


#----------------------
info = read.table(myinf2, sep="\t", header=T, row.names=1, stringsAsFactors=F)
#record in months
t.surv = as.numeric(info$Time.from.surgery.to.last.f)/12
e.surv = as.integer(info$Dead)
info = cbind(t.surv, e.surv,  info)
row.names(info) = gsub(" ", "", row.names(info))
xx = info$Stage
xx[grep("1", xx)] = 1
info$stage = xx

comxx = intersect(row.names(score), row.names(info))
score = score[comxx,]
info = info[comxx, ]


library(survival)
library(ggplot2)
library(ggfortify)
library(survminer)
library(scales)

# Scale method 1; used for univariate
scale.score = (score-mean(score))/sd(score)
xx= cbind(scale.score,info)

# Scale method 2 (0-100); used for multivariate analysis
# normalized.score =rescale(score, to = c(0, 100))
# xx = cbind(normalized.score, info)

mycox = coxph(Surv(t.surv, e.surv)~scale.score, data = xx) 
summary(mycox)

#KM with 3 groups
library(ggplot2)
library(ggfortify)
library(survminer)

xx= cbind(score,info)
xx$range =xx$score
q= quantile(score, probs =c(1/3,2/3))
se = which(xx$score<=q[1])
xx[se,]$range ="low-value"
se = which(xx$score>=q[2])
xx[se,]$range ="high-value"
se = which(xx$score<q[2]&xx$score>q[1])
xx[se,]$range ="intermediate-value"
model_fit <- survfit(Surv(t.surv, e.surv)~range, data = xx)

ggsurvplot(
  model_fit,
  data = xx,
  size = 1,                 # change line size
  pval = TRUE,       
  pval.coord = c(0, 0.03),
  #default color palette
  palette  = c("#F8766D","#7CAE00","#4d4d4d"),
  legend.labs =
    c("High-value", "Intermediate-Value" ,"Low-value"),    # Change legend labels
  xlab = "\n Survival Time (Years)",
  ylab ="Overall Survival",
)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#[2D] Use normalzied score for multivariate
# normalized.score =rescale(score, to = c(0, 1))
xx = cbind(scale.score,info)

se = which(xx$Stage%in%c(1,2,3,4))
yy = xx[se, ]
yy=within(yy, Histology <- relevel(as.factor(Histology), ref = "E"))
se= which(yy$Asbestos=="UNK")
yy=yy[-se,]

#Bott multivariate

# yy$normalized.score
library(gridExtra)
library(ggplot2)

yy <- within(yy, {
  "Stage"=stage
})
mycox = coxph(Surv(t.surv, e.surv)~scale.score+Age+stage+Sex+Histology+Asbestos, data = yy) 
mycox
summary(mycox)
print(ggforest(mycox,data=yy))

