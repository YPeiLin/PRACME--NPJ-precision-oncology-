#Shedden
# txt = "ca00182.HGU133A_EntrezCDF.MAS5.pcl"
# tsv= "ca00182.info.tsv"

# Okayama
# txt = "GSE31210.HGU133Plus2_EntrezCDF.MAS5.pcl.txt"
# tsv= "GSE31210.info.tsv"


################################################################################
rm(list=ls())
txt = #txt file name
tsv= #tsv file name
myinf1 = #path to gene expression
myinf2 = #path to info
myinf3 = #path to PRACME
myinf4 =  #path to ORACLE


data = read.table(myinf1, sep="\t", header=T, row.names=1, quote="")
xx = apply(data, 2, sum)
mymed = median(xx)
for(k in 1:ncol(data))
{
  data[,k] = data[,k]*(mymed/xx[k])
}
data = log2(data+1)

data <- na.omit(data)
mygen = read.table(myinf3, sep="\t", header=T, row.names=1, stringsAsFactors=F, quote="")

comxx = intersect(row.names(data), row.names(mygen))
data = data[comxx,]
mygen = mygen[comxx,]

score = t(data) %*% mygen /length(mygen)

library(survival)
myinf2 = #updated clinical info

info = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
t.surv = as.numeric(info$days.before.relapse.censor)
e.surv = ifelse(info$relapse ==" relapsed", 1, 0)
info = cbind(t.surv, e.surv, info)
tag = !is.na(info$t.surv)
info = info[tag==1,]
dim(info)

comxx = intersect(row.names(score), row.names(info))
score = score[comxx,]
info = info[comxx, ]
xx = cbind(score, info)

scale.score = (score -mean(score))/sd(score)
xx = cbind(scale.score, info)
mycox <-coxph(Surv(as.numeric(t.surv), as.numeric(e.surv))~scale.score, data = xx)
tmp =summary(mycox)

#kmplot
library(ggplot2)
library(ggfortify)
library(survminer)

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
#km plot
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


yy=cbind(t.surv,e.surv,normalized.score,info.2)
yy=xx

yy <- within(yy, {
  "Age" = age
  "Stage"=pathological.stage
  "Sex" =gender
  "Smoke" = smoking.status
})


mycox = coxph(Surv(t.surv, e.surv)~scale.score+Age+Stage+Sex+Smoke, id=Title, data = yy) 
summary(mycox)
print(ggforest(mycox,data=yy))


##################################################################################
#Shedden

rm(list=ls())
txt = #txt file name
tsv= #tsv file name
myinf1 = #path to gene expression
myinf2 = #path to info
myinf3 = #path to PRACME
myinf4 =  #path to ORACLE

data = read.table(myinf1, sep="\t", header=T, row.names=1, quote="")
xx = apply(data, 2, sum)
mymed = median(xx)
for(k in 1:ncol(data))
{
  data[,k] = data[,k]*(mymed/xx[k])
}
data = log2(data+1)

data <- na.omit(data)
mygen = read.table(myinf3, sep="\t", header=T, row.names=1, stringsAsFactors=F, quote="")

comxx = intersect(row.names(data), row.names(mygen))
data = data[comxx,]
mygen = mygen[comxx,]

score = t(data) %*% mygen /length(mygen)

library(survival)
info = read.table(myinf2, sep="\t", header=T, stringsAsFactors=F)
se  <- which(!is.na(info$OS_Time)&!is.na(info$OS_Status))
info = info[se,]
# 
t.surv = info$OS_Time/12
e.surv = info$OS_Status


info = cbind(t.surv, e.surv, info)
row.names(info) = info$Array

comxx = intersect(row.names(score), row.names(info))
score = score[comxx,]
info = info[comxx, ]
xx = cbind(score, info)

scale.score = (score -mean(score))/sd(score)
xx = cbind(scale.score, info)
mycox <-coxph(Surv(as.numeric(t.surv), as.numeric(e.surv))~scale.score, data = xx)
tmp =summary(mycox)

#kmplot
library(ggplot2)
library(ggfortify)
library(survminer)

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
#km plot
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

#Multivariate
myinf5= #path to clinical info

info.2 = read.table(myinf5,sep = "\t")
se  <- which(!is.na(info.2$days.before.death.censor)&!is.na(info.2$death))
info.2 = info.2[se,]
library(scales)
normalized.score =rescale(score, to = c(0, 1))
yy=cbind(normalized.score,info)


yy <- within(yy, {
  "Age" = age
  "Stage"=pathological.stage
  "Sex" =gender
  "Smoke" = smoking.status
})
mycox = coxph(Surv(t.surv, e.surv)~scale.score+Age+Stage+Sex+Smoke, data = yy) 

print(ggforest(mycox,data=yy))
summary(mycox)
