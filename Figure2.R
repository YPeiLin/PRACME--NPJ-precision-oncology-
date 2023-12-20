#  define gene signature
rm(list=ls())

myinf1 = #path to gexp
myinf2 = #path to prognostic genes
myinf3 = #path to calculated intra/inter ITH
myinf4 = #path to clinical info
myoutf1 = #output path

load(myinf1)
data = mydata
xx = apply(data, 1, mean)
se = which(xx>median(xx))
data = data[se,]
dim(data)		## S1: 10250 genes

#----------------------
#COX-regression result
comxx = intersect(row.names(data), row.names(prog))
data = data[comxx,]
prog = prog[comxx,]

se = which(prog$adj.pval<0.01)
data = data[se,]
dim(data)			## S2: 1667


##--
div.in = info$within.sam.diversity
div.ou = info$between.sam.diversity
cce = info$concordance.coeff
se = which(div.in<=mean(div.in, na.rm=T) & div.ou>=mean(div.ou, na.rm=T))

#----------------------
info = read.table(myinf3, sep="\t", header=T,  row.names=1, quote="")
q4.gene = row.names(info)[se]
se = which(row.names(data)%in%q4.gene)
data = data[se, ]
dim(data)			## S3: 323 
#---------
xx = info$concordance.coeff
hist(xx, br=20)
se = which(xx>=0.2)
xx = row.names(info)[se]
se = which(row.names(data)%in%xx)
length(se)

data = data[se, ]
dim(data)		## S4: 118 157

##########################################################################
#----------------------
#clinical data
info = read.table(myinf4, sep="\t", header=T, row.names=1, stringsAsFactors=F, quote="")
se = c("days_to_birth", "gender",  "stage_event.pathologic_stage", "vital_status","days_to_death", "days_to_last_followup", "history_asbestos_exposure")
info = info[,se]
colnames(info)[1:3] =c("age", "gender", "stage")
# ifelse(test, yes, no)
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

comxx = intersect(colnames(data), row.names(info))
data = data[,comxx]
info = info[comxx, ]

data = log2(data+1)
xx = cbind(info[, 1:5], t(data))


library(survival)
library(glmnet)
myx = xx[, 6:ncol(xx)]

fit <- glmnet(myx, Surv(xx$t.surv, xx$e.surv), family = "cox")
plot(fit)

plot(fit, "lambda", label=TRUE)
lamda = c(0.001*(1:100))
count = rep(0, length(lamda))


for(k in 1:length(lamda))
{
  tmp = coef(fit, s = lamda[k])
  count[k] = sum(tmp[,1]!=0)
  
}
data.frame(lamda, count)
tmp = coef(fit, s = 0.06)
se = which(tmp[,1]!=0)
tmp = tmp[se,]

# res
res = data.frame(gene=names(tmp), coeff= tmp)
row.names(res) = NULL
res


## S5: 29 genes 
write.table(res, myoutf1, sep="\t", row.names=F, quote=F)

#----------------------
# Other approaches
# myx = xx[, 6:ncol(xx)]
# 
# library(tidyverse)
# library(caret)
# library(leaps)
# 
# library(MASS)
# myx$survival = Surv(xx$t.surv, xx$e.surv)
# 
# full.model <- coxph(survival  ~., data = myx)
# # Stepwise regression model
# step.model <- stepAIC(full.model, direction = "both", 
#                       trace = FALSE)
# step.model.2 <- stepAIC(full.model, direction = "forward", 
#                         trace = FALSE)
# summary.both= summary(step.model)
# #50
# length(which(summary.both$coefficients[,5]<0.01))
# gene.both = colnames(myx[which(summary.both$coefficients[,5]<0.01),])
# #78
# summary.forward = summary(step.model.2)
# length(which(summary.forward$coefficients[,5]<0.01))
# gene.forward = colnames(myx[which(summary.forward$coefficients[,5]<0.01),])
# 
# # model <- lm(Surv(xx$t.surv, xx$e.surv)~myx )
# # ols_step_forward_p(model, details = TRUE)
# 
# 
# 
# # all 29 genes intersect with forward-selection and lasso regression
# length(intersect(rownames(mygen), gene.forward))
######################################################################
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

