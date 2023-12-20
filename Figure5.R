#Forest plot 5A 5B
################################################################################
# TCGA
# Univariate
# exp(coef) exp(-coef) lower .95 upper .95
# normalized.score     1.122     0.8913     1.091     1.153

#TCGA-PRACME_Multi
# exp(coef) exp(-coef) lower .95 upper .95
# normalized.score                1.1288     0.8859    1.0922     1.167
# gendermale                      0.7337     1.3629    0.3291     1.636
# age                             1.0291     0.9717    0.9961     1.063
# history_asbestos_exposureyes    1.0176     0.9827    0.4956     2.089
#Univariate ORACLE
# exp(coef) exp(-coef) lower .95 upper .95
# normalized.score     1.041     0.9608     1.027     1.055
##TCGA-ORACLE-Multi
# exp(coef) exp(-coef) lower .95 upper .95
# normalized.score                 1.042     0.9599    1.0254     1.058
# gendermale                       1.551     0.6447    0.7176     3.353
# age                              1.040     0.9612    1.0010     1.081
# history_asbestos_exposureyes     1.121     0.8924    0.5536     2.268
################################################################################
# Bueno
# after rescale from 1-100

# Uni-PRACME
# exp(coef) exp(-coef) lower .95 upper .95
# scale.score     1.427     0.7009     1.229     1.657
#Multi PRACME
# exp(coef) exp(-coef) lower .95 upper .95
# score                    1.367     0.7315    1.1477     1.628

#Uni- ORACLE
# exp(coef) exp(-coef) lower .95 upper .95
# scale.score      1.63     0.6136     1.391      1.91
#---------
# Multi
# exp(coef) exp(-coef) lower .95 upper .95
# exp(coef) exp(-coef) lower .95 upper .95
# score                   1.6061     0.6226    1.3119     1.966
################################################################################
# bott
#Univariate
# exp(coef) exp(-coef) lower .95 upper .95
# scale.score     1.429     0.6996     1.069     1.912
------
#Multivariate
#   exp(coef) exp(-coef) lower .95 upper .95
#   exp(coef) exp(-coef) lower .95 upper .95
# scale.score    1.6109     0.6208    1.1067     2.345
  
  
#ORACLE
#Univariate
# exp(coef) exp(-coef) lower .95 upper .95
# scale.score     1.513      0.661     1.125     2.034
#Multivariate
#   exp(coef) exp(-coef) lower .95 upper .95
# 2.077     0.4814    1.2612     3.421


#MESO data plot
datasrc =c(rep("TCGA-Univariable",2), rep("TCGA-Multivariable",2),rep("Bueno-Univariable",2),rep("Bueno-Multivariable",2),rep("Bott-Univariale",2),rep("Bott-Multivariale",2) )
sig = (rep(cbind("PRACME","ORACLE"),6))
# sig = (rep(cbind("PRACME-Uni","ORACLE-Uni", "PRACME-Multi", "PRACME_Multi"),3))
hr = c(1.12, 1.04, 1.1288 , 1.042,  1.43,1.63,1.37,1.6061 ,   1.43, 1.51,  1.6109 , 2.077 )
ll = c(1.091, 1.027 ,1.0922, 1.0254 , 1.229, 1.391, 1.1477 , 1.3119,     1.069, 1.125, 1.1067   ,1.2612)
ul = c(1.153, 1.055,1.167, 1.058,  1.657,1.942, 1.638, 1.966,      1.912, 2.034, 2.345, 3.421)
dat <- cbind(datasrc, as.numeric(hr), as.numeric(ll), as.numeric(ul), sig)
dat = as.data.frame(dat)
# , ymin = LowerLimit, ymax = UpperLimit
p = ggplot(data=dat,
           aes(x = sig,y = as.numeric(hr) , ymin = ll, ymax = ul))+
  geom_pointrange(aes(col=sig))+
  geom_hline(yintercept =1, linetype=2)+
  xlab('')+ ylab("Hazard Ratio (95% Confidence Interval)")+
  geom_errorbar(aes(col=sig, ymin = ll, ymax = ul),width=0.5,cex=1)+ 
  facet_wrap(~datasrc,strip.position="right",nrow=9,scales = "free_y") +
  scale_fill_manual(values=c("black","red"))+
  theme(plot.title=element_text(size=16,face="bold"),
        # panel.border =
        panel.background = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(size=6,hjust=0,vjust = 1,angle=360,face="bold"))+scale_y_continuous(breaks=c(0.8, 1, 1.5, 2, 2.5))+coord_flip()  
# +scale_y_discrete(breaks=1:2, labels=1:2)
# +scale_y_discrete(guide = guide_axis(check.overlap = TRUE))
p +scale_color_manual(values=c("#4d4d4d","#F8766D"))

#okayama
# PRACME
# exp(coef) exp(-coef) lower .95 upper .95
# normalized.score      1.0133     0.9869    0.9947     1.032
#
# exp(coef) exp(-coef) lower .95 upper .95
# normalized.score      1.0164     0.9838    0.9966     1.037

#Shedden
#Pracme
# exp(coef) exp(-coef) lower .95 upper .95
# normalized.score         9.984e-01  1.002e+00    0.9902     1.007
# Oracle
# exp(coef) exp(-coef) lower .95 upper .95
# normalized.score         1.002e+00  9.980e-01   0.99466     1.009

#LUNG data plot
datasrc =c(rep("Shedden-Uni",2), rep("Shedden-Multi",2),rep("Okayama-Uni",2),rep("Okayama-Multi",2),rep("Lee-Uni",2) )
sig = (rep(cbind("PRACME","ORACLE"),5))
hr = c(1.268, 1.2,0.9984,1.002,         1.654, 1.602, 1.33, 1.0164 , 1.385  ,	1.901  )
ll = c(1.109, 1.109, 0.9902,0.99466,    1.239, 1.209, 1.02,  0.9966,  1.05     , 1.4     )
ul = c(1.449,	1.364,  1.007, 1.009,     2.209, 2.123, 1.7,  1.037,1.827,	2.581 )
dat.2<- cbind(datasrc, as.numeric(hr), as.numeric(ll), as.numeric(ul), sig)
dat.2 = as.data.frame(dat.2)
# , ymin = LowerLimit, ymax = UpperLimit
p = ggplot(data=dat.2,
           aes(x = sig,y = hr , ymin = ll, ymax = ul))+
  geom_pointrange(aes(col=sig))+
  geom_hline(yintercept =1, linetype=2)+
  xlab('')+ ylab("Hazard Ratio (95% Confidence Interval)")+
  geom_errorbar(aes(col=sig, ymin = ll, ymax = ul),width=0.5,cex=1)+ 
  facet_wrap(~datasrc,strip.position="right",nrow=9,scales = "free_y") +
  scale_fill_manual(values=c("black","red"))+
  theme(plot.title=element_text(size=16,face="bold"),
        # panel.border =
        panel.background = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(size=6,hjust=0,vjust = 1,angle=360,face="bold"))+scale_y_continuous(breaks=c(0.8, 1, 1.5, 2, 2.5))+coord_flip()  

# +scale_y_discrete(breaks=1:2, labels=1:2)
p +scale_color_manual(values=c("#4d4d4d","#F8766D"))


################################################################################
# VennDiagram 5[C-E]
library(VennDiagram)

# move to new plotting page
grid.newpage()


totalpop<- 16286+20267
o.q4 = 1080+15206
m.q4 = 5495+35039
1-phyper(178,m.q4,totalpop - m.q4, o.q4)
1-phyper(178,5495,totalpop - 5495, 1080)
1-phyper(178,1080,totalpop - 1080, 5495)
phyper(178,1080,totalpop - 1080, 5495)


## how to calculate enrichement ration and p-value based on hypergenometric test (Fisher's exact test)
total = 17724
m.q4 = 2314
l.q4 = 1086
shared = 275
## Method 1
white = m.q4   
black = total - white	
take = l.q4   

## Method 1
## p-value is the probability of observing 1000 pr more white ball b chance
tmp=phyper(shared:take, white, black, take )
pval1 = sum(dhyper(shared:take, white, black, take ))
pval1
## Method 2:
tmp = c(shared, white-shared, take-shared, black-(take-shared))
tmp = matrix(tmp, 2, 2)
pval2 = fisher.test(tmp, alternative="g")$p.value
pval1
pval2
## note the two p-values are the same
################################################################################
# Figure 5F CNV correlation
#Hierarchical clustering in supplementary
rm(list=ls())
myinf3 = #CNV info
data = read.table(myinf3, sep="\t", header=T, row.names=1, quote="")
dim(data)

all_tmp = sapply(strsplit(row.names(data), "_"), function(x) x[[1]])
se=which(all_tmp=="MESO")
meso.df =data[se,]
xx = substr(row.names(meso.df), 33, 34)
se = which(xx=="01")
meso.df=meso.df[se,]
#172

amp.thr = log2(2.8/2)
del.thr = log2(1.4/2)



se=which(all_tmp=="LUAD")
luad.df = data[se,]
nchar(row.names(luad.df[1,]))
xx = substr(row.names(luad.df), 33, 34)
se = which(xx=="01")
luad.df=luad.df[se,]


amp.f = apply(luad.df>amp.thr, 2, sum)
del.f = apply(luad.df<del.thr, 2, sum)
# divided by number of samples is to calculate the frequency of loss/gain
amp.f = amp.f/nrow(luad.df)
del.f = del.f/nrow(luad.df)
amp1 = amp.f
del1 = del.f

amp.f = apply(meso.df>amp.thr, 2, sum)
del.f = apply(meso.df<del.thr, 2, sum)
amp.f = as.data.frame(amp.f/nrow(meso.df))
del.f = as.data.frame(del.f/nrow(meso.df))

amp2 = amp.f
del2 = del.f

data.amp <- data.frame(
  LUAD=amp1,
  MESO=amp2
)
colnames(data.amp)[1]="LUAD"
colnames(data.amp)[2]="MESO"

ggplot(data.amp) +
  geom_point( aes(x=LUAD, y=MESO), color="black", size=2 ) +
  labs(
    color="groups"
  )+
  theme_classic(
  ) +
  xlab("Frequency of chr-band amplified in LUAD") +
  ylab("Frequency of  chr-band amplified in MESO")
cor.test(data.amp$LUAD,data.amp$MESO)
# 0.7750603 p-value < 2.2e-16
cor.test(data.amp$LUAD,data.amp$MESO)
#Deletion Frequency
data.del <- data.frame(
  value1=del1,
  MESO=del2
)
colnames(data.del)[1]="LUAD"
colnames(data.del)[2]="MESO"
ggplot(data.del) +
  geom_point( aes(x=LUAD, y=MESO), color="black", size=2 ) +
  labs(
    color="groups"
  )+
  # theme_ipsum() +
  theme_classic(
  ) +
  xlab("Frequency of chr-band deletion in LUAD") +
  ylab("Frequency of chr-band deletion in MESO")

cor.test(data.del$LUAD,data.del$MESO)
# 0.4671548 p-value < 2.2e-16

###############################################################################
rm(list=ls())
myinf1 = # LUAD cnv
myinf2 = # MESO cnv
#LUAD
load(myinf1)
data = mydata

xx = substr(colnames(data), 14, 15)
se = which(xx=="01")
data = data[, se]
dim(data)

amp.thr = log2(2.8/2)
del.thr = log2(1.4/2)

amp.f = apply(data>amp.thr, 1, sum)
del.f = apply(data<del.thr, 1, sum)
amp.f = as.data.frame(amp.f/ncol(data))
del.f = as.data.frame(del.f/ncol(data))

amp1 = amp.f
del1 = del.f

#MESO
load(myinf2)
data = mydata
xx = substr(colnames(data), 14, 15)
se = which(xx=="01")
data = data[, se]
dim(data)
amp.thr = log2(2.8/2)
del.thr = log2(1.4/2)
amp.f = apply(data>amp.thr, 1, sum)
del.f = apply(data<del.thr, 1, sum)
# divided by number of samples is to calculate the frequency of loss/gain
amp.f = amp.f/ncol(data)
del.f = del.f/ncol(data)
amp2 = amp.f
del2 = del.f

# amp2
cor.test(amp1$`amp.f/ncol(data)`, amp2)
# 0.7705273,p-value < 2.2e-16

cor.test(del1$`del.f/ncol(data)`, del2)
# 0.4252737, p-value < 2.2e-16

library(ggplot2)
library(dplyr)


data.amp <- data.frame(
  LUAD=del1,
  MESO=del2
)
colnames(data.amp)[1]="LUAD"

ggplot(data.amp) +
  geom_point( aes(x=LUAD, y=MESO), color="black", size=2 ) +
  labs(
    color="groups"
  )+
  theme_classic(
  ) +
  xlab("Frequency of genes Amplified in LUAD") +
  ylab("Frequency of genes Amplified in MESO")

#Deletion Frequency
data.del <- data.frame(
  value1=amp1,
  MESO=amp2
)
colnames(data.del)[1]="LUAD"

ggplot(data.del) +
  geom_point( aes(x=LUAD, y=MESO), color="black", size=2 ) +
  labs(
    color="groups"
  )+
  # theme_ipsum() +
  theme_classic(
  ) +
  xlab("Frequency of genes Deleted in LUAD") +
  ylab("Frequency of genes Deleted in MESO")

#634genes amplified in both
index=intersect(which(data$LUAD>0.1),which(data$MESO>0.1))
data.del= data[index,]

cor.test(data.del$LUAD,data.del$MESO)
#amp t = 23.335, df = 632, p-value < 2.2e-16,       cor 0.6803105 
#del: t = 13.779, df = 24, p-value = 6.792e-13,       cor 0.9422201 
################################################################################
# VennDiagram 5[C-E]
rm(list=ls())

myinf1 = #path to bueno gexp
myinf2 = #path to bueno clinical info
myinf3 = #path to GSE31210 gexp
myinf4 = #path to LUAD prognostic test
myinf5 =  #path to MESO prognostic test

#prognostic
luad = read.table(myinf4, sep="\t", header=T,  row.names=1, quote="")
prog = read.table(myinf5, sep="\t", header=T,  row.names=1, quote="")
luad=luad[luad$adj.pval<0.01,]
prog=prog[prog$adj.pval<0.01,]

overlap=intersect(row.names(luad), row.names(prog))
library(VennDiagram)
grid.newpage()
draw.pairwise.venn(area1=dim(prog)[1], area2=dim(luad)[1],
                   cross.area=length(overlap),fill=c("Yellow","Red"),lab.cex=2)

total = 20501
# prognostic genes in mpm
m.q4 = dim(luad)[1]
# prognostic genes in luad
l.q4 = dim(prog)[1]
shared = length(overlap)
## Method 1
white = m.q4   
black = total - white	
take = l.q4    


tmp=phyper(shared:take, white, black, take )
pval1 = sum(dhyper(shared:take, white, black, take ))
pval1

#pval: 2.781295e-143

mpm = read.table("#path to 4 quadrants")
se=which(mpm$quadrant=="Q4")

# prognostic genes in okayama  
prog.luad= luad[which(row.names(luad)%in%row.names(ok.data)),]
prog.meso = prog[which(row.names(prog)%in%row.names(ok.data)), ]
overlap=intersect(row.names(prog.luad), row.names(prog.meso))
library(VennDiagram)
grid.newpage()
draw.pairwise.venn(area1=dim(prog.meso)[1], area2=dim(prog.meso)[1],
                   cross.area=length(overlap),fill=c("Red","Yellow"),lab.cex=2)
# pvalue for overlap
#space: all filtered genes in okayama
total = 17724
# prognostic genes in mpm
m.q4 = 2314
# prognostic genes in luad
l.q4 = 1086
shared = 275
white = m.q4   ## assume 2000 are white balls
black = total - white	## the number of black balls in a box
take = l.q4     ## now you take 3000 genes from the box

## Method 1
## p-value is the probability of observing 1000 pr more white ball b chance
tmp=phyper(shared:take, white, black, take )
pval1 = sum(dhyper(shared:take, white, black, take ))
pval1


# prognostic genes in bueno                   
prog.luad= luad[which(row.names(luad)%in%row.names(bueno.data)),]
prog.meso = prog[which(row.names(prog)%in%row.names(bueno.data)), ]
overlap= intersect(row.names(prog.luad), row.names(prog.meso))
grid.newpage()
draw.pairwise.venn(area1=dim(prog.luad)[1], area2=dim(prog.meso)[1],
                   cross.area=length(overlap),fill=c("Red","Yellow"),lab.cex=2)

total = 20267
# prognostic genes in mpm
m.q4 = 2259
# prognostic genes in luad
l.q4 = 1088
shared = 268
## Method 1
white = m.q4   ## assume 2000 are white balls
black = total - white	## the number of black balls in a box
take = l.q4     ## now you take 3000 genes from the box

tmp=phyper(shared:take, white, black, take )
pval1 = sum(dhyper(shared:take, white, black, take ))
pval1


#Separate them by hazard ratio (main/supplementary?)
upok.df.prog = prog.luad[prog.luad$HR1>1,]
#544
downok.df.prog = prog.luad[prog.luad$HR1<1,]
#542

upbueno.df.prog=prog.meso[prog.meso$HR1>1,]
#1390
downbueno.df.prog = prog.meso[prog.meso$HR1<1,]
#924

length(intersect(rownames(downok.df.prog),row.names(downbueno.df.prog)))
# 184
#70

grid.newpage()

#HR>1
# create pairwise Venn diagram
grid.newpage()
draw.pairwise.venn(area1=740+208, area2=208+264,cross.area=208,fill=c("Yellow","Red"),lab.cex=2)


## assume you have nn= 20000 genes; MPM has 2500 Q4 and LUAD has 3000 Q4 genes, shared = 1000 genes
total = 17724
m.q4 = 2314
l.q4 = 1086
shared = 275
white = m.q4   
black = total - white	
take = l.q4     

tmp=phyper(shared:take, white, black, take )
pval1 = sum(dhyper(shared:take, white, black, take ))
pval1

###########################################
# Bueno signatures
rm(list=ls())
myinf1 = #path to gexp
myinf3 = #path to ORACLE
myinf4 = #path to PRACME

oracle = read.table(myinf3, sep="\t", header=T, row.names=1, stringsAsFactors=F, quote="")
miracle = read.table(myinf4, sep="\t", header=T, row.names=1, stringsAsFactors=F, quote="")

#-------------
data = read.table(myinf1, sep="\t", header=T, row.names=1, quote="")
xx = apply(data, 2, sum)
mymed = median(xx)
for(k in 1:ncol(data))
{
  data[,k] = data[,k]*(mymed/xx[k])
}
data = log2(data+1)

mygen=oracle
comxx = intersect(row.names(data), row.names(mygen))
data = data[comxx,]
mygen = mygen[comxx,]
s1 = t(data) %*% mygen /length(mygen)


mygen = miracle
comxx = intersect(row.names(data), row.names(mygen))
data = data[comxx,]
mygen = mygen[comxx,]
s2 = t(data) %*% mygen /length(mygen)

cor.test(s1, s2)
# 0.6245028  p< 2.2e-16
df <- data.frame(
  PRACME=s2,
  ORACLE=s1
)


ggplot(df) +
  geom_point( aes(x=PRACME, y=ORACLE), color="black", size=2 ) +
  geom_smooth(aes(x=PRACME, y=ORACLE),method = "lm", se = FALSE,color="red")+
  labs(
    color="groups"
  )+
  theme_classic(
  ) +
  xlab("PRACME") +
  ylab("ORACLE")
 

#-------------
rm(list=ls())
myinf1 = #path to GSE31210 expression
myinf3 = #path to ORACLE
myinf4 = #path to MIRACLE

data = read.table(myinf1, sep="\t", header=T, row.names=1, quote="")
xx = apply(data, 2, sum)
mymed = median(xx)
for(k in 1:ncol(data))
{
  data[,k] = data[,k]*(mymed/xx[k])
}
data = log2(data+1)

data <- na.omit(data)

oracle = read.table(myinf3, sep="\t", header=T, row.names=1, stringsAsFactors=F, quote="")
miracle = read.table(myinf4, sep="\t", header=T, row.names=1, stringsAsFactors=F, quote="")

#-------------

mygen=oracle
xx = apply(data, 2, sum)
mymed = median(xx)
for(k in 1:ncol(data))
{
  data[,k] = data[,k]*(mymed/xx[k])
}
data = log2(data+1)
comxx = intersect(row.names(data), row.names(mygen))
data = data[comxx,]
mygen = mygen[comxx,]
s1 = t(data) %*% mygen /length(mygen)


mygen = miracle
xx = apply(data, 2, sum)
mymed = median(xx)
for(k in 1:ncol(data))
{
  data[,k] = data[,k]*(mymed/xx[k])
}
data = log2(data+1)
comxx = intersect(row.names(data), row.names(mygen))
data = data[comxx,]
mygen = mygen[comxx,]
s2 = t(data) %*% mygen /length(mygen)
cor.test(s1, s2)
# 0.6837352  p< 2.2e-16
df <- data.frame(
  PRACME=s2,
  ORACLE=s1
)


ggplot(df) +
  geom_point( aes(x=PRACME, y=ORACLE), color="black", size=2 ) +
  geom_smooth(aes(x=PRACME, y=ORACLE),method = "lm", se = FALSE,color="red")+
  labs(
    color="groups"
  )+
  theme_classic(
  ) +
  xlab("PRACME") +
  ylab("ORACLE")

################################################################################
rm(list=ls())
myinf2 = #path to MESO CNV
myinf3 = #path to ORACLE
myinf4 = #path to PRACME


oracle = read.table(myinf3, sep="\t", header=T, row.names=1, stringsAsFactors=F, quote="")
miracle = read.table(myinf4, sep="\t", header=T, row.names=1, stringsAsFactors=F, quote="")


#MESO
load(myinf2)
data = mydata
xx = substr(colnames(data), 14, 15)
se = which(xx=="01")
data = data[, se]
dim(data)

mygen=oracle
comxx = intersect(row.names(data), row.names(mygen))
data = data[comxx,]

amp.thr = log2(2.8/2)
del.thr = log2(1.4/2)
amp.f = apply(data>amp.thr, 1, sum)
del.f = apply(data<del.thr, 1, sum)
# divided by number of samples is to calculate the frequency of loss/gain
amp.f = amp.f/ncol(data)
del.f = del.f/ncol(data)

amp2 = amp.f
del2 = del.f
hist(amp2)
hist(del2)

df <- data.frame(
  amp=amp2,
  del=del2
)


ggplot(df) +
  geom_point( aes(x=del, y=amp), color="black", size=2 ) +
  geom_smooth(aes(x=del, y=amp),method = "lm", se = FALSE,color="red")+
  labs(
    color="groups"
  )+
  theme_classic(
  ) +
  xlab("Del") +
  ylab("Amp")
