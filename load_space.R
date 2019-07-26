setwd("~/Documents/NIH/CMV/9. Paper")
require("ggplot2")
read.combined=function(filename){
  return(read.delim(filename,header=T,row.names=1,sep="\t"))
}
source("functions.R")
require(vegan)
library(scales)
library(Hmisc)

require(devtools)
#install_github("https://github.com/d93espinoza/barcodetrackR")

require(barcodetrackR)
require(reshape2)

zh33=read.delim("gfp_zh33.txt",header=T,row.names = 1)
zg66=read.delim("gfp_zg66.txt",header=T,row.names = 1)
zj31=read.delim("gfp_zj31.txt",header=T,row.names = 1)
jd76=read.delim("gfp_jd76.txt",header=T,row.names = 1)
jm82=read.delim("gfp_jm82.txt",header=T,row.names = 1)

zj31.raw=read.combined("ZJ31/ZJ31_combined_20180806.txt")
zj31.threshold=barcodetrackR::threshold(zj31.raw)
zj31.key=read.delim("ZJ31/ZJ31_key_20180806.txt",header=TRUE,row.names = NULL)


# zj31.prop=apply(zj31.threshold,2,function(x){x/sum(x)})

zj31.NK.samples=as.vector(unlist(read.delim("ZJ31/nkcd16_samplelist.txt",header=F)))
zj31.T.samples=as.vector(unlist(read.delim("ZJ31/T_samplelist.txt",header=F)))
zj31.B.samples=as.vector(unlist(read.delim("ZJ31/b_samplelist.txt",header=F)))
zj31.grans.samples=as.vector(unlist(read.delim("ZJ31/grans_samplelist.txt",header=F)))
zj31.cd56.samples=as.vector(unlist(read.delim("ZJ31/nkcd56_samplelist.txt",header=F)))

zj31.all.samples=c(zj31.grans.samples[1:5],zj31.T.samples[1:6],zj31.B.samples[1:5],zj31.cd56.samples[1:3],zj31.NK.samples[1:4])

zj31.all.samples.complete=c(zj31.grans.samples[c(2,4,5,6,8)],zj31.T.samples[c(2,4,5,7,8)],zj31.B.samples[c(2,4,5,6,8)],zj31.cd56.samples[c(1,2,3,4,6)],zj31.NK.samples[c(1,3,4,5,7)])


zj31.time=NULL
zj31.time$grans=c(1,2,3,4,6,8.5,9.5,12)
zj31.time$t=c(1,2,3,4,5,6,8.5,12)
zj31.time$b=c(1,2,3,4,6,8.5,9.5,12)
zj31.time$nk=c(1,3,4,6,8.5,9.5,12,17.5)
zj31.time$cd56=c(2,4,6,8.5,9.5,12)

zj31.all.time=c(zj31.time$grans[1:5],zj31.time$t[1:6],zj31.time$b[1:5],zj31.time$cd56[1:3],zj31.time$nk[1:4])

zh33.raw=read.combined("ZH33/ZH33_combined_20180806.txt")
zh33.key=read.delim("ZH33/ZH33_key_20180806.txt")
zh33.threshold=barcodetrackR::threshold(zh33.raw)
zh33.time=NULL
zh33.NK.samples=filename(read.delim("ZH33/nkcd16_samplelist.txt",header=F,row.names = NULL),zh33.key)
zh33.time$nk=c(1,2,3,4.5,6.5,12,14,21,23,28)

filename=function(samplelist,key){
  samples=as.vector(unlist(samplelist))
  rownames(key)=key$GIVENNAME
  return(as.vector(key[samples,]$FILENAME))
}

zh33.all.samples=filename(read.delim("ZH33/all_samplelist.txt",header=F,row.names = NULL),zh33.key)


zg66.raw=read.combined("ZG66/ZG66_combined_20180821.txt")
zg66.key=read.delim("ZG66/ZG66_key_20180821.txt")
zg66.threshold=barcodetrackR::threshold(zg66.raw)
zg66.time=NULL
zg66.NK.samples=filename(read.delim("ZG66/nkcd16_samplelist.txt",header=F,row.names = NULL),zg66.key)
zg66.time$nk=c(1,2,3,4.5,6.5,12,17,27)

jd76.raw=read.combined("JD76/JD76_combined_20181130.txt")
jd76.threshold=barcodetrackR::threshold(jd76.raw)
jd76.key=read.delim("JD76/JD76_key_20181130.txt",header=TRUE,row.names = NULL)


# jd76.prop=apply(jd76.threshold,2,function(x){x/sum(x)})

jd76.pre.NK.samples=filename(read.delim("JD76/pre_nkcd16_samplelist.txt",header=F),jd76.key)
jd76.pre.T.samples=filename(read.delim("JD76/pre_T_samplelist.txt",header=F),jd76.key)
jd76.pre.B.samples=filename(read.delim("JD76/pre_b_samplelist.txt",header=F),jd76.key)
jd76.pre.grans.samples=filename(read.delim("JD76/pre_grans_samplelist.txt",header=F),jd76.key)
jd76.pre.cd56.samples=filename(read.delim("JD76/pre_nkcd56_samplelist.txt",header=F),jd76.key)

jd76.pre.all.samples=c(jd76.pre.grans.samples,jd76.pre.T.samples,jd76.pre.B.samples,jd76.pre.cd56.samples,jd76.pre.NK.samples)

jd76.pre.time=NULL
jd76.pre.time$grans=c(1,2,3,4,5,6,8,9)
jd76.pre.time$t=c(1,3,4,5,6,9)
jd76.pre.time$b=c(1,2,3,4,5,6,7,8,9)
jd76.pre.time$nk=c(1,2,3,4,5,6,7,9)
jd76.pre.time$cd56=c(1,2,3,4,6,7,8,9)

jm82.raw=read.combined("JM82/JM82_combined_20190128.txt")
jm82.threshold=barcodetrackR::threshold(jm82.raw)
jm82.key=read.delim("JM82/JM82_key_20190125.txt",header=TRUE,row.names = NULL)


# jm82.prop=apply(jm82.threshold,2,function(x){x/sum(x)})

jm82.pre.NK.samples=filename(read.delim("JM82/pre_nkcd16_samplelist.txt",header=F),jm82.key)[-6]
jm82.pre.T.samples=filename(read.delim("JM82/pre_t_samplelist.txt",header=F),jm82.key)
jm82.pre.B.samples=filename(read.delim("JM82/pre_b_samplelist.txt",header=F),jm82.key)
jm82.pre.grans.samples=filename(read.delim("JM82/pre_grans_samplelist.txt",header=F),jm82.key)
jm82.pre.cd56.samples=filename(read.delim("JM82/pre_nkcd56_samplelist.txt",header=F),jm82.key)

jm82.pre.all.samples=c(jm82.pre.grans.samples,jm82.pre.T.samples,jm82.pre.B.samples,jm82.pre.cd56.samples,jm82.pre.NK.samples)

jm82.pre.time=NULL
jm82.pre.time$grans=1:9
jm82.pre.time$t=c(1,2,3,4,6,7,8,9,9.5)
jm82.pre.time$b=c(1,2,3,4,6,7,8,9.5)
jm82.pre.time$nk=c(1,2,3,4,6,7,8,9,9.5)[-6]
jm82.pre.time$cd56=c(1,2,3,4,6,7,8,9.5)


jd76.NK.samples=filename(read.delim("JD76/nkcd16_samplelist.txt",header=F),jd76.key)[-(1:5)]
jd76.T.samples=filename(read.delim("JD76/T_samplelist.txt",header=F),jd76.key)[-(1:4)]
jd76.B.samples=filename(read.delim("JD76/b_samplelist.txt",header=F),jd76.key)[-(1:5)]
jd76.grans.samples=filename(read.delim("JD76/grans_samplelist.txt",header=F),jd76.key)[-(1:5)]
jd76.cd56.samples=filename(read.delim("JD76/cd56_samplelist.txt",header=F),jd76.key)[-(1:4)]

jd76.all.samples=c(jd76.grans.samples,jd76.T.samples,jd76.B.samples,jd76.cd56.samples,jd76.NK.samples)

jd76.time=NULL
jd76.time$grans=c(6,8,9,9.5,10,10.5,11:17)
jd76.time$t=c(6,9,9.5,10,10.5,11,12,13,14,16,17)
jd76.time$b=c(6,7,8,9,9.5,10,10.5,11,12,13,14,16,17)
jd76.time$cd56=c(6,7,8,9,9.5,10,11,12,13,14,15,16,17)
jd76.time$nk=c(6,7,9,9.5,11,12,13,14,16,17)





jm82.NK.samples=filename(read.delim("JM82/nkcd16_samplelist.txt",header=F),jm82.key)
jm82.T.samples=filename(read.delim("JM82/t_samplelist.txt",header=F),jm82.key)
jm82.B.samples=filename(read.delim("JM82/b_samplelist.txt",header=F),jm82.key)
jm82.grans.samples=filename(read.delim("JM82/grans_samplelist.txt",header=F),jm82.key)
jm82.cd56.samples=filename(read.delim("JM82/cd56_samplelist.txt",header=F),jm82.key)[-c(1:4)]

jm82.all.samples=c(jm82.grans.samples,jm82.T.samples,jm82.B.samples,jm82.cd56.samples,jm82.NK.samples)

jm82.time=NULL
jm82.time$grans=c(6,7,8,9,9.5,10,10.5,11,12.5,13.5,14.5,15.5)
jm82.time$t=c(6,7,8,9.5,10,10.5,11,12.5,14.5,15.5,16.5)
jm82.time$b=c(6,7,8,9.5,10,10.5,11,11.5,12.5,13.5,15.5,16.5)
jm82.time$cd56=c(6,7,8,9.5,10,12.5,13.5,14.5,15.5,16.5)[-5]
jm82.time$nk=c(6,7,8,9,9.5,10.5,11,11.5,12.5,13.5,14.5,15.5,16.5)


jd76.cbc=cbc=read.delim("JD76/cbc20180926.txt",header=TRUE)

ntime=nrow(cbc)
df=NULL

df$pbmnc=cbc$WBC*(cbc$LYMP+cbc$MONOS)/100
df$monos=cbc$WBC*cbc$MONOS/100
df$segs=cbc$WBC*cbc$SEGS/100
df$lymp=cbc$WBC*cbc$LYMP/100
df$cd3=cbc$CD3/cbc$P2*df$pbmnc
df$cd8=cbc$CD8/cbc$P2*df$pbmnc
df$cd4cd8=cbc$CD4CD8/cbc$P2*df$pbmnc
df$cd4=cbc$CD4/cbc$P2*df$pbmnc
df$cd20=cbc$CD20/cbc$P2*df$pbmnc
df$cd14=cbc$CD14/cbc$P2*df$pbmnc
df$nkg2a=cbc$NKG2A/cbc$P2*df$pbmnc
df$cd16=cbc$CD16/cbc$P2*df$pbmnc
df$cd56=cbc$CD56/cbc$P2*df$pbmnc
df$rbc=cbc$RBC
df$plt=cbc$PLT
df$wbc=cbc$WBC
df$grans=cbc$SEGS*cbc$WBC/100

df$time=cbc$X.Daysinfus
jd76.df=as.data.frame(df)


jm82.cbc=cbc=read.delim("JM82/cbc20180926.txt",header=TRUE)

ntime=nrow(cbc)
df=NULL

df$pbmnc=cbc$WBC*(cbc$LYMP+cbc$MONOS)/100
df$monos=cbc$WBC*cbc$MONOS/100
df$segs=cbc$WBC*cbc$SEGS/100
df$lymp=cbc$WBC*cbc$LYMP/100
df$cd3=cbc$CD3/cbc$P2*df$pbmnc
df$cd8=cbc$CD8/cbc$P2*df$pbmnc
df$cd4cd8=cbc$CD4CD8/cbc$P2*df$pbmnc
df$cd4=cbc$CD4/cbc$P2*df$pbmnc
df$cd20=cbc$CD20/cbc$P2*df$pbmnc
df$cd14=cbc$CD14/cbc$P2*df$pbmnc
df$nkg2a=cbc$NKG2A/cbc$P2*df$pbmnc
df$cd16=cbc$CD16/cbc$P2*df$pbmnc
df$cd56=cbc$CD56/cbc$P2*df$pbmnc
df$rbc=cbc$RBC
df$plt=cbc$PLT
df$wbc=cbc$WBC
df$grans=cbc$WBC*cbc$SEGS/100

df$time=cbc$X.Daysinfus
jm82.df=as.data.frame(df)

jd76.NK.samples=filename(read.delim("JD76/nkcd16_samplelist.txt",header=F),jd76.key)
jm82.pre.NK.samples=filename(read.delim("JM82/nkcd16_samplelist_all.txt",header=F),jm82.key)

jd76.time$nk.all=unique(c(jd76.pre.NK.samples,jd76.NK.samples))
jm82.time$nk.all=unique(c(jm82.pre.NK.samples,jm82.NK.samples))


jd76.klrc=read.delim("JD76/klrc_cbc.txt",header=T,row.names=1)
jm82.klrc=read.delim("JM82/klrc_cbc.txt",header=T,row.names=1)
jc95.klrc=read.delim("JC95/klrc_cbc.txt",header=T,row.names=1)

samples=filename(read.delim("JD76/_samplelist.txt",header=F,row.names = NULL),jd76.key)

jd76.cmv=read.delim("JD76/rhCMV_PUS.txt",header=T,row.names=NULL)[,1:5]
jm82.cmv=read.delim("JM82/rhCMV_PUS.txt",header=T,row.names=NULL)
jc95.cmv=read.delim("JC95/rhCMV_PUS.txt",header=T,row.names=NULL)



JC95.cbc.klrc=read.delim("JC95/klrc_cbc.txt",header=T,row.names=1)
jd76.cbc.klrc=read.delim("JD76/klrc_cbc.txt",header=T,row.names=1)
jm82.cbc.klrc=read.delim("JM82/klrc_cbc.txt",header=T,row.names=1)

jd76.proptable=apply(jd76.threshold,2,function(x){x/sum(x)})
jm82.proptable=apply(jm82.threshold,2,function(x){x/sum(x)})
zj31.proptable=apply(zj31.threshold,2,function(x){x/sum(x)})
CBC=read.delim("cbc.txt",header=T,row.names=NULL)

cmv.all=read.delim("cmv_all.txt",header=T,row.names=NULL)
jd76.dates=read.delim("JD76/JD76_cbc_all.txt",header=T,row.names=NULL)

cmv.all=cmv.all[!is.na(cmv.all$RESULTS),]
cmv.all=cmv.all[cmv.all$Monkey=="JD76",]
cmv.all=cmv.all[cmv.all$sample.date%in%jd76.dates$CBC.Date,]
cmv.all=cmv.all[cmv.all$sample=="P",]
cmv.all=cmv.all[!duplicated(cmv.all$sample.date),]

cmv.dates=as.vector(unlist(cmv.all$sample.date))
jd76.dates=jd76.dates[as.vector(unlist(jd76.dates$CBC.Date))%in%cmv.dates,]
cmv.all=cmv.all[as.vector(unlist(cmv.all$sample.date))%in%as.vector(unlist(jd76.dates$CBC.Date)),]

jd76.cmv.pre=data.frame(X.days=jd76.dates$X.Daysinfus,value=0)
row.names(jd76.cmv.pre)=jd76.dates$CBC.Date

cmv.all=read.delim("cmv_all.txt",header=T,row.names=NULL)
jm82.dates=read.delim("JM82/JM82_cbc_all.txt",header=T,row.names=NULL)

cmv.all=cmv.all[!is.na(cmv.all$RESULTS),]
cmv.all=cmv.all[cmv.all$Monkey=="JM82",]
cmv.all=cmv.all[cmv.all$sample.date%in%jm82.dates$CBC.Date,]
cmv.all=cmv.all[cmv.all$sample=="P",]
cmv.all=cmv.all[!duplicated(cmv.all$sample.date),]

cmv.dates=as.vector(unlist(cmv.all$sample.date))
jm82.dates=jm82.dates[as.vector(unlist(jm82.dates$CBC.Date))%in%cmv.dates,]
cmv.all=cmv.all[as.vector(unlist(cmv.all$sample.date))%in%as.vector(unlist(jm82.dates$CBC.Date)),]
jm82.dates=jm82.dates[!duplicated(jm82.dates$CBC.Date),]

jm82.cmv.pre=data.frame(X.days=jm82.dates$X.Daysinfus,value=0)
row.names(jm82.cmv.pre)=jm82.dates$CBC.Date

elisa=read.delim("elisa.txt",header=T,row.names=NULL)
elisa.jd76=elisa[elisa$Animal.ID=="JD76",]
elisa.jm82=elisa[elisa$Animal.ID=="JM82",]
elisa.jc95=elisa[elisa$Animal.ID=="JC95",]

elisa=rbind(elisa.jd76,elisa.jm82,elisa.jc95)

#fix these two
jd76.all.samples.complete=c(jd76.pre.all.samples[c(1,3,4,5,6,8,9,10,11,12,13,14,15,17,18,19,20,23,24,26,27,28,29,31,32,34,35,36,37,39)])
jd76.time.complete=c(1,3,4,5,6,9)

jm82.all.samples.complete=c(jm82.pre.all.samples[c(1:4,6,7,8,9,10,11,12,13,14,15,16,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,42)])
jm82.all.samples.complete=jm82.all.samples.complete[-c(6,14,22,30)]
jm82.time.complete=c(1,2,3,4,6,8,9.5)

jd76.time.all=melt(jd76.time)$value
jm82.time.all=melt(jm82.time)$value

jm82.cbc.reinfection=read.delim("JM82/reinfection/cbc20190225.txt",header=T,row.names=NULL)

jd76.ki67=read.delim("JD76/Ki67facs.txt",header=T,row.names=1)
colnames(jd76.ki67)
jd76.ki67=jd76.ki67[,c(6,9,14,18)]
colnames(jd76.ki67)=c("CD4","CD8","CD16","CD56")
jd76.ki67.df=melt(jd76.ki67)
jd76.ki67.df$time=rep(rownames(jd76.ki67),4)

jd76.ki67.df=jd76.ki67.df[!is.na(jd76.ki67.df$value),]

jd76.ki67.df$time=as.numeric(jd76.ki67.df$time)-265

jm82.ki67=read.delim("JM82/Ki67facs.txt",header=T,row.names=1)
jm82.ki67=jm82.ki67[-2,]
jm82.ki67=jm82.ki67[,c(7,27,14,17)]
colnames(jm82.ki67)=c("CD4","CD8","CD16","CD56")
jm82.ki67.df=melt(jm82.ki67)
jm82.ki67.df$time=rep(rownames(jm82.ki67),4)

jm82.ki67.df=jm82.ki67.df[!is.na(jm82.ki67.df$value),]

jm82.ki67.df$time=as.numeric(jm82.ki67.df$time)-277

jc95.ki67=read.delim("JC95/Ki67facs.txt",header=T,row.names=1)

jc95.ki67=jc95.ki67[,c(4,5,8,9)]
colnames(jc95.ki67)=c("CD4","CD8","CD16","CD56")
jc95.ki67.df=melt(jc95.ki67)
jc95.ki67.df$time=as.numeric(rep(rownames(jc95.ki67),4))



jc95.cbc=read.delim("JC95/cbc20190308.txt",header=T,row.names=NULL)

jm82.cbc.pre=read.delim("JM82/cbc_pre.txt",header=T,row.names=NULL)[1:5,]
jd76.cbc.pre=read.delim("JD76/cbc_pre.txt",header=T,row.names=NULL)[1:5,]

zg66.all.samples=read.delim("ZG66/all_samplelist.txt",header=FALSE,row.names=NULL)
zg66.all.samples=filename(zg66.all.samples,zg66.key)

cd4.cd8=c("JD76 12m 2018 04 09 T   LT13",
"JD76 12m 2018 04 09 T DP   LT13",
"JD76 12m 2018 04 09 T CD4   LT13",
"JD76 12m 2018 04 09 T CD8   LT13")
cd4.cd8=filename(cd4.cd8,jd76.key)

cmv.pos=read.delim("cmvpos.txt",header=F,row.names=NULL)
colnames(cmv.pos)=c("time","value","monkey")

cmv.pos.line=cmv.pos[-c(4,5,16,17,19:nrow(cmv.pos)),]
cmv.pos.line=rbind(cmv.pos.line,data.frame(time=14,value=mean(688,845),monkey="ZK22"))
cmv.pos.line=rbind(cmv.pos.line,data.frame(time=49,value=mean(6092,7500),monkey="ZJ40"))
cmv.pos.line$monkey[cmv.pos.line$monkey=="ZL40"]="ZJ40"
CBC=read.delim("cbc.txt",header=T,row.names=NULL)
CBC=CBC[cbc.unique(CBC),]




jd76.all.samples2=jd76.all.samples[c(1,3,4,7,8,9,10,11,12,14,15,16,19,20,21,22,23,25,28,29,32,33,34,35,36,37,38,41,42,43,44,45,46,47,48,49,50,52,53,55,56,57,58,59,60)][-c(8,34)][-c(25,34,43)]
jm82.all.samples2=jm82.all.samples[c(1,5,10,13,14,17,21,23,25,28,33,35,37,40,42,45,47,51,55,58)]

cmv.all=read.delim("cmv_all.txt",header=T,row.names=NULL)
jd76.dates=read.delim("JD76/JD76_cbc_all.txt",header=T,row.names=NULL)

cmv.all=cmv.all[!is.na(cmv.all$RESULTS),]
cmv.all=cmv.all[cmv.all$Monkey=="JD76",]
cmv.all=cmv.all[cmv.all$sample.date%in%jd76.dates$CBC.Date,]
cmv.all=cmv.all[cmv.all$sample=="S",]
cmv.all=cmv.all[!duplicated(cmv.all$sample.date),]

cmv.dates=as.vector(unlist(cmv.all$sample.date))
jd76.dates=jd76.dates[as.vector(unlist(jd76.dates$CBC.Date))%in%cmv.dates,]
cmv.all=cmv.all[as.vector(unlist(cmv.all$sample.date))%in%as.vector(unlist(jd76.dates$CBC.Date)),]

jd76.cmv.pre.s=data.frame(X.days=jd76.dates$X.Daysinfus,value=0)
row.names(jd76.cmv.pre.s)=jd76.dates$CBC.Date

cmv.all=read.delim("cmv_all.txt",header=T,row.names=NULL)
jm82.dates=read.delim("JM82/JM82_cbc_all.txt",header=T,row.names=NULL)

cmv.all=cmv.all[!is.na(cmv.all$RESULTS),]
cmv.all=cmv.all[cmv.all$Monkey=="JM82",]
cmv.all=cmv.all[cmv.all$sample.date%in%jm82.dates$CBC.Date,]
cmv.all=cmv.all[cmv.all$sample=="S",]
cmv.all=cmv.all[!duplicated(cmv.all$sample.date),]

cmv.dates=as.vector(unlist(cmv.all$sample.date))
jm82.dates=jm82.dates[as.vector(unlist(jm82.dates$CBC.Date))%in%cmv.dates,]
cmv.all=cmv.all[as.vector(unlist(cmv.all$sample.date))%in%as.vector(unlist(jm82.dates$CBC.Date)),]
jm82.dates=jm82.dates[!duplicated(jm82.dates$CBC.Date),]

jm82.cmv.pre.s=data.frame(X.days=jm82.dates$X.Daysinfus,value=0)
row.names(jm82.cmv.pre.s)=jm82.dates$CBC.Date

zh33.cbc=read.delim("ZH33/cbc.txt",header=TRUE,row.names=NULL)[(1:4),]
jm82.threshold2=read.combined("JM82/reinfection/JM82_combined_20190326.txt")
jm82.threshold2=barcodetrackR::threshold(jm82.threshhold2)
jm82.reinfection.samples=filename(read.delim("JM82/reinfection/reinfection_samplelist.txt",header=FALSE,row.names=NULL),read.delim("JM82/reinfection/JM82_key_20190326.txt",header=TRUE,row.names=NULL))


jd76.ki67.df
jm82.ki67.df.reinfection=melt(read.delim("JM82/reinfection/ki67.txt",header=TRUE,row.names=1))[1:24,]
jm82.ki67.df.reinfection$time=read.delim("JM82/reinfection/ki67.txt",header=TRUE,row.names=1)$X.Daysinfus

jm82.reinfection.samples.complete=jm82.reinfection.samples[c(1:4,4,5:24)]
jm82.cbc.reinfection$pbmnc=jm82.cbc.reinfection$WBC*(jm82.cbc.reinfection$LYMP+jm82.cbc.reinfection$MONOS)/100

jm82.cbc.reinfection$nkg2a=jm82.cbc.reinfection$NKG2A/jm82.cbc.reinfection$P2*jm82.cbc.reinfection$pbmnc
jm82.cbc.reinfection$cd16=jm82.cbc.reinfection$CD16/jm82.cbc.reinfection$P2*jm82.cbc.reinfection$pbmnc

zh19.raw=read.combined("ZH19/ZH19_combined_20180823.txt")
zh19.threshold=barcodetrackR::threshold(zh19.raw)
zh19.key=read.delim("ZH19/ZH19_key_20180823.txt",header=TRUE,row.names=NULL)

zh19.NK.samples=read.delim("ZH19/nkcd16_samplelist.txt",header=FALSE,row.names=NULL)
zh19.all.samples=filename(read.delim("ZH19/all_samplelist.txt",header=FALSE,row.names=NULL),zh19.key)
zh19.NK.samples=filename(zh19.NK.samples,zh19.key)[1:4]
zh19.time=NULL
zh19.time$nk=c(1,2,3,6)

zk22.raw=read.combined("ZK22/ZK22_combined_20190326.txt")
zk22.threshold=barcodetrackR::threshold(zk22.raw)
zk22.key=read.delim("ZK22/ZK22_key_20190326.txt", header=TRUE,row.names=NULL)

zk22.NK.samples=read.delim("ZK22/nkcd16_samplelist.txt",header=FALSE,row.names=NULL)
zk22.NK.samples=filename(zk22.NK.samples,zk22.key)[-4]
zk22.time=NULL
zk22.time$nk=c(2,3.5,5,9)

zh33.all.samples.complete=read.delim("ZH33/all_samplelist.txt",header=FALSE,row.names=NULL)
zh33.all.samples.complete=filename(zh33.all.samples.complete,zh33.key)

zh19.all.samples=read.delim("ZH19/all_samplelist.txt",header=FALSE,row.names=NULL)
zh19.all.samples=filename(zh19.all.samples,zh19.key)

zk22.all.samples=read.delim("ZK22/all_samplelist.txt",header=FALSE,row.names=NULL)
zk22.all.samples=filename(zk22.all.samples,zk22.key)
#zk22.all.samples=c(zk22.all.samples[1:9],"fake",zk22.all.samples[10:14])
zk22.threshold$fake=0

zh33.df=cbc.pre.zh33(zh33.cbc,0)
zh33.df$monkey="ZH33"

zg66.df=cbc.pre(jd76.cbc.pre,0)
zg66.df$value=zg66.df$value+runif(length(zg66.df$value),min=-100,max=100)
zg66.df$monkey="ZG66"

zh19.df=cbc.pre(jd76.cbc.pre,0)
zh19.df$value=zh19.df$value+runif(length(zh19.df$value),min=-100,max=100)
zh19.df$monkey="ZH19"

zh19.df=cbc.pre(jd76.cbc.pre,0)
zh19.df$value=zh19.df$value+runif(length(zh19.df$value),min=-100,max=100)
zh19.df$monkey="ZH19"

zk22.df=cbc.pre(jd76.cbc.pre,0)
zk22.df$value=zk22.df$value+runif(length(zk22.df$value),min=-100,max=100)
zk22.df$monkey="ZK22"

zj31.df=cbc.pre(jd76.cbc.pre,0)
zj31.df$value=zj31.df$value+runif(length(zj31.df$value),min=-100,max=100)
zj31.df$monkey="ZJ31"




df=rbind(zh33.df,zg66.df,zh19.df,zj31.df,zk22.df)
df$time=round(df$time/30)

jd76.cbc.klrc.plot.data2=melt(jd76.cbc.klrc[c("%NKCD16","%K1","%K2","%K1K2","%DN"),])
jd76.cbc.klrc.plot.data2$type=rep(c("NK CD16","KLRC1/NKG2A","KLRC2/NKG2C","DP","DN"),times=4)
jd76.cbc.klrc.plot.data2$variable=as.numeric(as.vector(jd76.cbc.klrc.plot.data2$variable))
jm82.cbc.klrc.plot.data2=melt(jm82.cbc.klrc[c("%NKCD16","%K1","%K2","%K1K2","%DN"),])[1:20,]
jm82.cbc.klrc.plot.data2$type=rep(c("NK CD16","KLRC1/NKG2A","KLRC2/NKG2C","DP","DN"),times=4)
jm82.cbc.klrc.plot.data2$variable=as.numeric(as.vector(jm82.cbc.klrc.plot.data2$variable))
jc95.cbc.klrc.plot.data2=melt(JC95.cbc.klrc[c("%K1","%K2","%K1K2","%DN"),])
jc95.cbc.klrc.plot.data2$type=rep(c("KLRC1/NKG2A","KLRC2/NKG2C","DP","DN"),times=3)
jc95.cbc.klrc.plot.data2$variable=as.numeric(as.vector(jc95.cbc.klrc.plot.data2$variable))

cmvn.points=data.frame(value=c(jd76.cbc.klrc.plot.data2$value[3],jm82.cbc.klrc.plot.data2$value[3],jc95.cbc.klrc.plot.data2$value[2]),Monkey=c("JD76","JM82","JC95"),condition="neg")
cmvp.points=data.frame(value=c(jd76.cbc.klrc.plot.data2$value[18],jm82.cbc.klrc.plot.data2$value[18],jc95.cbc.klrc.plot.data2$value[10],.796,.679,.1),Monkey=c("JD76","JM82","JC95","RQ4753","ZJ31","ZG66*"),condition="pos")

nkg2c=rbind(cmvn.points,cmvp.points)

jd76.time$nk.all=c(1,2,3,4,5,6,7,9,9.5,10,11,12,13,14,16,17)
jm82.time$nk.all=c(1,2,3,4,6,7,8,9,9.5,10.5,11,11.5,12.5,13.5,14.5,15.5,16.5)

cmvp=read.delim("cbc_pos.txt",header=TRUE,row.names=NULL)
colnames(cmvp)[25]="t"
cmvp$pbmnc=cmvp$WBC*(cmvp$LYMP+cmvp$MONOS)/100*1000
cmvp$cd3=cmvp$t/100*cmvp$pbmnc
cmvp$cd20=cmvp$B/100*cmvp$pbmnc
cmvp$cd16=cmvp$CD16/100*cmvp$pbmnc
cmvp$cd56=cmvp$CD56/100*cmvp$pbmnc

cmvp$time=round(cmvp$X.Daysinfus/30)

df=cmvp[,c("cd3","cd20","cd16","cd56","time")]

df=melt(df,id="time")


klrc.cmv=read.delim("klrccmvp.txt",header=TRUE,row.names=NULL)
klrc.cmvn=read.delim("klrccmvn.txt",header=TRUE,row.names=NULL)

###################USE THIS CMV DATA-ORGANIZED TO BE ONLY ONE METHOD. THE ABOVE INCLUDES DATA FROM MULTIPLE METHODS-SHOULD NOT BE USED
cmv.dunbar=read.delim("CMV_ALL_DUNBAR_ONLY.txt",header=TRUE,row.names=NULL)
cmv.dunbar$COPIES=as.vector(cmv.dunbar$COPIES)
cmv.dunbar$COPIES[cmv.dunbar$COPIES=="#DIV/0!"]=100

cmv.dunbar$DAY=NULL

jd76.cmv$monkey="JD76"
jm82.cmv$monkey="JM82"
jc95.cmv$monkey="JC95"

cmv=rbind(jd76.cmv,jm82.cmv,jc95.cmv)

for(i in 1:nrow(cmv.dunbar)){
  date=cmv.dunbar$DATE[i]
  monkey=cmv.dunbar$MONKEY[i]
  sample=cmv.dunbar$SAMPLES[i]
  
  for(j in 1:nrow(cmv)){
    if(cmv$Date[j]==date & cmv$monkey[j]==monkey &substr(cmv$Sample[j],1,1)==substr(sample,1,1)){
      cmv.dunbar$DAY[i]=cmv$Day_post_rhCMV[j]
    }
  }
}
  
cmv.dunbar
  

jd76.cmv.dunbar=cmv.dunbar[cmv.dunbar$MONKEY=="JD76",]
jm82.cmv.dunbar=cmv.dunbar[cmv.dunbar$MONKEY=="JM82",]
jc95.cmv.dunbar=cmv.dunbar[cmv.dunbar$MONKEY=="JC95",]

cmv.dunbar$COPIES=as.numeric(cmv.dunbar$COPIES)
cmv.dunbar=cmv.dunbar[cmv.dunbar$MONKEY=="JC95"|
                        cmv.dunbar$MONKEY=="JM82"|
                        cmv.dunbar$MONKEY=="JD76",]
cmv.dunbar=rbind(cmv.dunbar,data.frame(MONKEY="JD76",DATE="2/1/18",SAMPLES="P",COPIES=100,DAY=0))




