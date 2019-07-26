#source("load_space.R")
###############################################################################################
#################################Figure 1######################################################
###############################################################################################
dir.create("rhCMV")

pdf(file="rhCMV/positive.pdf")
ggplot(cmv.pos.line,aes(x=time,y=value,group=monkey,linetype=monkey))+
  geom_point(size=3)+
  geom_line(size=2)+
  scale_linetype_discrete("Monkey")+
  scale_x_continuous("Days post transplant")+
  scale_y_continuous("rhCMV DNA copy number per uL",trans = log10_trans(),breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)),limits=c(10^.5,10^4))+ theme_grey(base_size = 22)
dev.off()


jd76.cmv.pre$Monkey="JD76"
jm82.cmv.pre$Monkey="JM82"

cmv.pre=rbind(jd76.cmv.pre,jm82.cmv.pre)

pdf(file="rhCMV/negative.pdf")
ggplot(cmv.pre,aes(x=X.days,y=value,color=Monkey))+
  geom_point(size=3)+
  geom_line(size=2)+
  scale_color_manual(values=c("red","blue"))+ geom_line(linetype = 2)+
  scale_y_continuous("rhCMV DNA copy number per uL",trans = log10_trans(),breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)),limits=c(10^.5,10^4))+ theme_grey(base_size = 22)+
  scale_x_continuous("Days post transplant")
dev.off()

dir.create("cell_counts")
pdf("cell_counts/JM82_transplant.pdf")
cbc.pre(jm82.cbc.pre,0)
dev.off()

pdf("cell_counts/JD76_transplant.pdf")
cbc.pre(jd76.cbc.pre,0)
dev.off()

# pdf("cell_counts/ZH33_transplant.pdf")
# cbc.pre.zh33(zh33.cbc,0)
# dev.off()

pdf("cell_counts/CMVp.pdf")
df=cmvp[,c("cd3","cd20","cd16","cd56","time")]

df=melt(df,id="time")
ggplot()+geom_boxplot(df[df$variable=="cd20",],mapping=aes(time,value,group=time),fill=NA,color="black")+
  geom_boxplot(df[df$variable=="cd3",],mapping=aes(time,value,group=time),fill=NA,color="red")+
  geom_boxplot(df[df$variable=="cd16",],mapping=aes(time,value,group=time),fill=NA,color="green")+
  geom_boxplot(df[df$variable=="cd56",],mapping=aes(time,value,group=time),fill=NA,color="blue")+
  scale_y_continuous("Number of cells per uL PB",trans = log10_trans(),breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)),limits=c(10^0,10^4))+
  scale_x_continuous("Months post-transpant",breaks=seq(1,10,by=2))+theme_grey(base_size = 22)
dev.off()

###############################################################################################
#################################Figure 2######################################################
###############################################################################################
temp=barcode_ggheatmap_bar(jd76.threshold[,jd76.pre.all.samples],label_size = 10,names=c(paste(jd76.pre.time$grans,"m"),
                                                                                    paste(jd76.pre.time$t,"m"),
                                                                                    paste(jd76.pre.time$b,"m"),
                                                                                    paste(jd76.pre.time$cd56,"m"),
                                                                                    paste(jd76.pre.time$nk,"m")),grid=FALSE,cellnote_size = 3,printtable=TRUE)
dir.create("heatmaps")
pdf(file="heatmaps/JD76_pre.pdf")
temp=barcodetrackR::barcode_ggheatmap(jd76.threshold[,jd76.pre.all.samples],label_size = 10,names=c(paste(jd76.pre.time$grans,"m"),
                                                                                         paste(jd76.pre.time$t,"m"),
                                                                                         paste(jd76.pre.time$b,"m"),
                                                                                         paste(jd76.pre.time$cd56,"m"),
                                                                                         paste(jd76.pre.time$nk,"m")),grid=FALSE,cellnote_size = 3,printtable=TRUE)

input=as.data.frame(ifelse(temp[,ncol(temp)]>10*temp[,8],"red","black"))
rownames(input)=rownames(temp)
barcode_ggheatmap_bar(jd76.threshold[,jd76.pre.all.samples],label_size = 10,names=c("",paste(jd76.pre.time$grans,"m"),
                                                                                           paste(jd76.pre.time$t,"m"),
                                                                                           paste(jd76.pre.time$b,"m"),
                                                                                           paste(jd76.pre.time$cd56,"m"),
                                                                                           paste(jd76.pre.time$nk,"m")),grid=FALSE,cellnote_size = 3,input=input)
dev.off()

pdf(file="heatmaps/JM82_pre.pdf")

temp=barcodetrackR::barcode_ggheatmap(jm82.threshold[,jm82.pre.all.samples],label_size = 10,names=c(paste(jm82.pre.time$grans,"m"),
                                                                                               paste(jm82.pre.time$t,"m"),
                                                                                               paste(jm82.pre.time$b,"m"),
                                                                                               paste(jm82.pre.time$cd56,"m"),
                                                                                               paste(jm82.pre.time$nk,"m")),grid=FALSE,cellnote_size = 0,printtable=TRUE)
input=as.data.frame(ifelse(temp[,ncol(temp)]>10*temp[,9],"red","black"))
rownames(input)=rownames(temp)

barcode_ggheatmap_bar(jm82.threshold[,jm82.pre.all.samples],label_size = 10,COLORS=c("white","white","red","green"),names=c("",paste(jm82.pre.time$grans,"m"),
                                                                                               paste(jm82.pre.time$t,"m"),
                                                                                               paste(jm82.pre.time$b,"m"),
                                                                                               paste(jm82.pre.time$cd56,"m"),
                                                                                               paste(jm82.pre.time$nk,"m")),grid=FALSE,cellnote_size = 3,input=input)
dev.off()

pdf(file="heatmaps/ZJ31.pdf")
temp=barcodetrackR::barcode_ggheatmap(zj31.threshold[,zj31.all.samples[-c(8,16,24,30,37,38)]],label_size=10,names=c("",paste(zj31.time$grans,"m")[-8],
                                                                                                                    paste(zj31.time$t,"m")[-8],
                                                                                                                    paste(zj31.time$b,"m")[-8],
                                                                                                                    paste(zj31.time$cd56,"m")[-6],
                                                                                                                    paste(zj31.time$nk,"m")[-c(7,8)]),grid=FALSE,cellnote_size = 0,printtable=TRUE)

input=as.data.frame(ifelse(temp[,ncol(temp)]>10*temp[,7],"red","black"))
rownames(input)=rownames(temp)

barcode_ggheatmap_bar(zj31.threshold[,zj31.all.samples[-c(8,16,24,30,37,38)]],input=input,label_size=10,names=c("",paste(zj31.time$grans,"m")[-8],
                                                                                         paste(zj31.time$t,"m")[-8],
                                                                                         paste(zj31.time$b,"m")[-8],
                                                                                         paste(zj31.time$cd56,"m")[-6],
                                                                                         paste(zj31.time$nk,"m")[-c(7,8)]),grid=FALSE,cellnote_size = 3)
dev.off()

diversity=NULL
for(i in jd76.pre.all.samples){
  temp=jd76.threshold[,i]
  temp.df=data.frame(Month=i,Count=temp,Barcode=rownames(jd76.threshold))
  temp_diversity=vegan::diversity(jd76.threshold[,i])
  
  diversity=c(diversity,as.numeric(temp_diversity))
}

diversity=data.frame(Diversity=diversity,Cell=c(rep("Grans",length(jd76.pre.grans.samples)),rep("T",length(jd76.pre.T.samples)),rep("B",length(jd76.pre.B.samples)),rep("NK CD56",length(jd76.pre.cd56.samples)),rep("NK CD16",length(jd76.pre.NK.samples))),time=c(jd76.pre.time$grans,jd76.pre.time$t,jd76.pre.time$b,jd76.pre.time$cd56,jd76.pre.time$nk))

dir.create("diversity")
pdf(file="diversity/JD76_pre.pdf")
ggplot(diversity,aes(x=time,y=Diversity,color=Cell))+
  geom_line(size=2)+scale_y_continuous(limits=c(4,8.5))+
  geom_point(size=3)+theme_grey(base_size=22)+scale_x_continuous("Months post-transplant",breaks=seq(1,9.5,by=2),limits=c(1,10))+
  theme(legend.position = "none")
dev.off()

pdf(file="diversity/legend.pdf")
ggplot(diversity,aes(x=time,y=Diversity,color=Cell))+
  scale_color_discrete(labels=c("B","Grans","NK CD56-CD16+","NK CD56+CD16-","T"))+
  geom_line(size=2)+scale_y_continuous(limits=c(4,8.5))+
  geom_point(size=3)+theme_grey(base_size=22)+scale_x_continuous("Months post-transplant",breaks=seq(1,9.5,by=2),limits=c(1,10))
dev.off()

diversity=NULL
for(i in jm82.pre.all.samples){
  temp=jm82.threshold[,i]
  temp.df=data.frame(Month=i,Count=temp,Barcode=rownames(jm82.threshold))
  temp_diversity=vegan::diversity(jm82.threshold[,i])
  
  diversity=c(diversity,as.numeric(temp_diversity))
}

diversity=data.frame(Diversity=diversity,Cell=c(rep("Grans",length(jm82.pre.grans.samples)),rep("T",length(jm82.pre.T.samples)),rep("B",length(jm82.pre.B.samples)),rep("NK CD56",length(jm82.pre.cd56.samples)),rep("NK CD16",length(jm82.pre.NK.samples[c(1:5,7:9)]))),time=c(jm82.pre.time$grans,jm82.pre.time$t,jm82.pre.time$b,jm82.pre.time$cd56,jm82.pre.time$nk))

pdf(file="diversity/JM82_pre.pdf")
ggplot(diversity,aes(x=time,y=Diversity,color=Cell))+
  geom_line(size=2)+scale_y_continuous(limits=c(4,8.5))+
  geom_point(size=3)+theme_grey(base_size=22)+scale_x_continuous("Months post-transplant",breaks=seq(1,9.5,by=2),limits=c(1,10))+
  theme(legend.position = "none")
dev.off()

diversity=NULL
for(i in zj31.all.samples){
  temp=zj31.threshold[,i]
  temp.df=data.frame(Month=i,Count=temp,Barcode=rownames(zj31.threshold))
  temp_diversity=vegan::diversity(zj31.threshold[,i])
  
  diversity=c(diversity,as.numeric(temp_diversity))
}

diversity=data.frame(Diversity=diversity,Cell=c(rep("Grans",length(zj31.grans.samples)),rep("T",length(zj31.T.samples)),rep("B",length(zj31.B.samples)),rep("NK CD56",length(zj31.cd56.samples)),rep("NK CD16",length(zj31.NK.samples))),time=c(zj31.time$grans,zj31.time$t,zj31.time$b,zj31.time$cd56,zj31.time$nk))

pdf(file="diversity/ZJ31.pdf")
ggplot(diversity,aes(x=time,y=Diversity,color=Cell))+
  geom_line(size=2)+scale_y_continuous(limits=c(4,8.5))+
  geom_point(size=3)+theme_grey(base_size=22)+scale_x_continuous("Months post-transplant",breaks=seq(1,9.5,by=2),limits=c(1,10))+
  theme(legend.position = "none")
dev.off()


diversity=NULL
#samples=c(jd76.pre.NK.samples,jm82.pre.NK.samples,zj31.NK.samples,zh33.NK.samples,zg66.NK.samples)
for(i in jd76.pre.NK.samples){
  temp=jd76.threshold[,i]
  temp.df=data.frame(Month=i,Count=temp,Barcode=rownames(jd76.threshold))
  temp_diversity=vegan::diversity(jd76.threshold[,i])
  
  diversity=c(diversity,as.numeric(temp_diversity))
}
for(i in jm82.pre.NK.samples[1:9][-6]){
  temp=jm82.threshold[,i]
  temp.df=data.frame(Month=i,Count=temp,Barcode=rownames(jm82.threshold))
  temp_diversity=vegan::diversity(jm82.threshold[,i])
  
  diversity=c(diversity,as.numeric(temp_diversity))
}
for(i in zj31.NK.samples){
  temp=zj31.threshold[,i]
  temp.df=data.frame(Month=i,Count=temp,Barcode=rownames(zj31.threshold))
  temp_diversity=vegan::diversity(zj31.threshold[,i])
  
  diversity=c(diversity,as.numeric(temp_diversity))
}
for(i in zh33.NK.samples){
  temp=zh33.threshold[,i]
  temp.df=data.frame(Month=i,Count=temp,Barcode=rownames(zh33.threshold))
  temp_diversity=vegan::diversity(zh33.threshold[,i])
  
  diversity=c(diversity,as.numeric(temp_diversity))
}
for(i in zg66.NK.samples[-c(6:8)]){
  temp=zg66.threshold[,i]
  temp.df=data.frame(Month=i,Count=temp,Barcode=rownames(zg66.threshold))
  temp_diversity=vegan::diversity(zg66.threshold[,i])
  
  diversity=c(diversity,as.numeric(temp_diversity))
}
for(i in zh19.NK.samples){
  temp=zh19.threshold[,i]
  temp.df=data.frame(Month=i,Count=temp,Barcode=rownames(zh19.threshold))
  temp_diversity=vegan::diversity(zh19.threshold[,i])
  
  diversity=c(diversity,as.numeric(temp_diversity))
}
for(i in zk22.NK.samples){
  temp=zk22.threshold[,i]
  temp.df=data.frame(Month=i,Count=temp,Barcode=rownames(zk22.threshold))
  temp_diversity=vegan::diversity(zk22.threshold[,i])
  
  diversity=c(diversity,as.numeric(temp_diversity))
}
diversity=data.frame(Diversity=diversity,Cell=c(rep("JD76",length(jd76.pre.NK.samples)),rep("JM82",length(jm82.pre.NK.samples[1:9][-6])),rep("ZJ31",length(zj31.NK.samples)),rep("ZH33",length(zh33.NK.samples)),rep("ZG66",length(zg66.NK.samples[-c(6:8)])),rep("ZH19",length(zh19.NK.samples)),rep("ZK22",length(zk22.NK.samples))),
                     time=c(jd76.pre.time$nk,jm82.pre.time$nk,zj31.time$nk,zh33.time$nk,zg66.time$nk[-c(6:8)],zh19.time$nk,zk22.time$nk))
#diversity$shape=as.factor(diversity$shape)
colnames(diversity)=c("Diversity","Monkey","time")
pdf(file="diversity/all.pdf")
ggplot(diversity,aes(x=time,y=Diversity,color=Monkey,group=Monkey))+
  geom_line(size=1)+scale_y_continuous("NK CD16+ Diversity",limits=c(4,8.5))+scale_color_manual(values=c("red","blue",rep("black",5)))+
  geom_point(aes(shape=Monkey),size=4)+theme_grey(base_size=22)+scale_x_continuous("Months post-transplant",breaks=seq(1,9.5,by=2),limits=c(1,10))+
  scale_shape_manual(values=c(19,19,19,17,15,1,2,3))
dev.off()

dir.create("biased_barcodes")
pdf("biased_barcodes/JD76_pre.pdf")
biased_barcode_barplot(jd76.threshold[,jd76.all.samples.complete], length(jd76.time.complete), 5,jd76.time.complete, need_table = FALSE)+theme_grey(base_size = 22)
dev.off()

pdf("biased_barcodes/JM82_pre.pdf")
biased_barcode_barplot(jm82.threshold[,jm82.all.samples.complete], length(jm82.time.complete),5,jm82.time.complete, need_table = FALSE)+theme_grey(base_size = 22)
dev.off()

pdf("biased_barcodes/ZJ31.pdf")
biased_barcode_barplot(zj31.threshold[,zj31.all.samples.complete[-seq(5,25,by=5)]], length(zj31.time$cd56[-c(5,6)]), 5,zj31.time$cd56[-c(5,6)], need_table = FALSE)+theme_grey(base_size = 22)
dev.off()

pdf("biased_barcodes/all.pdf")
biased_barcode_summary()
dev.off()
###############################################################################################
#################################Figure 3######################################################
###############################################################################################

jd76.cmv$monkey="JD76"
jm82.cmv$monkey="JM82"
jc95.cmv$monkey="JC95"

cmv=rbind(jd76.cmv,jm82.cmv,jc95.cmv)
cmv$Mean=as.numeric(gsub("[^0-9\\.]", "", cmv$Mean) )

pdf(file="rhCMV/infection_plasma.pdf")
ggplot(cmv[cmv$Sample=="Plasma",],aes(x=Day_post_rhCMV,y=Mean,color=monkey))+
  geom_line(size=2)+
  geom_point(size=3)+
  scale_color_manual("Monkey",values=c("green","red","blue"))+
  scale_y_continuous("rhCMV DNA copy number per uL",trans = log10_trans(),breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)),limits=c(10^.5,10^4))+
  scale_x_continuous("Days post rhCMV infection",limits=c(-7,110),breaks=seq(0,1000,by=30))+theme_grey(base_size = 22)
dev.off()

pdf(file="rhCMV/infection_plasma_new.pdf")
ggplot(cmv.dunbar[cmv.dunbar$SAMPLE=="Plasma"|cmv.dunbar$SAMPLE=="P",],aes(x=DAY,y=COPIES,color=MONKEY))+
  geom_line(size=2)+
  geom_point(size=3)+
  scale_color_manual("Monkey",values=c("green","red","blue"))+
  scale_y_continuous("rhCMV DNA copy number per uL",trans = log10_trans(),breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)),limits=c(10^2,10^3))+
  scale_x_continuous("Days post rhCMV infection",limits=c(-7,110),breaks=seq(0,1000,by=30))+theme_grey(base_size = 22)
dev.off()

dir.create("IgG")
pdf(file="IgG/IgG.pdf",width=14)
ggplot(elisa,aes(x=Day.post.RhCMV,y=Average,color=Animal.ID))+
  geom_line(size=2)+geom_point(size=3)+
  scale_color_manual("Monkey",values=c("green","red","blue"))+
  scale_x_continuous("Days post rhCMV nfection",limits=c(-7,max(elisa$Day.post.RhCMV)),breaks=seq(0,1000,by=30))+theme_grey(base_size = 22)+
  scale_y_continuous("RhCMV Binding IgG")
dev.off()

pdf(file="cell_counts/legend.pdf")
cbc(jd76.cbc,265,legend=TRUE)
dev.off()
pdf(file="cell_counts/JD76.pdf")
cbc(jd76.cbc,265,legend=FALSE)
dev.off()

pdf(file="cell_counts/JM82.pdf")
cbc(jm82.cbc,277,legend=FALSE)
dev.off()

pdf(file="cell_counts/JC95.pdf")
cbc(jc95.cbc,0,legend=FALSE)
dev.off()

dir.create("Ki67")
pdf(file="Ki67/JD76.pdf")
ggplot(jd76.ki67.df,aes(x=time,y=value,color=variable,group=variable))+
  geom_line(na.rm=TRUE,size=2)+
  geom_point(size=3)+
  scale_color_manual(labels=c("T CD4+","T CD8+","NK CD16","NKCD56"),values=c("black","red","green","blue"))+
  scale_y_continuous("Percent Ki67+",limits=c(0,50))+
  scale_x_continuous("Days post rhCMV infection",limits=c(-30,225),breaks=seq(-30,1000,by=30))+theme_grey(base_size=22)+theme(legend.position="none")
dev.off()

pdf(file="Ki67/JM82.pdf")
ggplot(jm82.ki67.df,aes(x=time,y=value,color=variable,group=variable))+
  geom_line(na.rm=TRUE,size=2)+
  geom_point(size=3)+
  scale_color_manual(labels=c("T CD4+","T CD8+","NK CD16","NKCD56"),values=c("black","red","green","blue"))+
  scale_y_continuous("Percent Ki67+",limits=c(0,50))+
  scale_x_continuous("Days post rhCMV infection",limits=c(-30,225),breaks=seq(-30,1000,by=30))+theme_grey(base_size=22)+theme(legend.position="none")
dev.off()

pdf(file="Ki67/JC95.pdf")
ggplot(jc95.ki67.df,aes(x=time,y=value,color=variable,group=variable))+
  geom_line(na.rm=TRUE,size=2)+
  geom_point(size=3)+
  scale_color_manual(labels=c("T CD4+","T CD8+","NK CD16","NKCD56"),values=c("black","red","green","blue"))+
  scale_y_continuous("Percent Ki67+",limits=c(0,50))+
  scale_x_continuous("Days post rhCMV infection",limits=c(-30,225),breaks=seq(-30,1000,by=30))+theme_grey(base_size=22)+theme(legend.position="none")
dev.off()

###############################################################################################
#################################Figure 4######################################################
###############################################################################################
pdf("heatmaps/JD76_post_NKCD162.pdf")
temp=barcodetrackR::barcode_ggheatmap(jd76.threshold[,jd76.all.samples[c(15:24,41:60)]],names=paste(c(jd76.time.all[c(15:24,41:60)][1:24],10,jd76.time.all[c(15:24,41:60)][25:30]),"m")[-13],label_size = 7,grid=FALSE,cellnote_size = 0,printtable=TRUE)
input=as.data.frame(ifelse(temp[,ncol(temp)]>10*temp[,22],"red","black"))
rownames(input)=rownames(temp)
barcode_ggheatmap_bar(jd76.threshold[,jd76.all.samples[c(15:24,41:60)]],COLORS=c("white","white","green","black"),input=input,names=c("",paste(c(jd76.time.all[c(15:24,41:60)][1:24],10,jd76.time.all[c(15:24,41:60)][25:30]),"m")[-13]),label_size = 7,grid=FALSE,cellnote_size = 3)
dev.off()



pdf("heatmaps/JD76_post_NKCD161.pdf")
temp=barcodetrackR::barcode_ggheatmap(jd76.threshold[,jd76.all.samples[c(15:24,41:60)]],names=paste(c(jd76.time.all[c(15:24,41:60)][1:24],10,jd76.time.all[c(15:24,41:60)][25:30]),"m")[-13],label_size = 7,grid=FALSE,cellnote_size = 0,printtable=TRUE)
input=as.data.frame(ifelse(temp[,22]>10*temp[,ncol(temp)],"red","black"))
rownames(input)=rownames(temp)
barcode_ggheatmap_bar(jd76.threshold[,jd76.all.samples[c(15:24,41:60)]],COLORS=c("white","white","black","green"),input=input,names=c("",paste(c(jd76.time.all[c(15:24,41:60)][1:24],10,jd76.time.all[c(15:24,41:60)][25:30]),"m")[-13]),label_size = 7,grid=FALSE,cellnote_size = 3)
dev.off()

pdf("heatmaps/JD76_post_T.pdf")
temp=barcodetrackR::barcode_ggheatmap(jd76.threshold[,jd76.all.samples[c(15:24,41:60)]],names=paste(c(jd76.time.all[c(15:24,41:60)][1:24],10,jd76.time.all[c(15:24,41:60)][25:30]),"m")[-13],label_size = 7,grid=FALSE,cellnote_size = 0,printtable=TRUE)
input=as.data.frame(ifelse(temp[,6]>10*temp[,3],"red","black"))
rownames(input)=rownames(temp)
barcode_ggheatmap_bar(jd76.threshold[,jd76.all.samples[c(15:24,41:60)]],input=input,names=c("",paste(c(jd76.time.all[c(15:24,41:60)][1:24],10,jd76.time.all[c(15:24,41:60)][25:30]),"m")[-13]),label_size = 7,grid=FALSE,cellnote_size = 3)
dev.off()

pdf("heatmaps/JM82_post_T.pdf")
temp=barcodetrackR::barcode_ggheatmap(jm82.threshold[,jm82.all.samples[-c(9,41)][c(13:23,36:57)][-22]],names=paste(jm82.time.all[c(13:23,36:57)][-22],"m"),label_size = 7,grid=FALSE,cellnote_size = 0,printtable=TRUE)
input=as.data.frame(ifelse(temp[,9]>10*temp[,4],"red","black"))
rownames(input)=rownames(temp)
barcode_ggheatmap_bar(jm82.threshold[,jm82.all.samples[-c(9,41)][c(13:23,36:57)]][-22],input=input,names=c("",paste(jm82.time.all[c(13:23,36:57)][-22],"m")),label_size = 7,grid=FALSE,cellnote_size = 3)
dev.off()
pdf("heatmaps/JM82_post_NK1.pdf")
temp=barcodetrackR::barcode_ggheatmap(jm82.threshold[,jm82.all.samples[-c(9,41)][c(13:23,36:57)][-22]],names=paste(jm82.time.all[c(13:23,36:57)][-22],"m"),label_size = 7,grid=FALSE,cellnote_size = 0,printtable=TRUE)

input=as.data.frame(ifelse(temp[,24]>10*temp[,ncol(temp)],"red","black"))
rownames(input)=rownames(temp)
barcode_ggheatmap_bar(jm82.threshold[,jm82.all.samples[-c(9,41)][c(13:23,36:57)]][-22],COLORS=c("white","white","black","green"),input=input,names=c("",paste(jm82.time.all[c(13:23,36:57)][-22],"m")),label_size = 7,grid=FALSE,cellnote_size = 3)
dev.off()
pdf("heatmaps/JM82_post_NK2.pdf")
temp=barcodetrackR::barcode_ggheatmap(jm82.threshold[,jm82.all.samples[-c(9,41)][c(13:23,36:57)][-22]],names=paste(jm82.time.all[c(13:23,36:57)][-22],"m"),label_size = 7,grid=FALSE,cellnote_size = 0,printtable=TRUE)
input=as.data.frame(ifelse(temp[,ncol(temp)]>10*temp[,24],"red","black"))
rownames(input)=rownames(temp)
barcode_ggheatmap_bar(jm82.threshold[,jm82.all.samples[-c(9,41)][c(13:23,36:57)]][-22],COLORS=c("white","white","green","black"),input=input,names=c("",paste(jm82.time.all[c(13:23,36:57)][-22],"m")),label_size = 7,grid=FALSE,cellnote_size = 3)
dev.off()

jm82.all.samples[-c(9,41)][c(13:23,36:57)]

dir.create("normalized_heatmaps")
pdf("normalized_heatmaps/JD761.pdf")
NK.post.cmv=jd76.NK.samples[8:16]
jd76.time$nk=c(6,7,9,9.5,10,11,12,13,14,16,17)
temp=barcodetrackR::barcode_ggheatmap(jd76.threshold[,jd76.all.samples[c(15:24,41:60)]],names=paste(c(jd76.time.all[c(15:24,41:60)][1:24],10,jd76.time.all[c(15:24,41:60)][25:30]),"m")[-13],label_size = 7,grid=FALSE,cellnote_size = 0,printtable=TRUE)
input=as.data.frame(ifelse(temp[,22]>10*temp[,ncol(temp)],"black","white"))
rownames(input)=rownames(temp)
input=as.data.frame(input[rownames(barcode_ggheatmap(jd76.threshold[,NK.post.cmv],jd76.df$cd16,names=paste(jd76.time$nk[3:11],"m"),label_size=10,printtable=TRUE)),])
rownames(input)=rownames(barcode_ggheatmap(jd76.threshold[,NK.post.cmv],jd76.df$cd16,names=paste(jd76.time$nk[3:11],"m"),label_size=10,printtable=TRUE))
barcode_ggheatmapN_bar(jd76.threshold[,NK.post.cmv],jd76.df$cd16,COLORS=c("white","black","white","white"),input=input,names=c("",paste(jd76.time$nk[3:11],"m")),label_size=10)
dev.off()
pdf("normalized_heatmaps/JD762.pdf")
NK.post.cmv=jd76.NK.samples[8:16]
jd76.time$nk=c(6,7,9,9.5,10,11,12,13,14,16,17)
temp=barcodetrackR::barcode_ggheatmap(jd76.threshold[,jd76.all.samples[c(15:24,41:60)]],names=paste(c(jd76.time.all[c(15:24,41:60)][1:24],10,jd76.time.all[c(15:24,41:60)][25:30]),"m")[-13],label_size = 7,grid=FALSE,cellnote_size = 0,printtable=TRUE)
input=as.data.frame(ifelse(temp[,ncol(temp)]>10*temp[,22],"black","white"))
rownames(input)=rownames(temp)
input=as.data.frame(input[rownames(barcode_ggheatmap(jd76.threshold[,NK.post.cmv],jd76.df$cd16,names=paste(jd76.time$nk[3:11],"m"),label_size=10,printtable=TRUE)),])
rownames(input)=rownames(barcode_ggheatmap(jd76.threshold[,NK.post.cmv],jd76.df$cd16,names=paste(jd76.time$nk[3:11],"m"),label_size=10,printtable=TRUE))
barcode_ggheatmapN_bar(jd76.threshold[,NK.post.cmv],jd76.df$cd16,COLORS=c("white","green","white","white"),input=input,names=c("",paste(jd76.time$nk[3:11],"m")),label_size=10)
dev.off()

pdf("normalized_heatmaps/JM821.pdf")
NK.post.cmv=jm82.pre.NK.samples[5:13]
temp=barcodetrackR::barcode_ggheatmap(jm82.threshold[,jm82.all.samples[-c(9,41)][c(13:23,36:57)]],names=paste(jm82.time.all[c(13:23,36:57)],"m"),label_size = 7,grid=FALSE,cellnote_size = 0,printtable=TRUE)
input=as.data.frame(ifelse(temp[,ncol(temp)]>10*temp[,21],"red","black"))
rownames(input)=rownames(temp)
input=as.data.frame(input[rownames(barcode_ggheatmap(jm82.threshold[,NK.post.cmv],jm82.df$cd16,names=paste(jm82.time$nk[5:13],"m"),label_size=10,printtable=TRUE)),])
rownames(input)=rownames(barcode_ggheatmap(jm82.threshold[,NK.post.cmv],jm82.df$cd16,names=paste(jm82.time$nk[5:13],"m"),label_size=10,printtable=TRUE))
barcode_ggheatmapN_bar(jm82.threshold[,NK.post.cmv],jm82.df$cd16,COLORS=c("white","white","green","black"),input=input,names=c("",paste(jm82.time$nk[5:13],"m")),label_size=10)
dev.off()
pdf("normalized_heatmaps/JM822.pdf")
NK.post.cmv=jm82.pre.NK.samples[5:13]
temp=barcodetrackR::barcode_ggheatmap(jm82.threshold[,jm82.all.samples[-c(9,41)][c(13:23,36:57)]],names=paste(jm82.time.all[c(13:23,36:57)],"m"),label_size = 7,grid=FALSE,cellnote_size = 0,printtable=TRUE)
input=as.data.frame(ifelse(temp[,21]>10*temp[,ncol(temp)],"red","black"))
rownames(input)=rownames(temp)
input=as.data.frame(input[rownames(barcode_ggheatmap(jm82.threshold[,NK.post.cmv],jm82.df$cd16,names=paste(jm82.time$nk[5:13],"m"),label_size=10,printtable=TRUE)),])
rownames(input)=rownames(barcode_ggheatmap(jm82.threshold[,NK.post.cmv],jm82.df$cd16,names=paste(jm82.time$nk[5:13],"m"),label_size=10,printtable=TRUE))
barcode_ggheatmapN_bar(jm82.threshold[,NK.post.cmv],jm82.df$cd16,input=input,COLORS=c("white","white","black","green"),names=c("",paste(jm82.time$nk[5:13],"m")),label_size=10)
dev.off()

# plot(c(jd76.time$nk.all[-1]),(diag(cor(jd76.threshold[,jd76.NK.samples])[-1,-ncol(jd76.threshold[,jd76.NK.samples])])),type="b",ylab="Pair-wise Pearson Correlation",xlab="(Later) Month")
# abline(v=9.25)
# 
# plot(c(jm82.time$nk.all[-1]),(diag(cor(jm82.threshold[,jm82.pre.NK.samples])[-1,-ncol(jm82.threshold[,jm82.pre.NK.samples])])),type="b",ylab="Pair-wise Pearson Correlation",xlab="(Later) Month")
# abline(v=9.75)
# 
# plot(c(zh33.time$nk[-1]),(diag(cor(zh33.threshold[,zh33.NK.samples])[-1,-ncol(zh33.threshold[,zh33.NK.samples])])),type="b",ylab="Pair-wise Pearson Correlation",xlab="(Later) Month",ylim=c(0,1))
# 
# plot(c(zg66.time$nk[-1]),(diag(cor(zg66.threshold[,zg66.NK.samples])[-1,-ncol(zg66.threshold[,zg66.NK.samples])])),type="b",ylab="Pair-wise Pearson Correlation",xlab="(Later) Month",ylim=c(0,1))
# 
# pdf("biased_barcodes/JD76_post.pdf")
# jd76.time.complete2=c(6,9,9.5,11,12,13,14,16)
# biased_barcode_barplot(jd76.threshold[,jd76.all.samples2], length(jd76.time.complete2), 5,jd76.time.complete2, need_table = FALSE)+theme_grey(base_size = 22)
# dev.off()
# 
# pdf("biased_barcodes/JM82_post.pdf")
# jm82.time.complete2=c(6,9.5,12.5,15.5)
# biased_barcode_barplot(jm82.threshold[,jm82.all.samples2], length(jm82.time.complete2), 5,jm82.time.complete2, need_table = FALSE)+theme_grey(base_size = 22)
# dev.off()

jd76.time$nk.all=c(1,2,3,4,5,6,7,9,9.5,10,11,12,13,14,16,17)
jm82.time$nk.all=c(1,2,3,4,6,7,8,9,9.5,10.5,11,11.5,12.5,13.5,14.5,15.5,16.5)

dir.create("autocorrelation")
pdf("autocorrelation/all.pdf",width=10)
autocorrelation.cmvpos(list(jd76.threshold[,jd76.NK.samples],
                            jm82.threshold2[,c(jm82.pre.NK.samples[-c(2,6)],jm82.reinfection.samples[12])],
                            zj31.threshold[,zj31.NK.samples],
                            zg66.threshold[,zg66.NK.samples],
                            zh19.threshold[,c(zh19.NK.samples,"zh1910mCD16p.fastq","zh1915hmNKredoCD96R24.fastq")],
                            zh33.threshold[,zh33.NK.samples],
                            zk22.threshold[,c(zk22.NK.samples,"zk22_15hm_PB_NK_CD16sp_NKP80p_RNA_20170530_sampled_FX36_N721_S16_L008_R1_001.fastq")]),
                       list(jd76.time$nk.all,
                            c(jm82.time$nk.all[-c(2,6)],20),
                            c(zj31.time$nk,8.5,9.5,12,17),
                            c(zg66.time$nk,12,17,27),
                            c(zh19.time$nk,10,15.5),
                            c(zh33.time$nk,12,14,21,23,42),
                            c(zk22.time$nk,9,15.5)),
                       Monkey=c("JD76","JM82","ZJ31","ZG66","ZH19","ZH33","ZK22"))
dev.off()

###############################################################################################
#################################Figure 5######################################################
###############################################################################################


pdf("rhCMV/JM82_reinfection.pdf")
ggplot(cmv[(cmv$Sample=="Plasma"|cmv$Sample=="P")&cmv$monkey=="JM82"&cmv$Day_post_rhCMV<325,],aes(x=Day_post_rhCMV-298,y=Mean,color=monkey))+
  geom_line(size=2)+
  geom_point(size=3)+
  scale_color_manual("Monkey",values=c("blue"))+
  scale_y_continuous("rhCMV DNA copy number per uL",trans = log10_trans(),breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)),limits=c(10^.5,10^10))+
  scale_x_continuous("Days post re-rhCMV infection\n(Months post transplant)",limits=c(-3,30),breaks=c(seq(0,100,by=30)),labels=C(seq(0,100,by=30),paste("(",seq(19,22)," m)",sep="")))+theme_grey(base_size=22)
dev.off()

pdf("cell_counts/JM82_reinfection.pdf")
(jm82.cbc.reinfection[1:8,][-6,],575)
dev.off()


pdf("rhCMV/JM82_reinfection_new.pdf")
ggplot(cmv.dunbar[(cmv.dunbar$SAMPLE=="Plasma"|cmv.dunbar$SAMPLE=="P")&cmv.dunbar$MONKEY=="JM82"&cmv.dunbar$DAY<325,],aes(x=DAY-298,y=COPIES,color=MONKEY))+
  geom_line(size=2)+
  geom_point(size=3)+
  scale_color_manual("Monkey",values=c("blue"))+
  scale_y_continuous("rhCMV DNA copy number per uL",trans = log10_trans(),breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)),limits=c(10^2,10^4.5))+
  scale_x_continuous("Days post re-rhCMV infection\n(Months post transplant)",limits=c(-3,30),breaks=c(seq(0,100,by=30)),labels=C(seq(0,100,by=30),paste("(",seq(19,22)," m)",sep="")))+theme_grey(base_size=22)
dev.off()


jm82.reinfection.samples=c(jm82.reinfection.samples[10],
                           jm82.T.samples[7],
                           jm82.reinfection.samples[13],
                           jm82.reinfection.samples[14],
                           jm82.reinfection.samples[15],
                           jm82.cd56.samples[6],
                           jm82.reinfection.samples[18],
                           jm82.reinfection.samples[19],
                           jm82.reinfection.samples[20],
                           jm82.NK.samples[8],
                           jm82.reinfection.samples[23],
                           jm82.reinfection.samples[24])

pdf("heatmaps/JM82_reinfection.pdf")
temp=barcodetrackR::barcode_ggheatmap(jm82.threshold2[,jm82.reinfection.samples],label_size=10,names=c("",rep(c("18m","20m"),times=3)),cellnote_size = 6,grid=FALSE,printtable=TRUE)
input=as.data.frame(ifelse(temp[,ncol(temp)]>2*temp[,ncol(temp)-1],"red","white"))
rownames(input)=rownames(temp)
barcode_ggheatmap_bar(jm82.threshold2[,jm82.reinfection.samples],COLORS=c("white","red","white"),label_size=10,cellnote_size = 6,input=input,grid=FALSE,names=c("","9.5m","11m","18.5m","20m","9.5m","12.5m","18.5m","20m","9.5m","11.5m","17.5m","20m"))
dev.off()

pdf("Ki67/reinfection.pdf")
ggplot(jm82.ki67.df.reinfection,aes(x=time-575,y=value,color=variable,group=variable))+
  geom_line(na.rm=TRUE,size=2)+
  geom_point(size=3)+
  scale_color_manual("Cell type",labels=c("T","NK CD56+CD16-","NK CD56-CD16+"),values=c("red","blue","green"))+
  scale_y_continuous("Percent Ki67+",limits=c(0,100))+
  scale_x_continuous("Days post re-rhCMV infection\n(Months post transplant)",limits=c(-3,30*1),breaks=c(seq(0,100,by=30)),labels=C(seq(0,100,by=30),paste("(",seq(19,22)," m)",sep="")))+theme_grey(base_size=22)
dev.off()

pdf("autocorrelation/JM82_reinfection.pdf")
jm82.reinfection.samples.all=c(jm82.NK.samples,jm82.reinfection.samples[3:4])
autocorrelation(jm82.threshold2[,jm82.reinfection.samples.all],c(6,7,8,9,9.5,10.5,11,11.5,12.5,13.5,14.5,15.5,16.5,17.5,20))

dev.off()
# pdf("normalized_heatmaps/JM82_reinfection.pdf")
# NK.reinfection=jm82.reinfection.samples[20:24]
# 
# temp=rownames(barcodetrackR::barcode_ggheatmap(jm82.threshold2[,NK.reinfection[c(4,5)]],jm82.cbc.reinfection$cd16[c(1,8)],printtable=TRUE))
# input=as.data.frame(input[temp,])
# rownames(input)=temp
# barcode_ggheatmapN_bar(jm82.threshold2[,NK.reinfection[c(4,5)]],input=input,
#                   #rep(1,times=4),
#                   jm82.cbc.reinfection$cd16[c(1,8)],
#                   names=c("","18m","20m"),label_size=10,cellnote_size = 7)
# dev.off()
###############################################################################################
#################################Figure 6######################################################
###############################################################################################
colnames(jd76.cbc.klrc)=c(-8,9,20,82)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 5
cols = gg_color_hue(n)[1:4]

plot(1:n, pch = 16, cex = 2, col = cols)

cols2=gg_color_hue(n)[c(5,1:4)]
plot(1:n, pch = 16, cex = 2, col = cols2)


klrc=data.frame(variable=rep(c("ZJ31","RQ4753"),each=4),values=c(.004,.796,.142,.0575,.0037,.679,.236,.088),type=rep(c("KLRC1/NKG2A","KLRC2/NKG2C","DP","DN"),times=2))
pdf("barplots/CMVp_KLRC.pdf")

ggplot(klrc,aes(x=variable,y=values,fill=type,width=.5))+
  geom_bar(stat="identity",color="black")+
  scale_fill_manual(values=cols)+
  scale_y_continuous("Percent cell contribution",breaks=seq(0,1,by=.1),labels=paste(seq(0,100,by=10),"%"))+
  scale_x_discrete("",limits=c("ZJ31","RQ4753"))+
  theme_grey(base_size = 22)
dev.off()

jd76.cbc.klrc.plot.data2=melt(jd76.cbc.klrc[c("%NKCD16","%K1","%K2","%K1K2","%DN"),])
jd76.cbc.klrc.plot.data2$type=rep(c("NK CD16","KLRC1/NKG2A","KLRC2/NKG2C","DP","DN"),times=4)
jd76.cbc.klrc.plot.data2$variable=as.numeric(as.vector(jd76.cbc.klrc.plot.data2$variable))
dir.create("barplots")
pdf("barplots/JD76_KLRC.pdf")
ggplot(jd76.cbc.klrc.plot.data2,aes(x=variable,y=value,fill=type))+
  geom_bar(stat="identity",color="black")+
  scale_fill_manual(values=cols)+
  scale_y_continuous("Percent cell contribution",breaks=seq(0,1,by=.1),labels=paste(seq(0,100,by=10),"%"))+
  scale_x_continuous("Days post rhCMV infection",breaks=seq(0,500,by=30))+
  theme_grey(base_size = 22)
dev.off()

jd76.cbc.klrc.plot.data=melt(jd76.cbc.klrc[c("#NKCD16","#K1","#K2","#K1K2","#DN"),])
jd76.cbc.klrc.plot.data$type=rep(c("NK CD16","KLRC1/NKG2A","KLRC2/NKG2C","DP","DN"),times=4)
jd76.cbc.klrc.plot.data$variable=as.numeric(as.vector(jd76.cbc.klrc.plot.data$variable))
pdf("cell_counts/JD76_KLRC.pdf")
ggplot(jd76.cbc.klrc.plot.data,aes(x=variable,y=value*10^3,group=type,color=type))+
  geom_line(size=2)+
  scale_x_continuous("Days post rhCMV infection",breaks=seq(0,500,by=30))+theme_grey(base_size = 22)+
  scale_y_continuous("Number of cells per uL PB",limits=c(0,500))+
  geom_point(size=3)+scale_color_manual(values=cols2[c(2,3,4,5,1)])+theme_grey(base_size=22)
dev.off()

colnames(JC95.cbc.klrc)=c(-3,14,57)

jc95.cbc.klrc.plot.data2=melt(JC95.cbc.klrc[c("%K1","%K2","%K1K2","%DN"),])
jc95.cbc.klrc.plot.data2$type=rep(c("KLRC1/NKG2A","KLRC2/NKG2C","DP","DN"),times=3)
jc95.cbc.klrc.plot.data2$variable=as.numeric(as.vector(jc95.cbc.klrc.plot.data2$variable))

pdf("barplots/JC95_KLRC.pdf")
ggplot(jc95.cbc.klrc.plot.data2,aes(x=variable,y=value,fill=type))+
  geom_bar(stat="identity",color="black")+
  scale_fill_manual(values=cols)+
  scale_y_continuous("Percent cell contribution",breaks=seq(0,1,by=.1),labels=paste(seq(0,100,by=10),"%"))+
  scale_x_continuous("Days post rhCMV infection",breaks=seq(-30,500,by=30),limits=c(-15,30*2.5))+
  theme_grey(base_size = 22)
dev.off()


JC95.cbc.klrc.plot.data=melt(JC95.cbc.klrc[c("#NKCD16","#K1","#K2","#K1K2","#DN"),])
JC95.cbc.klrc.plot.data$variable=rep(c(-3,14,57),each=5)
JC95.cbc.klrc.plot.data$type=rep(c("NK CD16","KLRC1/NKG2A","KLRC2/NKG2C","DP","DN"),times=3)
JC95.cbc.klrc.plot.data$variable=as.numeric(as.vector(JC95.cbc.klrc.plot.data$variable))
pdf("cell_counts/JC95_KLRC.pdf")
ggplot(JC95.cbc.klrc.plot.data,aes(x=variable,y=value*10^3,group=type,color=type))+
  geom_line(size=2)+
  scale_x_continuous("Days post rhCMV infection",breaks=seq(0,500,by=30))+theme_grey(base_size = 22)+
  scale_y_continuous("Number of cells per uL PB",limits=c(0,500))+
  geom_point(size=3)+scale_color_manual(values=cols2[c(2,3,4,5,1)])+theme_grey(base_size=22)
dev.off()


colnames(jm82.cbc.klrc)=c(-4,14,36,183)

jm82.cbc.klrc.plot.data2=melt(jm82.cbc.klrc[c("%NKCD16","%K1","%K2","%K1K2","%DN"),])[1:20,]
jm82.cbc.klrc.plot.data2$type=rep(c("NK CD16","KLRC1/NKG2A","KLRC2/NKG2C","DP","DN"),times=4)
jm82.cbc.klrc.plot.data2$variable=as.numeric(as.vector(jm82.cbc.klrc.plot.data2$variable))
jm82.cbc.klrc.plot.data2$value[20]=.0477
pdf("barplots/JM82_KLRC.pdf")
ggplot(jm82.cbc.klrc.plot.data2,aes(x=variable,y=value,fill=type))+
  geom_bar(stat="identity",color="black")+
  scale_fill_manual(values=cols)+
  scale_y_continuous("Percent cell contribution",breaks=seq(0,1,by=.1),labels=paste(seq(0,100,by=10),"%"))+
  scale_x_continuous("Days post rhCMV nfection",breaks=seq(0,500,by=30))+
  theme_grey(base_size = 22)
dev.off()

jm82.cbc.klrc.plot.data=melt(jm82.cbc.klrc[c("#NKCD16","#K1","#K2","#K1K2","#DN"),])[1:20,]
jm82.cbc.klrc.plot.data$type=rep(c("NK CD16","KLRC1/NKG2A","KLRC2/NKG2C","DP","DN"),times=4)
jm82.cbc.klrc.plot.data$variable=as.numeric(as.vector(jm82.cbc.klrc.plot.data$variable))
pdf("cell_counts/JM82_KLRC.pdf")
ggplot(jm82.cbc.klrc.plot.data,aes(x=variable,y=value*10^3,group=type,color=type))+
  geom_line(size=2)+
  scale_x_continuous("Days post rhCMV infection",breaks=seq(0,500,by=30))+theme_grey(base_size = 22)+
  scale_y_continuous("Number of cells per uL PB",limits=c(0,500))+
  geom_point(size=3)+scale_color_manual(values=cols2[c(2,3,4,5,1)])+theme_grey(base_size=22)
dev.off()

nkg2c$condition=c("neg","neg","neg","pos","pos","pos","pos*","pos*","pos*")

pdf("barplots/nkg2c.pdf")
ggplot(nkg2c[-9,],aes(condition,value,group=Monkey,color=Monkey))+
  geom_line(size=2)+scale_x_discrete("rhCMV Status")+
  geom_point(size=5)+scale_y_continuous("Percent cell contribution",breaks=seq(0,1,by=.1),labels=paste(seq(0,100,by=10),"%"),limits=c(0,1))+
  scale_color_manual(values=c("green","red","blue","black","black"))+theme_grey(base_size=22)
dev.off()

pdf("barplots/kircmvp.pdf")
ggplot(klrc.cmv,aes(x=Monkey,y=Value,fill=Variable))+
  geom_bar(stat="identity",color="black")+
  scale_fill_manual(values=cols)+
  scale_y_continuous("Percent cell contribution",breaks=seq(0,100,by=10),labels=paste(seq(0,100,by=10),"%"))+
  theme_grey(base_size = 22)
dev.off()


pdf("barplots/kircmvn.pdf",width=10)
klrc.cmvn$label=paste(klrc.cmvn$Monkey,klrc.cmvn$KIR,klrc.cmvn$Pre.post)
ggplot(klrc.cmvn,aes(x=label,y=Value,fill=Variable))+
  geom_bar(stat="identity",color="black")+
  scale_fill_manual(values=cols)+
  scale_y_continuous("Percent cell contribution",breaks=seq(0,100,by=10),labels=paste(seq(0,100,by=10),"%"))+
  scale_x_discrete(limits=unique(klrc.cmvn$label))+
  theme_grey(base_size = 22)
dev.off()
###############################################################################################
#################################Supplemental##################################################
###############################################################################################

###############################################################################################
#################################SFigure 1#####################################################
###############################################################################################

jd76.cmv.pre.s$Monkey="JD76"
jm82.cmv.pre.s$Monkey="JM82"

cmv.pre.s=rbind(jd76.cmv.pre.s,jm82.cmv.pre.s)

pdf(file="rhCMV/negative_s.pdf")
ggplot(cmv.pre.s,aes(x=X.days,y=value,color=Monkey))+
  geom_line()+
  scale_color_manual(values=c("red","blue"))+ geom_line(linetype = 2)+
  scale_y_continuous("rhCMV DNA copy number per uL",trans = log10_trans(),breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)),limits=c(10^.5,10^4))+ theme_grey(base_size = 22)+
  scale_x_continuous("Days post transplant")
dev.off()

jd76=jd76[,c(6,1,2,3,5)]
rownames(jd76)=as.numeric(sub("h",".5",strsplit(rownames(jd76),"m")))

jm82=jm82[,c(6,1,2,3,5)]
rownames(jm82)=as.numeric(sub("h",".5",strsplit(rownames(jm82),"m")))

dir.create("GFP")
pdf("GFP/Monos.pdf")
gfp_plot(jd76,jm82,
         cell_type=1,
         monkeys=c("jd76","jm82"))+theme_grey(base_size=22)
dev.off()
pdf("GFP/Grans.pdf")
gfp_plot(jd76,jm82,
         cell_type=2,
         monkeys=c("jd76","jm82"))+theme_grey(base_size=22)
dev.off()
pdf("GFP/T.pdf")
gfp_plot(jd76,jm82,
         cell_type=3,
         monkeys=c("jd76","jm82"))+theme_grey(base_size=22)
dev.off()
pdf("GFP/B.pdf")
gfp_plot(jd76,jm82,
         cell_type=4,
         monkeys=c("jd76","jm82"))+theme_grey(base_size=22)
dev.off()
pdf("GFP/NK.pdf")
gfp_plot(jd76,jm82,
         cell_type=5,
         monkeys=c("jd76","jm82"))+theme_grey(base_size=22)
dev.off()

#CBC=melt(CBC,id.vars=c("Monkey","X.Daysinfus"))
pdf("cell_counts/PLTRBCSEGSJD76.pdf")
cbc.all(CBC[CBC$Monkey=="JD76",],0)+theme_grey(base_size=22)
dev.off()

pdf("cell_counts/PLTRBCSEGSJM82.pdf")
cbc.all(CBC[CBC$Monkey=="JM82",],0)+theme_grey(base_size=22)
dev.off()
###############################################################################################
#################################SFigure 2#####################################################
###############################################################################################

pdf("heatmaps/ZG66.pdf")
temp=barcodetrackR::barcode_ggheatmap(zg66.threshold[,zg66.all.samples],names=rep(paste(zg66.time$nk,"m"),times=4),label_size=10,printtable=TRUE)
input=as.data.frame(ifelse(temp[,ncol(temp)]>10*temp[,5],"red","black"))
rownames(input)=rownames(temp)
barcode_ggheatmap_bar(zg66.threshold[,zg66.all.samples],input=input,COLORS=c("white","white","red","white"),names=c("",rep(paste(zg66.time$nk[1:5],"m"),times=4)),label_size=10,grid=FALSE)
dev.off()
pdf("heatmaps/ZH33.pdf")
temp=barcodetrackR::barcode_ggheatmap(zh33.threshold[,zh33.all.samples[c(1:5,11:14,21:25,31:35)]],label_size = 10,names=c("",rep(paste(zh33.time$nk,"m"),5)[c(1:5,11:14,21:25,31:35)]),printtable=TRUE)
input=as.data.frame(ifelse(temp[,ncol(temp)]>10*temp[,5],"red","black"))
rownames(input)=rownames(temp)
barcode_ggheatmap_bar(zh33.threshold[,zh33.all.samples[c(1:5,11:14,21:25,31:35)]],input=input,COLORS=c("white","white","red","white"),label_size = 10,names=c("",rep(paste(zh33.time$nk,"m"),5)[c(1:5,11:14,21:25,31:35)]),grid=FALSE)
dev.off()
pdf("heatmaps/ZJ31.pdf")
temp=barcodetrackR::barcode_ggheatmap(zj31.threshold[,zj31.all.samples],names=c("",paste(zj31.all.time,"m")),label_size=10,printtable=TRUE)
input=as.data.frame(ifelse(temp[,ncol(temp)]>10*temp[,5],"red","black"))
rownames(input)=rownames(temp)
barcode_ggheatmap_bar(zj31.threshold[,zj31.all.samples],input=input,COLORS=c("white","white","red","white"),names=c("",paste(zj31.all.time,"m")),label_size=10,grid=FALSE)
dev.off()
pdf("heatmaps/ZH19.pdf")
temp=barcodetrackR::barcode_ggheatmap(zh19.threshold[,zh19.all.samples],names=c("",paste(rep(zh19.time$nk,4),"m")),label_size=10,printtable=TRUE)
input=as.data.frame(ifelse(temp[,ncol(temp)]>10*temp[,4],"red","black"))
rownames(input)=rownames(temp)
barcode_ggheatmap_bar(zh19.threshold[,zh19.all.samples],names=c("",paste(rep(zh19.time$nk,4),"m")),input=input,COLORS=c("white","white","red","white"),label_size=10,grid=FALSE)
dev.off()
pdf("heatmaps/ZK22.pdf")
temp=barcodetrackR::barcode_ggheatmap(zk22.threshold[,zk22.all.samples],names=c("",paste(rep(zk22.time$nk[1:3],5),"m")[-10]),label_size=10,printtable=TRUE)
input=as.data.frame(ifelse(temp[,ncol(temp)]>10*temp[,3],"red","black"))
rownames(input)=rownames(temp)
barcode_ggheatmap_bar(zk22.threshold[,zk22.all.samples],input=input,names=c("",paste(rep(zk22.time$nk[1:3],5),"m")[-10]),label_size=10,grid=FALSE)
dev.off()
library(ggplot2)
pdf("diversity/ZG66.pdf")
zg66.time$grans=zg66.time$t=zg66.time$b=zg66.time$nk=zg66.time$nk[1:5]
diversity(zg66.threshold[,zg66.all.samples],zg66.time)
dev.off()
pdf("diversity/ZH33.pdf")
zh33.time$grans=zh33.time$t=zh33.time$b=zh33.time$nk=zh33.time$nk[1:5]
diversity(zh33.threshold[,zh33.all.samples],zh33.time)
dev.off()
pdf("diversity/ZJ31.pdf")
zj31.time$grans=zj31.time$grans[1:5]
zj31.time$t=zj31.time$t[1:6]
zj31.time$b=zj31.time$b[1:5]
zj31.time$nk=zj31.time$nk[1:4]
zj31.time$cd56=zj31.time$cd56[1:3]
diversity2(zj31.threshold[,zj31.all.samples],zj31.time)
dev.off()
pdf("diversity/ZH19.pdf")
zh19.time$grans=zh19.time$t=zh19.time$b=zh19.time$nk
diversity(zh19.threshold[,zh19.all.samples],zh19.time)
dev.off()
pdf("diversity/ZK22.pdf")
zk22.time$grans=zk22.time$t=zk22.time$b=zk22.time$nk[1:3]
zk22.time$cd56=zk22.time$nk[2:3]
zk22.time$nk=zk22.time$nk[1:3]
diversity2(zk22.threshold[,zk22.all.samples],zk22.time)
dev.off()



#zg66.all.samples2=zg66.all.samples[c(seq(1,32,by=8),seq(2,32,by=8),seq(3,32,by=8),seq(4,32,by=8),seq(5,32,by=8),seq(6,32,by=8),seq(7,32,by=8),seq(8,32,by=8))]
pdf("biased_barcodes/ZG66.pdf")
biased_barcode_barplot(zg66.threshold[,zg66.all.samples], length(zg66.all.samples)/4, 4,zg66.time$nk[1:5], need_table = FALSE)
dev.off()
pdf("biased_barcodes/ZH33.pdf")
zh33.all.samples=zh33.all.samples[c(1:5,11:15,21:25,31:35)]
biased_barcode_barplot(zh33.threshold[,zh33.all.samples], length(zh33.all.samples)/4, 4,zh33.time$nk[1:5], need_table = FALSE)
dev.off()
pdf("biased_barcodes/ZJ31.pdf")
biased_barcode_barplot(zj31.threshold[,zj31.all.samples.complete[-seq(5,25,by=5)]], (length(zj31.all.samples.complete)-5)/5, 5,c(1,4,6,8.5), need_table = FALSE)
dev.off()
pdf("biased_barcodes/ZH19.pdf")
biased_barcode_barplot(zh19.threshold[,zh19.all.samples], (length(zh19.all.samples))/4, 4,zh19.time$nk, need_table = FALSE)
dev.off()
pdf("biased_barcodes/ZK22.pdf")
biased_barcode_barplot(zk22.threshold[,zk22.all.samples[-c(1,4,7,12)]], (length(zk22.all.samples)-4)/5, 5,zk22.time$nk[-1], need_table = FALSE)
dev.off()
###############################################################################################
#################################SFigure 3#####################################################
###############################################################################################
#bias.plot(jd76.proptable[,jd76.pre.all.samples[c(1,9,15,24,32)]])
#bias.plot(jd76.proptable[,jd76.pre.all.samples[c(3,10,17,26,34)]])
bias.plot(zj31.proptable[,zj31.all.samples.complete[c(2,7,12,17,22)]])
bias.plot(zj31.proptable[,zj31.all.samples.complete[c(3,8,13,18,23)]])
bias.plot(zj31.proptable[,zj31.all.samples.complete[c(4,9,14,19,24)]])
bias.plot(zj31.proptable[,zj31.all.samples.complete[c(4,5,10,15,20,25)]])

bias.plot(jd76.proptable[,jd76.pre.all.samples[c(4,11,18,27,35)]])
bias.plot(jd76.proptable[,jd76.pre.all.samples[c(5,12,19,29,36)]])
bias.plot(jd76.proptable[,jd76.pre.all.samples[c(6,13,20,28,37)]])
bias.plot(jd76.proptable[,jd76.pre.all.samples[c(8,14,23,31,39)]])

bias.plot(jm82.proptable[,jm82.pre.all.samples[c(4,13,22,30,38)]])
bias.plot(jm82.proptable[,jm82.pre.all.samples[c(6,14,23,31,39)]])
bias.plot(jm82.proptable[,jm82.pre.all.samples[c(8,16,25,33,40)]])
bias.plot(jm82.proptable[,jm82.pre.all.samples[c(9,18,26,34,42)]])

###############################################################################################
#################################SFigure 4#####################################################
###############################################################################################

pdf("cell_counts/supp_jd76.pdf")
cbc.all2(jd76.cbc,265)
dev.off()
pdf("cell_counts/supp_jm82.pdf")
cbc.all2(jm82.cbc,277)
dev.off()
pdf("cell_counts/supp_jc95.pdf")
cbc.all2(jc95.cbc,0)
dev.off()

barcodetrackR::barcode_ggheatmap(jd76.threshold[,cd4.cd8[c(1,3,4)]],label_size = 10,names=c("JD76 Bulk T (12m)","JD76 CD4 (12m)","JD76 Bulk CD8 (12m)"))

###############################################################################################
#################################SFigure 5#####################################################
###############################################################################################

bias.plot(jd76.proptable[,jd76.all.samples2[seq(3,40,by=8)]])
bias.plot(jd76.proptable[,jd76.all.samples2[seq(4,40,by=8)]])
bias.plot(jd76.proptable[,jd76.all.samples2[seq(5,40,by=8)]])
bias.plot(jd76.proptable[,jd76.all.samples2[seq(6,40,by=8)]])
bias.plot(jd76.proptable[,jd76.all.samples2[seq(7,40,by=8)]])
bias.plot(jd76.proptable[,jd76.all.samples2[seq(8,40,by=8)]])

bias.plot(jm82.proptable[,jm82.all.samples2[seq(2,20,by=4)]])
bias.plot(jm82.proptable[,jm82.all.samples2[seq(3,20,by=4)]])
bias.plot(jm82.proptable[,jm82.all.samples2[seq(4,20,by=4)]])

pdf("biased_barcodes/prereinfection.pdf")
jm82.proptable2=apply(jm82.threshold2,2,function(x){x/sum(x)})
bias.plot(jm82.proptable2[,jm82.reinfection.samples.complete[seq(4,25,by=5)]])
dev.off()

pdf("biased_barcodes/reinfection.pdf")
bias.plot(jm82.proptable2[,jm82.reinfection.samples.complete[seq(5,25,by=5)]])
dev.off()
###############################################################################################
#################################Other#########################################################
###############################################################################################

ggplot(cmv.dunbar[cmv.dunbar$SAMPLES=="Saliva"|cmv.dunbar$SAMPLES=="S",],aes(x=DAY,y=COPIES,color=MONKEY))+
  geom_line()+
  geom_point()+
  scale_color_manual(values=c("green","red","blue"))+
  scale_y_continuous("rhCMV DNA copy number per uL",trans = log10_trans(),breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)),limits=c(10^1,10^9))+
  scale_x_continuous("Days post rhCMV infection",limits=c(0,max(cmv$Day_post_rhCMV)))

ggplot(cmv.dunbar[cmv.dunbar$SAMPLE=="Urine"|cmv.dunbar$SAMPLES=="U",],aes(x=DAY,y=COPIES,color=MONKEY))+
  geom_line()+
  geom_point()+
  scale_color_manual(values=c("green","red","blue"))+
  scale_y_continuous("rhCMV DNA copy number per uL",trans = log10_trans(),breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)),limits=c(10^1,10^9))+
  scale_x_continuous("Days post rhCMV infection",limits=c(0,max(cmv$Day_post_rhCMV)))

dir.create("autocorrelation")
autocorrelation(jd76.threshold[,jd76.NK.samples],jd76.time$nk.all)
autocorrelation(jm82.threshold[,jm82.pre.NK.samples],jm82.time$nk.all)
autocorrelation(zj31.threshold[,zj31.NK.samples],zj31.time$nk)

pdf("autocorrelation/all.pdf",width=10)
autocorrelation.cmvpos(list(jd76.threshold[,jd76.NK.samples],
                            jm82.threshold2[,c(jm82.pre.NK.samples[-2])],
                            zj31.threshold[,zj31.NK.samples],
                            zg66.threshold[,zg66.NK.samples],
                            zh19.threshold[,c(zh19.NK.samples,"zh1910mCD16p.fastq","zh1915hmNKredoCD96R24.fastq")],
                            zh33.threshold[,zh33.NK.samples],
                            zk22.threshold[,c(zk22.NK.samples,"zk22_15hm_PB_NK_CD16sp_NKP80p_RNA_20170530_sampled_FX36_N721_S16_L008_R1_001.fastq")]),
                       list(jd76.time$nk.all,
                            c(jm82.time$nk.all[-2]),
                            c(zj31.time$nk,8.5,9.5,12,17.5),
                            zg66.time$nk,
                            c(zh19.time$nk,10,15.5),
                            zh33.time$nk,
                            c(.5,zk22.time$nk,15.5)),
                       Monkey=c("JD76","JM82","ZJ31","ZG66","ZH19","ZH33","ZK22"))
dev.off()


