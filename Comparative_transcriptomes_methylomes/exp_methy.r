library(ggplot2)
library(data.table)
library(DESeq2)
library(pheatmap)

###exp data processing
##load data
exp_all=data.frame(fread("gene_count_matrix.csv"))
row.names(exp_all)=gsub("(.*)\\.1","\\1",exp_all[,1])
exp_all=exp_all[,-1]
gene_pair=data.frame(fread("gene_pair.txt"))
gene_pair=gene_pair[,1:3]
names(gene_pair)=c("Homolog_S","Homolog_T","Similarity")

##sub genome
sub_genome=data.frame(fread("Sub_genome.csv"))
gene_position=data.frame(fread("gene_position.txt"))
gene_position=gene_position[,c(1,3,5,6)]
gene_position$Gene.ID=gsub("(.*)\\.1","\\1",gene_position$Gene.ID)
gene_position$Pos=(gene_position$Start+gene_position$End)/2
gene_position$Sub=""
gene_position=gene_position[gene_position$Reference%in%paste0("Chr",1:24),]
head(gene_position)
table(gene_position$Reference)
for(i in 1:nrow(gene_position)){
  sub=sub_genome[sub_genome$Chromsome==gene_position[i,"Reference"],]
  if(nrow(sub)==1){
    gene_position[i,"Sub"]=sub$Subgenome
  }else{
    gene_position[i,"Sub"]=sub[findInterval(gene_position[i,"Pos"],sub[,"Start"]),"Subgenome"]
  }
}
table(gene_position$Sub)
S_genome_gene=gene_position[gene_position$Sub=="S","Gene.ID"]
T_genome_gene=gene_position[gene_position$Sub=="T","Gene.ID"]

##Diff exp
#S vs Ssub
S_all=exp_all[S_genome_gene,c(1:3,7:9)]
condition=factor(rep(c("S","ZY300"),each=3))
colData <- data.frame(row.names=colnames(S_all), condition)
dds <- DESeqDataSetFromMatrix(S_all, DataFrame(condition), design= ~ condition)
dds <- DESeq(dds) 
sizeFactors(dds)
res <- results(dds)
summary(res)
S_Ssub <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)
#T vs Tsub
T_all=exp_all[T_genome_gene,c(4:6,7:9)]
condition=factor(rep(c("T","ZY300"),each=3))
colData <- data.frame(row.names=colnames(T_all), condition)
dds <- DESeqDataSetFromMatrix(T_all, DataFrame(condition), design= ~ condition)
dds <- DESeq(dds) 
sizeFactors(dds)
res <- results(dds)
summary(res)
T_Tsub <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)
#S vs T
a=gene_pair
a$num=paste0(gene_pair$Homolog_S,"_",gene_pair$Homolog_T)
b=a[duplicated(a$num),]
table(duplicated(a$num))
a=a[!duplicated(a$num),]
table(duplicated(a$Homolog_S))
table(duplicated(a$Homolog_T))
S_all=exp_all[unique(gene_pair$Homolog_S),c(1:3)]
T_all=exp_all[unique(gene_pair$Homolog_T),c(4:6)]
data=cbind(S_all[a$Homolog_S,],T_all[a$Homolog_T,])
row.names(data)=a$num
condition=factor(rep(c("S","T"),each=3))
colData <- data.frame(row.names=colnames(data), condition)
dds <- DESeqDataSetFromMatrix(data, DataFrame(condition), design= ~ condition)
dds <- DESeq(dds) 
sizeFactors(dds)
res <- results(dds)
summary(res)
S_T <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)
#Ssub vs Tsub
a=gene_pair
a$num=paste0(gene_pair$Homolog_S,"_",gene_pair$Homolog_T)
b=a[duplicated(a$num),]
table(duplicated(a$num))
a=a[!duplicated(a$num),]
table(duplicated(a$Homolog_S))
table(duplicated(a$Homolog_T))
S_all=exp_all[unique(gene_pair$Homolog_S),c(7:9)]
T_all=exp_all[unique(gene_pair$Homolog_T),c(7:9)]
data=cbind(S_all[a$Homolog_S,],T_all[a$Homolog_T,])
row.names(data)=a$num
names(data)=c("S1","S2","S3","T1","T2","T3")
rep=factor(rep(c("1","2","3"),2))
sub=factor(rep(c("S","T"),each=3))
dds <- DESeqDataSetFromMatrix(data, DataFrame(rep,sub), design= ~ rep+sub)
dds <- DESeq(dds) 
sizeFactors(dds)
res <- results(dds)
summary(res)
Ssub_Tsub <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)

### methy data processing
library(EnrichedHeatmap)
library(GenomicFeatures)
library(circlize)
library(ChIPseeker)
library(rtracklayer)
library(methylKit)
library(data.table)
library(RcppRoll)
## gene methylation
gene=data.frame(fread("/data/public/hanyu/Q7/exp_data/TSV/RNA-S-1.tsv"))
gene=gene[,c(1,3,5,6,4)]
gene[,1]=gsub("(.*)\\.1","\\1",gene[,1])
head(gene)
names(gene)=c("gene","chr","start","end","strand")
row.names(gene)=gene$gene
gene$start=gene$start-2000
gene$end=gene$end+2000
region=GRanges(seqnames=gene$chr,ranges=IRanges(start=gene$start,end=gene$end),strand=gene$strand)

CG_file=c()
CHG_file=c()
CHH_file=c()
for(i in list.dirs("extracted")[2:10]){
    CG_file=c(CG_file,paste0(i,"/",list.files(i)[grep("*CG_methykit.txt",list.files(i))]))
    CHG_file=c(CHG_file,paste0(i,"/",list.files(i)[grep("*CHG_methykit.txt",list.files(i))]))
    CHH_file=c(CHH_file,paste0(i,"/",list.files(i)[grep("*CHH_methykit.txt",list.files(i))]))
}
CG_all=methRead(as.list(CG_file),
    sample.id=list("S-1","S-2","S-3","T-1","T-2","T-3","ZY300-1","ZY300-2","ZY300-3"),
    treatment=c(1,1,1,1,1,1,1,1,1),
    assembly="ZY300",header=T, context="CpG", resolution="base", mincov = 3)

CHG_all=methRead(as.list(CHG_file),
    sample.id=list("S-1","S-2","S-3","T-1","T-2","T-3","ZY300-1","ZY300-2","ZY300-3"),
    treatment=c(1,1,1,1,1,1,1,1,1),
    assembly="ZY300",header=T, context="CHG", resolution="base", mincov = 3)

CHH_all=methRead(as.list(CHH_file),
    sample.id=list("S-1","S-2","S-3","T-1","T-2","T-3","ZY300-1","ZY300-2","ZY300-3"),
    treatment=c(1,1,1,1,1,1,1,1,1),
    assembly="ZY300",header=T, context="CHH", resolution="base", mincov = 3)

gene_CG_sub=regionCounts(CG_all,region,strand.aware=T)
gene_CHG_sub=regionCounts(CHG_all,region,strand.aware=T)
gene_CHH_sub=regionCounts(CHH_all,region,strand.aware=T)


gene_CG=gene
gene_CHG=gene
gene_CHH=gene
names=c("S-1","S-2","S-3","T-1","T-2","T-3","ZY300-1","ZY300-2","ZY300-3")
all_gene=gene$gene
names(all_gene)=paste0(gene[,2],"_",gene[,3],"_",gene[,4])
for(i in 1:9){
    print(i)
    data=getData(gene_CG_sub[[i]])
    row.names(data)=paste0(data[,1],"_",data[,2],"_",data[,3])
    data$rate=data[,6]/(data[,6]+data[,7])
    data=data[row.names(data)%in%names(all_gene),]
    name=names[i]
    gene_CG[,name]=NA
    gene_CG[all_gene[row.names(data)],name]=data$rate
    
    data=getData(gene_CHG_sub[[i]])
    row.names(data)=paste0(data[,1],"_",data[,2],"_",data[,3])
    data$rate=data[,6]/(data[,6]+data[,7])
    data=data[row.names(data)%in%names(all_gene),]
    name=names[i]
    gene_CHG[,name]=NA
    gene_CHG[all_gene[row.names(data)],name]=data$rate

    data=getData(gene_CHH_sub[[i]])
    row.names(data)=paste0(data[,1],"_",data[,2],"_",data[,3])
    data$rate=data[,6]/(data[,6]+data[,7])
    data=data[row.names(data)%in%names(all_gene),]
    name=names[i]
    gene_CHH[,name]=NA
    gene_CHH[all_gene[row.names(data)],name]=data$rate
}
S_CG=apply(gene_CG[S_genome_gene,6:8],1,function(x){mean(x,na.rm=T)})
Sinzy_CG=apply(gene_CG[S_genome_gene,13:14],1,function(x){mean(x,na.rm=T)})
S_CHG=apply(gene_CHG[S_genome_gene,6:8],1,function(x){mean(x,na.rm=T)})
Sinzy_CHG=apply(gene_CHG[S_genome_gene,13:14],1,function(x){mean(x,na.rm=T)})
S_CHH=apply(gene_CHH[S_genome_gene,6:8],1,function(x){mean(x,na.rm=T)})
Sinzy_CHH=apply(gene_CHH[S_genome_gene,13:14],1,function(x){mean(x,na.rm=T)})

## profile
gene_pair=data.frame(fread("homologous_info.txt"))
gene_pair=gene_pair[,c(1,2,3)]

gr <- import("ZY300.gff3")
gr$Name=gr$ID
txdb <- makeTxDbFromGRanges(gr)
genebody <- getBioRegion(TxDb = txdb,by = "gene",type = "body")

names(genebody)=gsub("(.*)\\.1","\\1",names(genebody))
genebody$gene_id=gsub("(.*)\\.1","\\1",genebody$gene_id)
genebody=genebody[unique(c(gene_pair$Homolog_S,gene_pair$Homolog_T))]

methy_profile=list()
for(i in c(1:9)){
    cat(paste(i,"is profiling \n"))
    list_name=c("S_1","S_2","S_3","T_1","T_2","T_3","ZY300_1","ZY300_2","ZY300_3")

    data=data.frame(CG_all[i])
    meth=data[,6]/data[,5]
    data=makeGRangesFromDataFrame(data[,c(1:4)])
    data$meth=meth
    mat = normalizeToMatrix(data,genebody, value_column = "meth", mean_mode = "absolute",extend = 2000, w = 50, background = NA)
    row.names(mat)=names(genebody)
    methy_profile[[paste0("CG","_",list_name[i])]]=mat

    data=data.frame(CHG_all[i])
    meth=data[,6]/data[,5]
    data=makeGRangesFromDataFrame(data[,c(1:4)])
    data$meth=meth
    mat = normalizeToMatrix(data,genebody, value_column = "meth", mean_mode = "absolute",extend = 2000, w = 50, background = NA)
    row.names(mat)=names(genebody)
    methy_profile[[paste0("CHG","_",list_name[i])]]=mat

    data=data.frame(CHH_all[i])
    meth=data[,6]/data[,5]
    data=makeGRangesFromDataFrame(data[,c(1:4)])
    data$meth=meth
    mat = normalizeToMatrix(data,genebody, value_column = "meth", mean_mode = "absolute",extend = 2000, w = 50, background = NA)
    row.names(mat)=names(genebody)
    methy_profile[[paste0("CHH","_",list_name[i])]]=mat
}

### analysis
## Sfig exp change
par(mfrow=c(2,2))
hist(S_Ssub$log2FoldChange,xlim=c(-10,10),
     breaks=200,col="tomato",border = "tomato",
     main="S vs Ssub",xlab=expression("Log"["2"]*"(Fold Change)"))
hist(T_Tsub$log2FoldChange,xlim=c(-10,10),
     breaks=200,col="tomato",border = "tomato",
     main="T vs Tsub",xlab=expression("Log"["2"]*"(Fold Change)"))
hist(S_T$log2FoldChange,xlim=c(-10,10),
     breaks=200,col="tomato",border = "tomato",
     main="S vs T",xlab=expression("Log"["2"]*"(Fold Change)"))
hist(Ssub_Tsub$log2FoldChange,xlim=c(-10,10),
     breaks=200,col="tomato",border = "tomato",
     main="Ssub vs Tsub",xlab=expression("Log"["2"]*"(Fold Change)"))
##2AB
diff_s=na.omit(S_Ssub$Row.names[S_Ssub$padj<0.05&abs(S_Ssub$log2FoldChange)>2])
diff_t=na.omit(T_Tsub$Row.names[T_Tsub$padj<0.05&abs(T_Tsub$log2FoldChange)>2])
sig_s=diff_s
sig_t=diff_t
up_s=sig_s[S_Ssub[sig_s,"log2FoldChange"]>0]
down_s=sig_s[S_Ssub[sig_s,"log2FoldChange"]<0]
up_t=sig_t[T_Tsub[sig_t,"log2FoldChange"]>0]
down_t=sig_t[T_Tsub[sig_t,"log2FoldChange"]<0]


volcan_fun=function(data,gene,main,...){
  plot(data$log2FoldChange,-log10(data$padj),col="gray",pch=20,main=main,...)
  points(data[gene,"log2FoldChange"],-log10(data[gene,"padj"]),col="tomato",pch=20)
  abline(h=-log10(0.05),col="black",lty="dashed")
  abline(v=c(-2,2),col="black",lty="dashed")
}

volcan_fun(data=S_Ssub,gene=diff_s,main="S vs Sub",xlim=c(-20,20),
    bty="l",xlab=expression("Log"["2"]*"(Fold Change)"),ylab=expression("-Log"["10"]*"(FDR)"))

volcan_fun(data=T_Tsub,gene=diff_t,main="T vs Tub",xlim=c(-20,20),
    bty="l",xlab=expression("Log"["2"]*"(Fold Change)"),ylab=expression("-Log"["10"]*"(FDR)"))

##2CDEF

genes=up_s
data=data.frame(data=c(S_CG[genes],Sinzy_CG[genes],S_CHG[genes],Sinzy_CHG[genes],S_CHH[genes],Sinzy_CHH[genes]),
    names=rep(c("CG","CHG","CHH"),each=length(genes)*2),
    type=rep(rep(c("S","S_sub"),each=length(genes)),3))
ggplot(data,aes(x=names,y=data,fill=type))+
  geom_boxplot(alpha=1,width=0.6,position=position_dodge(width=0.8),
               size=0.75,outlier.color = "white",outlier.size=0.01)+
  theme_bw()+ylim(c(0,1.1))+ ylab("Gene methylation")+xlab("")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_fill_manual(values=c(S="#FDE3E1",S_sub="#E4EECB"))+
    geom_signif(annotations = c(rep("**",3)),
        y_position = c(1.05,1.05,0.5),
        xmin = c(0.8,1.8,2.8),
        xmax = c(1.2,2.2,3.2),
        tip_length = c(0.03,0.03,0.03),
        color="black")


genes=down_s
data=data.frame(data=c(S_CG[genes],Sinzy_CG[genes],S_CHG[genes],Sinzy_CHG[genes],S_CHH[genes],Sinzy_CHH[genes]),
    names=rep(c("CG","CHG","CHH"),each=length(genes)*2),
    type=rep(rep(c("S","S_sub"),each=length(genes)),3))
ggplot(data,aes(x=names,y=data,fill=type))+
  geom_boxplot(alpha=1,width=0.6,position=position_dodge(width=0.8),
               size=0.75,outlier.color = "white",outlier.size=0.01)+
  theme_bw()+ylim(c(0,1.1))+ ylab("Gene methylation")+xlab("")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_fill_manual(values=c(S="#FDE3E1",S_sub="#E4EECB"))+
    geom_signif(annotations = c(rep("**",3)),
        y_position = c(1.05,1.05,0.5),
        xmin = c(0.8,1.8,2.8),
        xmax = c(1.2,2.2,3.2),
        tip_length = c(0.03,0.03,0.03),
        color="black")



genes=up_t
data=data.frame(data=c(T_CG[genes],Tinzy_CG[genes],T_CHG[genes],Tinzy_CHG[genes],T_CHH[genes],Tinzy_CHH[genes]),
    names=rep(c("CG","CHG","CHH"),each=length(genes)*2),
    type=rep(rep(c("T","T_sub"),each=length(genes)),3))
ggplot(data,aes(x=names,y=data,fill=type))+
  geom_boxplot(alpha=1,width=0.6,position=position_dodge(width=0.8),
               size=0.75,outlier.color = "white",outlier.size=0.01)+
  theme_bw()+ylim(c(0,1.1))+ ylab("Gene methylation")+xlab("")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_fill_manual(values=c(T="#CBF1F2",T_sub="#F3E4FE"))+
    geom_signif(annotations = c(rep("**",3)),
        y_position = c(1.05,1.05,0.5),
        xmin = c(0.8,1.8,2.8),
        xmax = c(1.2,2.2,3.2),
        tip_length = c(0.03,0.03,0.03),
        color="black")


genes=down_t
data=data.frame(data=c(T_CG[genes],Tinzy_CG[genes],T_CHG[genes],Tinzy_CHG[genes],T_CHH[genes],Tinzy_CHH[genes]),
    names=rep(c("CG","CHG","CHH"),each=length(genes)*2),
    type=rep(rep(c("T","T_sub"),each=length(genes)),3))
ggplot(data,aes(x=names,y=data,fill=type))+
  geom_boxplot(alpha=1,width=0.6,position=position_dodge(width=0.8),
               size=0.75,outlier.color = "white",outlier.size=0.01)+
  theme_bw()+ylim(c(0,1.1))+ ylab("Gene methylation")+xlab("")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_fill_manual(values=c(T="#CBF1F2",T_sub="#F3E4FE"))+
    geom_signif(annotations = c(rep("**",3)),
        y_position = c(1.05,1.05,0.5),
        xmin = c(0.8,1.8,2.8),
        xmax = c(1.2,2.2,3.2),
        tip_length = c(0.03,0.03,0.03),
        color="black")

##2GK
fold=2
pthr=0.05
sig_zy=na.omit(Ssub_Tsub[abs(Ssub_Tsub$log2FoldChange)>fold&Ssub_Tsub$padj<pthr,1])
sig_st=na.omit(S_T[abs(S_T$log2FoldChange)>fold&S_T$padj<pthr,1])
data2=list("Sig in ZY300"=sig_zy,"Sig in ST"=sig_st)
ggvenn(data2,set_name_size=4,fill_color=c("skyblue","tomato"))

cols=ifelse(abs(S_T$log2FoldChange)<2&abs(Ssub_Tsub$log2FoldChange)<2,"gray","tomato")
plot(S_T$log2FoldChange,Ssub_Tsub$log2FoldChange,pch=20,col=cols,bty="l",xlab="S vs T",ylab="Ssub vs Tsub",xlim=c(-21,21),ylim=c(-21,21),fg="black")
abline(h=c(-2,2),col="black",lty="dashed")
abline(v=c(-2,2),col="black",lty="dashed")
abline(a=0,b=1,col="black",lty="dashed")
text(x=c(-18,0,18,18,18,0,-18,-18),y=c(15,15,15,0,-15,-15,-15,0),labels=c(1:8),cex=2.5)

##2HJNO
homo_up_s=na.omit(S_Ssub$Row.names[S_Ssub$padj<0.05&S_Ssub$log2FoldChange>2])
homo_up_s=homo_up_s[homo_up_s%in%gene_pair$Homolog_S]
homo_down_s=na.omit(S_Ssub$Row.names[S_Ssub$padj<0.05&S_Ssub$log2FoldChange< -2])
homo_down_s=homo_down_s[homo_down_s%in%gene_pair$Homolog_S]

homo_up_t=na.omit(T_Tsub$Row.names[T_Tsub$padj<0.05&T_Tsub$log2FoldChange>2])
homo_up_t=homo_up_t[homo_up_t%in%gene_pair$Homolog_T]
homo_down_t=na.omit(T_Tsub$Row.names[T_Tsub$padj<0.05&T_Tsub$log2FoldChange< -2])
homo_down_t=homo_down_t[homo_down_t%in%gene_pair$Homolog_T]

length(homo_up_s);length(homo_down_s);length(homo_up_t);length(homo_down_t)
length(gene_pair[gene_pair[,1]%in%homo_up_s,4]);length(gene_pair[gene_pair[,1]%in%homo_down_s,4])
length(gene_pair[gene_pair[,2]%in%homo_up_t,4]);length(gene_pair[gene_pair[,2]%in%homo_down_t,4])
homo_up_s=gene_pair[gene_pair[,1]%in%homo_up_s,4]
homo_down_s=gene_pair[gene_pair[,1]%in%homo_down_s,4]
homo_up_t=gene_pair[gene_pair[,2]%in%homo_up_t,4]
homo_down_t=gene_pair[gene_pair[,2]%in%homo_down_t,4]

data2=list("Ssub up"=homo_up_s,"S up"=homo_down_s,"Tsub up"=homo_up_t,"T up"=homo_down_t)
ggvenn(data2,set_name_size=3,text_size=3,fill_color=c("tomato","skyblue","lightgreen","#e3eb87"))

gene=intersect(homo_up_s,homo_down_t)
par(mfrow=c(1,1))
par(mar=c(2,4,2,0))
boxplot(
  S_T[gene,"mean_1"],
  Ssub_Tsub[gene,"mean_1"],
  S_T[gene,"mean_2"],
  Ssub_Tsub[gene,"mean_2"],
  ylim=c(0,2000),main="",ylab="Expression",names=c("S","Ssub","T","Tsub"),boxwex = 0.3,staplewex=0,
  frame.plot="F",col=c("#FDE3E1","#E4EECB","#CBF1F2","#F3E4FE"),outline=F,at=c(1,1.7,2.4,3.1)
)

methy_plot_fun=function(methy_profile,test="two.sided",genes,type1,type2,main,...){
  if(length(grep("ZY300",type1))==1){
    data_1=paste0(type1,c("_1","_2"))
  }else{
    data_1=paste0(type1,c("_1","_2","_3"))
  }
  if(length(grep("ZY300",type2))==1){
    data_2=paste0(type2,c("_1","_2"))
  }else{
    data_2=paste0(type2,c("_1","_2","_3"))
  }
  
  type1_data=data.frame(matrix(ncol=133,nrow=length(data_1)))
  row.names(type1_data)=data_1
  for(i in data_1){
    sub=methy_profile[[i]]
    sub=sub[genes,]
    type1_data[i,]=apply(sub,2,function(x){mean(x,na.rm=T)})
  }
  
  type2_data=data.frame(matrix(ncol=133,nrow=length(data_2)))
  row.names(type2_data)=data_2
  for(i in data_2){
    sub=methy_profile[[i]]
    sub=sub[genes,]
    type2_data[i,]=apply(sub,2,function(x){mean(x,na.rm=T)})
  }
  
  type1_data=apply(type1_data,2,function(x){mean(x,na.rm=T)})
  type2_data=apply(type2_data,2,function(x){mean(x,na.rm=T)})
  plot(type1_data,type="l",col="gray",frame.plot=F,xaxt="n",xlab="",main=main,ylim=c(0,1),lwd=2,...)
  axis(1,c(1,40,93,133),c("-2k","TSS","TES","2k"),las=1,cex.axis=1)
  lines(type2_data,col="tomato",lwd=2)
  legend(x=2,y=1, c(type1,type2),col=c("gray","tomato"),lwd=3.5,border="white",box.lty=0)
  cat(wilcox.test(type1_data,type2_data,alternative=test)$p.value)
}

methy_plot_fun(methy_profile = methy_profile,test="greater",genes=str_split_fixed(gene,"_",2)[,1],
               type1="CHG_S",type2="CHG_ZY300",main="",ylab="CHG methylation")
methy_plot_fun(methy_profile = methy_profile,test="less",genes=str_split_fixed(gene,"_",2)[,2],
               type1="CHG_T",type2="CHG_ZY300",main="",ylab="CHG methylation")