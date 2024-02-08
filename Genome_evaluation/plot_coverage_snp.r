###ZY300 as example
library(data.table)
library(RcppRoll)
grep_fun=function(lai,class){
    head(lai)
    lai <- lai[!grepl("scaffold",lai$V1),]
    lai <- lai[!grepl("Mt",lai$V1),]
    lai <- lai[!grepl("Pt",lai$V1),]
    
    if(class=="ZY300"){
        chrlen <- fread("ZY300.fa.fai")
        chrlen <- chrlen[!grepl("scaffold",chrlen$V1),]
        chr <- gsub("Chr(.*)","\\1",chrlen$V1)
        chrlen <- chrlen[order(as.numeric(chr)),]
        chrlen <- chrlen[1:24,]
        len <- c(0)
        for (i in 2:24){
            len <- c(len,sum(chrlen$V2[1:(i-1)]))
        }
        names(len) <- chrlen$V1
        lai$linear <- lai$V2 + len[lai$V1]
        lai$V1 <- as.numeric(gsub("Chr(.*)","\\1",lai$V1))
        #lai <- lai[order(as.numeric(lai$V1),lai$V2),]
        names(lai)=c("chr","position","depth","linear")
        return(list(lai,len))  
    }
}

### coverage
data <- fread("zy300.depth.all.txt")
data_plot=grep_fun(lai=data,class="ZY300")
len=data_plot[[2]]
data_plot=data_plot[[1]]
tic <- aggregate(data_plot$linear,list(data_plot$chr),mean)
tic <- tic[order(as.numeric(tic$Group.1)),][1:24,]
deep=roll_mean(x=data_plot$depth,n=1000000,by=500000)
position=roll_mean(x=data_plot$linear,n=1000000,by=500000)
cols=ifelse(findInterval(position,len)%%2==1,"skyblue","tomato")
chrs=findInterval(position,len)
plot(x=c(0,max(position)),y=c(0,max(deep)),type="l",xaxt="n",frame.plot=F,ylim=c(0,100),
    col="white",ylab="Coverage",xlab="Chr",main="")
for(i in names(table(chrs))){
    sub=chrs==i
    lines(x=position[sub],y=deep[sub],col=cols[sub])
}
smooth_data=smooth.spline(x=position,y=deep)
lines(x=smooth_data$x,y=smooth_data$y,col="black",lwd=1)
abline(v=len,col="gray",lty="dashed",lwd=1)
axis(1,at=tic$x,labels = 1:24)

### snp
CF=data.frame(fread("ZY300.raw.vcf"))
names(VCF)[1:2]=c("V1","V2")
indel=VCF[grepl("INDEL",VCF[,8]),1:2]
hete_snp=VCF[!grepl("INDEL",VCF[,8])&grepl("0/1",VCF[,10]),1:2]
homo_snp=VCF[!grepl("INDEL",VCF[,8])&grepl("1/1",VCF[,10]),1:2]
indel$type="indel"
hete_snp$type="hete"
homo_snp$type="homo"
snps=rbind(indel,hete_snp,homo_snp)
ZY300=grep_fun(lai=snps,class="ZY300")
len=ZY300[[2]]
ZY300=ZY300[[1]]
ZY300$value=1
tic <- aggregate(ZY300$linear,list(ZY300$chr),mean)
tic <- tic[order(as.numeric(tic$Group.1)),][1:24,]
par(mfrow=c(4,1))
par(mar=c(0,2,0,0))
plot(ZY300[ZY300$depth=="indel","linear"],ZY300[ZY300$depth=="indel","value"],type="h",col=scales::alpha("tomato", alpha=0.01),
    ylim=c(0,1),xlim=c(0,genome_length),ylab="",xlab="",xaxt="n",yaxt="n",frame.plot=F)
axis(2,las=2,at=0.5,labels = "INDEL",lwd=0,cex.axis=0.6,line=-2.4)
abline(v=len[-1],col="gray",lty="dashed",lwd=0.5)
plot(ZY300[ZY300$depth=="homo","linear"],ZY300[ZY300$depth=="homo","value"],type="h",col=scales::alpha("skyblue", alpha=0.05),
    ylim=c(0,1),xlim=c(0,genome_length),ylab="",xlab="",xaxt="n",yaxt="n",frame.plot=F)
axis(2,las=2,at=0.5,labels = "Homo SNPs",lwd=0,cex.axis=0.6,line=-2.4)
abline(v=len[-1],col="gray",lty="dashed",lwd=0.5)
plot(ZY300[ZY300$depth=="hete","linear"],ZY300[ZY300$depth=="hete","value"],type="h",col=scales::alpha("#a8eea8", alpha=0.05),
    ylim=c(0,1),xlim=c(0,genome_length),ylab="",xlab="",xaxt="n",yaxt="n",frame.plot=F)
axis(2,las=2,at=0.5,labels = "Heter SNPs",lwd=0,cex.axis=0.6,line=-2.4)
abline(v=len[-1],col="gray",lty="dashed",lwd=0.5)
axis(1,at=tic$x,labels = 1:24,cex.axis=0.6,lwd=0.5)
