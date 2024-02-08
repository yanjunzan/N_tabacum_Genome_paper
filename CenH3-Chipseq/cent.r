blast=data.frame(fread("quartet_blast"))
blast=blast[,c(1,2,3,11)]
blast$chr=substring(blast$V1,1,5)
blast$TR=paste0(str_split_fixed(blast[,1],"_",3)[,1],"_",str_split_fixed(blast[,1],"_",3)[,2])
table(blast$chr,blast[,2])
table(blast$chr)

genome=data.frame(fread("genome.fasta.fai"))
genome=genome[1:24,1:2]
sum(genome[1:24,2])
Cent=data.frame(fread("cent_region.txt"))
Cent=Cent[,1:3]

best=data.frame(fread("TR_region.txt"))
best=best[,2:5]
best$chr=str_split_fixed(best[,1],"@",2)[,1]
head(best)

chip_1=data.frame(fread(paste0("deep/",1,".","5k","bin.rev.bedgraph")))
chip_2=data.frame(fread(paste0("deep/",2,".","5k","bin.rev.bedgraph")))

par(mfrow=c(6,4))
par(mar=c(1,2,2,1))
mark=1
for(i in genome[,1]){
    print(i)
    TR=data.frame(fread(paste0("T.",i,".tr.blast.gff3")))   
    TR$TR=sub("ID=","",str_split_fixed(TR[,9],";",2)[,1])
    
    sub_blast=best[best$chr==i,]
    get_TR=unique(sub_blast[,1])
    sub=TR[TR$TR%in%get_TR,]
    chr=seq(0,genome[genome[,1]==i,2],100000)
    height=rep(0,length(chr))
    pos=findInterval(sub[,4],chr)
    height[as.numeric(names(table(pos)))]=table(pos)
    
    sub_blast=blast[blast$chr==i,]
    get_TR=unique(sub_blast$TR)
    sub=TR[TR$TR%in%get_TR,]
    chr2=seq(0,genome[genome[,1]==i,2],100000)
    height2=rep(0,length(chr2))
    pos=findInterval(sub[,4],chr2)
    height2[as.numeric(names(table(pos)))]=table(pos)

    sub=chip_1[chip_1[,1]==i,]
    sub=sub[order(sub[,2]),]
    pos1=roll_mean(sub[,2],n=200,by=200)
    sub1=roll_mean(sub[,4],n=200,by=200)
    sub=chip_2[chip_2[,1]==i,]
    sub=sub[order(sub[,2]),]
    pos3=roll_mean(sub[,2],n=200,by=200)
    sub3=roll_mean(sub[,4],n=200,by=200)
    
    wd=1.5
    plot(x=pos1,y=sub1,col="skyblue",type="l",xaxt="n",frame.plot=F,ylim=c(-(abs(min(sub1,sub3))+max(sub1,sub3))*2,max(sub1,sub3)),
        ylab="",xlab="",main="",xlim=c(0,genome[genome[,1]==i,2]),yaxt="n",lwd=1.5)
    lines(x=pos3,y=sub3,col="tomato",lwd=1.5)
    par(new = TRUE)
    plot(x=chr,y=height,xlim=c(0,genome[genome[,1]==i,2]),ylab="",xlab="",
        ,main="",type = "l", xaxt = "n",ylim=c(-max(height),max(height)*2),col="black",frame.plot=F,yaxt="n",lwd=0.3)
    if(nrow(sub_blast)!=0){
        par(new = TRUE)
        lines(chr2,height2,col="red",lwd=0.3)
    }
    par(new = TRUE)
    plot(x=c(0,genome[genome[,1]==i,2]),y=c(0,3),main="",ylab="",xlab="",xaxt="n",yaxt="n",col="white",frame.plot=T)
    axis(2,at=c(0.5,1.5,2.5),labels=c("Predicted","TR","CENH3"),lwd=0,line=-1,cex.axis=0.7)
    rect(0,0.3,genome[genome[,1]==i,2],0.6,col=scales::alpha("black", alpha=0.05),angle = c(90,90, 30, -30),lwd=0.3)   
    rect(Cent[Cent[,1]==i,2],0.3,Cent[Cent[,1]==i,3],0.6,col=scales::alpha("red", alpha=0.2),lty="dashed",lwd=0.3)
    mtext(paste("Chromosome",mark), side = 1, line = 0,cex=0.7)
    mark=mark+1
}
