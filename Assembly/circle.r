library(circlize)
library(data.table)
library(RcppRoll)
library(RColorBrewer)

### data preparing ----------------------------------------------------------------
samtools faidx genome.fasta 
bedtools makewindows -g genome.fasta.fai -w 100000 > 100k_win
grep -w "gene" final.gene.gff3 |awk '{print $1"\t"$4"\t"$5}'|uniq > gene.pos
bedtools intersect -a 100k_win -b gene.pos -c >gene_density.100k.txt
grep -w "ID" genome.fasta.mod.EDTA.TEanno.gff3 |awk '{print $1"\t"$4"\t"$5}'|uniq > TE.pos
bedtools intersect -a 100k_win -b TE.pos -c >TE_density.100k.txt
seqkit sliding -s 100000 -W 100000 genome.fasta |seqkit fx2tab -n -g >gc.txt

###plot ----------------------------------------------------------------
## ZY300 ----------------------------------------------------------------
setwd("ZY300")
circos.clear()
chrlen <- read.table("genome.fa.fai",sep="\t")
chrlen <- chrlen[1:24,]
ref <- data.frame("chr" = chrlen$V1, "start" = 0,"end"=chrlen$V2)
rownames(ref)<- ref$chr
chr_order <- paste0("Chr",c(21,22,7,8,1,5,20,16,6,10,11,18,3,17,9,13,2,4,12,15,24,23,19,14))
ref <- ref[chr_order,]

load("circos-data.RData")
load(file="avg_S_T.RData")
out2 <- read.table("linkfile.txt",sep="\t",stringsAsFactors = F)
# panel ----------------------------------------------------------------

pdf(file = "ZY300.pdf",width = 4,height = 4)
circos.clear()
circos.par(gap.degree= c(rep(1,4),8,rep(1,19)))
circos.initialize(factors=chr_order,xlim =cbind(0,ref$end/1e6) )
circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
  chr=gsub("Chr(.*)","\\1",CELL_META$sector.index)
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
  circos.text(mean(xlim), mean(ylim), chr, cex=0.4, col="black", 
              facing="bending.inside", niceFacing=TRUE)
}, bg.col="grey95", bg.border=F, track.height=0.05)

l.c <- 0.44
for (i in 1:nrow(ref)){
  mt <-seq(0,ref$end[i]/1e6,by=40)
  mt=c(mt,max(mt)+40)
  circos.axis(h = "top",sector.index = ref$chr[i],track.index = 1,labels.cex=l.c,major.at = mt,
    labels.facing = "clockwise",direction="outside",minor.ticks = 1,labels=mt,major.tick.length=0.8,lwd=0.4)
}
# Gene density ----------------------------------------------------------------
gd=data.frame(fread("ZY300/gene_density.100k.txt"))
gd <- gd[!grepl("scaffold",gd[,1]),]
colnames(gd) <- c("chr","start","end","value1")
gd$start <- as.numeric(gd$start )/1e6
gd$end <- as.numeric(gd$end )/1e6
head(gd)

gd_10win=data.frame()
window=20
for(i in 1:nrow(ref)){
  gd_now <-gd[gd$chr == ref$chr[i],]
  sub=data.frame("start"=roll_mean(x=gd_now$start,n=window,by=window),
    "end"=roll_mean(x=gd_now$end,n=window,by=window),
    "value1"=roll_sum(x=gd_now$value1,n=window,by=window))
  sub$start=sub$start-(window*7/100)/2
  sub$end=sub$end+(window*7/100)/2
  sub=cbind(gd_now[1,1],sub)
  gd_10win=rbind(gd_10win,sub)
}
head(gd_10win)
color_assign <- colorRamp2(breaks = c(min(gd_10win[,4]), mean(gd_10win[,4]), max(gd_10win[,4])), 
                           col =c("#d4eaf3ec","#5ac2eb","blue"))
circos.genomicTrackPlotRegion(
  gd_10win, track.height = 0.08, stack = TRUE, bg.border = NA,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)
  } )

# TE density ----------------------------------------------------------------
gc=data.frame(fread("ZY300/TE_density.100k.txt"))
gc <- gc[!grepl("scaffold",gc[,1]),]
colnames(gc) <- c("chr","start","end","value1")
gc$start <- as.numeric(gc$start )/1e6
gc$end <- as.numeric(gc$end )/1e6
head(gc)
gc_10win=data.frame()
window=20
for(i in 1:nrow(ref)){
  gc_now <-gc[gc$chr == ref$chr[i],]
  sub=data.frame("start"=roll_mean(x=gc_now$start,n=window,by=window),
    "end"=roll_mean(x=gc_now$end,n=window,by=window),
    "value1"=roll_sum(x=gc_now$value1,n=window,by=window))
  sub$start=sub$start-(window*7/100)/2
  sub$end=sub$end+(window*7/100)/2
  sub=cbind(gc_now[1,1],sub)
  gc_10win=rbind(gc_10win,sub)
}
head(gc_10win)
color_assign <- colorRamp2(breaks = c(min(gc_10win[,4]), mean(gc_10win[,4]), max(gc_10win[,4])*0.8,max(gc_10win[,4])), 
                           col =c("#fdd3d3","#ff7861","red","#ff0000c5"))
circos.genomicTrackPlotRegion(
  gc_10win, track.height = 0.08, stack = TRUE, bg.border = NA,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)
  } )
# GC --------------------------------
library(stringr)
GC=data.frame(fread("ZY300/gc.txt"))
GC=cbind(str_split_fixed(GC[,1],":",2),GC[,2])
GC[,1]=gsub("_sliding","",GC[,1])
start=str_split_fixed(GC[,2],"-",2)[,1]
end=str_split_fixed(GC[,2],"-",2)[,2]
head(GC)
GC=cbind.data.frame(as.character(GC[,1]),as.numeric(start),as.numeric(end),as.numeric(GC[,3]))
GC <- GC[!grepl("scaffold",GC[,1]),]
GC <- GC[!grepl("Mt",GC[,1]),]
GC <- GC[!grepl("Pt",GC[,1]),]
colnames(GC) <- c("chr","start","end","value1")
gc=GC
gc$start <- as.numeric(gc$start )/1e6
gc$end <- as.numeric(gc$end )/1e6
head(gc)

gc_10win=data.frame()
window=20
for(i in 1:nrow(ref)){
  gc_now <-gc[gc$chr == ref$chr[i],]
  sub=data.frame("start"=roll_mean(x=gc_now$start,n=window,by=window),
    "end"=roll_mean(x=gc_now$end,n=window,by=window),
    "value1"=roll_mean(x=gc_now$value1,n=window,by=window))
  sub$start=sub$start-(window*7/100)/2
  sub$end=sub$end+(window*7/100)/2
  sub=cbind(gc_now[1,1],sub)
  gc_10win=rbind(gc_10win,sub)
}
head(gc_10win)
color_assign <- colorRamp2(breaks = c(min(gc_10win[,4]), mean(gc_10win[,4]), max(gc_10win[,4])*0.95,max(gc_10win[,4])), 
                           col =c("#d7f6db","#b2f3bb","#2fff00c4","#2fff00c5"))
circos.genomicTrackPlotRegion(
  gc_10win, track.height = 0.08, stack = TRUE, bg.border = NA,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)
  } )

# illumina mapped coverage ----------------------------------------------------------------
###S genome

circos.track(ylim=c(0, 20), panel.fun=function(x, y) {
  chr=gsub("Chr(.*)","\\1",CELL_META$sector.index)
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
}, bg.col="white", bg.border=F, track.height=0.06)

l.c <- 0.5
for (i in 1:nrow(ref)){
  Gs_now <-Gs[Gs$chr == ref$chr[i],]
  circos.barplot(Gs_now$value,Gs_now$start, col="grey70", lwd=0.6,
                 track.index = CELL_META$track.index,
                 sector.index= ref$chr[i],border="grey70")

}
### T gneome 
circos.track(ylim=c(0, 15), panel.fun=function(x, y) {
  chr=gsub("Chr(.*)","\\1",CELL_META$sector.index)
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
}, bg.col="white", bg.border=F, track.height=0.06)
l.c <- 0.5
for (i in 1:nrow(ref)){
  #i=1
  Gt_now <-Gt[Gt$chr == ref$chr[i],]
  circos.barplot(Gt_now$value,Gt_now$start, col="grey70", lwd=0.6,
                 track.index = CELL_META$track.index,
                 sector.index= ref$chr[i],border="grey70")
}


#colinear ----------------------------------------------------------------
#col.track <- brewer.pal(10,"Set3")
#rcols <- scales::alpha(sample(x = col.track,size = 38,replace=TRUE), alpha=0.7)
rcols=c("#BEBADAB2", "#BC80BDB2", "#D9D9D9B2", "#B3DE69B2", "#FFFFB3B2", "#BC80BDB2",
 "#BC80BDB2", "#80B1D3B2", "#FDB462B2", "#FB8072B2", "#BC80BDB2", "#80B1D3B2",
 "#FDB462B2", "#BC80BDB2", "#80B1D3B2", "#BEBADAB2", "#FDB462B2", "#D9D9D9B2",
 "#8DD3C7B2", "#B3DE69B2", "#FDB462B2", "#FB8072B2", "#FDB462B2", "#BEBADAB2",
 "#FFFFB3B2", "#D9D9D9B2", "#FDB462B2", "#BC80BDB2", "#FDB462B2", "#FB8072B2",
 "#8DD3C7B2", "#FDB462B2", "#D9D9D9B2", "#BEBADAB2", "#BEBADAB2", "#FDB462B2",
 "#80B1D3B2", "#FFFFB3B2")
circos.genomicLink(out2[1:38,1:3], out2[1:38,4:6], col=rcols, border="black",lwd=0.2)
dev.off()

##T genome ----------------------------------------------------------------
setwd("T")
circos.clear()
chrlen <- read.table("genome.fasta.fai",sep="\t")
chrlen <- chrlen[1:12,]
ref <- data.frame("chr" = chrlen$V1, "start" = 0,"end"=chrlen$V2)
rownames(ref)<- ref$chr

# panel ----------------------------------------------------------------
pdf(file = "T.pdf",width = 2,height = 2)
circos.clear()
circos.par(gap.degree= c(1,8,rep(1,10)),circle.margin=c(0.05,0.05,0.05,0.05))
circos.initialize(factors=ref$chr,xlim =cbind(0,ref$end/1e6) )
circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
  chr=gsub("chr(.*)","\\1",CELL_META$sector.index)
  chr[chr=="22"]="21"
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
  circos.text(mean(xlim), mean(ylim), chr, cex=0.4, col="black", 
              facing="bending.inside", niceFacing=TRUE)
}, bg.col="grey95", bg.border=F, track.height=0.1)

l.c <- 0.22
for (i in 1:nrow(ref)){
  mt <-seq(0,ref$end[i]/1e6,by=40)
  mt=c(mt,max(mt)+40)
  #pos=c(1,round(length(mt)/4,0),round(length(mt)/2,0),round(length(mt)*3/4,0),length(mt))
  circos.axis(h = "top",sector.index = ref$chr[i],track.index = 1,labels.cex=l.c,major.at = mt,labels.facing = "clockwise",direction="outside",minor.ticks = 1,labels=mt, major.tick.length = 0.4,lwd=0.2)
}
# Gene density ----------------------------------------------------------------
gd=data.frame(fread("T/gene_density.100k.txt"))
gd <- gd[!grepl("scaffold",gd[,1]),]
gd <- gd[!grepl("Mt",gd[,1]),]
gd <- gd[!grepl("Pt",gd[,1]),]
colnames(gd) <- c("chr","start","end","value1")
gd$start <- as.numeric(gd$start )/1e6
gd$end <- as.numeric(gd$end )/1e6
head(gd)

gd_10win=data.frame()
window=20
for(i in 1:nrow(ref)){
  gd_now <-gd[gd$chr == ref$chr[i],]
  sub=data.frame("start"=roll_mean(x=gd_now$start,n=window,by=window),
    "end"=roll_mean(x=gd_now$end,n=window,by=window),
    "value1"=roll_sum(x=gd_now$value1,n=window,by=window))
  sub$start=sub$start-(window*7/100)/2
  sub$end=sub$end+(window*7/100)/2
  sub=cbind(gd_now[1,1],sub)
  gd_10win=rbind(gd_10win,sub)
}
head(gd_10win)
color_assign <- colorRamp2(breaks = c(min(gd_10win[,4]), mean(gd_10win[,4]), max(gd_10win[,4])), 
                           col =c("#d4eaf3ec","#5ac2eb","blue"))
circos.genomicTrackPlotRegion(
  gd_10win, track.height = 0.2, stack = TRUE, bg.border = NA,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)
  } )

# TE density ----------------------------------------------------------------
gc=data.frame(fread("T/TE_density.100k.txt"))
gc <- gc[!grepl("scaffold",gc[,1]),]
colnames(gc) <- c("chr","start","end","value1")
gc$start <- as.numeric(gc$start )/1e6
gc$end <- as.numeric(gc$end )/1e6
head(gc)
gc_10win=data.frame()
window=20
for(i in 1:nrow(ref)){
  gc_now <-gc[gc$chr == ref$chr[i],]
  sub=data.frame("start"=roll_mean(x=gc_now$start,n=window,by=window),
    "end"=roll_mean(x=gc_now$end,n=window,by=window),
    "value1"=roll_sum(x=gc_now$value1,n=window,by=window))
  sub$start=sub$start-(window*7/100)/2
  sub$end=sub$end+(window*7/100)/2
  sub=cbind(gc_now[1,1],sub)
  gc_10win=rbind(gc_10win,sub)
}
head(gc_10win)
color_assign <- colorRamp2(breaks = c(min(gc_10win[,4]), mean(gc_10win[,4]), max(gc_10win[,4])*0.8,max(gc_10win[,4])), 
                           col =c("#fdd3d3","#ff7861","red","#ff0000c5"))
circos.genomicTrackPlotRegion(
  gc_10win, track.height = 0.2, stack = TRUE, bg.border = NA,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)
  } )
library(stringr)
GC=data.frame(fread("T/gc.txt"))
GC=cbind(str_split_fixed(GC[,1],":",2),GC[,2])
GC[,1]=gsub("_sliding","",GC[,1])
start=str_split_fixed(GC[,2],"-",2)[,1]
end=str_split_fixed(GC[,2],"-",2)[,2]
head(GC)
GC=cbind.data.frame(as.character(GC[,1]),as.numeric(start),as.numeric(end),as.numeric(GC[,3]))
GC <- GC[!grepl("scaffold",GC[,1]),]
GC <- GC[!grepl("Mt",GC[,1]),]
GC <- GC[!grepl("Pt",GC[,1]),]
colnames(GC) <- c("chr","start","end","value1")
gc=GC
gc$start <- as.numeric(gc$start )/1e6
gc$end <- as.numeric(gc$end )/1e6
head(gc)

gc_10win=data.frame()
window=20
for(i in 1:nrow(ref)){
  gc_now <-gc[gc$chr == ref$chr[i],]
  sub=data.frame("start"=roll_mean(x=gc_now$start,n=window,by=window),
    "end"=roll_mean(x=gc_now$end,n=window,by=window),
    "value1"=roll_mean(x=gc_now$value1,n=window,by=window))
  sub$start=sub$start-(window*7/100)/2
  sub$end=sub$end+(window*7/100)/2
  sub=cbind(gc_now[1,1],sub)
  gc_10win=rbind(gc_10win,sub)
}
head(gc_10win)
color_assign <- colorRamp2(breaks = c(min(gc_10win[,4]), mean(gc_10win[,4]), max(gc_10win[,4])*0.90,max(gc_10win[,4])), 
                           col =c("white","#91f1a0","#2fff00","#2fff00c5"))
circos.genomicTrackPlotRegion(
  gc_10win, track.height = 0.2, stack = TRUE, bg.border = NA,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)
  } )

dev.off()
## S genome ----------------------------------------------------------------
setwd("S")
circos.clear()
chrlen <- read.table("genome.fasta.fai",sep="\t")
chrlen <- chrlen[1:12,]
ref <- data.frame("chr" = chrlen$V1, "start" = 0,"end"=chrlen$V2)
rownames(ref)<- ref$chr

# panel ----------------------------------------------------------------
pdf(file = "S.pdf",width = 2,height = 2)
circos.clear()
circos.par(gap.degree= c(1,1,8,rep(1,9)),circle.margin=c(0.05,0.05,0.05,0.05))
circos.initialize(factors=ref$chr,xlim =cbind(0,ref$end/1e6) )
circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
  chr=gsub("chr(.*)","\\1",CELL_META$sector.index)
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
  circos.text(mean(xlim), mean(ylim), chr, cex=0.4, col="black", 
              facing="bending.inside", niceFacing=TRUE)
}, bg.col="grey95", bg.border=F, track.height=0.1)

l.c <- 0.22
for (i in 1:nrow(ref)){
  mt <-seq(0,ref$end[i]/1e6,by=40)
  mt=c(mt,max(mt)+40)
  #pos=c(1,round(length(mt)/4,0),round(length(mt)/2,0),round(length(mt)*3/4,0),length(mt))
  circos.axis(h = "top",sector.index = ref$chr[i],track.index = 1,labels.cex=l.c,major.at = mt,
    labels.facing = "clockwise",direction="outside",minor.ticks = 1,labels=mt, major.tick.length = 0.4,lwd=0.2)
}

# Gene density ----------------------------------------------------------------
gd=data.frame(fread("S/gene_density.100k.txt"))
gd <- gd[!grepl("scaffold",gd[,1]),]
gd <- gd[!grepl("Mt",gd[,1]),]
gd <- gd[!grepl("Pt",gd[,1]),]
colnames(gd) <- c("chr","start","end","value1")
gd$start <- as.numeric(gd$start )/1e6
gd$end <- as.numeric(gd$end )/1e6
head(gd)

gd_10win=data.frame()
window=20
for(i in 1:nrow(ref)){
  gd_now <-gd[gd$chr == ref$chr[i],]
  sub=data.frame("start"=roll_mean(x=gd_now$start,n=window,by=window),
    "end"=roll_mean(x=gd_now$end,n=window,by=window),
    "value1"=roll_sum(x=gd_now$value1,n=window,by=window))
  sub$start=sub$start-(window*7/100)/2
  sub$end=sub$end+(window*7/100)/2
  sub=cbind(gd_now[1,1],sub)
  gd_10win=rbind(gd_10win,sub)
}
head(gd_10win)
color_assign <- colorRamp2(breaks = c(min(gd_10win[,4]), mean(gd_10win[,4]), max(gd_10win[,4])), 
                           col =c("#d4eaf3ec","#5ac2eb","blue"))
circos.genomicTrackPlotRegion(
  gd_10win, track.height = 0.2, stack = TRUE, bg.border = NA,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)
  } )

# TE density ----------------------------------------------------------------
gc=data.frame(fread("S/TE_density.100k.txt"))
gc <- gc[!grepl("scaffold",gc[,1]),]
colnames(gc) <- c("chr","start","end","value1")
gc$start <- as.numeric(gc$start )/1e6
gc$end <- as.numeric(gc$end )/1e6
head(gc)
gc_10win=data.frame()
window=20
for(i in 1:nrow(ref)){
  gc_now <-gc[gc$chr == ref$chr[i],]
  sub=data.frame("start"=roll_mean(x=gc_now$start,n=window,by=window),
    "end"=roll_mean(x=gc_now$end,n=window,by=window),
    "value1"=roll_sum(x=gc_now$value1,n=window,by=window))
  sub$start=sub$start-(window*7/100)/2
  sub$end=sub$end+(window*7/100)/2
  sub=cbind(gc_now[1,1],sub)
  gc_10win=rbind(gc_10win,sub)
}
head(gc_10win)
color_assign <- colorRamp2(breaks = c(min(gc_10win[,4]), mean(gc_10win[,4]), max(gc_10win[,4])*0.8,max(gc_10win[,4])), 
                           col =c("#fdd3d3","#ff7861","red","#ff0000c5"))
circos.genomicTrackPlotRegion(
  gc_10win, track.height = 0.2, stack = TRUE, bg.border = NA,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)
  } )

# GC ------------------------------------
GC=data.frame(fread("S/gc.txt"))
GC=cbind(str_split_fixed(GC[,1],":",2),GC[,2])
GC[,1]=gsub("_sliding","",GC[,1])
start=str_split_fixed(GC[,2],"-",2)[,1]
end=str_split_fixed(GC[,2],"-",2)[,2]
head(GC)
GC=cbind.data.frame(as.character(GC[,1]),as.numeric(start),as.numeric(end),as.numeric(GC[,3]))
GC <- GC[!grepl("scaffold",GC[,1]),]
GC <- GC[!grepl("Mt",GC[,1]),]
GC <- GC[!grepl("Pt",GC[,1]),]
colnames(GC) <- c("chr","start","end","value1")
gc=GC
gc$start <- as.numeric(gc$start )/1e6
gc$end <- as.numeric(gc$end )/1e6
head(gc)

gc_10win=data.frame()
window=20
for(i in 1:nrow(ref)){
  gc_now <-gc[gc$chr == ref$chr[i],]
  sub=data.frame("start"=roll_mean(x=gc_now$start,n=window,by=window),
    "end"=roll_mean(x=gc_now$end,n=window,by=window),
    "value1"=roll_mean(x=gc_now$value1,n=window,by=window))
  sub$start=sub$start-(window*7/100)/2
  sub$end=sub$end+(window*7/100)/2
  sub=cbind(gc_now[1,1],sub)
  gc_10win=rbind(gc_10win,sub)
}
head(gc_10win)
color_assign <- colorRamp2(breaks = c(min(gc_10win[,4]), mean(gc_10win[,4]), max(gc_10win[,4])*0.90,max(gc_10win[,4])), 
                           col =c("white","#91f1a0","#2fff00","#2fff00c5"))
circos.genomicTrackPlotRegion(
  gc_10win, track.height = 0.2, stack = TRUE, bg.border = NA,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, col = color_assign(value[[1]]), border = NA, ...)
  } )
dev.off()



