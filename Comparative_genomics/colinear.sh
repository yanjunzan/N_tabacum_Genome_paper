### mapping
minimap2 -x asm5 -t 15 zy300.genome.fasta T.genome.fasta -o TtoZY300.paf & 
minimap2 -x asm5 -t 15 zy300.genome.fasta S.genome.fasta -o StoZY300.paf & 
awk -F'\t' '{if($12>=30) print $0}'  StoZY300.paf > StoZY300_q30.paf &
awk -F'\t' '{if($12>=30) print $0}'  TtoZy300.paf > TtoZY300_q30.paf &
###plot
perl  GetTwoGenomeSyn.pl  Paf2Link   StoZY300_q30.paf    1000000   StoZY300_NGsyn.link 
perl  GetTwoGenomeSyn.pl  Paf2Link   TtoZY300_q30.paf    1000000   TtoZY300_NGsyn.link 
awk '{print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3}' StoZY300_NGsyn.link > ZY300toS_NGsyn.link
awk '{print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3}' TtoZY300_NGsyn.link > ZY300toT_NGsyn.link
/data/public/hanyu/Q1/NGenomeSyn-1.41/bin/NGenomeSyn  -InConf   my_config.conf   -OutPut    three_genome

### dotplot run in R
for( len in c("1000","10000","100000","1000000")){
    print(len)
    data_S=data.frame(fread(paste0("StoZY300_NGsyn_",len,".link")))
    data_T=data.frame(fread(paste0("TtoZY300_NGsyn_",len,".link")))
    zy_len=data.frame(fread("ZY300.chr.len"))
    S_len=data.frame(fread("S.chr.len"))
    T_len=data.frame(fread("T.chr.len"))
    ZY_genome_position=c(0)
    ZY_mark=c(zy_len[1,3])
    for(i in 2:24){
        ZY_genome_position=c(ZY_genome_position,sum(zy_len[1:(i-1),3]))
        ZY_mark=c(ZY_mark,sum(zy_len[1:i,3]))
    }
    S_genome_position=c(0)
    S_mark=c(S_len[1,3])
    for(i in 2:12){
        S_genome_position=c(S_genome_position,sum(S_len[1:(i-1),3]))
        S_mark=c(S_mark,sum(S_len[1:i,3]))
    }
    T_genome_position=c(0)
    T_mark=c(T_len[1,3])
    for(i in 2:12){
        T_genome_position=c(T_genome_position,sum(T_len[1:(i-1),3]))
        T_mark=c(T_mark,sum(T_len[1:i,3]))
    }
    names(ZY_genome_position)=zy_len[,1]
    names(S_genome_position)=S_len[,1]
    names(T_genome_position)=T_len[,1]
    new_S=data_S
    new_T=data_T

    new_S[,2]=new_S[,2]+S_genome_position[new_S[,1]]
    new_S[,3]=new_S[,3]+S_genome_position[new_S[,1]]
    new_S[,5]=new_S[,5]+ZY_genome_position[new_S[,4]]
    new_S[,6]=new_S[,6]+ZY_genome_position[new_S[,4]]
    new_T[,2]=new_T[,2]+T_genome_position[new_T[,1]]
    new_T[,3]=new_T[,3]+T_genome_position[new_T[,1]]
    new_T[,5]=new_T[,5]+ZY_genome_position[new_T[,4]]
    new_T[,6]=new_T[,6]+ZY_genome_position[new_T[,4]]
    pdf(paste0("All_colinear_",len,".pdf"),width=11,height=15)
    plot(x=c(0,sum(zy_len[,3])),y=c(0,sum(c(S_len[,3],T_len[,3]))),col="white",frame.plot=F,
        xaxt="n",yaxt="n",xlab="ZY300",ylab="ST",main="")
    mid_zy=(ZY_genome_position+ZY_mark)/2
    mid_s=(S_genome_position+S_mark)/2
    mid_t=(T_genome_position+T_mark)/2
    axis(1,at=mid_zy,labels = paste0("Chr",1:24))
    axis(2,at=c(mid_s,(mid_t)+max(S_mark)),labels = paste0("Chr",c(c(1,3,5,6,7,8,10,11,16,18,20,22),c(2,4,9,12,13,14,15,17,19,21,23,24))))
    abline(v=c(0,ZY_mark),lty="dashed",col="grey")
    abline(h=c(0,S_mark,T_mark+max(S_mark)),lty="dashed",col="grey")
    abline(h=max(S_mark),col="black")
    for(i in 1:nrow(new_S)){
        lines(x=new_S[i,c(5,6)],y=new_S[i,c(2,3)],col="tomato",lwd=4)
    }
    for(i in 1:nrow(new_T)){
        lines(x=new_T[i,c(5,6)],y=sum(S_len[,3])+new_T[i,c(2,3)],col="skyblue",lwd=4)
    }
    dev.off()
}