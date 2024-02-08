### mapping. Run time 2d
for file in BS* ; do
    fq1=$file/*_1.clean.fq.gz
    fq2=$file/*_2.clean.fq.gz
    if [[ "${fq1:62:5}" != "zy300" ]]; then
        log=${fq1:62:3}.log.txt
    else
        log=${fq1:62:7}.log.txt
    fi  
    nohup bismark --parallel 3 --bowtie2 --non_directional -o "final_aligned/"$out -N 1 --genome $ref -1 $fq1 -2 $fq2 > $log &
done


### remove duplication
for file in final_aligned/*.bam ; do
    echo $file
    if [[ "${file:14:5}" != "zy300" ]]; then
        log=${file:14:3}.log.txt
    else
        log=${file:14:7}.log.txt
    fi
    echo $log
    nohup deduplicate_bismark -p $file --output_dir deduplicated > $log &
done  


#extraction
cd extracted
for file in deduplicated/*.bam ; do
    if [[ "${file:16:5}" != "zy300" ]]; then
        log=${file:16:3}/${file:16:3}.log.txt
        #mkdir ${file:16:3}
        nohup bismark_methylation_extractor --bedGraph --parallel 10 --CX_context --cytosine_report --no_overlap \
        --genome_folder bismark/\
        $file -o extracted/${file:16:3} > $log &
    else
        log=${file:16:7}/${file:16:7}.log.txt
        #mkdir ${file:16:7}
        nohup bismark_methylation_extractor --bedGraph --parallel 10 --CX_context --cytosine_report --no_overlap \
        --genome_folder /data/home/yanjun/project/tobacco/data/genome/bismark/\
        $file -o /data/public/hanyu/Q7/extracted/${file:16:7} > $log &
    fi
done  

### convert format for methylkit analysis
for path in $(ls extracted/); do
    file=$(ls $path/*deduplicated.CX_report.txt)
    perl bismarkCXmethykit.pl $file &
done

