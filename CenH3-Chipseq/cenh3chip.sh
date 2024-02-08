### mapping
for name in {1,2,Input-1,Input-2}; do
    file1="$name"_R1.fq.gz"
    file2="$name"_R2.fq.gz"
    fastp -i $file1 -I $file2 \
        -o "yancao_"$name"_R1.fq.gz" -O "yancao_"$name"_R2.fq.gz" \
        -h "yancao_"$name".html" -w 8 -F 5 -f 5 &&
    bwa mem -t 8 ZY300.fa "yancao_"$name"_R1.fq.gz" "yancao_"$name"_R2.fq.gz" | samtools view -Sb -@ 8 > "yancao_"$name".bam" &&
    samtools view -@ 5 -q 20 -F 256 -F 4 -b "yancao_"$name".bam" > "yancao_"$name".filter.bam" &&
    sambamba markdup -r -t 5 "yancao_"$name".filter.bam" "yancao_"$name".filter.rm.bam" &&
	sambamba sort -t 5 "yancao_"$name".filter.rm.bam" -o "yancao_"$name".filter.rm.sorted.bam" &
done

### call depth
for name in {1,2};do
    bamCompare -b1 "mapped/yancao_"$name".filter.rm.sorted.bam" \
        -b2 "mapped/yancao_Input-"$name".filter.rm.sorted.bam" \
        -o "deep/"$name".5kbin.bedgraph"  --binSize 5000 --operation log2 --outFileFormat bedgraph -p 5 &
done

### qc 
multiBamSummary bins --bamfiles mapped/*sorted.bam -p max --smartLabels -o readCounts.npz
plotCorrelation \
    -in readCounts.npz \
    --corMethod spearman --skipZeros \
    --plotTitle "Spearman Correlation of Read Counts" \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o heatmap_SpearmanCorr_readCounts.png   \
    --outFileCorMatrix SpearmanCorr_readCounts.tab &
plotPCA -in readCounts.npz \
    -o PCA_readCounts.png \
    -T "PCA of read counts" &
### coverage 
for name in {1,2,Input-1,Input-2}; do
    bamCoverage --bam "mapped/yancao_"$name".filter.rm.sorted.bam" -o "deep/"$name".5kbin.coverage.bedgraph" --outFileFormat bedgraph -p 10 --binSize 5000 --normalizeUsing RPKM &
done


