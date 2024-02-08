### illumina coverage 
#fq1 fq2 illumina sequencing data
#genome.fasta corresponding genome fasta file
fastp -i  fq1 -I  fq1 -o trimmed_fq1 -O trimmed_fq1 -Q --thread=15

bwa mem -aM -t 20 genome.fasta trimmed_fq1 trimmed_fq2 | samtools view -Sb -@ 20 -F 256 -F 4 > mapped.bam 
sambamba sort -t 10 mapped.bam -o sorted.bam
samtools view -q 30 -b -@ 10 -o sorted.Q30.bam sorted.bam
samtools depth -a sorted.Q30.bam >Q30.depth.txt 

### pacbio coverage
#pacbio.bam  pacbio sequencing data
bam2fastq -o fq pacbio.bam 
minimap2 -ax map-hifi -t 20 --split-prefix temp genome.fasta fq | samtools view -Sb -@ 20 -F 256 -F 4 > mapped.bam 
sambamba sort -t 10 mapped.bam -o sorted.bam
samtools view -q 30 -b -@ 10 -o sorted.Q30.bam sorted.bam
samtools depth -a sorted.Q30.bam >Q30.depth.txt 


### SNP calling 
#sorted.Q30.bam mapped bam file from illumina data
bcftools mpileup --threads 15 sorted.Q30.bam --fasta-ref genome.fasta | bcftools call -mv -o raw.vcf 
bcftools +counts raw.vcf&&
echo "heterozygous snp: "&&
grep -v "#" raw.vcf | grep -v "INDEL" | grep "0/1"| wc -l&&
echo "homozygous snp: "&&
grep -v "#" raw.vcf | grep -v "INDEL" | grep "1/1"| wc -l
grep -v "#" raw.vcf | grep  "INDEL"|awk '{print $10}'|cut -b 1-3|sort|uniq -c
