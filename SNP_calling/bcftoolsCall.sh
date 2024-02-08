

### SNP calling 
bcftools mpileup --threads 15 -b all.bam.list.txt --fasta-ref genome.fasta | bcftools call -mv -o raw.vcf 
bcftools +counts raw.vcf&&
