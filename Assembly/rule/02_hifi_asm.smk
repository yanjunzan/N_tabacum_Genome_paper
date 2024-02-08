###################   hifi asm   ######################
rule hifi_asm:
	input:
		fa = "rawdata/02_ccs/{sample}.ccs.fasta.gz",
		R1 = "rawdata/03_hic/{sample}_1.fq.gz",
		R2 = "rawdata/03_hic/{sample}_2.fq.gz"
	output:
		pgfa = "results/02_assembly/02_hifi_asm/{sample}_hifiasm.hic.p_ctg.gfa"
	params:
		dir = "results/02_assembly/02_hifi_asm/{sample}_hifiasm"
	threads: 64
	shell:
		"""
		singularity exec GenomeAssemblyContainer_v0.2 \
			hifiasm -t {threads} -o {params.dir} -l {purge_level} --h1 {input.R1} --h2 {input.R2} {input.fa}
		"""

rule gfa2fa:
	input:
		primary_gfa = rules.hifi_asm.output
	output:
		primary_fa = "results/02_assembly/02_hifi_asm/{sample}_hifiasm.hic.p_ctg.fa"
	threads: 1
	shell:
		"""
		singularity exec GenomeAssemblyContainer_v0.2 \
			awk '/^S/{{print ">"$2;print $3}}' {input.primary_gfa} > {output.primary_fa}
		"""

rule seqkit_stats:
	input:
		primary_fa = rules.gfa2fa.output.primary_fa,
		hap1_fa = "results/02_assembly/02_hifi_asm/{sample}_hifiasm.hic.hap1.p_ctg.fa",
		hap2_fa = "results/02_assembly/02_hifi_asm/{sample}_hifiasm.hic.hap2.p_ctg.fa"
	output:
		primary_fa_stats = "results/03_assembly_assement/01_seqkit_stats/{sample}_hifiasm.hic.p_ctg",
		hap1_fa_stats = "results/03_assembly_assement/01_seqkit_stats/{sample}_hifiasm.hic.hap1.p_ctg",
		hap2_fa_stats = "results/03_assembly_assement/01_seqkit_stats/{sample}_hifiasm.hic.hap2.p_ctg"
	shell:
		"""
		singularity exec GenomeAssemblyContainer_v0.2 \
		seqkit stats -a -j 64 -t dna -o {output.primary_fa_stats} {input.primary_fa}
		singularity exec GenomeAssemblyContainer_v0.2 \
		seqkit stats -a -j 64 -t dna -o {output.hap1_fa_stats} {input.hap1_fa}
		singularity exec GenomeAssemblyContainer_v0.2 \
		seqkit stats -a -j 64 -t dna -o {output.hap2_fa_stats} {input.hap2_fa}
		"""

rule stats_n10_n90:
	input:
		primary_fa = rules.gfa2fa.output.primary_fa,
		hap1_fa = "results/02_assembly/02_hifi_asm/{sample}_hifiasm.hic.hap1.p_ctg.fa",
		hap2_fa = "results/02_assembly/02_hifi_asm/{sample}_hifiasm.hic.hap2.p_ctg.fa"
	output:
		primary_fa_stats = "results/03_assembly_assement/01_seqkit_stats/{sample}_hifiasm.hic.p_ctg_n10_n90",
		hap1_fa_stats = "results/03_assembly_assement/01_seqkit_stats/{sample}_hifiasm.hic.hap1.p_ctg_n10_n90",
		hap2_fa_stats = "results/03_assembly_assement/01_seqkit_stats/{sample}_hifiasm.hic.hap2.p_ctg_n10_n90"
	shell:
		"""
		singularity exec GenomeAssemblyContainer_v0.2 \
		/home/software/seq_n50.pl {input.primary_fa} > {output.primary_fa_stats}
		singularity exec GenomeAssemblyContainer_v0.2 \
		/home/software/seq_n50.pl {input.hap1_fa} > {output.hap1_fa_stats}
		singularity exec GenomeAssemblyContainer_v0.2 \
		/home/software/seq_n50.pl {input.hap2_fa} > {output.hap2_fa_stats}
		"""


rule busco:
	input:
		primary_fa = rules.gfa2fa.output.primary_fa,
		hap1_fa = "results/02_assembly/02_hifi_asm/{sample}_hifiasm.hic.hap1.p_ctg.fa",
		hap2_fa = "results/02_assembly/02_hifi_asm/{sample}_hifiasm.hic.hap2.p_ctg.fa"
	output:
		primary_short_summary = "results/03_assembly_assement/02_busco/{sample}_hifiasm.hic.p_ctg/short_summary.specific.solanales_odb10.{sample}_hifiasm.hic.p_ctg_busco.txt",
		hap1_short_summary = "results/03_assembly_assement/02_busco/{sample}_hifiasm.hic.hap1.p_ctg/short_summary.specific.solanales_odb10.{sample}_hifiasm.hic.hap1.p_ctg_busco.txt",
		hap2_short_summary = "results/03_assembly_assement/02_busco/{sample}_hifiasm.hic.hap2.p_ctg/short_summary.specific.solanales_odb10.{sample}_hifiasm.hic.hap2.p_ctg_busco.txt"
	shell:
		"""
		cd results/03_assembly_assement/02_busco

		singularity exec ../../../busco.sif busco -m genome \
					-i ../../../{input.primary_fa} \
					-o {Psample}_hifiasm.hic.p_ctg \
					-l ../../../../solanales_odb10 \
					-c 64 -f \
					--offline
		singularity exec ../../../busco.sif busco -m genome \
					-i ../../../{input.hap1_fa} \
					-o {Psample}_hifiasm.hic.hap1.p_ctg \
					-l ../../../../solanales_odb10 \
					-c 64 -f \
					--offline
		singularity exec ../../../busco.sif busco -m genome \
					-i ../../../{input.hap2_fa} \
					-o {Psample}_hifiasm.hic.hap2.p_ctg \
					-l ../../../../solanales_odb10 \
					-c 64 -f \
					--offline

		cd ../../../
		"""















