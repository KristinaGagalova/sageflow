# Author: Kristina Gagalova
# Submodule: assembly
# to be run with the main-snakefile

from snakemake.utils import min_version
min_version("5.24.0")
from pathlib import Path

#get path current dir
MY_PATH = os.getcwd()

#config file
configfile: "../assembly-vars.yaml"

SAMPLES, = glob_wildcards("../" + config['long_reads'] + "/" + "{sample}_pc.fastq.gz")
print(SAMPLES)

rule all:
	input:
		".make_assembly.touch"

rule assemble_lr:
	input:
		fasta = "../" + config["long_reads"] + "/" + "{sample}_pc.fastq.gz",
	output:
		"01-assemblies/{sample}/assembly.fasta"
	log:
		"logs/{sample}/01-assemblies/{sample}.log"
	threads: 48
	params:
		out_dir = "01-assemblies/{sample}",
		gen_size = '11.6m'
	shell: 
		"/usr/bin/time -v flye --threads {threads} \
		--asm-coverage 100 \
		-o {params.out_dir}\
		--genome-size {params.gen_size}\
		--nano-hq {input} 2> {log}"


rule stats_lr_assembly:
	input:
		rules.assemble_lr.output
	output:
		"stats/{sample}/01-assembly/{sample}_01.stats"
	shell:
		"abyss-fac {input} > {output}"


rule prepare_ls:
	input:
		assembly = rules.assemble_lr.output,
		fastq = "../" + config["long_reads"] + "/" + "{sample}.fastq.gz"
	output:
		assembly_out = "{sample}.fa",
		fastq_out =  "{sample}.fq.gz"
	shell:
		"ln -s {input.assembly} {output.assembly_out} && ln -s {input.fastq} {output.fastq_out}"


rule assembly_corr:
	input:
		assembly = rules.prepare_ls.output.assembly_out,
		fastq = rules.prepare_ls.output.fastq_out
	output:
		"02-assembly_corr/{sample}/{sample}.scaffolds.fa"
	log:
		"logs/{sample}/02-assembly_corr/{sample}.log"
	params:
		cur_dir = MY_PATH
	shell:
		"""
		mkdir -p 02-assembly_corr/{wildcards.sample} && \
		/home/kgagalova/src/LongStitch/longstitch make_links -C {params.cur_dir}/02-assembly_corr/{wildcards.sample} reads_path={params.cur_dir}/{input.fastq} draft_path={params.cur_dir}/{input.assembly} && \
		/home/kgagalova/src/LongStitch/longstitch run -C {params.cur_dir}/02-assembly_corr/{wildcards.sample} out_prefix={wildcards.sample} draft={wildcards.sample} reads={wildcards.sample} G=12e6 1>{log} 2>&1 && \
		rm ./{wildcards.sample}*
		"""

rule stats_lr_ls:
        input:
                rules.assembly_corr.output
        output:
                "stats/{sample}/02-assembly_corr/{sample}_02.stats"
        shell:
                "abyss-fac {input} > {output}"


rule align_lr:
	input:
		fastq = "../" + config["long_reads"] + "/" + "{sample}.fastq.gz",
		assembly = rules.assembly_corr.output
	output:
		temp("03-polish_lr/{sample}.sam")
	threads:
		24
	log:
		"logs/{sample}/03-polish_lr/{sample}_minimap.log"
	shell:
		"minimap2 -t {threads} -ax map-ont {input.assembly} {input.fastq} -o {output} 1>{log} 2>&1"  



rule polish_lr:
	input:
		fastq = "../" + config["long_reads"] + "/" + "{sample}.fastq.gz",
		assembly = rules.assembly_corr.output,
		sam = rules.align_lr.output
	output:
		"03-polish_lr/{sample}_racon.fasta"
	log:
		"logs/{sample}/03-polish_lr/{sample}_racon.log"
	threads:
		24
	shell:
		"racon -t {threads} -u {input.fastq} {input.sam} {input.assembly} | cut -d' ' -f1 | cut -d':' -f1 1> {output} 2> {log}"

#medaka requirements
#Program    Version    Required  Local_installation
#bcftools   1.11-2-g67c95ee  1.11	
#bgzip      1.9        1.11		*
#minimap2   2.22       2.11
#samtools   1.9        1.11	* 
#tabix      1.9        1.11	*

rule run_consensus:
	input:
		rc = rules.polish_lr.output,
		fastq = "../" + config["long_reads"] + "/" + "{sample}.fastq.gz",
	output:
		"04-consensus/{sample}/consensus.fasta"
	log:
		"logs/{sample}/04-consensus/{sample}_medaka.log"
	threads:
		24
	params:
		out = "04-consensus/{sample}",
		source_path = "../" + config["source_medaka"]
	shell:
		"""
		source {params.source_path}
		medaka_consensus -t {threads} -i {input.fastq} -d {input.rc} -o {params.out} -m r941_min_high_g303 2> {log}
		"""

rule bwa_index:
	input:
		rules.run_consensus.output
	output:
		"04-consensus/{sample}/consensus.fasta.amb",
		"04-consensus/{sample}/consensus.fasta.ann",
	log:
		"logs/{sample}/04-consensus/{sample}_bwa_index.log"
	wrapper:
		"0.79.0/bio/bwa/index"

rule align_sr:
	input:
		fasta = rules.run_consensus.output,
		fq1 = "../" + config["sr_dir"] + "/{sample}_R1.fastq.gz",
		fq2 = "../" + config["sr_dir"] + "/{sample}_R2.fastq.gz",
		index = rules.bwa_index.output
	output:
		"05-polish_sr/{sample}/{sample}.posSorted.bam"
	log:
		"logs/{sample}/05-polish_sr/{sample}_align.log"
	threads:
		24
	shell:
		"bwa mem -t {threads} {input.fasta} {input.fq1} {input.fq2} | samtools view -Sb - | samtools sort -o {output} 2>{log}"


rule index_bams:
	input:
		rules.align_sr.output
	output:
		"05-polish_sr/{sample}/{sample}.posSorted.bam.bai"
	wrapper:
		"0.66.0/bio/samtools/index"


rule polish_sr:
	input:
		fasta = rules.run_consensus.output,
		bam = rules.align_sr.output,
		bai = rules.index_bams.output,
	output:
		"05-polish_sr/{sample}/assembly.fasta"
	log:
		"logs/{sample}/05-polish_sr/{sample}_pilon.log"
	threads: 
		48
	params:
		out_dir = "05-polish_sr/{sample}"
	shell:
		"pilon --genome {input.fasta} \
		--frags {input.bam} \
		--output assembly \
		--outdir {params.out_dir} \
		--fix all \
		--threads {threads} \
		--verbose \
		--vcf \
		--changes \
		--vcfqe \
		--tracks &> {log}"

rule stats_polish_sr:
	input:
		rules.polish_sr.output,
	output:
		"stats/{sample}/03-assembly_polish_sr/{sample}_03.stats"
	shell:
		"abyss-fac {input} > {output}"


rule list_rename:
	input:
		rules.polish_sr.output
	output:
		"06-final_out/names/{sample}.tab"
	shell:
		"""
		seqkit sort -lr {input} | grep ">" | sed 's/>//' | awk 'BEGIN {{FS=OFS="\t"}} {{print  $1,NR}}' | awk 'BEGIN {{FS=OFS="\t"}} {{$2 = sprintf("%02d", $2)}}1' | 
		awk -v s={wildcards.sample} 'BEGIN {{FS=OFS="\t"}} {{print $1,s"_scaf"$2}}' > {output}
		"""

rule rename_fasta:
	input:
		fasta = rules.polish_sr.output,
		tab = rules.list_rename.output
	output:
		temp("06-final_out/{sample}_rename.fasta")
	log:
		"logs/{sample}/06-final_out/{sample}_rename.log"
	shell:
		"seqkit replace -p '(.+)$' -r '{{kv}}' -k {input.tab} {input.fasta} 1> {output} 2> {log}"


rule get_mt:
	input:
		rules.rename_fasta.output
	output:
		"06-final_out/names/{sample}_MT.tab"
	log:
		"logs/{sample}/06-final_out/{sample}_renameMinimap.log"
	params:
		mt = config["mt_yeast"]
	shell:
		"""
		minimap2 -x asm5 {input} {params.mt} |\
		awk '{{print $6,$6,"MT"}}' | sed 's/ /\t/' | sort | uniq 1> {output} 2> {log}
		"""

rule rename_mt:
	input:
		mt = rules.get_mt.output,
		fasta = rules.rename_fasta.output
	output:
		"06-final_out/{sample}.fasta"
	log:
		"logs/{sample}/06-final_out/{sample}_renameMT.log"
	shell:
		"seqkit replace --keep-key -p '(.+)$' -r '{{kv}}' -k {input.mt} {input.fasta} 1> {output} 2> {log} && \
		samtools faidx {output}"



rule contamination_assembly:
	input:
		rules.rename_mt.output
	output:
		"contamination_post/{sample}/.{sample}_cont.out"
	params:
		out_dir = "contamination_post/{sample}",
		tax_db = config['tax_db']
	log:
		log_db = "logs/{sample}/contamination/{sample}_db.log",
		log_tax = "logs/{sample}/contamination/{sample}_cont.log"
	shell:
		"""
		mkdir -p {params.out_dir}/query && 
		mmseqs createdb {input} {params.out_dir}/query/{wildcards.sample} > {log.log_db} && 
		mmseqs taxonomy {params.out_dir}/query/{wildcards.sample} {params.tax_db} {params.out_dir}/{wildcards.sample}_taxonomyResult.out tmp_{wildcards.sample} > {log.log_tax} &&
		touch {output} && rm -r tmp_{wildcards.sample}
		"""

rule contamination_report:
	input:
		rules.contamination_assembly.output
	output:
		"contamination_post/{sample}/{sample}_taxonomyResult.out"
	params: 
		out_dir = "contamination_post/{sample}",
		tax_db = config['tax_db']
	log:
		"logs/{sample}/contamination/{sample}_report.log"
	shell:
		"""
		mmseqs createtsv {params.out_dir}/query/{wildcards.sample} {params.out_dir}/{wildcards.sample}_taxonomyResult.out {params.out_dir}/{wildcards.sample}_taxonomyResult.out > {log}
		"""



rule time_stamp:
	input:
		expand("stats/{sample}/01-assembly/{sample}_01.stats",sample=SAMPLES),
		expand("stats/{sample}/02-assembly_corr/{sample}_02.stats",sample=SAMPLES),
		expand("stats/{sample}/03-assembly_polish_sr/{sample}_03.stats",sample=SAMPLES),
		expand("05-polish_sr/{sample}/assembly.fasta",sample=SAMPLES),
		expand("06-final_out/{sample}.fasta",sample=SAMPLES),
		expand("contamination_post/{sample}/.{sample}_cont.out",sample=SAMPLES),
		expand("contamination_post/{sample}/{sample}_taxonomyResult.out",sample=SAMPLES)
	output:
		touch(".make_assembly.touch")
