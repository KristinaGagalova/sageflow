# Author: Kristina Gagalova
# decomntamination for short reads

from snakemake.utils import min_version
min_version("5.24.0")
from pathlib import Path
import os, errno

#-----------------------------------------
#parse index file

configfile: "decontamination-vars.yaml"
import pandas as pd

sample_file = config["index_file"]
sample_df = pd.read_table(sample_file, sep="\t", header=0)
sampleID = sample_df.nam_cont.tolist()
print(sampleID)

rule all:
	input:
		expand("cont_lr/{sample}/{sample}-bac-reads.txt",sample=sampleID),
		expand("cont_lr/{sample}/{sample}-fun-reads.txt",sample=sampleID),
		expand("cont_lr/{sample}/{sample}_clean.fastq.gz",sample=sampleID),
		expand("cont_sr/{sample}/{sample}.kraken",sample=sampleID),
		expand("cont_sr/{sample}/{sample}_clean_R1.fastq.gz",sample=sampleID),
		expand("cont_sr/{sample}/{sample}_cont_R1.fastq.gz",sample=sampleID),
		expand("cont_lr/{sample}/{sample}_cont.fastq.gz",sample=sampleID)
		


rule bacterial_lr:
	input:
		reads_lr = lambda w: sample_df[sample_df.nam_cont == w.sample].lr1.tolist(),
	output:
		report = "cont_lr/{sample}/{sample}-bac-report.out",
		reads ="cont_lr/{sample}/{sample}-bac-reads.txt"
	log:	
		"log/{sample}/01-bact_lr/{sample}_centrifuge.log"
	threads:
		24
	params:
		bact_lr = config['bacteria_lr']		
	shell:
		"""
		mkdir -p cont_lr/{wildcards.sample} &&
		centrifuge -x {params.bact_lr} \
			-U {input.reads_lr} \
			--report-file {output.report} \
			-S {output.reads} \
			-p {threads} 2> {log}
		"""

rule fun_lr:
	input:
		reads_lr = lambda w: sample_df[sample_df.nam_cont == w.sample].lr1.tolist(),
	output:
		report = "cont_lr/{sample}/{sample}-fun-report.out",
		reads ="cont_lr/{sample}/{sample}-fun-reads.txt"
	log:
		"log/{sample}/01-fungi_lr/{sample}_centrifuge.log"
	threads: 24
	params:
		fun_lr = config['fungi_lr']
	shell:
		"""
		mkdir -p cont_lr/{wildcards.sample} &&
		centrifuge -x {params.fun_lr} \
			-U {input.reads_lr} \
			--report-file {output.report} \
			-S {output.reads} \
			-p {threads} 2> {log}
		"""

rule remove_reads_lr:
	input:
		reads_lr = lambda w: sample_df[sample_df.nam_cont == w.sample].lr1.tolist(),
		cont_bac = rules.bacterial_lr.output.reads,
		cont_fun = rules.fun_lr.output.reads
	output:
		reads_clean ="cont_lr/{sample}/{sample}_clean.fastq.gz",
		reads_cont = "cont_lr/{sample}/{sample}_cont.fastq.gz"
	shell:
		"""
		join -a 1 -a 2 -1 1 -2 1 <(sort -k1,1 {input.cont_fun}) <(sort -k1,1 {input.cont_bac}) | tr ' ' '\t' |\
			 awk '$4 < $11 {{print $0}}' > cont_lr/{wildcards.sample}/{wildcards.sample}_remove_reads.in &&
		grep -vFf <(awk '{{print $1}}' cont_lr/{wildcards.sample}/{wildcards.sample}_remove_reads.in) <(zcat {input.reads_lr} | paste - - - - ) | tr '\t' '\n' | gzip > {output.reads_clean} &&
		grep -Ff <(awk '{{print $1}}' cont_lr/{wildcards.sample}/{wildcards.sample}_remove_reads.in) <(zcat {input.reads_lr} | paste - - - - ) | tr '\t' '\n' | gzip > {output.reads_cont}
		"""

rule bacterial_sr:
	input:
		fastq1 = lambda w: sample_df[sample_df.nam_cont == w.sample].sr1.tolist(),
		fastq2 = lambda w: sample_df[sample_df.nam_cont == w.sample].sr2.tolist(),
	output:
		report = "cont_sr/{sample}/{sample}.report.txt",
		kraken_out = "cont_sr/{sample}/{sample}.kraken"
	log:
		"log/{sample}/02-bact_sr/{sample}_kraken.log"
	params:
		con_db = config['bacterial_sr']
	threads:
		48
	shell:
		"""
		kraken2 -db {params.con_db} \
		--report {output.report} \
		--gzip-compressed \
		--threads {threads} \
		--paired {input.fastq1} {input.fastq2} > {output.kraken_out} 2>{log}
		"""


rule remove_reads_sr:
	input:
		fastq1 = lambda w: sample_df[sample_df.nam_cont == w.sample].sr1.tolist(),
		fastq2 = lambda w: sample_df[sample_df.nam_cont == w.sample].sr2.tolist(),
		cont = rules.bacterial_sr.output.kraken_out
	output:
		reads_clean_r1 = "cont_sr/{sample}/{sample}_clean_R1.fastq.gz",
		reads_clean_r2 = "cont_sr/{sample}/{sample}_clean_R2.fastq.gz",
		reads_cont_r1 = "cont_sr/{sample}/{sample}_cont_R1.fastq.gz",
		reads_cont_r2 = "cont_sr/{sample}/{sample}_cont_R2.fastq.gz"
	shell:
		"""
		awk '$1 == "C" {{print $2}}' {input.cont} > cont_sr/{wildcards.sample}/{wildcards.sample}_remove_reads.in &&
		grep -vFf cont_sr/{wildcards.sample}/{wildcards.sample}_remove_reads.in <(zcat {input.fastq1} | paste - - - - ) | tr '\t' '\n' | gzip > {output.reads_clean_r1} &&
		grep -vFf cont_sr/{wildcards.sample}/{wildcards.sample}_remove_reads.in <(zcat {input.fastq2} | paste - - - - ) | tr '\t' '\n' | gzip > {output.reads_clean_r2} &&
		grep -Ff cont_sr/{wildcards.sample}/{wildcards.sample}_remove_reads.in <(zcat {input.fastq1} | paste - - - - ) | tr '\t' '\n' | gzip > {output.reads_cont_r1} &&
                grep -Ff cont_sr/{wildcards.sample}/{wildcards.sample}_remove_reads.in <(zcat {input.fastq2} | paste - - - - ) | tr '\t' '\n' | gzip > {output.reads_cont_r2}
		"""
