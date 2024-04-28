# Author: Kristina Gagalova
# Submodule: pre-assembly
# to be run with the main-snakefile

from snakemake.utils import min_version
min_version("5.24.0")
from pathlib import Path
import os, errno

#-----------------------------------------
#parse index file
configfile: "../assembly-vars.yaml"

#parse out dir name
OUT_DIR = config["pre_dir"]

import pandas as pd

sample_file = config["index_file"]
sample_df = pd.read_table(sample_file, sep="\t", header=0)
sampleID = sample_df.nam

rule all:
	input:
		".make_preassembly.touch"

rule concat:
	input:
		dir1 = lambda w: sample_df[sample_df.nam == w.sample].lr1.tolist(),
		dir2 = lambda w: sample_df[sample_df.nam == w.sample].lr2.tolist() if "lr2" in list(sample_df.columns) else [],
		dir3 = lambda w: sample_df[sample_df.nam == w.sample].lr3.tolist() if "lr3" in list(sample_df.columns) else []
	output:
		"01-long_reads/{sample}.fastq.gz"
	version: "1.0"
	run:
		if "lr3" in list(sample_df.columns):
			shell("zcat {input.dir1}/*fastq.gz {input.dir2}/*fastq.gz {input.dir3}/*fastq.gz | gzip > {output}")
		elif "lr2" in list(sample_df.columns):
			shell("zcat {input.dir1}/*fastq.gz {input.dir2}/*fastq.gz | gzip > {output}")
		else:
			shell("zcat {input.dir1}/*fastq.gz | gzip > {output}")


rule porechop:
	input:
		rules.concat.output
	output:
		"01-long_reads/{sample}_pc.fastq.gz"
	log:	
		"log/{sample}/02-porechop/{sample}_porechop.log"
	threads:
		24
	version:
		"1.0"
	shell:
		"porechop -i {input} -o {output} --format fastq.gz -v 2 -t {threads} &> {log}"
	


rule longQC_post:
	input:
		rules.porechop.output
	output:
		directory("02-longQC/{sample}")
	log:
		"log/{sample}/02-longQC/{sample}_longQC.log"
	threads:
		24
	version:
		"1.0"
	shell:
		"python /home/kgagalova/src/LongQC/longQC.py sampleqc -p {threads} -x ont-rapid -o {output} {input} 2> {log}"



rule qc_short:
	input:
		fastq1 = "../" + config['sr_dir'] + "/" + "{sample}_R1.fastq.gz",
		fastq2 = "../" + config['sr_dir'] + "/" + "{sample}_R2.fastq.gz",
	output:
		"02-fastpQC/{sample}.html"
	log:
		"log/{sample}/02-shortQC/{sample}_fastpQC.log"
	threads:
		24
	shell:
		"fastp -i {input.fastq1} -I {input.fastq2} -w {threads} -p -h {output} 2> {log}"


rule contamination_short:
	input:
		fastq1 = "../" + config['sr_dir'] + "/" + "{sample}_R1.fastq.gz",
		fastq2 = "../" + config['sr_dir'] + "/" + "{sample}_R2.fastq.gz",
	output:
		report = "03-contamination/kraken2/{sample}.report.txt",
		kraken_out = "03-contamination/kraken2/{sample}.kraken"
	log:
		"log/{sample}/03-contamination/{sample}_kraken.log"
	params:
		con_db = config['db_kraken']
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


rule plot_krona:
	input:
		rules.contamination_short.output.kraken_out
	output:
		"03-contamination/kraken2/{sample}.kraken.krona.html"
	log:
		"log/{sample}/03-contamination/kraken2/{sample}_krona.log"
	shell:
		"""
		cat {input} | cut -f 2,3 > {wildcards.sample}.kraken.krona && 
		ktImportTaxonomy {wildcards.sample}.kraken.krona -o {output} 2> {log} && 
		rm {wildcards.sample}.kraken.krona
		"""


rule time_stamp:
	input:
		expand("01-long_reads/{sample}.fastq.gz",sample=sampleID),
		expand("02-longQC/{sample}",sample=sampleID),
		expand("02-fastpQC/{sample}.html",sample=sampleID),
		expand("03-contamination/kraken2/{sample}.kraken",sample=sampleID),
		expand("03-contamination/kraken2/{sample}.kraken.krona.html",sample=sampleID),
	output:
		touch(".make_preassembly.touch")
