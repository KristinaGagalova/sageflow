# Author: Kristina Gagalova
# Submodule: mapping reads
# to be run with the main-snakefile

from snakemake.utils import min_version
min_version("5.24.0")
from pathlib import Path

#config file
configfile: "../assembly-vars.yaml"

SAMPLES, = glob_wildcards("../" + config['long_reads'] + "/" + "{sample}_pc.fastq.gz")
print(SAMPLES)

rule all:
	input:
		".make_mapping.touch"


rule map_long:
	input:
		assembl = "../" + config["assemblies"] + "/" + "{sample}.fasta",
		fastq = "../" + config["long_reads"] + "/" + "{sample}_pc.fastq.gz"
	output:
		"mapping/self/long/{sample}.posSorted_lrSelf.bam"
	log:
		"logs/{sample}/mapping/{sample}_Selflong.log"
	threads:
		24
	shell:
		"minimap2 -t {threads} -ax map-ont {input.assembl} {input.fastq} -o {wildcards.sample}.sam 1>{log} 2>&1 && \
		samtools view -Sb {wildcards.sample}.sam | samtools sort -o {output} && \
		rm {wildcards.sample}.sam"

rule bwa_index:
	input:
		"../" + config["assemblies"] + "/" + "{sample}.fasta"
	output:
		"../" + config["assemblies"] + "/" + "{sample}.fasta.amb",
		"../" + config["assemblies"] + "/" + "{sample}.fasta.ann"
	log:
		"logs/{sample}/mapping/{sample}_bwa_index.log"
	wrapper:
		"0.79.0/bio/bwa/index"

rule align_sr:
	input:
		assembly = "../" + config["assemblies"] + "/" + "{sample}.fasta",
		fq1 = "../" + config["sr_dir"] + "/{sample}_R1.fastq.gz",
		fq2 = "../" + config["sr_dir"] + "/{sample}_R2.fastq.gz",
		index = rules.bwa_index.output
	output:
		"mapping/self/short/{sample}.posSorted_srSelf.bam"
	log:
		"logs/{sample}/mapping/{sample}_Selfshort.log"
	threads: 24
	shell:
		"bwa mem -t {threads} {input.assembly} {input.fq1} {input.fq2} | samtools view -Sb - | samtools sort -o {output} 2> {log}"

rule index_short:
	input:
		rules.align_sr.output
	output:
		"mapping/self/short/{sample}.posSorted_srSelf.bam.bai"
	wrapper:
		"0.66.0/bio/samtools/index"

rule index_long:
	input:
		rules.map_long.output
	output:
		"mapping/self/long/{sample}.posSorted_lrSelf.bam.bai"
	wrapper:
		"0.66.0/bio/samtools/index"	
	


rule map_longRef:
	input:
		ref = config["ref_cenpk"],
		fastq = "../" + config["long_reads"] + "/" + "{sample}_pc.fastq.gz"
	output:
		"mapping/reference/long/{sample}.posSorted_lrRef.bam"
	log:
		"logs/{sample}/mapping/{sample}_Reflong.log"
	threads:
		24
	shell:
		"minimap2 -t {threads} -ax map-ont {input.ref} {input.fastq} -o {wildcards.sample}ref.sam 1>{log} 2>&1 && \
		samtools view -Sb {wildcards.sample}ref.sam | samtools sort -o {output} && \
		rm {wildcards.sample}ref.sam"


rule index_longRef:
	input:
		rules.map_longRef.output
	output:
		"mapping/reference/long/{sample}.posSorted_lrRef.bam.bai"
	wrapper:
		"0.66.0/bio/samtools/index"



rule align_srRef:
	input:
		assembly = config["ref_cenpk"],
		fq1 = "../" + config["sr_dir"] + "/{sample}_R1.fastq.gz",
		fq2 = "../" + config["sr_dir"] + "/{sample}_R2.fastq.gz"
	output:
		"mapping/reference/short/{sample}.posSorted_srRef.bam"
	log:
		"logs/{sample}/mapping/{sample}_Refshort.log"
	threads:
		24
	shell:
		"bwa mem -t {threads} {input.assembly} {input.fq1} {input.fq2} | samtools view -Sb - | samtools sort -o {output} 2> {log}"
	

rule index_shortRef:
	input:
		rules.align_srRef.output
	output:
		"mapping/reference/short/{sample}.posSorted_srRef.bam.bai"
	wrapper:
		"0.66.0/bio/samtools/index"


rule map_assemblyRef:
	input:
		ref = config["ref_cenpk"],
		assembl = "../" + config["assemblies"] + "/" + "{sample}.fasta",
	output:
		"mapping/reference/assembly/{sample}.Ref.paf"
	log:
		"logs/{sample}/mapping/{sample}_Ref.log"
	shell:
		"minimap2 -cx asm5 {input.ref} {input.assembl} -o {output} 1>{log} 2>&1"




rule time_stamp:
	input:
		expand("mapping/self/long/{sample}.posSorted_lrSelf.bam",sample=SAMPLES),
		expand("mapping/self/short/{sample}.posSorted_srSelf.bam",sample=SAMPLES),
		expand("mapping/self/short/{sample}.posSorted_srSelf.bam.bai",sample=SAMPLES),
		expand("mapping/self/long/{sample}.posSorted_lrSelf.bam.bai",sample=SAMPLES),
		expand("mapping/reference/long/{sample}.posSorted_lrRef.bam",sample=SAMPLES),
		expand("mapping/reference/long/{sample}.posSorted_lrRef.bam.bai",sample=SAMPLES),
		expand("mapping/reference/short/{sample}.posSorted_srRef.bam",sample=SAMPLES),
		expand("mapping/reference/short/{sample}.posSorted_srRef.bam.bai",sample=SAMPLES),
		expand("mapping/reference/assembly/{sample}.Ref.paf",sample=SAMPLES)
	output:
		touch(".make_mapping.touch")
