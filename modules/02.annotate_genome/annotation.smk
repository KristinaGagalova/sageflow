# Author: Kristina Gagalova
# Submodule: annotation
# to be run with the main-snakefile

from snakemake.utils import min_version
min_version("5.24.0")
from pathlib import Path

#get path current dir
MY_PATH = os.getcwd()
#print(MY_PATH)

#create directories
#new_directory = Path('02-annotation')
#new_directory.mkdir(exist_ok=True)

#config file
configfile: "../assembly-vars.yaml"

SAMPLES, = glob_wildcards("../" + config['assemblies'] + "/" + "{sample}.fasta")
print(SAMPLES)

rule all:
	input:
		".make_annotation.touch"

rule prepare_anno:
	input:
		assembly = "../" + config["assemblies"] + "/" + "{sample}.fasta",
	output:
		temp("00-prerun/{sample}_cut.fasta")
	shell:
		"""
		cut -d" " -f1 {input} > {output}
		"""

rule prepare_funanno:
	input:
		rules.prepare_anno.output
	output:
		"00-prerun/{sample}sort.fasta"
	log:
		"logs/{sample}/00-prerun/{sample}_prepare.log"
	shell:
		"""
		seqkit sort -lr {input} > {output} 2> {log}  
		"""

rule funanno_clean:
	input:
		rules.prepare_funanno.output
	output:
		"00-prerun/{sample}clean.fasta"
	log:
		"logs/{sample}/00-prerun/{sample}_clean.log"
	shell:
		"funannotate clean -i {input} \
		-o {output} \
		--exhaustive 2> {log}" 

rule mask_reps:
	input:
		rules.prepare_funanno.output
	output:
		"00-prerun/{sample}mask.fasta"
	log:
		"logs/{sample}/00-prerun/{sample}_reps.log"
	shell:
		"""
		funannotate mask -i {input} \
		-s 'Saccharomyces cerevisiae' \
		-o {output} 2> {log}
		"""

rule fuannotate_predict:
	input:
		rules.mask_reps.output
	output:
		"01-predict/status/{sample}_predict.done"
	log:
		"logs/{sample}/01-prediction/{sample}_fun_predict.log"
	params:
		source_fa = "../" + config["source_funanno"],
		prot_sc = config["prot_sc"],
		prot_uniprot = config["prot_uni"]
	threads:
		12
	shell:	
		"source {params.source_fa} && \
		funannotate predict -i {input} -o 01-predict/{wildcards.sample} \
		-s 'Saccaromyces cerevisiae' \
		--name {wildcards.sample} --force \
		--protein_evidence {params.prot_sc} {params.prot_uniprot} \
		--cpus {threads} \
		--augustus_species saccharomyces 2> {log} && \
		touch {output}"


rule funannotate_anno:
	input:
		rules.fuannotate_predict.output
	output:
		status = "01-predict/status/{sample}_annotate.done",
		gff = "01-predict/{sample}/annotate_results/Saccharomyces_cerevisiae_{sample}.gff3",
		gbk = "01-predict/{sample}/annotate_results/Saccharomyces_cerevisiae_{sample}.gbk",
		mrna = "01-predict/{sample}/annotate_results/Saccharomyces_cerevisiae_{sample}.mrna-transcripts.fa",
		cds = "01-predict/{sample}/annotate_results/Saccharomyces_cerevisiae_{sample}.cds-transcripts.fa",
		prot = "01-predict/{sample}/annotate_results/Saccharomyces_cerevisiae_{sample}.proteins.fa",
		tab = "01-predict/{sample}/annotate_results/Saccharomyces_cerevisiae_{sample}.annotations.txt"
	log:
		"logs/{sample}/02-annotation/{sample}_fun_anno.log"
	threads:
		12
	params:
		source_fa = "../" + config["source_funanno"],
	shell:
		"source {params.source_fa} && \
		funannotate annotate -i 01-predict/{wildcards.sample}/predict_results \
		-o {output.status} \
		--cpus {threads} \
		--force \
		--strain {wildcards.sample} \
		--busco_db saccharomycetales 2> {log} && \
		touch {output}"
		
rule sort_gff:
	input:
		"01-predict/{sample}/annotate_results/Saccharomyces_cerevisiae_{sample}.gff3"
	output:
		"02-annotation_final/{sample}/{sample}.gff3.gz"
	log:
		"logs/{sample}/reformat/{sample}_gff_rename.log"
	shell:
		"../src/gff3sort.pl {input} > {wildcards.sample}s.gff3 && mv {wildcards.sample}s.gff3 02-annotation_final/{wildcards.sample}/{wildcards.sample}.gff3 && \
		bgzip -c 02-annotation_final/{wildcards.sample}/{wildcards.sample}.gff3 > {output} && \
		tabix {output}"

rule copy_final:
	input:
		gbk = "01-predict/{sample}/annotate_results/Saccharomyces_cerevisiae_{sample}.gbk",
		mrna = "01-predict/{sample}/annotate_results/Saccharomyces_cerevisiae_{sample}.mrna-transcripts.fa", 
		cds = "01-predict/{sample}/annotate_results/Saccharomyces_cerevisiae_{sample}.cds-transcripts.fa",
		prot = "01-predict/{sample}/annotate_results/Saccharomyces_cerevisiae_{sample}.proteins.fa",
		tab = "01-predict/{sample}/annotate_results/Saccharomyces_cerevisiae_{sample}.annotations.txt"
	output:
		gbk = "02-annotation_final/{sample}/{sample}.gbk",
		mrna = "02-annotation_final/{sample}/{sample}.mrna.fa",
		cds = "02-annotation_final/{sample}/{sample}.cds.fa",
		prot = "02-annotation_final/{sample}/{sample}.proteins.fa",
		tab = "02-annotation_final/{sample}/{sample}.annotations.txt",
	shell:
		"cp {input.gbk} {output.gbk} && \
		cp {input.mrna} {output.mrna} && \
		cp {input.cds} {output.cds} && \
		cp {input.prot} {output.prot} && \
		cp {input.tab} {output.tab}"


#------------------------------------------------------------------------------------

rule annotate_insert:
	input:
		assembly = "../" + config["assemblies"] + "/" + "{sample}.fasta",
		ins = config["insert"] + "/" + "{sample}_integrated.fasta"
	output:
		"03-insert/inserts/{sample}_insert.gff3.gz"
	log:
		"logs/{sample}/insert/{sample}_insert.log"
	params:
		con = config["con_insert"]
	shell:
		"../src/blast_to_gffv2.sh -q {input.ins} -d {input.assembly} -c {params.con} -p {wildcards.sample}_insert -k true -o 03-insert/inserts > {log} "


rule report_insert:
	input:
		rules.annotate_insert.output
	output:
		summary = "04-reports/inserts/{sample}.tsv",
		summary_counts = "04-reports/inserts/{sample}_c.tsv",
	shell:
		"../src/Rscript/insert_coord/insert_coord.R {input} {output.summary} {output.summary_counts}"

rule annotate_plasmid:
	input:
		assembly = "../" + config["assemblies"] + "/" + "{sample}.fasta",
		plasmid = config["plasmid"]
	output:
		"03-insert/plasmid/{sample}_plasmid.gff3.gz"
	log:
		"logs/{sample}/insert/{sample}_plasmid.log"
	shell:
		"../src/blast_to_gff_plasmid.sh -q {input.plasmid} -d {input.assembly} -p {wildcards.sample}_plasmid -k true -o 03-insert/plasmid > {log}"


rule time_stamp_annotation:
	input:
		expand("03-insert/inserts/{sample}_insert.gff3.gz",sample=SAMPLES),
		expand("03-insert/plasmid/{sample}_plasmid.gff3.gz",sample=SAMPLES),
		expand("01-predict/status/{sample}_annotate.done",sample=SAMPLES),
		expand("01-predict/status/{sample}_predict.done",sample=SAMPLES),
		expand("02-annotation_final/{sample}/{sample}.gff3.gz",sample=SAMPLES),
		expand("02-annotation_final/{sample}/{sample}.annotations.txt",sample=SAMPLES),
		expand("04-reports/inserts/{sample}_c.tsv",sample=SAMPLES)
	output:
		touch(".make_annotation.touch")
