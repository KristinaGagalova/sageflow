#!/usr/bin/bash

HELP="""
Wrapper for producing gff output after BLASTN.
By: Kristina Gagalova

Usage: blast_to_gff_wrapper.sh -q <query file> -d <database file> [OPTIONS]

Options:
	-h	:	Help. What you are reading now.
	-q	:	Query. Put the name/path to your query fasta here. [required]
	-d	:	Database. Put the name/path to your blast database here. [required]
	-p	:	Prefix for all intermediate files, Default: 'blast' [optional]
	-o	:	Output. Put the name to your output gff file, default 'blast' [optional]
	-t	:	Threads. Number of threads/processors you want the blast analysis to run on. (Default: 1)
	-k	:	Keep. This will keep the intermediate blast output. Otherwise only gff is saved. Default: 'true' [true,false]
"""

OUTPUT=0
QUERY=0
DATABASE=0
THREADS=1
KEEP=true
PREFIX=blast
OUTPUT=./blast

while getopts :q:d:o:p:t:k:h: opt; do
  case $opt in
	q)
		echo "-q (query) was input as $OPTARG" >&2
		QUERY=$OPTARG
	;;
	d)
		echo "-d (database) was input as $OPTARG" >&2
		DATABASE=$OPTARG
	;;
	o)
		echo "-o (output) for gff3 was input as $OPTARG" >&2
		if [ $OPTARG != 0 ]; then
			OUTPUT=$OPTARG
		fi
	;;
	p)
		echo "-p (output) for blastdb was input as $OPTARG" >&2
		if [ $OPTARG != 0 ]; then
			PREFIX=$OPTARG
		fi
	;;
	t)
		echo "-t (threads) was input as $OPTARG" >&2
		THREADS=$OPTARG
	;;
	k)
		echo "-k was triggered as $OPTARG" >&2
		KEEP=$OPTARG
	;;
	h)
		echo "$HELP"
		exit 1
	;;
	\?)
		echo "Invalid option: -$OPTARG" >&2
		echo "Type $0 -h for usage"
		exit 1
	;;
  esac
done


if [ $QUERY == 0 ] ; then
	echo "No query detected."
	echo "$HELP"
	exit 1

elif [ $DATABASE == 0 ] ; then
        echo "No database detected."
        echo "$HELP"
        exit 1


elif [ ! "$(command -v makeblastdb)" ] || [ ! "$(command -v blastn)" ]; then
	echo "ERROR: BLAST is not found - exiting program"
	exit 1

elif [ ! "$(command -v gff3sort.pl)" ] || [ ! "$(command -v bgzip)" ] || [ ! "$(command -v tabix)" ]; then
	echo "ERROR: gff3sort.pl or annexed files are not found - exiting program"
	exit 1

else
	mkdir -p ${OUTPUT}/blastdb

	makeblastdb -in $DATABASE -out ${OUTPUT}/blastdb/${PREFIX} -dbtype 'nucl' -hash_index

	#rm previous files
	[ -e ${OUTPUT}/$PREFIX.out ] && rm ${OUTPUT}/$PREFIX.out

	blastn -db ${OUTPUT}/blastdb/$PREFIX \
	-query $QUERY \
	-outfmt "6 qseqid sseqid pident mismatch gapopen qstart qend sstart send evalue length qcovs qcovhsp sstrand" \
	-num_threads $THREADS > ${OUTPUT}/${PREFIX}.out

	#remove previous gff3
	[ -e ${OUTPUT}/${PREFIX}.gff3 ] && rm ${OUPUT}/$PREFIX.gff3

	cat <(echo "##gff-version 3")\
                <(awk -v OFS='\t' '{if ($14 == "plus") print $2,"plasmid","mRNA",$8,$9,$10,"+","0","ID="NR";Name="$1";id="$3";cov="$12";covhsp="$13";len_align="$11;
                else print $2,"plasmid","mRNA",$9,$8,$10,"-","0","ID="NR";Name="$1";id="$3";cov="$12";covhsp="$13";len_align="$11;}' ${OUTPUT}/$PREFIX.out) > ${OUTPUT}/${PREFIX}.gff3

	#DIST=50000
	bedtools merge -d 50000 -i <(sort -k1,1 -k4,4n -k5,5n ${OUTPUT}/${PREFIX}.gff3) | awk -F '\t' -v OFS='\t'  '{print $0,NR}' > ${OUTPUT}/${PREFIX}.tmp
	bedtools intersect -a <(sort -k1,1 -k4,4n -k5,5n ${OUTPUT}/${PREFIX}.gff3) -b ${OUTPUT}/${PREFIX}.tmp -wao | awk -F'\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9";group="$13}' > ${OUTPUT}/${PREFIX}2.gff3
	rm ${OUTPUT}/${PREFIX}.tmp

	paste <(sort -k2,2 -k8,9n ${OUTPUT}/${PREFIX}.out) \
		<(sort -k1,1 -k4,5n ${OUTPUT}/${PREFIX}2.gff3 | awk -F';' '{print $5}') > ${OUTPUT}/${PREFIX}_group.out


	#sort/index
	gff3sort.pl ${OUTPUT}/${PREFIX}2.gff3 > ${OUTPUT}/${PREFIX}.gff3 && rm ${OUTPUT}/${PREFIX}2.gff3
	bgzip ${OUTPUT}/${PREFIX}.gff3 && tabix ${OUTPUT}/${PREFIX}.gff3.gz
fi


if [ $KEEP != true ] ; then
	echo $KEEP
	wait
	rm -r blastdb $PREFIX.out
fi

echo "Exited at: " ; date
exit
