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
OUTPUT=blast

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
	mkdir -p blastdb

	makeblastdb -in $DATABASE -out blastdb/$PREFIX -dbtype 'nucl' -hash_index

	#rm previous files
	[ -e $PREFIX.out ] && rm $PREFIX.out

	blastn -db blastdb/$PREFIX \
	-query $QUERY \
	-outfmt "6 qseqid sseqid pident mismatch gapopen qstart qend sstart send evalue qcovs sstrand" \
	-num_threads $THREADS > $PREFIX.out

	#remove previous gff3
	[ -e $PREFIX.gff3 ] && rm $PREFIX.gff3

	cat <(echo "##gff-version 3")\
		<(awk -v OFS='\t' '{if ($12 == "plus") print $2,"integration","mRNA",$8,$9,$10,"+","0","ID="NR";Name="$1";id="$3";cov="$11;
		else print $2,"integration","mRNA",$9,$8,$10,"-","0","ID="NR";Name="$1";id="$3";cov="$11;}' $PREFIX.out) > $OUTPUT.gff3
	#sort/index
	gff3sort.pl $OUTPUT.gff3 > ${OUTPUT}s.gff3 && mv ${OUTPUT}s.gff3 $OUTPUT.gff3
	bgzip $OUTPUT.gff3 && tabix $OUTPUT.gff3.gz
fi


if [ $KEEP != true ] ; then
	echo $KEEP
	wait
	rm -r blastdb $PREFIX.out
fi

echo "Exited at: " ; date
exit
