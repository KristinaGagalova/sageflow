## Install dependencies: software and libraries

Python libraries list in requirements_py     
```
pip install -r requirements_py.txt
```

Software from list - environment.yaml, to be installed with mamba to env as    
```
mamba env update -n <your-env> --file environment.yml
```

The following are other smaller scripts and tools used in the pipeline.

Wrappers from snakemake (stable releses and used as index)
```
0.79.0/bio/bwa/index
0.66.0/bio/samtools/index
```

Software from ```/home/kgagalova/src``` (called with absolute path in the pipeline)      
```
/home/kgagalova/src/LongStitch/longstitch #-> LongStitch v1.0.1
/home/kgagalova/src/LongQC/longQC.py #-> LongQC 1.2.0c
#called with *python longQC.py*
```

Software available in the repo ```src``` directory, called with relative path in the snakemake.      
```
../src/blast_to_gffv2.sh
../src/gff3sort.pl*
```
