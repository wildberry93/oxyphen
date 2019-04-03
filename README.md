# oxyphen
Program for oxygen phenotype assignment based on proteome sequences

## Requirements
- Python2.7
- scikitlearn
- Biopython
- numpy
- pandas

## How to run oxyphen
Specify input, Blast, and number of Blast threads (>1 for faster execution) in SETTINGS file, eg. 

INPUT_FILE=UP000245539_1247513.fasta
BLAST_PATH=/usr/local/ncbi/blast/bin
NUM_THREADS=2

Run oxyphen with a proteome of choice as an input:
python oxyphen.py
