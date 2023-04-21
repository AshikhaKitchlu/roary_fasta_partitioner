# roary_fasta_partitioner
To get fasta files for core, soft-core, shell and cloud genes from Roary output files.

While performing pan-genome analysis using Roary, a single fasta file containing all the gene sequences named "pan-genome_reference.fasta" is generated. Using this python program roary_fasta_partitioner.py, 4 fasta files i.e., one for core, one for soft-core, one for shell and one for cloud genes can be generated from the single fasta file.  


### Input:

Of the many output files generated as a result of Roary, the following files should be given as input for the python program roary_fasta_partitioner.py: the tabular file containing summary statistics, the CSV file containing gene presence absence and the fasta file containing pan-genome reference sequences. Along with these input files, the number of genomes submitted for pan-genome analysis has to be given as input.


### Pre-requirement:

Python

Pandas

csv

argparse

os

SeqIO


### Usage:

`python roary_fasta_partitioner.py -tf <location of tabular file containing summary statistics> -ff <location of fasta file containing pan-genome reference> -cf <location of csv file containing gene presence absence> -n <number of isolates/genomes submitted to Roary for pan-genome analysis>` 

