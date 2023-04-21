#!/usr/bin/env python

__author__="Ashikha"

#To move the fasta sequences of core, soft-core, shell and cloud genes into four separate files

import csv, argparse, os

#Using argparse to collect the arguments for file locations of Roary tabular file, fasta file and csv file
parser=argparse.ArgumentParser(description="Description: Gives core, soft-core, shell and cloud gene sequences as four separate output fasta files when the single fasta file containing all the sequences and gene_presence_absence csv file are given as input")
parser.add_argument('-tf', '--tabular_file', help="Give location of the tabular_file containing the summary statistics", required=True)
parser.add_argument('-ff', '--fasta_file', help="Give location of the fasta file containing all the sequences i.e., core, soft core, shell and cloud sequences together in one file", required=True)
parser.add_argument('-cf', '--csv_file', help="Give location of the csv_file containing the gene presence or absence data", required=True)
parser.add_argument('-n', '--no_of_isolates', help="Give the total number of isolates submitted to Roary", required=True)
args=parser.parse_args()

#convert tabular file to csv file
import pandas as pd
sum_stat=pd.read_csv(args.tabular_file, delimiter="\t")
sum_stat.to_csv("summary_statistics.csv", index=None)

#create an empty list to which parameters for filtering % of isolates are added
list_of_filter_no=[]
with open('summary_statistics.csv','r+') as stat:
    for row in csv.reader(stat, delimiter=','):
        row_split=[i for i in row[1].split("%", 1)]
        number=row_split[0].split('(')
        no_list=number[1]
        list_of_filter_no.append(no_list)
trunc_percent_list=list_of_filter_no[0:3]

#Calculating number of isolates from % of isolates

no_of_isolates_list=[float((int(i)/100)*int(args.no_of_isolates)) for i in trunc_percent_list]

#open the CSV file and filter them based on the elements in the no_of_isolates_list 

list_of_core_genes=[]
list_of_softcore_genes=[]
list_of_shell_genes=[]
list_of_cloud_genes=[]

with open(args.csv_file,'r+') as csvf:
    next(csvf) #skip the header row
    for row in csv.reader(csvf, delimiter=','):
        if int(row[3])>=no_of_isolates_list[0]:
            genes=row[0]
            list_of_core_genes.append(genes+"\n")
        elif int(row[3])>=no_of_isolates_list[1] and int(row[3])<no_of_isolates_list[0]:
            genes=row[0]
            list_of_softcore_genes.append(genes+"\n")
        elif int(row[3])>=no_of_isolates_list[2] and int(row[3])<no_of_isolates_list[1]:
            genes=row[0]
            list_of_shell_genes.append(genes+"\n")
        elif int(row[3])<no_of_isolates_list[2]:
            genes=row[0]
            list_of_cloud_genes.append(genes+"\n")

#For creating an empty output folder and getting its path
os.mkdir('output')
os.chdir('output')
out_dir=os.getcwd()

#create and call function to write the resulting list of genes after filtering the gene_presence_absence csv file based on the elements in the no_of_isolates_list
def writeafile(filename, content):
    file1_loc=os.path.join(out_dir, filename)
    file1=open(file1_loc, "w")
    file1.writelines(content)
    file1.close()

writeafile("list_of_core_genes.txt", list_of_core_genes)
writeafile("list_of_softcore_genes.txt", list_of_softcore_genes)
writeafile("list_of_shell_genes.txt", list_of_shell_genes)
writeafile("list_of_cloud_genes.txt", list_of_cloud_genes)

#create a list of listnames which will be used for subsequent list comprehension
listnames=["core_genes", "softcore_genes", "shell_genes", "cloud_genes"]

from Bio import SeqIO
for seq_record in SeqIO.parse(args.fasta_file,"fasta"): #parsing the fasta file which contains all the sequences
    desc=seq_record.description #getting the parsed description data
    genename=desc.split(" ") # splitting the description data
    for list in listnames:
        out_loc=os.path.join(out_dir, list+'_seq.fasta') #create the location for the output fasta file
        out_file=open(out_loc,'a') #open the output fasta in append mode
        in_loc=os.path.join(out_dir,'list_of_'+list+'.txt') #get the location of the output list_of_genes text file
        gene_list=open(in_loc, 'r').read().splitlines() #read the lines in the output list_of_genes text file
        for gene in gene_list:
            if (str(gene)==genename[1]): #if a gene from the list_of_genes is equal to the gene mentioned in the fasta file, then its fasta sequence is collected
                seqs=str(">"+desc+"\n"+seq_record.seq+"\n")
                out_file.write(seqs)
        out_file.close()