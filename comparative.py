#!/usr/bin/env python3

# script to perform comparative genomics at whole genome, gene, and SNP level
# command ./comparative.py -i /path/to/directory/having/assembled_contig_sequence/sub_directories -ok /path/to/output/directory -c core[optional; by default = 1] -ism path/of/the/raw/unassembled/sequence/directory

import subprocess
import sys
import re
import argparse
import os


parser = argparse.ArgumentParser(description="Comparative Genomics Analysis. Type comparative.py -help for options")
parser.add_argument("-i", "--inputpath", help="path to the directory that has all the assembled genomes as sub-directories", required=True)
parser.add_argument("-ok", "--outputpath", help="path to the directory that will store the results", required=True)
parser.add_argument("-c", "--cores", help="number of cores; by default = 1", type=str, default="1")
#parser.add_argument("-osm", "--output_stringMLST", help="path to the directory that will store the results", required=True)
parser.add_argument("-ism", "--input_stringMLST", help="path to the directory that has the data", required=True)


args = parser.parse_args()

def run_ksnp(inputpath, outputpath, cores):
    # Creating a kSNP input file contaning the list of genome sequence paths tab separated with genomeID
    ksnp_input_file_name = "in_list.txt"
    in_list_file_path = os.path.join(outputpath, ksnp_input_file_name)
    with open(in_list_file_path, "w") as fp: # creates a file named in_list.txt at the args.outputpath location
        for genome_dir in os.listdir(inputpath): # to read each genome sequence directory
            fp.write(inputpath+"/"+genome_dir+"/"+"contigs.fasta")
            fp.write("\t")
            fp.write(genome_dir) # writing the name of genome_dir in the list file as the contigs file have the same name for all sequences
            fp.write("\n")


    # Creating a combined fasta file to be used as input for MakeFasta utility
    # MakeFasta combines all fasta files into one file to be used as input for next step (Kchooser)
    combine_fasta_loc = outputpath+"/"+"comb_fasta.fasta"
    subprocess.call(["MakeFasta", in_list_file_path, combine_fasta_loc])

    # Finding optimal k-mer using Kchooser utility
    subprocess.call(["Kchooser", combine_fasta_loc])

    #Reading Kchooser.report file to find the optimum K-mer value
    with open("Kchooser.report", "r") as fp:
        lines = fp.readlines() # entire file read as a list named lines

    opt_k = lines[4][26:28] # finding the optimum k which is at line 5 of Kchooser.report
    fck = lines[14][6:11] # finding FCK value from Kchooser.report file

    #Remove combine_fasta file. Path is combine_fasta_loc
    subprocess.call(["rm", combine_fasta_loc])

    # Running kSNP
    subprocess.call(["kSNP3", "-in", in_list_file_path, "-outdir", outputpath, "-k", opt_k, "-ML", "-NJ", "-CPU", cores])

    print ("FCK value =", fck)

def makelist_stringMLST(wd, td): 
    fn = wd+'/list.txt'
    os.chdir(td)
    samples = os.listdir(".")

    #Find all the paired samples and seperate them with tab,then write into the file
    for sample in samples:
        os.chdir(sample)
        #os.system("gunzip *")
        with open(fn, 'a') as f:
            f.write(td+sample+"/"+sample+"_1.fq\t"+td+sample+"/"+sample+"_2.fq\n")
        os.chdir("..")
    f.close()

def run_stringMLST(wd, td): # results are stored as a stringMLST_result.txt in the current working directory
# 1. Database
    # 1.1 Download datasets:
    # subprocess.call(["stringMLST.py", "--getMLST", "-P", "datasets/", "--species", "all"])
    # Or Download the individual Neisseria spp. scheme
    os.system("stringMLST.py --getMLST -P Salmonella/se --species Salmonella -k 133")

    # 1.2 Database building:
    # can also create database on your own
    #os.system("stringMLST.py --buildDB -k 133 -P Salmonella/se")
# 2. Prediction
    makelist_stringMLST(wd, td)
    os.chdir(wd)
    # Usage
    os.system("stringMLST.py --predict -l list.txt -P Salmonella/se -k 133 -o stringMLST_result.txt")
    # deleting the created database
    os.system("rm -r Salmonella")

def movefile_chewbbaca(inputpath, workingpath):
    dirs = os.listdir(inputpath)	
    for dire in dirs:
        os.chdir(inputpath+"/"+dire)
        print(os.getcwd())
        fs = os.listdir(".")
        for f in fs:
            if f[-5:] == "fasta":
                os.system("cp "+f+" "+workingpath+"/fasta_files/"+dire+".fasta")
        #os.chdir(inputpath)

def chewbbaca(inputpath, workingpath):
    #0: fetch contigs into a dir
    os.system("mkdir "+"fasta_files")
    movefile_chewbbaca(inputpath, workingpath)
    os.chdir(workingpath)

    #1: download database Salmonella enterica
    os.system("chewBBACA.py DownloadSchema -sp 8 -sc 1 -o SE")
    
    os.system("mv SE/* SE_db")
    
    os.system("rm -r SE")

	#2: Allele calling
     
    os.system("chewBBACA.py AlleleCall -i fasta_files -g SE_db -o AC --cpu 4 --ptf SE_db/Salmonella_enterica.trn")
    
    os.system("mv AC/* ACresults")
    
    os.system("rm -r AC")
	
	#3: quality check
    os.system("chewBBACA.py TestGenomeQuality -i ACresults/results_alleles.tsv -n 1 -t 100 -s 5 -o QC")

	#4: extract loci in 95% matrix
    os.system("chewBBACA.py ExtractCgMLST -i ACresults/results_alleles.tsv --r ACresults/RepeatedLoci.txt --g QC/removedGenomes.txt -o output --t 0")
	
	#5: tree drawing
    os.system("grapetree --profile output/cgMLST.tsv --method MSTreeV2 > result.nwk")
    os.system("plottree result.nwk -l 8.4 -o chewbbaca_tree")



if __name__ == "__main__":
    print("running kSNP")
    run_ksnp(args.inputpath, args.outputpath, args.cores)
    wd = os.getcwd()
    td = args.input_stringMLST
    print("running stringMLST")
    run_stringMLST(wd, td)
    print("running chewbbaca")
    workingpath = os.getcwd()
    chewbbaca(args.inputpath, workingpath)
    print("Done! Exiting Now")
    
