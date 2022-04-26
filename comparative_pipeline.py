#!/usr/bin/env python3

# script to perform comparative genomics at whole genome, gene, and SNP level
# command ./comparative_pipeline.py -i /path/to/directory/having/assembled_contig_sequence/sub_directories -ok /path/to/output/directory -c core[optional; by default = 1] -ism path/of/the/raw/unassembled/sequence/directory

import subprocess
import sys
import re
import argparse
import os
from multiprocessing import Pool

parser = argparse.ArgumentParser(description="Comparative Genomics Analysis. Type comparative_pipeline.py -help for options")
parser.add_argument('-f','--fastani', help = "To run fastANI tool", type = str)
parser.add_argument('-k','--ksnp', help = "To run kSNP3", type = str)
parser.add_argument('-ch','--chewbbaca', help = "To run chewBBACA", type = str)
parser.add_argument('-s','--stringmlst', help = "To run stringMLST", type = str)
parser.add_argument('-a','--amrfinder', help = "To run AMRFinderPlus", type = str)
parser.add_argument('-m','--mummer', help = "To run mummer4", type = str)
parser.add_argument("-i", "--inputpath", help="path to the directory that has all the assembled genomes as sub-directories", required=True)
parser.add_argument("-ok", "--outputpath", help="path to the directory that will store the results", required=True)
parser.add_argument("-o","--output_file", type = str, help = "to specify output file")
parser.add_argument("-c", "--cores", help="number of cores; by default = 1", type=str, default="1")
parser.add_argument("-ism", "--input_stringMLST", help="path to the directory that has the data", required=True)


args = parser.parse_args()

for files in os.listdir(args.i):
	f = os.path.join(arg.i, files) + "/contigs.fasta"
	print(f, file = open('inputpath_file.txt', 'a'))

# Function for fastANI
def run_fastani(inputpath_file, output):
	# Enter the paths of all the isolates you want to run with the fastANI tool
	os.system(f"fastANI --ql {inputpath_file} --rl {inputpath_file} -o {output}_fastani.out")
	os.system(f"Rscript fastANI2tree.R {output} {output}.nwk")
	# return output



# Function for AMRFinder:
def run_amr(fasta, output):
	output = output + "_amr.txt"
	with open(fasta, 'r') as f:
		files = [file.strip().split("\n") for file in f.readlines()]
	# print(files)
	for file in files:
		for x in file:
			os.system(f"amrfinder --nucleotide {x} --plus --organism Salmonella > {output}")



# Function for MUMmer4
def run_mummer(input_dir, output_dir, output_file, thread):

    if os.path.exists(output_dir) is False:
        os.makedirs(output_dir)
    name_list = os.listdir(input_dir)
    name_list.sort()
    CWD = os.getcwd()

    list1 = []
    for i in range(len(name_list)):
        list1.append(os.path.join(input_dir, name_list[i], 'contigs.fasta'))
    parameter_list = []
    count_list = []

    # creates a list of each unique combination of fasta files
    for i in range(len(list1)):
        for j in range(len(list1)):
            if i != j and [j, i] not in count_list:
                parameter_list.append(['output_{i}_{j}'.format(i=i+1, j=j+1), list1[i], list1[j], i+1, j+1])
                count_list.append([i, j])

    # creates a subprocess list set for each unique combination of fasta files
    subprocess_list = []
    for i in range(len(parameter_list)):
        subprocess_list.append(['dnadiff', '-p', parameter_list[i][0], parameter_list[i][1], parameter_list[i][2]])

    # matrix to store results
    matrix = []
    for i in range(len(list1) + 1):
        matrix.append([])
        for j in range(len(list1) + 1):
            matrix[i].append('')

    # gives the first row and first column their genomeX.fasta names
    for i in range(len(matrix) - 1):
        matrix[i + 1][0] = name_list[i]
        matrix[0][i + 1] = name_list[i]

    # substitutes the 1:1 combinations with 100 since they'd be a perfect match
    for i in range(len(matrix) - 1):
        matrix[i+1][i+1] = '100.0'

    # runs the generated subprocess lists in parallel depending on the number of threads given in -c
    pool = Pool(int(thread))
    pool.map(subprocess.call, subprocess_list)
    pool.close()

    # gets the 19th line of each .report files, then enters the last number into the matrix
    # deletes all files made using dnadiff
    file_list = os.listdir(CWD)
    file_list.sort()
    for i in range(len(parameter_list)):
        filename = open('{}.report'.format(parameter_list[i][0]), 'r')
        content = filename.readlines()
        list3 = []
        piece = content[18].rstrip()
        list3.append(piece.split(' '))
        matrix[parameter_list[i][3]][parameter_list[i][4]] = list3[0][-1]
        matrix[parameter_list[i][4]][parameter_list[i][3]] = list3[0][-1]
        filename.close()
        if '{}.1coords'.format(parameter_list[i][0]) in file_list:
            subprocess.call(['rm', '{}.1coords'.format(parameter_list[i][0])])
            file_list.remove('{}.1coords'.format(parameter_list[i][0]))
        if '{}.1delta'.format(parameter_list[i][0]) in file_list:
            subprocess.call(['rm', '{}.1delta'.format(parameter_list[i][0])])
            file_list.remove('{}.1delta'.format(parameter_list[i][0]))
        if '{}.delta'.format(parameter_list[i][0]) in file_list:
            subprocess.call(['rm', '{}.delta'.format(parameter_list[i][0])])
            file_list.remove('{}.delta'.format(parameter_list[i][0]))
        if '{}.mcoords'.format(parameter_list[i][0]) in file_list:
            subprocess.call(['rm', '{}.mcoords'.format(parameter_list[i][0])])
            file_list.remove('{}.mcoords'.format(parameter_list[i][0]))
        if '{}.mdelta'.format(parameter_list[i][0]) in file_list:
            subprocess.call(['rm', '{}.mdelta'.format(parameter_list[i][0])])
            file_list.remove('{}.mdelta'.format(parameter_list[i][0]))
        if '{}.qdiff'.format(parameter_list[i][0]) in file_list:
            subprocess.call(['rm', '{}.qdiff'.format(parameter_list[i][0])])
            file_list.remove('{}.qdiff'.format(parameter_list[i][0]))
        if '{}.rdiff'.format(parameter_list[i][0]) in file_list:
            subprocess.call(['rm', '{}.rdiff'.format(parameter_list[i][0])])
            file_list.remove('{}.rdiff'.format(parameter_list[i][0]))
        if '{}.report'.format(parameter_list[i][0]) in file_list:
            subprocess.call(['rm', '{}.report'.format(parameter_list[i][0])])
            file_list.remove('{}.report'.format(parameter_list[i][0]))
        if '{}.snps'.format(parameter_list[i][0]) in file_list:
            subprocess.call(['rm', '{}.snps'.format(parameter_list[i][0])])
            file_list.remove('{}.snps'.format(parameter_list[i][0]))
        if '{}.unqry'.format(parameter_list[i][0]) in file_list:
            subprocess.call(['rm', '{}.unqry'.format(parameter_list[i][0])])
            file_list.remove('{}.unqry'.format(parameter_list[i][0]))
        if '{}.unref'.format(parameter_list[i][0]) in file_list:
            subprocess.call(['rm', '{}.unref'.format(parameter_list[i][0])])
            file_list.remove('{}.unref'.format(parameter_list[i][0]))


    # opens the output file and writes the matrix into it
    final_file = open(str(output_file), 'w')
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if j != len(matrix) - 1:
                final_file.write(str(matrix[i][j])+'\t')
            else:
                final_file.write(str(matrix[i][j])+'\n')
    final_file.close()

    subprocess.call(['mv', output_file, output_dir])


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
	if args.k:
		print("running kSNP")
		run_ksnp(args.inputpath, args.outputpath, args.cores)
	if args.s:
		wd = os.getcwd()
		td = args.input_stringMLST
		print("running stringMLST")
		run_stringMLST(wd, td)
	if args.ch:
		print("running chewbbaca")
		workingpath = os.getcwd()
		chewbbaca(args.inputpath, workingpath)
	if args.f:
		print("Running fastANI")
		run_fastani('inputpath_file.txt', args.o)
	if args.a:
		print("Running AMRFinderPlus")
		run_amr('inputpath_file.txt', output)
	if args.m:
		print("Running MUMmer4")
		run_mummer(args.i, args.ok , args.o, args.c)

	print("Done! Exiting Now")
    
    
