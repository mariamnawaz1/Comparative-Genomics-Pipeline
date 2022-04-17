# Team1-ComparativeGenomics
## Members
Kaiqin Bian, Rakin Choudhury, Upaasana Krishnan, Mariam Nawaz, Geetha Priyanka Yerradoddi, Wang Zun
FastANI:
FastANI is a tool that calculates the fast alignment-free computation of whole-genome Average Nucleotide Identity (ANI). ANI is defined as the mean nucleotide identity of orthologous gene pairs shared between two microbial genomes. FastANI supports pairwise comparison of both complete and draft genome assemblies. 
Input - Assembled fasta sequences
Output - ANI matrix, Newark file, Calculated Identity file

Download:
An easy option to install fastANI in the provided server environment (T1G4_CG2) is conda install.
CODE: conda install -c bioconda fastANI

Usage:
To know the available command-line options, version, and software usage.
CODE: ./fastANI -h
There are many options to run the command on Linux terminal. For more information, you can check the fastANI developers' GitHub page here.
In our project, we would like to run multiple assembled fasta files against multiple reference genomes. Many to Many. 
CODE: fastANI --ql [queryfasta_list.txt] --rl [reference_list.txt] -o [output_file]
Output_file = <output_file.out>

Visualization:
For visualization purposes, we used two R scripts, fastANI2tree.R and fastANI_heatmap.R.
To run the scripts:
For the phylogenic tree visualization, we need nwk format from the output file of fastANI run.
CODE: Rscript fastANI2tree.R <output_file.out> tree.nwk
	We can visualize the tree using MEGA software. Output using MEGA software is 
	attached below.
	Another way to visualize the tree is using an online GUI website called ETI Toolkit,              
	you can open the website here. Output from this tool is attached below:
For the heatmap visualization, we used the output file from the first fastANI run. 
CODE: Rscript fastANI_heatmap.R 
	Heatmap output from the Rscript is attached below:

