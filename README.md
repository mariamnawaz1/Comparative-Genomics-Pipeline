# Team1-ComparativeGenomics
## Members
Kaiqin Bian, Rakin Choudhury, Upaasana Krishnan, Mariam Nawaz, Geetha Priyanka Yerradoddi, Wang Zun
## FastANI:
FastANI is a tool that calculates the fast alignment-free computation of whole-genome Average Nucleotide Identity (ANI). ANI is defined as the mean nucleotide identity of orthologous gene pairs shared between two microbial genomes. FastANI supports pairwise comparison of both complete and draft genome assemblies. 
Input - Assembled fasta sequences
Output - ANI matrix, Newark file, Calculated Identity file

### Download:
An easy option to install fastANI in the provided server environment (T1G4_CG2) is conda install.
``` 
conda install -c bioconda fastANI 
```

### Usage:
To know the available command-line options, version, and software usage.
```
./fastANI -h
```
There are many options to run the command on Linux terminal. For more information, you can check the fastANI developers' GitHub page here.
In our project, we would like to run multiple assembled fasta files against multiple reference genomes. Many to Many. 
```
fastANI --ql [queryfasta_list.txt] --rl [reference_list.txt] -o [output_file]
```
Output_file = <output_file.out>

### Visualization:
For visualization purposes, we used two R scripts, fastANI2tree.R and fastANI_heatmap.R.
To run the scripts:
For the phylogenic tree visualization, we need nwk format from the output file of fastANI run.
```
Rscript fastANI2tree.R <output_file.out> tree.nwk
```
We can visualize the tree using MEGA software. Output using MEGA software is attached below:
![fastanimegatree](https://github.gatech.edu/storage/user/57475/files/d9e6708e-7cb3-4f25-9b9a-20e8f7733f7c)

Another way to visualize the tree is using an online GUI website called ETI Toolkit, you can open the website here. Output from this tool is attached below:
![fastani_etetool](https://github.gatech.edu/storage/user/57475/files/25d1028c-b51b-44c7-a883-1238c0208f46)

For the heatmap visualization, we used the output file from the first fastANI run. 
```
Rscript fastANI_heatmap.R 
```
Heatmap output from the Rscript is attached below:
![fastani_heatmap](https://github.gatech.edu/storage/user/57475/files/5eac5793-dad5-4c4e-83e7-96908b042528)

## kSNP3.1
kSNP is a SNP discovery and annotation tool which identifies pan-genome SNPs in a set of genome sequences and estimates phylogenetic trees based on those SNPs. It is based on k-mer analysis and doesn't require a reference genome or multiple sequence alignment. Therefore, it can take 100's of microbial genomes as input and those could be finished or unfinished genomes in assembled or unassembled reads. Finished and unfinished genomes can be analyzed together. For annotation, kSNP automatically download the Genbank files of the finished genomes and incorporate information from those files.

### Installation of kSNP
```
wget https://sourceforge.net/projects/ksnp/files/kSNP3.1_Linux_package.zip
unzip kSNP3.1_Linux_package.zip
vi ~/.bashrc
export PATH=$PATH:/path/to/kSNP3.1_Linux_package/kSNP3
source ~/.bashrc
which tcsh
copy the tcsh path
cd /path/to/kSNP3.1_Linux_package/kSNP3
vi kSNP3
set the shebang line to: #!/copied/path/of/tcsh/starting/from/home
set kSNP=/path/to/kSNP3.1_Linux_package/kSNP3
close the terminal and re-open
```

### Usage
**Step 1**: Create a text file containing path and GenomeID of all the genome sequence files. There are several ways to do that are given in the [kSNP User Manual](https://sourceforge.net/projects/ksnp/files/kSNP3.1.2%20User%20Guide%20.pdf/download). Format should be:
```
Path/to/the/1st/sequence/file	1-tabspace	 GenomeID
Path/to/the/2nd/sequence/file	1-tabspace	 GenomeID
```

**Step 2**: Create a combined Fasta file which will be input for the next step of determing optimal k-mer
```
MakeFasta <in_list.txt> <outfile_name.fasta>
```
where in_list.txt is the file from step 1.

**Step 3**: Finding optimal k-mer size
To identify a SNP, kSNP matches k-mers between different genomes and if they are identical except for the central base, it is counted as a SNP. The utility program is Kchooser.
```
Kchooser outfile_name.fasta
```

**Step 4**: Open the Kchoose.report and find the optimum value for k-mer.

**Step 5**: Run kSNP3
```
kSNP3 -in <in_list.txt> -outdir <outdir> -k <optimal_k-mer> -ML -NJ -CPU <number_of_CPUs> 
```

### Phylogenetic Trees
kSNP generates different form of phylogenetic tree files using Parsimony, Maximum likelihood (ML), and Neighbor joining (NJ) methods. Below are the tipAlleleCounts ML and NJ trees from the kSNP (step 5 above) visualized using [dendroscope](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/dendroscope/).

Phylogenetic tree using the Maximum likelihood method. The internal nodel labels show the number of SNP alleles that are present in all descendants of that node and nowhere else. Allele counts for branch tips are not shown. Strain names have been modified to show the strain specific allele counts:
![tree_tipAlleleCounts.ML](https://github.gatech.edu/computationalgenomics2022/Team1-ComparativeGenomics/blob/main/results/tree_tipAlleleCounts.ML.png)

Phylogenetic tree using the Neighbor joining method:
![tree_tipAlleleCounts.NJ](https://github.gatech.edu/computationalgenomics2022/Team1-ComparativeGenomics/blob/main/results/tree_tipAlleleCounts.NJ.png)
