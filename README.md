# ESILs
The bash script is used for alignement of fastq files and preparing the list of Single Nucleotide Polymorphisms (SNPs) in a short genetic interval
for a studied Arabidopsis cross. This is applied to sequenced reads, coming from a long-range PCR reaction, however it might be used for whole-genome-sequencing 
data. Firstly, the paired-end reads are pooled and aligned to a reference sequence. Resulting bam file is sorted and indexed. Later samtools mpileup and bcftools 
are usedto prepare a .txt extension list of all the SNPs found in a studied DNA sequence. Similarly, individual paired-end fastq files are aligned to the sequence
and later a list of SNPs with number of reads per each SNP is produced.
The R script is used to filter out the SNPs found in the individuals, assign a genotype and visualize the data.
Script was prepared by Wojciech Dzięgielewski, Laboratory of Genome Biology, Institute of Molecular Biology and Biotechnology, Adam Mickiewicz University, Poznan.
