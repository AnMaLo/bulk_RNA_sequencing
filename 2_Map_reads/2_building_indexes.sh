#!/bin/bash
#SBATCH --mail-type=fail
#SBATCH --job-name="indexing_hisat2"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
#SBATCH --mem=8G

#make a directory to put all hisat files 
mkdir ../hisat2_files 
cd ../hisat2_files 

#call the module for hisat2 
module add UHTS/Aligner/hisat/2.2.1

#Download GTF and make exon, splicesite file.
hisat2_extract_splice_sites.py ../ref_genome/genome.gtf > genome.ss
hisat2_extract_exons.py ../ref_genome/genome.gtf > genome.exon

#Download SNP
#If you want to build HFM index, you can skip this step.
#wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/snp144Common.txt.gz
#gzip -d snp144Common.txt.gz
#Convert chromosome names of UCSC Database to Ensembl Annotation
#awk 'BEGIN{OFS="\t"} {if($2 ~ /^chr/) {$2 = substr($2, 4)}; if($2 == "M") {$2 = "MT"} print}' snp144Common.txt > snp144Common.txt.ensembl
#make SNPs and haplotype file
#hisat2_extract_snps_haplotypes_UCSC.py ../ref_genome/ref_genome.fa snp144Common.txt.ensembl genome_index

#make the index without introns and extrons
hisat2-build -p 16 ../ref_genome/ref_genome.fa genome_dx




