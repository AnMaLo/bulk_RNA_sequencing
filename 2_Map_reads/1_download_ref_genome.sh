#!/bin/bash
#SBATCH --mail-type=fail
#SBATCH --job-name="reference_genomes"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=10:00:00
#SBATCH --mem=25G

#make a directory for the reference genome and annotation 
mkdir ../ref_genome
cd ../ref_genome

#Download on Ensembl ftp site reference genome and Checksum
#sequence: file named species.assembly.dna.primary_assembly.fa.gz (under DNA(FASTA)). 
wget ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

#Download the annotation: named species.assembly.build.gtf.gz (gtf (under gene sets))
wget ftp://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz

#download checksums and print the relevant checksum into the specific file
wget ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/CHECKSUMS
awk '$3 =="Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"' CHECKSUMS > sum_reference_sequence.txt 
rm CHECKSUMS 

wget ftp://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/CHECKSUMS
awk '$3 =="Homo_sapiens.GRCh38.104.gtf.gz"' CHECKSUMS > sum_annotation.txt 
rm CHECKSUMS

#run sum to make sure the files are complete 
sum Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz >> sum_reference_sequence.txt 
sum Homo_sapiens.GRCh38.104.gtf.gz >> sum_annotation.txt

#unzip the files and rename them 
gzip -d Homo_sapiens.GRCh38.104.gtf.gz
gzip -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa ref_genome.fa
mv Homo_sapiens.GRCh38.104.gtf genome.gtf
