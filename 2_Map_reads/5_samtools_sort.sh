#!/bin/bash 
#SBATCH --mail-type=fail
#SBATCH --job-name="sort_samtools"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=25G
#SBATCH --dependency=afterok:7368470

#add module  
module add UHTS/Analysis/samtools/1.10;

#move to folder
cd ../hisat2_files/mapping/

#sort the bam files 
samtools sort -@ 4 -m 2G -o NonTNBC1.sorted.bam -T temp NonTNBC1.bam
samtools sort -@ 4 -m 2G -o NonTNBC2.sorted.bam -T temp NonTNBC2.bam 
samtools sort -@ 4 -m 2G -o NonTNBC3.sorted.bam -T temp NonTNBC3.bam

samtools sort -@ 4 -m 2G -o TNBC1.sorted.bam -T temp TNBC1.bam
samtools sort -@ 4 -m 2G -o TNBC2.sorted.bam -T temp TNBC2.bam
samtools sort -@ 4 -m 2G -o TNBC3.sorted.bam -T temp TNBC3.bam

samtools sort -@ 4 -m 2G -o Normal1.sorted.bam -T temp Normal1.bam 
samtools sort -@ 4 -m 2G -o Normal2.sorted.bam -T temp Normal2.bam
samtools sort -@ 4 -m 2G -o Normal3.sorted.bam -T temp Normal3.bam

samtools sort -@ 4 -m 2G -o HER21.sorted.bam -T temp HER21.bam
samtools sort -@ 4 -m 2G -o HER22.sorted.bam -T temp HER22.bam
samtools sort -@ 4 -m 2G -o HER23.sorted.bam -T temp HER23.bam
