#!/bin/bash
#SBATCH --mail-type=fail
#SBATCH --job-name="index_samtools"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --mem=4G
#SBATCH --dependency=afterok:7363425

#add module
module add UHTS/Analysis/samtools/1.10;

#move to folder
cd ../hisat2_files/mapping/

#sort the bam files
samtools index NonTNBC1.sorted.bam
samtools index NonTNBC2.sorted.bam
samtools index NonTNBC3.sorted.bam

samtools index TNBC1.sorted.bam
samtools index TNBC2.sorted.bam
samtools index TNBC3.sorted.bam

samtools index Normal1.sorted.bam
samtools index Normal2.sorted.bam
samtools index Normal3.sorted.bam

samtools index HER21.sorted.bam
samtools index HER22.sorted.bam
samtools index HER23.sorted.bam                                 
