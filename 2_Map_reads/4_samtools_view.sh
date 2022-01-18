#!/bin/bash
#SBATCH --mail-type=fail
#SBATCH --job-name="view_samtools"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --mem=8G

#add the samtools module 
module add UHTS/Analysis/samtools/1.10;

#move to folder
cd ../hisat2_files/mapping/

#convert sam files to bam 
samtools view -hbS NonTNBC1.sam > NonTNBC1.bam
samtools view -hbS NonTNBC2.sam > NonTNBC2.bam
samtools view -hbS NonTNBC3.sam > NonTNBC3.bam

samtools view -hbS TNBC1.sam > TNBC1.bam
samtools view -hbS TNBC2.sam > TNBC2.bam
samtools view -hbS TNBC3.sam > TNBC3.bam

samtools view -hbS Normal1.sam > Normal1.bam
samtools view -hbS Normal2.sam > Normal2.bam
samtools view -hbS Normal3.sam > Normal3.bam

samtools view -hbS HER21.sam > HER21.bam
samtools view -hbS HER22.sam > HER22.bam
samtools view -hbS HER23.sam > HER23.bam

