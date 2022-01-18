#!/bin/bash
#SBATCH --mail-type=fail
#SBATCH --job-name="mapping_Normal2"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=13:00:00
#SBATCH --mem=8G

#load the module 
module add UHTS/Aligner/hisat/2.2.1

#map each file
hisat2 -p 4 -q  --dta -x ../../genome_dx -1 ../../../../reads/Normal2_R1.fastq.gz -2 ../../../../reads/Normal2_R2.fastq.gz -S Normal2.sam
