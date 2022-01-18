#!/bin/bash
#SBATCH --mail-type=fail
#SBATCH --job-name="feature_counts"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=4G

#add module 
module add UHTS/Analysis/subread/2.0.1;

#move to folder make a new directory 
cd ../hisat2_files/mapping
mkdir feature_counts
cd feature_counts 

#run feature counts 
featureCounts -a ../../../ref_genome/genome.gtf -o feature_counts.txt ../*.sorted.bam
