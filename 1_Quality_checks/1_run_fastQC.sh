#!/bin/bash
#SBATCH --mail-type=fail
#SBATCH --job-name="FastQC"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00:00
#SBATCH --mem=1G

#load module
module add UHTS/Quality_control/fastqc/0.11.7

#create and go to the QC directory
mkdir /data/courses/rnaseq/breastcancer_de/ananda_workspace/fastQC
cd /data/courses/rnaseq/breastcancer_de/ananda_workspace/fastQC

#copy reads into own folder
cp -r /data/courses/rnaseq/breastcancer_de/reads/*.fastq.gz .

#check quality of your data with fastqc
fastqc -t 2 *.fastq.gz
