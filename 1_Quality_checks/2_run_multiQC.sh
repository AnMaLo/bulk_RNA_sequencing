#!/bin/bash
#SBATCH --mail-type=fail
#SBATCH --job-name="MultiQC"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00:00
#SBATCH --mem=1G

#load modules
#source /data/courses/rnaseq/breastcancer_de/ananda_workspace/scripts/fastqc_modul.sh
module use /software/module/
module add UHTS/Analysis/MultiQC/1.8

# go to the QC directory
cd /data/courses/rnaseq/breastcancer_de/ananda_workspace/fastQC

#do multiQC
multiqc .
