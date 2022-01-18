#!/bin/bash
# this script sbatches all mapping jobs 

sbatch HER21.sh
sbatch HER22.sh
sbatch HER23.sh
sbatch NonTNBC1.sh
sbatch NonTNBC2.sh
sbatch NonTNBC3.sh
sbatch TNBC1.sh
sbatch TNBC2.sh
sbatch TNBC3.sh
sbatch Normal1.sh
sbatch Normal2.sh
sbatch Normal3.sh

