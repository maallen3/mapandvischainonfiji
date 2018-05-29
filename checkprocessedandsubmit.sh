#!/bin/bash
#SBATCH --job-name=nascentsubmit  # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu# Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=1# Number of CPU (processer cores i.e. tasks) In this example I use 32 for bowtie2. I only need one, since none of the commands I run are parallelized.
#SBATCH --time=240:00:00 # Time limit hrs:min:sec
#SBATCH -p long
#SBATCH --mem=10gb # Memory limit
#SBATCH --output=/scratch/Shares/public/nascentdb/processedv2.0/e_and_o/nascentsubmit.%j.out
#SBATCH --error=/scratch/Shares/public/nascentdb/processedv2.0/e_and_o/nascentsubmit.%j.err




module load python/2.7.14/pandas
python checkprocessedandsubmit.py
