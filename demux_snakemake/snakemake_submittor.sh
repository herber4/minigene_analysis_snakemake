#!/bin/bash
#SBATCH --job-name=demux
#SBATCH --ntasks=1
#SBATCH -n 1
#SBATCH --partition=compute
#SBATCH --time=10-00:00:00
#SBATCH --mem=2gb
#SBATCH --output=/scratch2/herber4/output.%j.out
#SBATCH --error=/scratch2/herber4/error.%j.err
#SBATCH --mail-type=all
#SBATCH --mail-user=herber4@clemson.edu

cd /scratch2/herber4
#mkdir -p ./{log,logs_slurm}

source /opt/ohpc/pub/Software/anaconda3/etc/profile.d/conda.sh
conda activate snakemake

#--dag | display | dot
#-p -n \

snakemake \
-s Snakefile \
--latency-wait 120 \
--profile slurm \
