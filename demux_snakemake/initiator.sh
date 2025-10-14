# initiator.sh

#!/bin/bash
mkdir -p ./{log,logs_slurm} | sbatch ./snakemake_submittor.sh
