#!/bin/bash
#SBATCH --job-name=scanfold
#SBATCH -n 8
#SBATCH --partition=gen-mk-compute-1
#SBATCH --time=100:00:00
#SBATCH --mem=32G
#SBATCH --mail-type=all
#SBATCH --output=n/MAP_NovaSeq/structure_similarityscanfold.%j.out
#SBATCH --error=n/MAP_NovaSeq/structure_similarity/scanfold.%j.err
#SBATCH --mail-user=herber4@clemson.edu
# Load change dir
cd n/MAP_NovaSeq/structure_similarity/scanfold

source /opt/ohpc/pub/Software/anaconda3/etc/profile.d/conda.sh

conda activate ago

for f in *.fasta; do
	N=$(basename ${f} .fasta);
	python /data2/lackey_lab/austin/gblock/bin/ScanFold/ScanFold.py $f --react ${N}.shape --out_name $N ;
done
