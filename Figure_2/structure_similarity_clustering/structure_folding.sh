#!/bin/bash
#SBATCH --job-name=txsome
#SBATCH -n 32
#SBATCH --partition=compute
#SBATCH --time=100:00:00
#SBATCH --mem=64G
#SBATCH --mail-type=all
#SBATCH --output=/DownloadedSequenceData/austin/MAP_NovaSeq/nia_txsome.%j.out
#SBATCH --error=/DownloadedSequenceData/austin/MAP_NovaSeq/nia_txsome.%j.err

source /opt/ohpc/pub/Software/anaconda3/etc/profile.d/conda.sh
conda activate shapemapper_txsome

ml shapemapper/2.2.0
ml samtools/1.10
cd /DownloadedSequenceData/austin/MAP_NovaSeq
/DownloadedSequenceData/austin/bin/shapemapper-txome/shapemapper-txome --unpaired --modified up_nia_nodemux/ --untreated up_dmso_nodemux/ --target test.fasta --out nia_txsome_no_demux --nproc 32 --max-files-per-folder 5000

### cd in to txsome results shapemaper folder

mkdir shape/

cp */*.shape shape/

### make folder with all fasta files in it

### from shape folder run

ml rnastructure/6.5

for f in *.shape; do N=$(basename $f .shape); Fold ../fasta/${N}.fasta ${N}.ct -mfe -sh $f ; done


for f in *.ct ; do N=$(basename $f .ct); ct2dot $f 1 ${N}.dot ; done
