#!/bin/bash
#SBATCH --job-name=oh_oh_five
#SBATCH -n 24
#SBATCH --partition=compute
#SBATCH --time=100:00:00
#SBATCH --mem=64G
#SBATCH --mail-type=all
#SBATCH --output=/MAP_NovaSeq/bootstraps/o_o_o_five.%j.out
#SBATCH --error=/MAP_NovaSeq/bootstraps/o_o_o_five.%j.err

source /opt/ohpc/pub/Software/anaconda3/etc/profile.d/conda.sh
conda activate umi_tools
cd /MAP_NovaSeq
seqkit sample -j 24 -p 0.0005 in_vitro_5NIA_2_S11_R1_001.fastq.gz > o_o_o_five_in_vitro_5NIA_2_R1.sampled.fastq

seqkit sample -j 24 -p 0.001 in_vitro_5NIA_2_S11_R2_001.fastq.gz > o_o_o_five_in_vitro_5NIA_2_R2.sampled.fastq

seqkit sample -j 24 -p 0.001 in_vitro_DMSO_2_S12_R1_001.fastq.gz > o_o_o_five_in_vitro_DMSO_2_R1.sampled.fastq

seqkit sample -j 24 -p 0.001 in_vitro_DMSO_2_S12_R2_001.fastq.gz > o_o_o_five_in_vitro_DMSO_2_R2.sampled.fastq

cp o_o_o_five_in_vitro_5NIA_2_R1.sampled.fastq bootstraps/o_o_o_five_nia/
cp o_o_o_five_in_vitro_5NIA_2_R2.sampled.fastq bootstraps/o_o_o_five_nia/
cp o_o_o_five_in_vitro_DMSO_2_R1.sampled.fastq bootstraps/o_o_o_five_dmso/
cp o_o_o_five_in_vitro_DMSO_2_R2.sampled.fastq bootstraps/o_o_o_five_dmso/

cd /MAP_NovaSeq/bootstraps

ml shapemapper/2.2.0

shapemapper --nproc 24 --name One_127_o_o_o_five --min-depth 1000 --out One_127_o_o_o_five --target One_127.fasta --modified --unpaired-folder o_o_o_five_nia/ --untreated --unpaired-folder o_o_o_five_dmso/ --overwrite
shapemapper --nproc 24 --name One_169_o_o_o_five --min-depth 1000 --out One_169_o_o_o_five --target One_169.fasta --modified --unpaired-folder o_o_o_five_nia/ --untreated --unpaired-folder o_o_o_five_dmso/ --overwrite
shapemapper --nproc 24 --name One_213_o_o_o_five --min-depth 1000 --out One_213_o_o_o_five --target One_213.fasta --modified --unpaired-folder o_o_o_five_nia/ --untreated --unpaired-folder o_o_o_five_dmso/ --overwrite
shapemapper --nproc 24 --name One_229_o_o_o_five --min-depth 1000 --out One_229_o_o_o_five --target One_229.fasta --modified --unpaired-folder o_o_o_five_nia/ --untreated --unpaired-folder o_o_o_five_dmso/ --overwrite
shapemapper --nproc 24 --name One_85_o_o_o_five --min-depth 1000 --out One_85_o_o_o_five --target One_85.fasta --modified --unpaired-folder o_o_o_five_nia/ --untreated --unpaired-folder o_o_o_five_dmso/ --overwrite
