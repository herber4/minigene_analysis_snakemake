#!/bin/bash
#SBATCH --job-name=boots
#SBATCH -n 24
#SBATCH --partition=compute
#SBATCH --time=100:00:00
#SBATCH --mem=64G
#SBATCH --mail-type=all
#SBATCH --output=/MAP_NovaSeq/bootstraps/correct_bootstraps.%j.out
#SBATCH --error=/MAP_NovaSeq/bootstraps/correct_bootstraps.%j.err

source /opt/ohpc/pub/Software/anaconda3/etc/profile.d/conda.sh
conda activate umi_tools
cd /MAP_NovaSeq
seqkit sample -j 24 -p 0.0001 in_vitro_5NIA_2_S11_R1_001.fastq.gz > o_o_o_five_in_vitro_5NIA_2_R1.sampled.fastq

seqkit sample -j 24 -p 0.0001 in_vitro_5NIA_2_S11_R2_001.fastq.gz > o_o_o_five_in_vitro_5NIA_2_R2.sampled.fastq

seqkit sample -j 24 -p 0.0001 in_vitro_DMSO_2_S12_R1_001.fastq.gz > o_o_o_five_in_vitro_DMSO_2_R1.sampled.fastq

seqkit sample -j 24 -p 0.0001 in_vitro_DMSO_2_S12_R2_001.fastq.gz > o_o_o_five_in_vitro_DMSO_2_R2.sampled.fastq

cp o_o_o_five_in_vitro_5NIA_2_R1.sampled.fastq bootstraps/o_o_o_five_nia/
cp o_o_o_five_in_vitro_5NIA_2_R2.sampled.fastq bootstraps/o_o_o_five_nia/
cp o_o_o_five_in_vitro_DMSO_2_R1.sampled.fastq bootstraps/o_o_o_five_dmso/
cp o_o_o_five_in_vitro_DMSO_2_R2.sampled.fastq bootstraps/o_o_o_five_dmso/

cd /MAP_NovaSeq/bootstraps

ml shapemapper/2.2.0

shapemapper --nproc 24 --name One_205_one --min-depth 1000 --out One_205_one --target One_205.fasta --modified --unpaired-folder o_o_o_five_nia/ --untreated --unpaired-folder o_o_o_five_dmso/ --overwrite

cd ..

seqkit sample -s 127 -j 24 -p 0.0001 in_vitro_5NIA_2_S11_R1_001.fastq.gz > o_o_o_five_in_vitro_5NIA_2_R1.sampled.fastq

seqkit sample -s 127 -j 24 -p 0.0001 in_vitro_5NIA_2_S11_R2_001.fastq.gz > o_o_o_five_in_vitro_5NIA_2_R2.sampled.fastq

seqkit sample -s 127 -j 24 -p 0.0001 in_vitro_DMSO_2_S12_R1_001.fastq.gz > o_o_o_five_in_vitro_DMSO_2_R1.sampled.fastq

seqkit sample -s 127 -j 24 -p 0.0001 in_vitro_DMSO_2_S12_R2_001.fastq.gz > o_o_o_five_in_vitro_DMSO_2_R2.sampled.fastq

cp o_o_o_five_in_vitro_5NIA_2_R1.sampled.fastq bootstraps/o_o_o_five_nia/
cp o_o_o_five_in_vitro_5NIA_2_R2.sampled.fastq bootstraps/o_o_o_five_nia/
cp o_o_o_five_in_vitro_DMSO_2_R1.sampled.fastq bootstraps/o_o_o_five_dmso/
cp o_o_o_five_in_vitro_DMSO_2_R2.sampled.fastq bootstraps/o_o_o_five_dmso/

cd /MAP_NovaSeq/bootstraps

ml shapemapper/2.2.0

shapemapper --nproc 24 --name One_205_two --min-depth 1000 --out One_205_two --target One_205.fasta --modified --unpaired-folder o_o_o_five_nia/ --untreated --unpaired-folder o_o_o_five_dmso/ --overwrite

cd ..

seqkit sample -s 1521 -j 24 -p 0.0001 in_vitro_5NIA_2_S11_R1_001.fastq.gz > o_o_o_five_in_vitro_5NIA_2_R1.sampled.fastq

seqkit sample -s 1521 -j 24 -p 0.0001 in_vitro_5NIA_2_S11_R2_001.fastq.gz > o_o_o_five_in_vitro_5NIA_2_R2.sampled.fastq

seqkit sample -s 1521 -j 24 -p 0.0001 in_vitro_DMSO_2_S12_R1_001.fastq.gz > o_o_o_five_in_vitro_DMSO_2_R1.sampled.fastq

seqkit sample -s 1521 -j 24 -p 0.0001 in_vitro_DMSO_2_S12_R2_001.fastq.gz > o_o_o_five_in_vitro_DMSO_2_R2.sampled.fastq

cp o_o_o_five_in_vitro_5NIA_2_R1.sampled.fastq bootstraps/o_o_o_five_nia/
cp o_o_o_five_in_vitro_5NIA_2_R2.sampled.fastq bootstraps/o_o_o_five_nia/
cp o_o_o_five_in_vitro_DMSO_2_R1.sampled.fastq bootstraps/o_o_o_five_dmso/
cp o_o_o_five_in_vitro_DMSO_2_R2.sampled.fastq bootstraps/o_o_o_five_dmso/

cd /MAP_NovaSeq/bootstraps

ml shapemapper/2.2.0

shapemapper --nproc 24 --name One_205_three --min-depth 1000 --out One_205_three --target One_205.fasta --modified --unpaired-folder o_o_o_five_nia/ --untreated --unpaired-folder o_o_o_five_dmso/ --overwrite

cd ..

seqkit sample -s 2648 -j 24 -p 0.0001 in_vitro_5NIA_2_S11_R1_001.fastq.gz > o_o_o_five_in_vitro_5NIA_2_R1.sampled.fastq

seqkit sample -s 2648 -j 24 -p 0.0001 in_vitro_5NIA_2_S11_R2_001.fastq.gz > o_o_o_five_in_vitro_5NIA_2_R2.sampled.fastq

seqkit sample -s 2648 -j 24 -p 0.0001 in_vitro_DMSO_2_S12_R1_001.fastq.gz > o_o_o_five_in_vitro_DMSO_2_R1.sampled.fastq

seqkit sample -s 2648 -j 24 -p 0.0001 in_vitro_DMSO_2_S12_R2_001.fastq.gz > o_o_o_five_in_vitro_DMSO_2_R2.sampled.fastq

cp o_o_o_five_in_vitro_5NIA_2_R1.sampled.fastq bootstraps/o_o_o_five_nia/
cp o_o_o_five_in_vitro_5NIA_2_R2.sampled.fastq bootstraps/o_o_o_five_nia/
cp o_o_o_five_in_vitro_DMSO_2_R1.sampled.fastq bootstraps/o_o_o_five_dmso/
cp o_o_o_five_in_vitro_DMSO_2_R2.sampled.fastq bootstraps/o_o_o_five_dmso/

cd /MAP_NovaSeq/bootstraps

ml shapemapper/2.2.0

shapemapper --nproc 24 --name One_205_four --min-depth 1000 --out One_205_four --target One_205.fasta --modified --unpaired-folder o_o_o_five_nia/ --untreated --unpaired-folder o_o_o_five_dmso/ --overwrite

cd ..

seqkit sample -s 654 -j 24 -p 0.0001 in_vitro_5NIA_2_S11_R1_001.fastq.gz > o_o_o_five_in_vitro_5NIA_2_R1.sampled.fastq

seqkit sample -s 654 -j 24 -p 0.0001 in_vitro_5NIA_2_S11_R2_001.fastq.gz > o_o_o_five_in_vitro_5NIA_2_R2.sampled.fastq

seqkit sample -s 654 -j 24 -p 0.0001 in_vitro_DMSO_2_S12_R1_001.fastq.gz > o_o_o_five_in_vitro_DMSO_2_R1.sampled.fastq

seqkit sample -s 654 -j 24 -p 0.0001 in_vitro_DMSO_2_S12_R2_001.fastq.gz > o_o_o_five_in_vitro_DMSO_2_R2.sampled.fastq

cp o_o_o_five_in_vitro_5NIA_2_R1.sampled.fastq bootstraps/o_o_o_five_nia/
cp o_o_o_five_in_vitro_5NIA_2_R2.sampled.fastq bootstraps/o_o_o_five_nia/
cp o_o_o_five_in_vitro_DMSO_2_R1.sampled.fastq bootstraps/o_o_o_five_dmso/
cp o_o_o_five_in_vitro_DMSO_2_R2.sampled.fastq bootstraps/o_o_o_five_dmso/

cd /MAP_NovaSeq/bootstraps

ml shapemapper/2.2.0

shapemapper --nproc 24 --name One_205_five --min-depth 1000 --out One_205_five --target One_205.fasta --modified --unpaired-folder o_o_o_five_nia/ --untreated --unpaired-folder o_o_o_five_dmso/ --overwrite
