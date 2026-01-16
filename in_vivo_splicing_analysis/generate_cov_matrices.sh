in vivo analysis:

1. shapemapper-txsome for read demuxxing
- One thing i might need to do is run shapemapper txsome with only the reverse read


#!/bin/bash
#SBATCH --job-name=txsome
#SBATCH -n 32
#SBATCH --partition=compute
#SBATCH --time=100:00:00
#SBATCH --mem=128G
#SBATCH --mail-type=all
#SBATCH --output=/DownloadedSequenceData/austin/MAP_NovaSeq/in_vivo/rev_txsome_mapping.%j.out
#SBATCH --error=/DownloadedSequenceData/austin/MAP_NovaSeq/in_vivo/rev_txsome_mapping.%j.err

source /opt/ohpc/pub/Software/anaconda3/etc/profile.d/conda.sh
conda activate shapemapper_txsome

ml shapemapper/2.2.0
ml samtools/1.10
cd /DownloadedSequenceData/austin/MAP_NovaSeq/in_vivo/
/DownloadedSequenceData/austin/bin/shapemapper-txome/shapemapper-txome --unpaired --modified in_vivo_K700E_1_S7_R2_001.fastq.gz --untreated in_vivo_MAP_NTC_1_S5_R2_001.fastq.gz --target in_vivo_oligos.fasta --out k700e_vs_NTC_txsome_rep_one_rev/ --nproc 32 --max-files-per-folder 5000

/DownloadedSequenceData/austin/bin/shapemapper-txome/shapemapper-txome --unpaired --modified in_vivo_K700E_2_S8_R2_001.fastq.gz --untreated in_vivo_MAP_NTC_2_S6_R2_001.fastq.gz --target in_vivo_oligos.fasta --out k700e_vs_NTC_txsome_rep_two_rev/ --nproc 32 --max-files-per-folder 5000

/DownloadedSequenceData/austin/bin/shapemapper-txome/shapemapper-txome --unpaired --modified in_vivo_K700E_3_S13_R2_001.fastq.gz --untreated in_vivo_MAP_NTC_3_S14_R2_001.fastq.gz --target in_vivo_oligos.fasta --out k700e_vs_NTC_txsome_rep_three_rev/ --nproc 32 --max-files-per-folder 5000

/DownloadedSequenceData/austin/bin/shapemapper-txome/shapemapper-txome --unpaired --modified in_vivo_K700E_1_S7_R2_001.fastq.gz --untreated in_vivo_MAP_NTC_1_S5_R2_001.fastq.gz --target in_vivo_oligos.fasta --out txsome_rep_one_rev/ --nproc 32 --max-files-per-folder 5000

2. STAR alignment for split read mapping of txsome fastqs

### build STAR genome index

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir . --genomeFastaFiles in_vivo_oligos.fasta


### run star

STAR --genomeDir star_index/ --runThreadN 8 --readFilesIn txsome_rep_one_rev/fastq_by_target/modified/A/One_205.fastq --outFileNamePrefix One_205_txsome_rev_only  --outSAMtype BAM SortedByCoordinate  --outSAMunmapped None  --outSAMattributes Standard

samtools index One_205_txsome_rev_onlyAligned.sortedByCoord.out.bam

3. blah blah blah


#!/bin/bash
#SBATCH --job-name=star
#SBATCH -n 16
#SBATCH --partition=compute
#SBATCH --time=100:00:00
#SBATCH --mem=128G
#SBATCH --mail-type=all
#SBATCH --output=/DownloadedSequenceData/austin/MAP_NovaSeq/in_vivo/samtools_depth_counts/star_mapping.%j.out
#SBATCH --error=/DownloadedSequenceData/austin/MAP_NovaSeq/in_vivo/samtools_depth_counts/star_mapping.%j.err

source /opt/ohpc/pub/Software/anaconda3/etc/profile.d/conda.sh
conda activate shapemapper_txsome

ml shapemapper/2.2.0
ml samtools/1.10
ml star/2.7.10a


cd /DownloadedSequenceData/austin/MAP_NovaSeq/in_vivo/txsome_rep_one_rev/fastq_by_target/modified/A

for f in *.fastq; do
	N=$(basename $f .fastq) ;
	STAR --genomeDir /DownloadedSequenceData/austin/MAP_NovaSeq/in_vivo/star_index/ --runThreadN 16 --readFilesIn $f --outFilterMultimapNmax 1 --outFileNamePrefix $N --outSAMtype BAM SortedByCoordinate --outSAMunmapped None --outSAMattributes Standard ;
done

for b in *.bam; do
	samtools index $b ;
done

samtools depth -H -aa *.bam -d 0 -o rep_one_k700e_samtools.depth.counts.txt

cp rep_one_k700e_samtools.depth.counts.txt /DownloadedSequenceData/austin/MAP_NovaSeq/in_vivo/samtools_depth_counts

cd /DownloadedSequenceData/austin/MAP_NovaSeq/in_vivo/txsome_rep_one_rev/fastq_by_target/untreated/A

for f in *.fastq; do
	N=$(basename $f .fastq) ;
	STAR --genomeDir /DownloadedSequenceData/austin/MAP_NovaSeq/in_vivo/star_index/ --runThreadN 16 --readFilesIn $f --outFilterMultimapNmax 1 --outFileNamePrefix $N --outSAMtype BAM SortedByCoordinate --outSAMunmapped None --outSAMattributes Standard ;
done

for b in *.bam; do
	samtools index $b ;
done

samtools depth -H -aa *.bam -d 0 -o rep_one_WT_samtools.depth.counts.txt

cp rep_one_WT_samtools.depth.counts.txt /DownloadedSequenceData/austin/MAP_NovaSeq/in_vivo/samtools_depth_counts

cd /DownloadedSequenceData/austin/MAP_NovaSeq/in_vivo/txsome_rep_two_rev/fastq_by_target/modified/A

for f in *.fastq; do
	N=$(basename $f .fastq) ;
	STAR --genomeDir /DownloadedSequenceData/austin/MAP_NovaSeq/in_vivo/star_index/ --runThreadN 16 --readFilesIn $f --outFilterMultimapNmax 1 --outFileNamePrefix $N --outSAMtype BAM SortedByCoordinate --outSAMunmapped None --outSAMattributes Standard ;
done

for b in *.bam; do
	samtools index $b ;
done

samtools depth -H -aa *.bam -d 0 -o rep_two_k700e_samtools.depth.counts.txt

cp rep_two_k700e_samtools.depth.counts.txt /DownloadedSequenceData/austin/MAP_NovaSeq/in_vivo/samtools_depth_counts

cd /DownloadedSequenceData/austin/MAP_NovaSeq/in_vivo/txsome_rep_two_rev/fastq_by_target/untreated/A

for f in *.fastq; do
	N=$(basename $f .fastq) ;
	STAR --genomeDir /DownloadedSequenceData/austin/MAP_NovaSeq/in_vivo/star_index/ --runThreadN 16 --readFilesIn $f --outFilterMultimapNmax 1 --outFileNamePrefix $N --outSAMtype BAM SortedByCoordinate --outSAMunmapped None --outSAMattributes Standard ;
done

for b in *.bam; do
	samtools index $b ;
done

samtools depth -H -aa *.bam -d 0 -o rep_two_WT_samtools.depth.counts.txt

cp rep_two_WT_samtools.depth.counts.txt /DownloadedSequenceData/austin/MAP_NovaSeq/in_vivo/samtools_depth_counts

cd /DownloadedSequenceData/austin/MAP_NovaSeq/in_vivo/txsome_rep_three_rev/fastq_by_target/modified/A

for f in *.fastq; do
	N=$(basename $f .fastq) ;
	STAR --genomeDir /DownloadedSequenceData/austin/MAP_NovaSeq/in_vivo/star_index/ --runThreadN 16 --readFilesIn $f --outFilterMultimapNmax 1 --outFileNamePrefix $N --outSAMtype BAM SortedByCoordinate --outSAMunmapped None --outSAMattributes Standard ;
done

for b in *.bam; do
	samtools index $b ;
done

samtools depth -H -aa *.bam -d 0 -o rep_three_k700e_samtools.depth.counts.txt

cp rep_three_k700e_samtools.depth.counts.txt /DownloadedSequenceData/austin/MAP_NovaSeq/in_vivo/samtools_depth_counts

cd /DownloadedSequenceData/austin/MAP_NovaSeq/in_vivo/txsome_rep_three_rev/fastq_by_target/untreated/A

for f in *.fastq; do
	N=$(basename $f .fastq) ;
	STAR --genomeDir /DownloadedSequenceData/austin/MAP_NovaSeq/in_vivo/star_index/ --runThreadN 16 --readFilesIn $f --outFilterMultimapNmax 1 --outFileNamePrefix $N --outSAMtype BAM SortedByCoordinate --outSAMunmapped None --outSAMattributes Standard ;
done

for b in *.bam; do
	samtools index $b ;
done

samtools depth -H -aa *.bam -d 0 -o rep_three_WT_samtools.depth.counts.txt

cp rep_three_WT_samtools.depth.counts.txt /DownloadedSequenceData/austin/MAP_NovaSeq/in_vivo/samtools_depth_counts
