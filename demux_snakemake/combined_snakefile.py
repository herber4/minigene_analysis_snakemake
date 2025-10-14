# Snakefile

import glob
import os

# Directories
TXT_DIR = "/data2/lackey_lab/DownloadedSequenceData/austin/MAP_NovaSeq/fwd/txt_files"
FASTQ_DIR = "/data2/lackey_lab/DownloadedSequenceData/austin/MAP_NovaSeq/fwd"
OUTPUT = "/scratch2/herber4/demux_out"

REV_DIR = "/data2/lackey_lab/DownloadedSequenceData/austin/MAP_NovaSeq/rev/txt_files"
REV_FASTQ_DIR = "/data2/lackey_lab/DownloadedSequenceData/austin/MAP_NovaSeq/rev"

FASTQ_FILE = glob_wildcards("/data2/lackey_lab/DownloadedSequenceData/austin/MAP_NovaSeq/fwd/{fname}.fastq.gz").fname
TXT_FILES = glob_wildcards("/data2/lackey_lab/DownloadedSequenceData/austin/MAP_NovaSeq/fwd/txt_files/{fname}.txt").fname
rev_fq = glob_wildcards("/data2/lackey_lab/DownloadedSequenceData/austin/MAP_NovaSeq/rev/{fname}.fastq.gz").fname
rev_txt = glob_wildcards("/data2/lackey_lab/DownloadedSequenceData/austin/MAP_NovaSeq/rev/txt_files/{fname}.txt").fname


rule all:
    input:
        expand("/scratch2/herber4/demux_out/{iBARCODE}.{iFASTQ}.fastq.gz", iBARCODE = TXT_FILES, iFASTQ = FASTQ_FILE),
        expand("/scratch2/herber4/demux_out/{REVTXT}.{REVFQ}.fastq.gz", REVTXT = rev_txt, REVFQ = rev_fq)

rule grep_seq:
    input:
        txt = os.path.join(TXT_DIR, "{iBARCODE}.txt"),
        fqs = os.path.join(FASTQ_DIR, "{iFASTQ}.fastq.gz")
    output:
        joined = os.path.join(OUTPUT,"{iBARCODE}.{iFASTQ}.fastq.gz")
    resources: cpus=4, mem_mb=16000, time_min=1440
    shell: "/data2/lackey_lab/DownloadedSequenceData/austin/bin/seqkit grep -j 4 -s -m 1 -f {input.txt} {input.fqs} -o {output.joined}"

rule grep_seq_rev:
    input:
        txt = os.path.join(REV_DIR, "{REVTXT}.txt"),
        fqs = os.path.join(REV_FASTQ_DIR, "{REVFQ}.fastq.gz")
    output:
        joined = os.path.join(OUTPUT,"{REVTXT}.{REVFQ}.fastq.gz")
    resources: cpus=4, mem_mb=16000, time_min=1440
    shell: "/data2/lackey_lab/DownloadedSequenceData/austin/bin/seqkit grep -j 4 -s -m 1 -f {input.txt} {input.fqs} -o {output.joined}"
