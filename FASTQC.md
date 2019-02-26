# FASTQC for different DATA

## ATACseq
> Check if the format of the file and change it to mordern format if the format is old.

<pre><code>module load BBMap
testformat.sh in=A4_ED_2_ATGCATG-CTCCTTAC_4R009_L1_P004_R1.fq.gz
ls *1.fq.gz | sed 's/1.fq.gz//' > ATACseq.prefixes.txt
cd ..
mkdir fastqc_before
cd ./fastqc_before
mv ../renamedATACseq/ATACseq.prefixes.txt ./
</code></pre>

> As the format is sangers we can directly do the fastqc and the script is as follows 
{(wc -l ATACseq.prefixes.txt) will give you how many parallel jobs you need to run}:

<pre><code>#!/bin/bash
#$ -N fastqc_all_ATAC
#$ -q epyc,bio
#$ -pe openmp 8
#$ -R y
#$ -t 1-24

module load fastqc/0.11.2
prefix=`head - n $SGE_TASK_ID ATACSeq.prefixes.txt | tail -n 1`
fastqc --noextract -o ./ ../renamedATACseq/${prefix}1.fq.gz ../renamedATACseq/${prefix}2.fq.gz
</code></pre>

## RNAseq
> Check if the format of the file and change it to mordern format if the format is old.

<pre><code>module load BBMap
testformat.sh in=100_CGTTCTT_L003_R1_001.fastq.gz
ls *1_001.fastq.gz | sed 's/1_001.fastq.gz//' > RNAseq.prefixes.txt
cd ../fastqc_before
mv ../renamedRNAseq/RNAseq.prefixes.txt ./
</code></pre>

> As the format is sangers we can directly do the fastqc and the script is as follows 
{(wc -l RNAseq.prefixes.txt) will give you how many parallel jobs you need to run}:

<pre><code>#!/bin/bash
#$ -N fastqc_all_RNA
#$ -q epyc,bio
#$ -pe openmp 8
#$ -R y
#$ -t 1-384

module load fastqc/0.11.2
prefix=`head - n $SGE_TASK_ID RNASeq.prefixes.txt | tail -n 1`
fastqc --noextract -o ./ ../renamedRNAseq/${prefix}1_001.fastq.gz ../renamedRNAseq/${prefix}2_001.fastq.gz
</code></pre>

## DNAseq
> Check if the format of the file and change it to mordern format if the format is old.

<pre><code>module load BBMap
testformat.sh in=A4_ADL06_1_1.fq.gz
ls *_1.fq.gz | sed 's/_1.fq.gz//' > DNAseq.prefixes.txt
cd ../fastqc_before
mv ../renamedDNAseq/DNAseq.prefixes.txt ./
</code></pre>

> As the format is Illumina we have to converted the format to do the fastqc and the script is as follows 
{(wc -l DNAseq.prefixes.txt) will give you how many parallel jobs you need to run}, We can use bwa -I parameter
for illumina reads:
 
<pre><code>./seqtk/seqtk seq -Q64 -V  A4_ADL06_1_1.fq.gz | gzip -c > A4_ADL06_1_converted_1.fq.gz #do for all reads
</code></pre>
<pre><code>#!/bin/bash
#$ -N fastqc_all_DNA
#$ -q epyc,bio
#$ -pe openmp 8
#$ -R y
#$ -t 1-12

module load fastqc/0.11.2
prefix=`head - n $SGE_TASK_ID DNASeq.prefixes.txt | tail -n 1`
fastqc --noextract -o ./ ../renamedDNAseq/${prefix}_1.fq.gz ../renamedDNAseq/${prefix}_2.fq.gz
</code></pre>
