# Pipeline for RNA-seq Analysis

## Using STAR

### STEP 1: FASTQC (Raw FASTQ)
<pre><code>#!/bin/bash
#$ -N fastqc_bef
#$ -q epyc,bio
#$ -pe openmp 8
#$ -R y
#$ -t 1-24

module load fastqc

prefix=`head -n $SGE_TASK_ID fastq.prefixes.txt | tail -n 1` # {ls *1.fastq.gz | sed 's/_R1.fastq.gz//' > fastq.prefixes.txt}

gunzip ${prefix}*
fastqc -t 64 -o ./fastqc_bef --noextract ${prefix}_R1.fastq ${prefix}_R2.fastq #mkdir fastqc_bef
</code></pre>

### STEP 2: rRNA Removal
<pre><code>#!/bin/bash
#$ -N rRNA_removal
#$ -q epyc,bio
#$ -pe openmp 8
#$ -R y
#$ -t 1-24

prefix=`head -n $SGE_TASK_ID fastq.prefixes.txt | tail -n 1`

./sortmerna-2.1b/scripts/merge-paired-reads.sh ${prefix}_R1.fastq ${prefix}_R2.fastq ${prefix}-interleaved.fastq

cd ./sortmerna-2.1b

time ./sortmerna --ref ./rRNA_databases/silva-bac-16s-id90.fasta,./index/silva-bac-16s-db:./rRNA_databases/silva-bac-23s-id98.fasta,./index/silva-bac-23s-db:./rRNA_databases/silva-arc-16s-id95.fasta,./index/silva-arc-16s-db:./rRNA_databases/silva-arc-23s-id98.fasta,./index/silva-arc-23s-db:./rRNA_databases/silva-euk-18s-id95.fasta,./index/silva-euk-18s-db:./rRNA_databases/silva-euk-28s-id98.fasta,./index/silva-euk-28s:./rRNA_databases/rfam-5s-database-id98.fasta,./index/rfam-5s-db:./rRNA_databases/rfam-5.8s-database-id98.fasta,./index/rfam-5.8s-db --reads ../${prefix}-interleaved.fastq --paired_in --fastx --other ../${prefix}-sortmerna --log sample.${prefix}.log -a8 -v --aligned ../${prefix}-rRNA 
 
cd ..

./sortmerna-2.1b/scripts/unmerge-paired-reads.sh ${prefix}-sortmerna.fastq ${prefix}-sortmerna_1.fastq ${prefix}-sortmerna_2.fastq
</code></pre>

### STEP 3 : Trim Adapters (Over Represented Sequences)
<pre><code>#!/bin/bash
#$ -N adapter_trimming
#$ -q epyc,bio
#$ -pe openmp 8
#$ -R y
#$ -t 1-24

module load trimmomatic
module load BBMap

prefix=`head -n $SGE_TASK_ID fastq.prefixes.txt | tail -n 1`

bbmerge.sh in1=${prefix}_R1.fastq in2=${prefix}_R2.fastq outa=${prefix}_adapters.fa

java -jar /data/apps/trimmomatic/0.35/trimmomatic-0.35.jar PE -threads 8 -phred33 ${prefix}-sortmerna_1.fastq ${prefix}-sortmerna_2.fastq ${prefix}-sortmerna-trimmomatic_1.fq ${prefix}-sortmerna-unpaired_1.fq ${prefix}-sortmerna-trimmomatic_2.fq ${prefix}-sortmerna-unpaired_2.fq ILLUMINACLIP:"${prefix}_adapters.fa":2:30:10 SLIDINGWINDOW:5:20 MINLEN:22
</code></pre>

### STEP 4: FASTQC (Processed FASTQ)
<pre><code>#!/bin/bash
#$ -N fastqc_aft
#$ -q epyc,bio
#$ -pe openmp 8
#$ -R y
#$ -t 1-24

module load fastqc

prefix=`head -n $SGE_TASK_ID fastq.prefixes.txt |tail -n 1`

fastqc -t 8 -o ./fastqc_aft --noextract ${prefix}-sortmerna-trimmomatic_1.fq ${prefix}-sortmerna-trimmomatic_2.fq #mkdir fastqc_aft
</code></pre>

### STEP 5a : STAR Alignment (Genome Generation)
<pre><code>#!/bin/bash
#$ -N star_genome_generation
#$ -q epyc,bio
#$ -pe openmp 8

module load STAR/2.5.2a

STAR --runMode genomeGenerate --genomeDir ../star_genome/indices/genome --genomeFastaFiles ../star_genome/genome.fa --runThreadN 8 --sjdbOverhang 99 --sjdbGTFfile ../star_genome/genome.gff3 
</code></pre>

### STEP 5b : STAR Alignment (Alignment)
<pre><code>#!/bin/bash
#$ -N star_jobs
#$ -q epyc,bio
#$ -pe openmp 8
#$ -R y
#$ -t 1-24

module load STAR/2.5.2a

prefix=`head -n $SGE_TASK_ID fastq.prefixes.txt | tail -n 1`

STAR --genomeDir ../star_genome --readFilesIn ../fastq_files/${prefix}-sortmerna-trimmomatic_1.fq ../fastq_files/${prefix}-sortmerna-trimmomatic_2.fq --runThreadN 24 --outFileNamePrefix ${prefix}-STAR  --outSAMtype BAM Unsorted SortedByCoordinate
</code></pre>

## Using HISAT2

### STEP 1: FASTQC (Raw FASTQ)
<pre><code>#!/bin/bash
#$ -N fastqc_bef
#$ -q epyc,bio
#$ -pe openmp 8
#$ -R y
#$ -t 1-24

module load fastqc

prefix=`head -n $SGE_TASK_ID fastq.prefixes.txt | tail -n 1` # {ls *1.fq.gz | sed 's/_R1.fq.gz//' > fastq.prefixes.txt}

fastqc -t 64 -o ./fastqc_bef --noextract ${prefix}_R1.fastq ${prefix}_R2.fastq #mkdir fastqc_bef
</code></pre>

### STEP 2: rRNA Removal
<pre><code>#!/bin/bash
#$ -N rRNA_removal
#$ -q epyc,bio
#$ -pe openmp 8
#$ -R y
#$ -t 1-24

prefix=`head -n $SGE_TASK_ID fastq.prefixes.txt | tail -n 1`

./sortmerna-2.1b/scripts/merge-paired-reads.sh ${prefix}_R1.fastq ${prefix}_R2.fastq ${prefix}-interleaved.fastq

cd ./sortmerna-2.1b

time ./sortmerna --ref ./rRNA_databases/silva-bac-16s-id90.fasta,./index/silva-bac-16s-db:./rRNA_databases/silva-bac-23s-id98.fasta,./index/silva-bac-23s-db:./rRNA_databases/silva-arc-16s-id95.fasta,./index/silva-arc-16s-db:./rRNA_databases/silva-arc-23s-id98.fasta,./index/silva-arc-23s-db:./rRNA_databases/silva-euk-18s-id95.fasta,./index/silva-euk-18s-db:./rRNA_databases/silva-euk-28s-id98.fasta,./index/silva-euk-28s:./rRNA_databases/rfam-5s-database-id98.fasta,./index/rfam-5s-db:./rRNA_databases/rfam-5.8s-database-id98.fasta,./index/rfam-5.8s-db --reads ../${prefix}-interleaved.fastq --paired_in --fastx --other ../${prefix}-sortmerna --log sample.${prefix}.log -a8 -v --aligned ../${prefix}-rRNA 
 
cd ..

./sortmerna-2.1b/scripts/unmerge-paired-reads.sh ${prefix}-sortmerna.fastq ${prefix}-sortmerna_1.fastq ${prefix}-sortmerna_2.fastq
</code></pre>

### STEP 3 : Trim Adapters (Over Represented Sequences)
<pre><code>#!/bin/bash
#$ -N adapter_trimming
#$ -q epyc,bio
#$ -pe openmp 8
#$ -R y
#$ -t 1-24

module load trimmomatic
module load BBMap

prefix=`head -n $SGE_TASK_ID fastq.prefixes.txt | tail -n 1`

bbmerge.sh in1=${prefix}_R1.fastq in2=${prefix}_R2.fastq outa=${prefix}_adapters.fa

java -jar /data/apps/trimmomatic/0.35/trimmomatic-0.35.jar PE -threads 8 -phred33 ${prefix}-sortmerna_1.fastq ${prefix}-sortmerna_2.fastq ${prefix}-sortmerna-trimmomatic_1.fq ${prefix}-sortmerna-unpaired_1.fq ${prefix}-sortmerna-trimmomatic_2.fq ${prefix}-sortmerna-unpaired_2.fq ILLUMINACLIP:"${prefix}_adapters.fa":2:30:10 SLIDINGWINDOW:5:20 MINLEN:22
</code></pre>

### STEP 4: FASTQC (Processed FASTQ)
<pre><code>#!/bin/bash
#$ -N fastqc_aft
#$ -q epyc,bio
#$ -pe openmp 8
#$ -R y
#$ -t 1-24

module load fastqc

prefix=`head -n $SGE_TASK_ID fastq.prefixes.txt |tail -n 1`

fastqc -t 8 -o ./fastqc_aft --noextract ${prefix}-sortmerna-trimmomatic_1.fq ${prefix}-sortmerna-trimmomatic_2.fq #mkdir fastqc_aft
</code></pre>

### STEP 5a : HISAT Alignment (Genome Index)
<pre><code>#!/bin/bash
#$ -N hisat_genome_index
#$ -q epyc,bio
#$ -pe openmp 8

module load hisat2/2.1.0
module load python/2.7.15

hisat2_extract_splice_sites.py -v ./referenceData/annotations/Mus_musculus.GRCm38.80.gtf > ./referenceData/hisat2_index/splice_sites.txt

hisat2_extract_exons.py -v ./referenceData/annotations/Mus_musculus.GRCm38.80.gtf > ./referenceData/hisat2_index/exons.txt

hisat2-build ./referenceData/fasta/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa --ss ./referenceData/hisat2_index/splice_sites.txt --exon ./referenceData/hisat2_index/exons.txt ./referenceData/hisat2_index/GRCm38.hisat2
</code></pre>

### STEP 5b : STAR Alignment (Alignment)
<pre><code>#!/bin/bash
#$ -N star_jobs
#$ -q epyc,bio
#$ -pe openmp 8
#$ -R y
#$ -t 1-24

module load STAR/2.5.2a

prefix=`head -n $SGE_TASK_ID fastq.prefixes.txt | tail -n 1`

STAR --genomeDir ../star_genome --readFilesIn ../fastq_files/${prefix}-sortmerna-trimmomatic_1.fq ../fastq_files/${prefix}-sortmerna-trimmomatic_2.fq --runThreadN 24 --outFileNamePrefix ${prefix}-STAR  --outSAMtype BAM Unsorted SortedByCoordinate
</code></pre>
