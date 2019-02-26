# Indexing Reference Genome
> Indexing the reference genome for all the modules you would require later on to align the data.

<pre><code>#!/bin/bash
#$ -N ref_index
#$ -q epyc,bio
#$ -pe openmp 8
#$ -R y
#$ -t 1

module load bwa/0.7.8
module load samtools/1.3
module load bcftools/1.3
module load enthought_python/7.3.2
module load gatk/2.4-7
module load picard-tools/1.87
module load java/1.7
module load tophat/2.1.0
module load bowtie2/2.2.7

ref="dmel-all-chromosome-r6.13.fasta"

bwa index ${ref}
samtools faidx ${ref}
java -d64 -Xmx128g -jar /data/apps/picard-tools/1.87/CreateSequenceDictionary.jar R=${ref} O=dmel-all-chromosome-r6.13.dict
bowtie2-build ${ref} ${ref}.out
</code></pre>

# Alignment

## ATACseq
<pre><code>#!/bin/bash
#$ -N align_ATAC
#$ -q epyc,bio
#$ -pe openmp 8
#$ -R y
#$ -t 1-24

module load bwa/0.7.8
module load samtools/1.3
module load bcftools/1.3
module load enthought_python/7.3.2
module load gatk/2.4-7
module load picard-tools/1.87
module load java/1.7
module load tophat/2.1.0
module load bowtie2/2.2.7

ref="../ref/dmel-all-chromosome-r6.13.fasta"

prefix=`head -n $SGE_TASK_ID ATACseq.prefixes.txt | tail -n 1`

dict="../ref/dmel-all-chromosome-r6.13.dict"

bwa mem -t 8 -M ${ref} ${prefix}1.fq.gz ${prefix}2.fq.gz | samtools view -bS - > ./alignATAC/${prefix}.bam

samtools sort ./alignATAC/${prefix}.bam -o ./alignATAC/${prefix}.sort.bam

java -Xmx20g -jar /data/apps/picard-tools/1.87/AddOrReplaceReadGroups.jar I=./alignATAC/${prefix}.sort.bam O=./alignATAC/${prefix}.RG.bam SORT_ORDER=coordinate RGPL=Sangers RGPU=D109LACXX RGLB=Lib1 RGID=${prefix} RGSM=${prefix} VALIDATION_STRINGENCY=LENIENT

samtools index ./alignATAC/${prefix}.RG.bam
</code></pre>

## RNAseq
<pre><code>#!/bin/bash
#$ -N align_RNA
#$ -q epyc,bio
#$ -pe openmp 8
#$ -R y
#$ -t 1

module load bwa/0.7.8
module load samtools/1.3
module load bcftools/1.3
module load enthought_python/7.3.2
module load gatk/2.4-7
module load picard-tools/1.87
module load java/1.7
module load tophat/2.1.0
module load bowtie2/2.2.7

giff="../ref/dmel-all-r6.13.gtf"
ref="../ref/dmel-all-chromosome-r6.13.fasta.out"

prefix=`head -n $SGE_TASK_ID RNAseq.prefixes.txt | tail -n 1`

mkdir ${prefix}_alignRNA

tophat -p 8 -G ${giff} -o ./${prefix}_alignRNA ${ref} ${prefix}1_001.fastq.gz ${prefix}2_001.fastq.gz

samtools sort ./${prefix}_alignRNA/accepted_hits.bam -o ./${prefix}_alignRNA/${prefix}_accepted_hits.sort.bam

samtools index ./${prefix}_alignRNA/${prefix}_accepted_hits.sort.bam
</code></pre>

## DNAseq
<pre><code>#!/bin/bash
#$ -N align_DNA
#$ -q epyc,bio
#$ -pe openmp 8
#$ -R y
#$ -t 1-12

module load bwa/0.7.8
module load samtools/1.3
module load bcftools/1.3
module load enthought_python/7.3.2
module load gatk/2.4-7
module load picard-tools/1.87
module load java/1.7
module load tophat/2.1.0
module load bowtie2/2.2.7

ref="../ref/dmel-all-chromosome-r6.13.fasta"

prefix=`head -n $SGE_TASK_ID DNAseq.prefixes.txt | tail -n 1`

dict="../ref/dmel-all-chromosome-r6.13.dict"

bwa mem -t 8 -M ${ref} ${prefix}_1.fq.gz ${prefix}_2.fq.gz | samtools view -bS - > ./alignDNA/${prefix}.bam

samtools sort ./alignDNA/${prefix}.bam -o ./alignDNA/${prefix}.sort.bam

java -Xmx20g -jar /data/apps/picard-tools/1.87/AddOrReplaceReadGroups.jar I=./alignDNA/${prefix}.sort.bam O=./alignDNA/${prefix}.RG.bam SORT_ORDER=coordinate RGPL=illumina RGPU=D109LACXX RGLB=Lib1 RGID=${prefix} RGSM=${prefix} VALIDATION_STRINGENCY=LENIENT

samtools index ./alignDNA/${prefix}.RG.bam
</code></pre>
