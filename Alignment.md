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

ref="dmel-all-chromosome-r6.13.fasta"

bwa index ${ref}
samtools faidx ${ref}
java -d64 -Xmx128g -jar /data/apps/picard-tools/1.87/CreateSequenceDictionary.j\
ar R=${ref} O=dmel-all-chromosome-r6.13.dict
bowtie2-build ${ref}
</code></pre>

