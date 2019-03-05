<pre><code>#!/bin/bash
#$ -N coveragebeds_all
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
module load jje/ken/2014.02.19

ref="../ref/dmel-all-chromosome-r6.13.fasta"
prefix=`head -n $SGE_TASK_ID ATACseq.prefixes.txt | tail -n 1`
Nreads=`samtools view -c -F 4 ./alignATAC/${prefix}.RG.bam`
Scale=`echo "1.0/($Nreads/1000000)" | bc -l`

samtools view -b ./alignATAC/${prefix}.RG.bam | genomeCoverageBed -ibam -g ${ref} -bg -scale ${Scale} > ./alignATAC/${prefix}.coverage

bedGraphToBigWig ./alignATAC/${prefix}.coverage ${ref}.fai ./alignATAC/${prefix}.bw
</code></pre>
