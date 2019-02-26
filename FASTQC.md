# FASTQC for different DATA

## ATACseq
> Check if the format of the file and change it to mordern format if the format is old.

<pre><code>module load BBMap
testformat.sh in=A4_ED_2_ATGCATG-CTCCTTAC_4R009_L1_P004_R1.fq.gz
ls *1.fq.gz | sed 's/1.fq.gz//' > ATACseq.prefixes.txt
</code></pre>

> As the format is sangers we can directly do the fastqc and the script is as follows 
{(wc -l ATACseq.prefixes.txt) will give you how many parallel jobs you need to run}:

<pre><code>#!/bin/bash

</code></pre>

