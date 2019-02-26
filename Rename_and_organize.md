# Renaming and Organizing Files
Before making any changes to the actual file names, I just copied the files to a new directory as follows:
<pre><code>cp .../Bioinformatics_Course/ATACseq/* .../Bioinformatics_Course/renameATACseq
cp .../Bioinformatics_Course/DNAseq/* .../Bioinformatics_Course/renameDNAseq
cp **/**/*.gz ../../renamedRNAseq/
</code></pre>
## ATACseq

<pre><code>for file in Sample*;
do
    mv "$file" "${file#Sample_}"
done
</code></pre>

<pre><code> for i in  $(ls *.gz); 
do 
    mv ${i} $(grep $(echo "$i" | awk -F'[__]' '{print $4}') README.ATACseq.txt | awk '{print $2"_"$3"_"$4}')_$i 
done 
</code></pre>

## DNAseq

<pre><code>for i in ADL06*; 
do     
    mv "$i" "A4_${i}"
done
</code></pre>
<pre><code>for i in ADL09*; 
do     
    mv "$i" "A5_${i}"
done
</code></pre>
<pre><code>for i in ADL10*; 
do     
    mv "$i" "A6_${i}"
done
</code></pre>
<pre><code>for i in ADL14*; 
do     
    mv "$i" "A7_${i}"
done
</code></pre>

## RNAseq
<pre><code>for i in  $(ls *.gz);  
do      
    mv ${i} $(grep $(echo "$i" | awk -F'[__]' '{print $1}') RNAseq384_SampleCoding.txt | awk '{print $12"_"$8}')_$i
done
</code></pre>