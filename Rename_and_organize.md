


<pre><code>for file in Sample*;
do
    mv "$file" "${file#Sample_}"
done
</code></pre>

<pre><code> for i in  $(ls *.fq.gz); 
do 
    mv ${i} $(grep $(echo "$i" | awk -F'[__]' '{print $4}') README.ATACseq.txt | awk '{print $2"_"$3"_"$4}')_$i 
done 
</code></pre>
