
<pre><code>
for file in Sample*;
do
    mv "$file" "${file#Sample_}"
done
</code></pre>
