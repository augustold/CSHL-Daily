# On Helix!

conda activate salmon

```
cd /projects/augustold/CSHL/salmon
```

## building an index
```
salmon index -t /projects/augustold/CSHL/Saccharum_genome_refs/SP803280/sc.mlc.cns.sgl.utg.scga7.CDS.codingseq.fasta -i scga7_index

mkdir quants 
```

## Quantify!
```
# vi salmon_quant.sh

#!/bin/bash
for i in $(ls /projects/augustold/fastq/SP* | sed s/_[12].fq.gz// | sort -u);
do
samp=`basename ${i}`
echo "Processing sample ${samp}"
salmon quant -i scga7_index -l A \
         -1 ${i}_1.fq.gz \
         -2 ${i}_2.fq.gz \
         -p 30 --validateMappings -o quants/${samp}_quant
done
```
