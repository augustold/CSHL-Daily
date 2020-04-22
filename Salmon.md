# On Helix!

conda activate salmon

```
cd /projects/augustold/CSHL/salmon
conda create -n salmon salmon
```

## building an index
```
conda activate salmon # version salmon-1.2.0

salmon index -t /projects/augustold/CSHL/Saccharum_genome_refs/SP803280/sc.mlc.cns.sgl.utg.scga7.CDS.codingseq.fasta -i scga7_index

mkdir quants 
```

## Quantify!
```
conda activate salmon # version salmon-1.2.0

mkdir quants

cd /projects/augustold/fastq/
for fn in $(ls SPL1* | sed s/_[12].fq.gz// | sort -u);
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i /projects/augustold/salmon_energy/scga7.and.Sspon_index -l A \
         -1 ${samp}_1.fq.gz \
         -2 ${samp}_2.fq.gz \
         -p 60 --validateMappings -o /projects/augustold/salmon_energy/quants/${samp}_quant
done 

```
