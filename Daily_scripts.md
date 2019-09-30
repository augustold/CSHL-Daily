# Sep 30 2019

## Running Fastqc

```
module load FastQC/0.11.8-Java-1.8
```

```
cat > Fastqc.sh
```

Add the following to the file
```
for f in /sonas-hs/ware/hpc/home/diniz/rnaseq/*.fq.gz; do fastqc --outdir  /sonas-hs/ware/hpc/home/diniz/Fastqc_Out -t 8 $f  ; done
```

Running on CSHL cluster
```
qsub -cwd -pe threads 8 -l m_mem_free=5G Fastqc.sh
```

## Running Trinity

```
mkdir Trininty
cd Trinity/
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Trinity/2.8.4
```

Creating script
```
cat > Trinity.sh
```

Add the following to the file


```
# Combining all reads data sets

cat /sonas-hs/ware/hpc/home/diniz/rnaseq/*_1.fq.gz > ALL.LEFT.fq.gz
cat /sonas-hs/ware/hpc/home/diniz/rnaseq/*_2.fq.gz > ALL.RIGHT.fq.gz

#Running Trinity for de novo assembly
Trinity --seqType fq --SS_lib_type RF --left ALL.LEFT.fq --right ALL.RIGHT.fq --CPU 8 --M 5G
```

Running on CSHL cluster
```
qsub -cwd -pe threads 8 -l m_mem_free=5G Trinity.sh
```

# Sep 26 2019

## Running STAR
### Using the Sorghum-genome-and-annotation-projects from Xiaofei

###Star genome index

diniz@bnbdev1:~/Saccharum_genome_refs/SP803280$
```
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load STAR/2.7.0d
```
```
cat > STARmap.sh
```
Add the following to the file
```
cd /sonas-hs/ware/hpc/home/diniz/Saccharum_genome_refs/SP803280

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /sonas-hs/ware/hpc/home/diniz/Saccharum_genome_refs/SP803280/star_index --genomeFastaFiles /sonas-hs/ware/hpc/home/diniz/Saccharum_genome_refs/SP803280/sc.mlc.cns.sgl.utg.scga7.importdb.fa --sjdbGTFfile /sonas-hs/ware/hpc/home/diniz/Saccharum_genome_refs/SP803280/sc.mlc.cns.sgl.utg_scga7.sort.gff3  --genomeChrBinNbits 13
```
###  Notes
--genomeChrBinNbits = log2(GenomeLength/NumberOfReferences)
log2((4.26*1000000000)/450609) = 13.20669

Running on CSHL cluster

```
qsub -cwd -pe threads 8 -l m_mem_free=5G STARmap.sh
```


# Sep 26 2019

## Prepare the genomic region for alignment by GMAP

```
gmap_build -d genome -D . -k 13 genome.fa
```

## sugarcane orfeome vs SP80 genome

```
/projects/augustold/SP80_genome

gmap -D . -d scga7 -f samse -t 20 -n 1 sugarcane.fulllength.analysis.all.fasta > mapping.sam 2> mapping.sam.log
```

## sorghum full-length mRNA vs SP80 genome

```
/projects/augustold/CSHL

gmap -D . -d scga7 -f samse -t 20 -n 1 sorghum.combined11.collapsed.rep.fa > mapping.sam 2> mapping.sam.log
```

## count how many transcripts didn’t map to the reference

```
grep ‘No paths found’ mapping.sam.log | wc -l
```
## Result:

Regarding truncated gene models in the sugarcane assembly, I have mapped the full-length transcripts from sorghum and sugarcane (obtained from 454 sequencing in 2014).
Sorghum has 95,380 full-length transcripts and 2,738 (2.87%) didn't map to the sugarcane reference. Sugarcane has 40,587 predicted transcripts and 809 (1.99%) didn't map to the sugarcane reference. Maybe we can say that we have only a small fraction of truncated genes gene models in our reference (~2%).
