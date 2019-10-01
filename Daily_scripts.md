# Oct 1 2019

## Running GMAP: 80% coverage and identity

```
ssh diniz@bnbdev1
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load GMAP-GSNAP/2019-03-15
module load SAMtools/1.9
```

## Prepare the genomic region for alignment by GMAP

```
# cat > GMAP_build.sh

cd /sonas-hs/ware/hpc/home/diniz/Saccharum_genome_refs/SP803280/
gmap_build -d gmap_index -D . -k 13 sc.mlc.cns.sgl.utg.scga7.importdb.fa
```

Running on CSHL cluster
```
qsub -cwd -pe threads 8 -l m_mem_free=5G GMAP_build.sh
```

## Sugarcane orfeome *vs* SP80 genome: 80% cov and ident

```
# cat > GMAP_run.sh

gmap -D /sonas-hs/ware/hpc/home/diniz/Saccharum_genome_refs/SP803280/ -d gmap_index -f samse -t 20 -n 1 --min-trimmed-coverage=0.80 --min-identity=0.80 sugarcane.fulllength.analysis.all.fasta > mapping.sam 2> mapping.sam.log
```

Running on CSHL cluster
```
qsub -cwd -pe threads 20 -l m_mem_free=5G GMAP_run.sh
```

Out of the 40,587 predicted transcripts for sugarcane, 34547 (85.11%) map to the sugarcane reference considering a threshold of 80% coverage and identidy.
On the other hand, and considering the same threshold, 67707 (71%) out 95,380 sorghum full-length transcripts map to the sugarcane reference. 

## Sorghum orfeome *vs* SP80 genome: 80% cov and ident

```
# cat > GMAP_sorgh_run.sh

gmap -D /sonas-hs/ware/hpc/home/diniz/Saccharum_genome_refs/SP803280/ -d gmap_index -f samse -t 20 -n 1 --min-trimmed-coverage=0.80 --min-identity=0.80 sorghum.combined11.collapsed.rep.fa > mapping_sorgh.sam 2> mapping_sorgh.sam.log
```

Running on CSHL cluster
```
qsub -cwd -pe threads 20 -l m_mem_free=5G GMAP_sorgh_run.sh
```

## Count how many transcripts didn’t map to the reference

```
grep 'No paths found' mapping.sam.log | wc -l

#Alternative to grep
awk '/No paths found/ { count++ } END { print count }' mapping.sam.log
```

# Sep 30 2019

## Running Fastqc

```
ssh diniz@bnbdev1
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
ssh diniz@bnbdev1
mkdir Trininty
cd Trinity/
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Trinity/2.8.4
module load Jellyfish/2.2.10
module load SAMtools/1.9
module load Salmon
```

Creating file samples.txt: https://docs.google.com/spreadsheets/d/1gQn80BJNgY8LMA3DNF5LuaA6OgwkgGWJ8w-ijBq6n6w/edit#gid=901858591

Creating script
```
cat > Trinity.sh
```
Add the following to the file
```
Trinity --seqType fq --samples_file samples.txt --CPU 20 --max_memory 100G
```

Running on CSHL cluster
```
qsub -cwd -pe threads 20 -l m_mem_free=5G Trinity.sh
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

## sugarcane orfeome *vs* SP80 genome

```
/projects/augustold/SP80_genome

gmap -D . -d scga7 -f samse -t 20 -n 1 sugarcane.fulllength.analysis.all.fasta > mapping.sam 2> mapping.sam.log
```

## sorghum full-length mRNA *vs* SP80 genome

```
/projects/augustold/CSHL

gmap -D . -d scga7 -f samse -t 20 -n 1 sorghum.combined11.collapsed.rep.fa > mapping.sam 2> mapping.sam.log

# 80% coverage and identity
gmap -D . -d scga7 -f samse -t 20 -n 1 --min-trimmed-coverage=0.80 --min-identity=0.80 sorghum.combined11.collapsed.rep.fa > mapping_80.sam 2> mapping_80.sam.log
```

## count how many transcripts didn’t map to the reference

```
grep ‘No paths found’ mapping.sam.log | wc -l

#If grep doesn't work
awk '/No paths found/ { count++ } END { print count }' mapping.sam.log
```
## Result:

Regarding truncated gene models in the sugarcane assembly, I have mapped the full-length transcripts from sorghum and sugarcane (obtained from 454 sequencing in 2014).
Sorghum has 95,380 full-length transcripts and 2,738 (2.87%) didn't map to the sugarcane reference. Sugarcane has 40,587 predicted transcripts and 809 (1.99%) didn't map to the sugarcane reference. Maybe we can say that we have only a small fraction of truncated genes gene models in our reference (~2%).
