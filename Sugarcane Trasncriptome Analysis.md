## Augusto L. Diniz / USP - CSHL

# Transcriptome Assembly

## *De novo* Transcriptome

Building the transcriptome using Trinity (v2.8.5) in a *de novo* approach for SP80-3280, *S. spontaneum* and *S. officinarun*. For this, I will use samples collected by Maryke from leaf +1 (L1), upper internode (I1), young internode (I5) and mature internode (I9) for each genotype. Each tissue has three biological replicates. 

**Table 1**. Number of samples (RNA-seq libraries)

| Genotype         | L1 | I1 | I5 | I9 |  Total |
|------------------|----|----|----|----|--------|
| SP80-3280        | 3  | 3  | 3  | 3  |  12    |
| *S. spontaneum*  | 3  | 3  | 3  | 3  |  12    |
| *S. officinarum* | 3  | 3  | 3  | 3  |  12    |
| **Total**        |    |    |    |    |**36**  |

```
# Result files in
augustold@helix:/projects/augustold/CSHL/Trinity/Trinity_denovo/
officinarum  SP80-3280  spontaneum

# FastQ file links in
augustold@helix:/projects/augustold/fastq
```

## SP80-3280

```
cd /projects/augustold/CSHL/Trinity/Trinity_denovo/SP80-3280/
# create the samples.txt file
# cat > samples.txt
# add only sample ids from SP80-3280 (SP) files

# run Trinity
screen
/projects/augustold/CSHL/Trinity/trinityrnaseq-Trinity-v2.8.5/Trinity --seqType fq --max_memory 50G --CPU 10 --samples_file samples.txt --SS_lib_type RF 
```

## *S. spontaneum*

```
cd /projects/augustold/CSHL/Trinity/Trinity_denovo/spontaneum/
# create the samples.txt file
# cat > samples.txt
# add only sample ids from S. spontaneum (IN) files

# run Trinity
screen
/projects/augustold/CSHL/Trinity/trinityrnaseq-Trinity-v2.8.5/Trinity --seqType fq --max_memory 50G --CPU 10 --samples_file samples.txt --SS_lib_type RF 
```

## *S. officinarum*

```
cd /projects/augustold/CSHL/Trinity/Trinity_denovo/officinarum/
# create the samples.txt file
# cat > samples.txt
# add only sample ids from S. spontaneum (CB) files

# run Trinity
screen
/projects/augustold/CSHL/Trinity/trinityrnaseq-Trinity-v2.8.5/Trinity --seqType fq --max_memory 50G --CPU 10 --samples_file samples.txt --SS_lib_type RF 
```

## Genome-guided Trinity Transcriptome Assembly

Users must provide read alignments to Trinity as a coordinate-sorted bam file. In order to do that, I'll use STAR.

### STAR genome index
```
# SP80-3280
cd /projects/augustold/CSHL/Saccharum_genome_refs/SP803280/
mkdir star_index
screen

/projects/augustold/CSHL/STAR-2.7.2c/bin/Linux_x86_64/STAR --runThreadN 20 --runMode genomeGenerate --genomeDir /projects/augustold/CSHL/Saccharum_genome_refs/SP803280/star_index --genomeFastaFiles /projects/augustold/CSHL/Saccharum_genome_refs/SP803280/sc.mlc.cns.sgl.utg.scga7.importdb.fa --sjdbGTFfile /projects/augustold/CSHL/Saccharum_genome_refs/SP803280/sc.mlc.cns.utg.sgl_v2_without_mobile-element_prtn_coding_prim_transcrpt.gff --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 150 --limitGenomeGenerateRAM 350000000000

# specify '--limitGenomeGenerateRAM' as 350000000000 due to error message:
# EXITING because of FATAL PARAMETER ERROR: limitGenomeGenerateRAM=200000000000is too small for your genome
# SOLUTION: please specify --limitGenomeGenerateRAM not less than 315007636181 and make that much RAM available

#########################################################

# S. spontaneum
cd projects/augustold/CSHL/Saccharum_genome_refs/Sspon/
mkdir star_index
screen

/projects/augustold/CSHL/STAR-2.7.2c/bin/Linux_x86_64/STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /projects/augustold/CSHL/Saccharum_genome_refs/Sspon/star_index --genomeFastaFiles /projects/augustold/CSHL/Saccharum_genome_refs/Sspon/Sspon.HiC_chr_asm.fasta --sjdbGTFfile /projects/augustold/CSHL/Saccharum_genome_refs/Sspon/Sspon.v20190103.gff3 --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 150
```

### Align reads

```
# Result files in
augustold@helix:/projects/augustold/CSHL/Trinity/Trinity_guided/
officinarum  SP80-3280  spontaneum
```

```
# SP80-3280
cd /projects/augustold/CSHL/Trinity/Trinity_guided/SP80-3280/
mkdir aligned

for i in $(ls fastq | sed s/_[12].fq.gz// | sort -u)
do
    /projects/augustold/CSHL/STAR-2.7.2c/bin/Linux_x86_64/STAR --genomeDir /projects/augustold/CSHL/Saccharum_genome_refs/SP803280/star_index --readFilesIn fastq/${i}_1.fq.gz,fastq/${i}_2.fq.gz --runThreadN 10 --outFileNamePrefix aligned/${i}. --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outSAMunmapped Within --outFilterType BySJout --outSAMattributes All --chimSegmentMin 12 --chimSegmentReadGapMax 3 --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterIntronMotifs RemoveNoncanonical --clip5pNbases 0 --seedSearchStartLmax 50 --genomeLoad LoadAndKeep --limitBAMsortRAM 35000000000
done
```

```
# S. spontaneum
cd /projects/augustold/CSHL/Trinity/Trinity_guided/spontaneum/
mkdir aligned

for i in $(ls fastq | sed s/_[12].fq.gz// | sort -u)
do
    /projects/augustold/CSHL/STAR-2.7.2c/bin/Linux_x86_64/STAR --genomeDir /projects/augustold/CSHL/Saccharum_genome_refs/Sspon/star_index --readFilesIn fastq/${i}_1.fq.gz,fastq/${i}_2.fq.gz --runThreadN 10 --outFileNamePrefix aligned/${i}. --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outSAMunmapped Within --outFilterType BySJout --outSAMattributes All --chimSegmentMin 12 --chimSegmentReadGapMax 3 --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterIntronMotifs RemoveNoncanonical --clip5pNbases 0 --seedSearchStartLmax 50 --genomeLoad LoadAndKeep --limitBAMsortRAM 35000000000
done
```

### Merge BAM files and build Genome-guided Trinity

```
# SP80-3280
screen
cd /projects/augustold/CSHL/Trinity/Trinity_guided/SP80-3280/

samtools merge atlas_merged.bam aligned/*.bam
samtools index atlas_merged.bam

# New run 30 Oct 2019
/projects/augustold/CSHL/Trinity/trinityrnaseq-Trinity-v2.8.5/Trinity --genome_guided_bam atlas_merged.bam --genome_guided_max_intron 10000 --max_memory 100G --CPU 40 --SS_lib_type RF --output /projects/augustold/CSHL/Trinity/Trinity_guided/SP80-3280/trinity_out_dir

# S. spontaneum
screen
cd /projects/augustold/CSHL/Trinity/Trinity_guided/spontaneum/

samtools merge atlas_merged.bam aligned/*.bam
samtools index atlas_merged.bam
/projects/augustold/CSHL/Trinity/trinityrnaseq-Trinity-v2.8.5/Trinity --genome_guided_bam atlas_merged.bam --genome_guided_max_intron 10000 --max_memory 50G --CPU 10
```

### CLASS2: takes too long! try psyCLASS

```
#SP80
cd /projects/augustold/CSHL/Class/
screen
perl /projects/augustold/CSHL/Class/CLASS-2.1.7/run_class.pl -a /projects/augustold/CSHL/Trinity/Trinity_guided/SP80-3280/atlas_merged.bam -o /projects/augustold/CSHL/Class/SP80-3280_atlas_cl.gtf -p 16 --wd /projects/augustold/CSHL/Class

#SP80
cd /projects/augustold/CSHL/Class/
screen
perl /projects/augustold/CSHL/Class/CLASS-2.1.7/run_class.pl -a /projects/augustold/CSHL/Trinity/Trinity_guided/spontaneum/atlas_merged.bam -o /projects/augustold/CSHL/Class/spontaneum_atlas_cl.gtf -p 16 --wd /projects/augustold/CSHL/Class
```

### psyCLASS

```
cd /projects/augustold/CSHL/psiclass/

./psiclass \
-b /projects/augustold/CSHL/Trinity/Trinity_guided/SP80-3280/atlas_merged.bam \
-p 20
```

### StringTie

```
cd /projects/augustold/CSHL/StringTie/
screen
/projects/augustold/CSHL/StringTie/stringtie/stringtie \
/projects/augustold/CSHL/Trinity/Trinity_guided/SP80-3280/atlas_merged.bam \
-G /projects/augustold/CSHL/Saccharum_genome_refs/SP803280/sc.mlc.cns.sgl.utg_scga7.sort.gff3 \
-A /projects/augustold/CSHL/StringTie/gene_abund.tab \
-C /projects/augustold/CSHL/StringTie/cov_refs.gtf --rf \
-B -e -p 16 \
-o /projects/augustold/CSHL/StringTie/SP80-3280_atlas_ST.gtf
```

### Cufflinks - to run!

Run on BNB!

```
cat > cufflinks.sh

#Paste this
cd /sonas-hs/ware/hpc/home/diniz/cufflinks/SP80-3280/

/sonas-hs/ware/hpc/home/bwang/software/cufflinks-2.1.1.Linux_x86_64/cufflinks \
-p 16 \
-o /sonas-hs/ware/hpc/home/diniz/cufflinks/SP80-3280/ \
-g /sonas-hs/ware/hpc/home/diniz/Saccharum_genome_refs/SP803280/sc.mlc.cns.sgl.utg_scga7.sort.gff3 \
--library-type fr-firststrand \
-u atlas_merged.bam
```

```
qsub -cwd -pe threads 16 -l m_mem_free=5G cufflinks.sh
```

Run on Helix

```
cd /projects/augustold/CSHL/Cufflinks/

cufflinks \
-p 20 \
-o /projects/augustold/CSHL/Cufflinks/ \
-g /projects/augustold/CSHL/Saccharum_genome_refs/SP803280/sc.mlc.cns.sgl.utg_scga7.sort.gff3 \
--library-type fr-firststrand \
-u /projects/augustold/CSHL/Trinity/Trinity_guided/SP80-3280/atlas_merged.bam
```
