## Augusto L. Diniz / USP - CSHL

# Transcriptome Assembly

## *De novo* Transcriptome

Building the transcriptome using Trinity (v2.8.5) in a *de novo* approach for SP80-3280, *S. spontaneum* and *S. officinarun*. For this, I will use samples collected by Mryke from leaf +1 (L1), upper internode (I1), young internode (I5) and mature internode (I9) for each genotype. Each tissue has three biological replicates. 

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

### Merge BAM files

```
# SP80-3280
cd /projects/augustold/CSHL/Trinity/Trinity_guided/SP80-3280/

samtools merge atlas_merged.bam aligned/*.bam
samtools index atlas_merged.bam

# S. spontaneum
cd /projects/augustold/CSHL/Trinity/Trinity_guided/spontaneum/

samtools merge atlas_merged.bam aligned/*.bam
samtools index atlas_merged.bam
```

### Genome-guided Trinity

```
# SP80-3280
cd /projects/augustold/CSHL/Trinity/Trinity_guided/SP80-3280/

/projects/augustold/CSHL/Trinity/trinityrnaseq-Trinity-v2.8.5/TrinityTrinity --genome_guided_bam atlas_merged.bam --genome_guided_max_intron 10000 --max_memory 50G --CPU 10 --SS_lib_type RF

# S. spontaneum
cd /projects/augustold/CSHL/Trinity/Trinity_guided/spontaneum/

/projects/augustold/CSHL/Trinity/trinityrnaseq-Trinity-v2.8.5/TrinityTrinity --genome_guided_bam atlas_merged.bam --genome_guided_max_intron 10000 --max_memory 50G --CPU 10 --SS_lib_type RF
```
