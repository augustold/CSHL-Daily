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

## Genome-guided Trinity *De novo* Transcriptome Assembly

###

```
STAR --genomeDir /sonas-hs/ware/hpc_norepl/data/xwang/NAM/Sorg_newWF/ref/star_index_2ndPass --readFilesIn ${ARRAYFILE}_1.fastq.gz ${ARRAYFILE}_2.fastq.gz --outSAMunmapped Within --outFilterType BySJout \
 --outSAMattributes All --chimSegmentMin 12 --chimSegmentReadGapMax 3 --outFilterMultimapNmax 20   --outFilterMismatchNmax 999   \
 --outFilterMismatchNoverReadLmax 0.04   --alignIntronMin 20   --alignIntronMax 1000000   --alignMatesGapMax 1000000   \
 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outSAMtype BAM SortedByCoordinate  --readFilesCommand zcat \
 --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix atlas_${ARRAYFILE}_2ndPass_ \
 --clip5pNbases 0 --seedSearchStartLmax 50 --runThreadN 12 --genomeLoad LoadAndKeep --limitBAMsortRAM 30000000000
 ```
