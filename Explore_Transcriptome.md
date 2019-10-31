# Transcriptome Assembly Quality Assessment
https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Assembly-Quality-Assessment

```
# Transcriptome results folder
/projects/augustold/CSHL/Trinity/Trinity_denovo/SP80-3280/trinity_out_dir
```

## Assessing the Read Content of the Transcriptome Assembly

First, build a bowtie2 index for the transcriptome:
```
bowtie2-build --threads 20 Trinity.fasta Trinity.fasta
```

## Then perform the alignment to just capture the read alignment statistics:
```
# Example using one sample
# Need to write a for loop to parse multiple fastq files

bowtie2 -p 40 -q --no-unal -k 20 -x Trinity.fasta \
-1 /projects/augustold/fastq/SPI11_1.fq.gz \
-2 /projects/augustold/fastq/SPI11_2.fq.gz \
2>align_stats.txt| samtools view -@10 -Sb -o bowtie2.bam 
```

## Visualize statistics:
```
cat 2>&1 align_stats.txt
```
