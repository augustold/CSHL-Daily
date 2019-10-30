# Transcriptome Assembly Quality Assessment
https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Assembly-Quality-Assessment

```
# Transcriptome results folder
/projects/augustold/CSHL/Trinity/Trinity_denovo/SP80-3280/trinity_out_dir
```

## Assessing the Read Content of the Transcriptome Assembly

First, build a bowtie2 index for the transcriptome:
```
bowtie2-build Trinity.fasta Trinity.fasta
```

## Then perform the alignment to just capture the read alignment statistics
```
# Example using one sample
# Need to write a for loop to parse multiple fastq files



```



