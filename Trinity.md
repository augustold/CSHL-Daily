## Oct 2 2019
## Augusto L. Diniz
## USP - CSHL

# Sugarcane SP80-3280 Transcriptome Assembly

## Trinity - *de novo*

```
# FastQ file links in
augustold@helix:/projects/augustold/fastq

# Analysis will be parsed in
augustold@helix:/projects/augustold/CSHL/Trinity/Trinity_denovo
```

### Running trinity
```
cd /projects/augustold/CSHL/Trinity/Trinity_denovo
/projects/augustold/CSHL/Trinity/trinityrnaseq-Trinity-v2.8.5/Trinity --seqType fq --max_memory 100G --CPU 40 --samples_file samples.txt --SS_lib_type RF 
```
