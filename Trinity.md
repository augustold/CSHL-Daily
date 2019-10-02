## Oct 2 2019
## Augusto L. Diniz
## USP - CSHL

# Sugarcane SP80-3280 Transcriptome Assembly

## Trinity - *de novo*

```
# FastQ file links in
augustold@helix:/projects/augustold/fastq
```

### Running trinity
```
Trinity --seqType fq --max_memory 100G --CPU 40 \
         --left reads_1.fq.gz  --right reads_2.fq.gz --CPU 6

```
