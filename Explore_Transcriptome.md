# Transcriptome Assembly Quality Assessment
https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Assembly-Quality-Assessment

```
# Transcriptome results folder
cd /projects/augustold/CSHL/Trinity/Trinity_denovo/SP80-3280/trinity_out_dir
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

## Full-length transcript analysis for model and non-model organisms using BLAST+

Protein databases on Helix:

```
#SwissProt
cd /projects/augustold/CSHL/Trinity/Trinity_denovo
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz 
makeblastdb -in uniprot_sprot.fasta -dbtype prot
```

Running BLAST

```
cd /projects/augustold/CSHL/Trinity/Trinity_denovo/SP80-3280/trinity_out_dir

blastx \
-query Trinity.fasta \
-db /projects/augustold/CSHL/Trinity/Trinity_denovo/uniprot_sprot.fasta \
-out blastx.outfmt6 \
-evalue 1e-20 -num_threads 40 -max_target_seqs 1 -outfmt 6
```

Next, examine the percent of the target being aligned to by the best matching Trinity transcript, like so:

```
cd /projects/augustold/CSHL/Trinity/Trinity_denovo/SP80-3280/trinity_out_dir

/projects/augustold/CSHL/Trinity/trinityrnaseq-Trinity-v2.8.5/util/analyze_blastPlus_topHit_coverage.pl \
blastx.outfmt6 \
Trinity.fasta \
/projects/augustold/CSHL/Trinity/Trinity_denovo/uniprot_sprot.fasta
```

```
#hit_pct_cov_bin        count_in_bin    >bin_below
100     4285    4285
90      1540    5825
80      1088    6913
70      858     7771
60      758     8529
50      796     9325
40      823     10148
30      774     10922
20      624     11546
10      184     11730
```

Statements we can make based on this table include:

- There are 1540 proteins that each match a Trinity transcript by >80% and <= 90% of their protein lengths.
- There are 5825 proteins that are represented by nearly full-length transcripts, having >80% alignment coverage.
- There are 4285 proteins that are covered by more than 90% of their protein lengths.

