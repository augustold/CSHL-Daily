## Install

```
conda create -n circlator
conda activate circlator

conda install -c bioconda circlator
conda install -c bioconda canu=1.6-1 #circlator only works with canu version 1.6
```

## Run
```
conda activate circlator
cd /projects/augustold/circlator

circlator all --assembler canu \
/projects/yutaka/Analysis/Genome/AlignmentSugarcaneGenomevsChloroplastMitochondria/Mitochondrial.contigs.LC107875.1.fa \
/projects/augustold/SSPACE/SSPACE-LongRead_v1-1/SP80_filtered_subreads.fastq \
Mitochondrial_1_directory \
--threads 20 \
--verbose
```
