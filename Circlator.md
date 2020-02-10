## Install

```
conda create -n circlator
conda activate circlator

conda install -c bioconda circlator
conda install -c bioconda canu=1.6-1 #circlator only works with canu version 1.6
```

## Get pacbio scaffolds for plastid genomes
```
#Run on BNB
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/pacbio_plastid

/sonas-hs/ware/hpc_norepl/data/programs/maker_v3_updated/maker/bin/fasta_tool --select Chloroplast_scaffold_unique_IDs.txt scaffolds.fasta > Chloroplast_scaffolds.fasta

/sonas-hs/ware/hpc_norepl/data/programs/maker_v3_updated/maker/bin/fasta_tool --select Mitochondrial_1_scaffold_unique_IDs.txt scaffolds.fasta > Mitochondrial_1_scaffolds.fasta

/sonas-hs/ware/hpc_norepl/data/programs/maker_v3_updated/maker/bin/fasta_tool --select Mitochondrial_2_scaffold_unique_IDs.txt scaffolds.fasta > Mitochondrial_2_scaffolds.fasta

```

## Run
```
conda activate circlator
cd /projects/augustold/circlator

#Chloroplast
circlator all --assembler canu \
Chloroplast_scaffolds.fasta \
/projects/yutaka/Analysis/Genome/AlignmentSugarcaneGenomevsChloroplastMitochondria/Chloroplast.contigs.NC_005878.2.fa \
Chloroplast_directory \
--threads 20 \
--verbose

#Mitochondrial_1
circlator all --assembler canu \
/projects/yutaka/Analysis/Genome/AlignmentSugarcaneGenomevsChloroplastMitochondria/LC107874.1.fa \
/projects/yutaka/Analysis/Genome/AlignmentSugarcaneGenomevsChloroplastMitochondria/Mitochondrial.contigs.LC107874.1.fa \
Mitochondrial_1_directory \
--threads 20 \
--verbose

#Mitochondrial_2
circlator all --assembler canu \
Mitochondrial_2_scaffolds.fasta \
/projects/yutaka/Analysis/Genome/AlignmentSugarcaneGenomevsChloroplastMitochondria/Mitochondrial.contigs.LC107875.1.fa \
Mitochondrial_2_directory \
--threads 20 \
--verbose
```
