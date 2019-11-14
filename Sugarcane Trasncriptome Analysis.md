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

### STAR genome index for 1-pass
```
# SP80-3280
cd /projects/augustold/CSHL/Saccharum_genome_refs/SP803280/
mkdir star_index
mkdir star_index_2ndPass
screen

# generate genome for 1-pass

/projects/augustold/CSHL/STAR-2.7.2c/bin/Linux_x86_64/STAR --runThreadN 40 --runMode genomeGenerate --genomeDir /projects/augustold/CSHL/Saccharum_genome_refs/SP803280/star_index --genomeFastaFiles /projects/augustold/CSHL/Saccharum_genome_refs/SP803280/sc.mlc.cns.sgl.utg.scga7.importdb.fa --sjdbGTFfile /projects/augustold/CSHL/Saccharum_genome_refs/SP803280/sc.mlc.cns.sgl.utg_scga7.sort.gff3 --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 150 --limitGenomeGenerateRAM 350000000000

# specify '--limitGenomeGenerateRAM' as 350000000000 due to error message:
# EXITING because of FATAL PARAMETER ERROR: limitGenomeGenerateRAM=200000000000is too small for your genome
# SOLUTION: please specify --limitGenomeGenerateRAM not less than 315007636181 and make that much RAM available
```

### Align reads for 1-pass

```
# Result files in
augustold@helix:/projects/augustold/CSHL/Trinity/Trinity_guided/
officinarum  SP80-3280  spontaneum

# SP80-3280
cd /projects/augustold/CSHL/Trinity/Trinity_guided/SP80-3280/
mkdir aligned_1

for i in $(ls fastq | sed s/_[12].fq.gz// | sort -u)
do
    /projects/augustold/CSHL/STAR-2.7.2c/bin/Linux_x86_64/STAR --genomeDir /projects/augustold/CSHL/Saccharum_genome_refs/SP803280/star_index --readFilesIn fastq/${i}_1.fq.gz,fastq/${i}_2.fq.gz --runThreadN 15 --outFileNamePrefix aligned_1/${i}. --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outSAMunmapped Within --outFilterType BySJout --outSAMattributes All --chimSegmentMin 12 --chimSegmentReadGapMax 3 --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterIntronMotifs RemoveNoncanonical --clip5pNbases 0 --seedSearchStartLmax 50 --genomeLoad LoadAndKeep --limitBAMsortRAM 35000000000
done
```

### Re-generate genome index for 2nd-pass

```
cd /projects/augustold/CSHL/Trinity/Trinity_guided/SP80-3280/aligned_1/

/projects/augustold/CSHL/STAR-2.7.2c/bin/Linux_x86_64/STAR --runThreadN 40 --runMode genomeGenerate --genomeDir /projects/augustold/CSHL/Saccharum_genome_refs/SP803280/star_index_2ndPass --genomeFastaFiles /projects/augustold/CSHL/Saccharum_genome_refs/SP803280/sc.mlc.cns.sgl.utg.scga7.importdb.fa --sjdbGTFfile /projects/augustold/CSHL/Saccharum_genome_refs/SP803280/sc.mlc.cns.sgl.utg_scga7.sort.gff3 --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 150 --limitGenomeGenerateRAM 350000000000 --limitSjdbInsertNsj 2000000 --sjdbFileChrStartEnd SPI11.SJ.out.tab SPI12.SJ.out.tab SPI13.SJ.out.tab SPI51.SJ.out.tab SPI52.SJ.out.tab SPI53.SJ.out.tab SPI91.SJ.out.tab SPI92.SJ.out.tab SPI93.SJ.out.tab SPL11.SJ.out.tab SPL12.SJ.out.tab SPL13.SJ.out.tab

# Fatal LIMIT error: the number of junctions to be inserted on the fly =1924868 is larger than the limitSjdbInsertNsj=1000000
# SOLUTION: re-run with at least --limitSjdbInsertNsj 1924868
# added: --limitSjdbInsertNsj 2000000

#########################################################

# S. spontaneum
cd projects/augustold/CSHL/Saccharum_genome_refs/Sspon/
mkdir star_index
screen

/projects/augustold/CSHL/STAR-2.7.2c/bin/Linux_x86_64/STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /projects/augustold/CSHL/Saccharum_genome_refs/Sspon/star_index --genomeFastaFiles /projects/augustold/CSHL/Saccharum_genome_refs/Sspon/Sspon.HiC_chr_asm.fasta --sjdbGTFfile /projects/augustold/CSHL/Saccharum_genome_refs/Sspon/Sspon.v20190103.gff3 --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 150
```


### Align reads for 2nd-pass

```
# RUNING it again (nov 11th 2019)
# unsing star_index_2ndPass

cd /projects/augustold/CSHL/Trinity/Trinity_guided/SP80-3280/

# 2nd pass
for i in $(ls fastq | sed s/_[12].fq.gz// | sort -u)
do
    /projects/augustold/CSHL/STAR-2.7.2c/bin/Linux_x86_64/STAR --genomeDir /projects/augustold/CSHL/Saccharum_genome_refs/SP803280/star_index_2ndPass --readFilesIn fastq/${i}_1.fq.gz,fastq/${i}_2.fq.gz --outSAMunmapped Within --outFilterType BySJout --outSAMattributes All --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix aligned_2/${i}. --clip5pNbases 0 --seedSearchStartLmax 50 --runThreadN 15 --genomeLoad LoadAndKeep --limitBAMsortRAM 35000000000
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

samtools merge atlas_merged.bam aligned_2/*.bam
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

### Align Trinity transcripts back to the genome using gmap

Running on CSHL cluster
```
# cat > GMAP_run.sh

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load GMAP-GSNAP/2019-03-15
module load SAMtools/1.9

gmap -D /sonas-hs/ware/hpc/home/diniz/Saccharum_genome_refs/SP803280/ -d gmap_index -f gff3_gene -t 16 -n 1 --min-trimmed-coverage=0.80 --min-identity=0.80 Trinity-GG.fasta > trinity.gff3

# Running on CSHL cluster
qsub -cwd -pe threads 16 -l m_mem_free=5G GMAP_run.sh
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


### PORTCULLIS

Installing following Peter's instructions:

1) download and install miniconda. In a shell, run these commands:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh  #(you'll need to answer "yes" to the license agreement and accept some default parameters)
```

2) Log out and then log back in so miniconda becomes your default python

3) run these commands:
```
conda create -n portcullis #just once
source activate portcullis

conda install portcullis --channel=bioconda
conda install pandas=0.24.2
```

Testing Portcullis step by step

```
#Prepare

portcullis prep \
   --threads 50 \
   --verbose \
   /projects/augustold/CSHL/Saccharum_genome_refs/SP803280/sc.mlc.cns.sgl.utg.scga7.importdb.fa /projects/augustold/CSHL/Trinity/Trinity_guided/SP80-3280/aligned_2/SPI12.Aligned.sortedByCoord.out.bam

#Junction

portcullis junc \
   --threads 50 \
   --verbose \
   --orientation RF \
   --strandedness firststrand \
   portcullis_prep

#Junction Filtering

portcullis filter \
   --threads 50 \
   --verbose \
   portcullis_prep portcullis_prep/portcullis.junctions.tab
```

Note: The filter step needs python 3.7!

Run on Helix
```
screen
source activate portcullis

cd /projects/augustold/CSHL/mikado_SP80/

portcullis full \
   --threads 50 \
   --verbose \
   --output portcullis_out \
   /projects/augustold/CSHL/Saccharum_genome_refs/SP803280/sc.mlc.cns.sgl.utg.scga7.importdb.fa /projects/augustold/CSHL/Trinity/Trinity_guided/SP80-3280/atlas_merged.bam

```

### MIKADO

Run on BNB!
```
/sonas-hs/ware/hpc/home/diniz/mikado_SP80
```

```
vi list.txt

#Paste the following
class.gtf   cl      True
cufflinks.gtf       cuff    True
stringtie.gtf       st      True    1
trinity.gff3        tr      False   -0.5

#save and quit!
```

Make a 'mikado_subSmp.sh' file. You can use 'vi'.
```
cd /sonas-hs/ware/hpc/home/diniz/mikado_SP80/

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6  
module load TransDecoder/5.5.0-Perl-5.28.0

#Creating the configuration file for Mikado
mikado configure \
--list list.txt \
--reference /sonas-hs/ware/hpc/home/diniz/Saccharum_genome_refs/SP803280/sc.mlc.cns.sgl.utg.scga7.importdb.fa \
-t 10 \
--mode permissive \
--scoring plant.yaml  \
--copy-scoring plant.yaml \
--junctions portcullis_all.junctions.bed \
-bt uniprot_sprot_plants.fasta \
configuration.yaml

#Mikado prepare
mikado prepare --json-conf configuration.yaml

#BLAST of the candidate transcripts
makeblastdb -in uniprot_sprot_plants.fasta -dbtype prot -parse_seqids > blast_prepare.log

blastx -max_target_seqs 5 -num_threads 10 -query mikado_prepared.fasta -outfmt 5 -db uniprot_sprot_plants.fasta -evalue 0.000001 2> blast.log | sed '/^$/d' | gzip -c - > mikado.blast.xml.gz

TransDecoder.LongOrfs -t mikado_prepared.fasta
TransDecoder.Predict -t mikado_prepared.fasta

#Mikado serialise
mikado serialise --json-conf configuration.yaml --xml mikado.blast.xml.gz --orfs mikado_prepared.fasta.transdecoder.bed --blast_targets uniprot_sprot_plants.fasta --transcripts mikado_prepared.fasta

#Mikado pick
mikado pick --json-conf configuration.yaml --subloci-out mikado.subloci.gff3 --procs 10
```

```
qsub -cwd -pe threads 10 -l m_mem_free=5G mikado_subSmp.sh
```
