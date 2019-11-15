## Genome-guided Trinity Transcriptome Assembly

Users must provide read alignments to Trinity as a coordinate-sorted bam file. In order to do that, I'll use STAR.

### STAR genome index for 1-pass
```
# Sontaneum
cd /projects/augustold/CSHL/Saccharum_genome_refs/Sspon/
mkdir star_index
mkdir star_index_2ndPass
screen

# generate genome for 1-pass

/projects/augustold/CSHL/STAR-2.7.2c/bin/Linux_x86_64/STAR --runThreadN 40 --runMode genomeGenerate --genomeDir /projects/augustold/CSHL/Saccharum_genome_refs/Sspon/star_index --genomeFastaFiles /projects/augustold/CSHL/Saccharum_genome_refs/Sspon/Sspon.HiC_chr_asm.fasta --sjdbGTFfile /projects/augustold/CSHL/Saccharum_genome_refs/Sspon/Sspon.v20190103.gff3 --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 150 --limitGenomeGenerateRAM 350000000000
```

### Align reads for 1-pass

```
# Result files in
augustold@helix:/projects/augustold/CSHL/Trinity/Trinity_guided/
officinarum  SP80-3280  spontaneum

# SP80-3280
cd /projects/augustold/CSHL/Trinity/Trinity_guided/spontaneum/
mkdir aligned_1

for i in $(ls fastq | sed s/_[12].fq.gz// | sort -u)
do
    /projects/augustold/CSHL/STAR-2.7.2c/bin/Linux_x86_64/STAR --genomeDir /projects/augustold/CSHL/Saccharum_genome_refs/Sspon/star_index --readFilesIn fastq/${i}_1.fq.gz,fastq/${i}_2.fq.gz --runThreadN 15 --outFileNamePrefix aligned_1/${i}. --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outSAMunmapped Within --outFilterType BySJout --outSAMattributes All --chimSegmentMin 12 --chimSegmentReadGapMax 3 --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterIntronMotifs RemoveNoncanonical --clip5pNbases 0 --seedSearchStartLmax 50 --genomeLoad LoadAndKeep --limitBAMsortRAM 35000000000
done
```

### Re-generate genome index for 2nd-pass

```
cd /projects/augustold/CSHL/Trinity/Trinity_guided/spontaneum/aligned_1/

/projects/augustold/CSHL/STAR-2.7.2c/bin/Linux_x86_64/STAR --runThreadN 40 --runMode genomeGenerate --genomeDir /projects/augustold/CSHL/Saccharum_genome_refs/Sspon/star_index_2ndPass --genomeFastaFiles /projects/augustold/CSHL/Saccharum_genome_refs/Sspon/Sspon.HiC_chr_asm.fasta --sjdbGTFfile /projects/augustold/CSHL/Saccharum_genome_refs/Sspon/Sspon.v20190103.gff3 --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 150 --limitGenomeGenerateRAM 350000000000 --limitSjdbInsertNsj 2000000 --sjdbFileChrStartEnd INI11.SJ.out.tab INI12.SJ.out.tab INI13.SJ.out.tab INI51.SJ.out.tab INI52.SJ.out.tab INI53.SJ.out.tab INI91.SJ.out.tab INI92.SJ.out.tab INI93.SJ.out.tab INL11.SJ.out.tab INL12.SJ.out.tab INL13.SJ.out.tab
```

### Align reads for 2nd-pass

```
cd /projects/augustold/CSHL/Trinity/Trinity_guided/spontaneum/

# 2nd pass
for i in $(ls fastq | sed s/_[12].fq.gz// | sort -u)
do
    /projects/augustold/CSHL/STAR-2.7.2c/bin/Linux_x86_64/STAR --genomeDir /projects/augustold/CSHL/Saccharum_genome_refs/Sspon/star_index_2ndPass --readFilesIn fastq/${i}_1.fq.gz,fastq/${i}_2.fq.gz --outSAMunmapped Within --outFilterType BySJout --outSAMattributes All --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix aligned_2/${i}. --clip5pNbases 0 --seedSearchStartLmax 50 --runThreadN 15 --genomeLoad LoadAndKeep --limitBAMsortRAM 35000000000
done

```

### Merge BAM files and build Genome-guided Trinity

```
# S. spontaneum
screen
cd /projects/augustold/CSHL/Trinity/Trinity_guided/spontaneum/

samtools merge atlas_merged.bam aligned_2/*.bam
samtools index atlas_merged.bam

# Running Trinity
/projects/augustold/CSHL/Trinity/trinityrnaseq-Trinity-v2.8.5/Trinity --genome_guided_bam atlas_merged.bam --genome_guided_max_intron 10000 --max_memory 100G --CPU 40 --SS_lib_type RF --output /projects/augustold/CSHL/Trinity/Trinity_guided/SP80-3280/trinity_out_dir
```
