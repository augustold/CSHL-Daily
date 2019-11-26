### STAR genome index for 1-pass

script: star_1stpass.sh

```
cd /sonas-hs/ware/hpc/home/diniz/AP85-441/
mkdir star_index

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load STAR/2.7.0d

STAR --runThreadN 16 \
--runMode genomeGenerate \
--genomeDir /sonas-hs/ware/hpc/home/diniz/AP85-441/star_index \
--genomeFastaFiles /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/AP85-441/Sspon.HiC_chr_asm.fasta \
--sjdbGTFfile /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/AP85-441/Sspon.v20190103.gff3 \
--genomeChrBinNbits 13 \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbOverhang 150 \
--limitGenomeGenerateRAM 350000000000
```
```
qsub -cwd -pe threads 16 -l m_mem_free=5G star_1stpass.sh
```

### Align reads for 1-pass

script: star_align_1stpass.sh

```
cd /sonas-hs/ware/hpc/home/diniz/AP85-441/
mkdir aligned_1

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load STAR/2.7.0d

for i in $(ls /sonas-hs/ware/hpc_norepl/data/diniz/fastq/spont | sed s/_[12].fq.gz// | sort -u)
do
    STAR \
    --genomeDir /sonas-hs/ware/hpc/home/diniz/AP85-441/star_index \
    --readFilesIn /sonas-hs/ware/hpc_norepl/data/diniz/fastq/spont/${i}_1.fq.gz,/sonas-hs/ware/hpc_norepl/data/diniz/fastq/spont/${i}_2.fq.gz \
    --runThreadN 16 \
    --outFileNamePrefix aligned_1/${i}. \
    --outSAMtype BAM SortedByCoordinate \
    --readFilesCommand zcat \
    --outSAMunmapped Within \
    --outFilterType BySJout \
    --outSAMattributes All \
    --chimSegmentMin 12 \
    --chimSegmentReadGapMax 3 \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterIntronMotifs RemoveNoncanonical \
    --clip5pNbases 0 \
    --seedSearchStartLmax 50 \
    --genomeLoad LoadAndKeep \
    --limitBAMsortRAM 35000000000
done
```
```
qsub -cwd -pe threads 16 -l m_mem_free=5G star_align_1stpass.sh
```

### STAR genome index for 2-pass

script: star_2ndpass.sh

```
cd /sonas-hs/ware/hpc/home/diniz/AP85-441/
mkdir star_index_2ndPass
cd /sonas-hs/ware/hpc/home/diniz/AP85-441/aligned_1/

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load STAR/2.7.0d

STAR --runThreadN 16 \
--runMode genomeGenerate \
--genomeDir /sonas-hs/ware/hpc/home/diniz/AP85-441/star_index_2ndPass \
--genomeFastaFiles /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/AP85-441/Sspon.HiC_chr_asm.fasta \
--sjdbGTFfile /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/AP85-441/Sspon.v20190103.gff3 \
--genomeChrBinNbits 13 \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbOverhang 150 \
--limitGenomeGenerateRAM 350000000000 \
--limitSjdbInsertNsj 2000000 \
--sjdbFileChrStartEnd INI11.SJ.out.tab INI12.SJ.out.tab INI13.SJ.out.tab INI51.SJ.out.tab INI52.SJ.out.tab INI53.SJ.out.tab INI91.SJ.out.tab INI92.SJ.out.tab INI93.SJ.out.tab INL11.SJ.out.tab INL12.SJ.out.tab INL13.SJ.out.tab
```
```
qsub -cwd -pe threads 16 -l m_mem_free=5G star_2ndpass.sh
```
