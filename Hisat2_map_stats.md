```
#conda create -n hisat2
conda activate hisat2
#conda install hisat2

cd /projects/augustold/CSHL/Hisat2_map_stats
mkdir R570 SP80-3280 AP85-441

#Build genome index

#AP85-441
screen
conda activate hisat2
cd /projects/augustold/CSHL/Saccharum_genome_refs/Sspon
mkdir hisat2_index
cd hisat2_index
hisat2-build -p 40 /projects/augustold/CSHL/Saccharum_genome_refs/Sspon/Sspon.HiC_chr_asm.fasta hisat2_index

#SP80-3280
screen
conda activate hisat2
cd /projects/augustold/CSHL/Saccharum_genome_refs/R570
mkdir hisat2_index
cd hisat2_index
hisat2-build -p 20 /projects/augustold/CSHL/Saccharum_genome_refs/SP803280/sc.mlc.cns.sgl.utg.scga7.importdb.fa hisat2_index

#R570
screen
conda activate hisat2
cd /projects/augustold/CSHL/Saccharum_genome_refs/R570
mkdir hisat2_index
cd hisat2_index
hisat2-build -p 20 /projects/augustold/CSHL/Saccharum_genome_refs/R570/single_tiling_path_assembly.fna hisat2_index


#Map RNAseq samples

#AP85-441
screen
cd /projects/augustold/CSHL/Hisat2_map_stats/AP85-441
conda activate hisat2

for SAMPLE in $(ls /projects/augustold/fastq/ | sed s/_[12].fq.gz// | sort -u); do
    hisat2 \
    -p 20 \
    -x /projects/augustold/CSHL/Saccharum_genome_refs/Sspon/hisat2_index/hisat2_index \
    -1 /projects/augustold/fastq/${SAMPLE}_1.fq.gz \
    -2 /projects/augustold/fastq/${SAMPLE}_2.fq.gz \
    -S ${SAMPLE}.sam \
    2> summary_${SAMPLE}.txt

    rm -rf ${SAMPLE}.sam
done

#SP80-3280
screen
cd /projects/augustold/CSHL/Hisat2_map_stats/SP80-3280
conda activate hisat2

for SAMPLE in $(ls /projects/augustold/fastq/ | sed s/_[12].fq.gz// | sort -u); do
    hisat2 \
    -p 20 \
    -x /projects/augustold/CSHL/Saccharum_genome_refs/SP803280/hisat2_index/hisat2_index \
    -1 /projects/augustold/fastq/${SAMPLE}_1.fq.gz \
    -2 /projects/augustold/fastq/${SAMPLE}_2.fq.gz \
    -S ${SAMPLE}.sam \
    2> summary_${SAMPLE}.txt

    rm -rf ${SAMPLE}.sam
done

#R570
screen
cd /projects/augustold/CSHL/Hisat2_map_stats/R570
conda activate hisat2

for SAMPLE in $(ls /projects/augustold/fastq/ | sed s/_[12].fq.gz// | sort -u); do
    hisat2 \
    -p 20 \
    -x /projects/augustold/CSHL/Saccharum_genome_refs/R570/hisat2_index/hisat2_index \
    -1 /projects/augustold/fastq/${SAMPLE}_1.fq.gz \
    -2 /projects/augustold/fastq/${SAMPLE}_2.fq.gz \
    -S ${SAMPLE}.sam \
    2> summary_${SAMPLE}.txt

    rm -rf ${SAMPLE}.sam
done

# Get summary tables

#R570
cd /projects/augustold/CSHL/Hisat2_map_stats/R570

python3 /projects/aliemelo/resources/scripts/make_hisat2_table.py -l summary_list.txt -o R570_map_summary.csv
```
