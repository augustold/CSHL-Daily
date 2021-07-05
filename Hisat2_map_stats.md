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
```

# Get summary tables
```
#R570
cd /projects/augustold/CSHL/Hisat2_map_stats/R570

ls /projects/augustold/CSHL/Hisat2_map_stats/R570/summary* | xargs -n 1 basename | cut -d "." -f1 | cut -d "_" -f2 | sort -u > sample_lst_A.txt
ls /projects/augustold/CSHL/Hisat2_map_stats/R570/summary* | xargs -n 1 basename | sort -u > sample_lst_B.txt
paste -d'\t' sample_lst_A.txt sample_lst_B.txt > sample_lst.txt
rm -rf sample_lst_A.txt sample_lst_B.txt

python3 /projects/aliemelo/resources/scripts/make_hisat2_table.py -l sample_lst.txt -o R570_map_summary.csv

```
# Get expression tables
```
# AP85-441
#screen
#conda activate RNAseq_SP803280_2021
#cd /projects/augustold/CSHL/Hisat2_map_stats/AP85-441/

for i in $(ls /projects/augustold/CSHL/Hisat2_map_stats/AP85-441/*bam | xargs -n 1 basename | cut -d "." -f1 | sort -u)
do
   stringtie ${i}.sorted.bam \
   -G /projects/augustold/CSHL/Saccharum_genome_refs/Sspon/Sspon.v20190103.gff3 \
   -e --rf -p 20 \
   -o stringtie/${i}.gtf \
   -A stringtie/${i}_gene_abund.tab
done


# SP80-3280
#screen
#conda activate RNAseq_SP803280_2021
#cd /projects/augustold/CSHL/Hisat2_map_stats/SP80-3280

for i in $(ls /projects/augustold/CSHL/Hisat2_map_stats/SP80-3280/*bam | xargs -n 1 basename | cut -d "." -f1 | sort -u)
do
   stringtie ${i}.sorted.bam \
   -G /projects/augustold/CSHL/Saccharum_genome_refs/SP803280/scga7.v2.gff3 \
   -e --rf -p 20 \
   -o stringtie_scga7_v2/${i}.gtf \
   -A stringtie_scga7_v2/${i}_gene_abund.tab
done
```

