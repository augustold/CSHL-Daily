# SSPACE

```
cd /projects/augustold

mkdir SSPACE; cd SSPACE

wget https://www.baseclear.com/wp-content/uploads/SSPACE-STANDARD-v.-3.0-linux-x86_64.tar.gz
wget https://www.baseclear.com/wp-content/uploads/SSPACE-longread-v.-1-1.tar.gz
wget https://www.baseclear.com/wp-content/uploads/GapFiller-v.1-10-linux-x86_64.tar.gz
```

### SSPACE-longread
```
tar xvzf  SSPACE-longread-v.-1-1.tar.gz
cd SSPACE-LongRead_v1-1/
```

### PacBio long reads FASTQ to FASTA
```
cd /projects/augustold/SSPACE/SSPACE-LongRead_v1-1
sed -n '1~4s/^@/>/p;2~4p' /projects/augustold/PacBio_Microsoft_2015_09_08/017021_filtered_subreads.fastq > 017021_filtered_subreads.fasta 
sed -n '1~4s/^@/>/p;2~4p' /projects/augustold/PacBio_Microsoft_2015_09_08/17136_filtered_subreads.fastq > 17136_filtered_subreads.fasta 
```

### PacBio long reads FASTQ to FASTA
```
cd /projects/augustold/SSPACE/SSPACE-LongRead_v1-1
```

### Running
```
cd /projects/augustold/SSPACE/SSPACE-LongRead_v1-1 
perl SSPACE-LongRead.pl \
-c /projects/augustold/SP80_genome/sc.mlc.cns.sgl.utg.scga7.importdb.fa \
-p SP80_filtered_subreads.fastq \
-t 50
```

### Scaffold Stats

```
pip install assembly_stats

cd /projects/augustold/SSPACE/SSPACE-LongRead_v1-1/PacBio_scaffolder_results

assembly_stats scaffolds.fasta
```

```
{
  "Contig Stats": {
    "L10": 1478,
    "L20": 3951,
    "L30": 7175,
    "L40": 11095,
    "L50": 15693,
    "N10": 75345,
    "N20": 54499,
    "N30": 43541,
    "N40": 36628,
    "N50": 31578,
    "gc_content": 43.68660190391076,
    "longest": 468011,
    "mean": 21532.10471666851,
    "median": 16240.5,
    "sequence_count": 72424,
    "shortest": 1427,
    "total_bps": 1559441152
  },
  "Scaffold Stats": {
    "L10": 927,
    "L20": 2275,
    "L30": 3932,
    "L40": 5884,
    "L50": 8156,
    "N10": 140957,
    "N20": 109479,
    "N30": 91658,
    "N40": 78812,
    "N50": 67675,
    "gc_content": 43.68660190391076,
    "longest": 589528,
    "mean": 58781.60309626105,
    "median": 49154.0,
    "sequence_count": 28163,
    "shortest": 22753,
    "total_bps": 1655466288
  }
}
```

## BUSCO

### Instaling
```
## Instaling
conda create -n busco
conda activate busco

conda install -c bioconda -c conda-forge busco=4.0.2
```

### Running BUSCO
```
conda activate busco
cd /projects/augustold/SSPACE/SSPACE-LongRead_v1-1/PacBio_scaffolder_results

date # Tue Feb  4 15:13:52 -02 2020
busco -m genome -i scaffolds.fasta -o busco_results --auto-lineage-euk --cpu 20
date #
```

### Results

```
INFO:   

        --------------------------------------------------
        |Results from generic domain eukaryota_odb10      |
        --------------------------------------------------
        |C:78.8%[S:29.0%,D:49.8%],F:3.5%,M:17.7%,n:255    |
        |201    Complete BUSCOs (C)                       |
        |74     Complete and single-copy BUSCOs (S)       |
        |127    Complete and duplicated BUSCOs (D)        |
        |9      Fragmented BUSCOs (F)                     |
        |45     Missing BUSCOs (M)                        |
        |255    Total BUSCO groups searched               |
        --------------------------------------------------

        --------------------------------------------------
        |Results from dataset eudicots_odb10              |
        --------------------------------------------------
        |C:64.3%[S:34.4%,D:29.9%],F:3.7%,M:32.0%,n:2326   |
        |1495   Complete BUSCOs (C)                       |
        |800    Complete and single-copy BUSCOs (S)       |
        |695    Complete and duplicated BUSCOs (D)        |
        |85     Fragmented BUSCOs (F)                     |
        |746    Missing BUSCOs (M)                        |
        |2326   Total BUSCO groups searched               |
        --------------------------------------------------
INFO:   BUSCO analysis done. Total running time: 56254 seconds
```

## BLAST mitochondria and chloroplast

```
# Subfiles in
/projects/yutaka/Analysis/Genome/AlignmentSugarcaneGenomevsChloroplastMitochondria/
- Chloroplast.contigs.NC_005878.2.fa  
- Mitochondrial.contigs.LC107874.1.fa
- Mitochondrial.contigs.LC107875.1.fa
```

### BLAST conda install

```
mkdir -p ~/blast
cd ~/blast
conda create -n blast
conda activate blast
conda install -y blast
```

### Make data base: pacbio scaffold

```
conda activate blast
cd /projects/augustold/SSPACE/SSPACE-LongRead_v1-1/PacBio_scaffolder_results
makeblastdb -in scaffolds.fasta -dbtype nucl -parse_seqids > blast_prepare.log
```

### Run BLASTx

```
conda activate blast
#Chloroplast
blastn -query /projects/yutaka/Analysis/Genome/AlignmentSugarcaneGenomevsChloroplastMitochondria/Chloroplast.contigs.NC_005878.2.fa -db scaffolds.fasta -outfmt 7 -num_threads 12 -out Chloroplast_blast_output.txt

cat Chloroplast_blast_output.txt |awk '/hits found/{getline;print}' | grep -v "#" > Chloroplast_top_hits.txt

#Mitochondrial_1
blastn -query /projects/yutaka/Analysis/Genome/AlignmentSugarcaneGenomevsChloroplastMitochondria/Mitochondrial.contigs.LC107874.1.fa -db scaffolds.fasta -outfmt 7 -num_threads 12 -out Mitochondrial_1_blast_output.txt

cat Mitochondrial_1_blast_output.txt |awk '/hits found/{getline;print}' | grep -v "#" > Mitochondrial_1_top_hits.txt

Mitochondrial_2
blastn -query /projects/yutaka/Analysis/Genome/AlignmentSugarcaneGenomevsChloroplastMitochondria/Mitochondrial.contigs.LC107875.1.fa -db scaffolds.fasta -outfmt 7 -num_threads 12 -out Mitochondrial_2_blast_output.txt

cat Mitochondrial_2_blast_output.txt |awk '/hits found/{getline;print}' | grep -v "#" > Mitochondrial_2_top_hits.txt

```

### Get unique scaffold IDs using R

```
#Expample for Mitochondrial_1

R
data = read.table("Mitochondrial_1_top_hits.txt")
data2 = as.data.frame(table(data$V2))
data3 = data2$Var1
write.table(data3, file = "Mitochondrial_1_scaffold_unique_IDs.txt", row.names = F, col.names = F)
quit()

#remove "!
sed -i 's/"//g' Mitochondrial_1_scaffold_unique_IDs.txt 
```
