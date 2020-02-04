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
