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
-c /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg.scga7.importdb.fa \
-p SP80_filtered_subreads.fastq \
-t 50
```
