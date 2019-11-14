## Nov 14th 2019

```
# GFF3 to fasta - TRINITY OK
/sonas-hs/ware/hpc/home/mcampbel/applications/cufflinks-2.2.1.Linux_x86_64/gffread trinity.gff3 \
-g /sonas-hs/ware/hpc/home/diniz/Saccharum_genome_refs/SP803280/sc.mlc.cns.sgl.utg.scga7.importdb.fa \
-x trinity.cds.fasta \
-y trinity.protein.fasta

# GFF3 to fasta - Stringtie FAIL
/sonas-hs/ware/hpc/home/mcampbel/applications/cufflinks-2.2.1.Linux_x86_64/gffread stringtie.gff3 \
-g /sonas-hs/ware/hpc/home/diniz/Saccharum_genome_refs/SP803280/sc.mlc.cns.sgl.utg.scga7.importdb.fa \
-x stringtie.cds.fasta \
-y stringtie.protein.fasta


# Removing redundant transcripts - TRINITY
/sonas-hs/ware/hpc/home/mcampbel/applications/cdhit/cd-hit-est -i  trinity.cds.fasta -o trinity.cds.fasta_cdhit_95.fasta -c 0.95 -n 10 -d 0 -M 3000 -t 10

grep '>' trinity.cds.fasta | wc -l # 1096846 seqs
grep '>' trinity.cds.fasta_cdhit_95.fasta | wc -l # 377049

# Align filtered Trinity transcripts back to the genome using gmap

# cat > GMAP_filt_trinity_run.sh

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load GMAP-GSNAP/2019-03-15
module load SAMtools/1.9

gmap -D /sonas-hs/ware/hpc/home/diniz/Saccharum_genome_refs/SP803280/ -d gmap_index -f gff3_gene -t 16 -n 1 --min-trimmed-coverage=0.80 --min-identity=0.80 trinity.cds.fasta_cdhit_95.fasta > trinity.cdhit_95.gff3

# Running on CSHL cluster
qsub -cwd -pe threads 16 -l m_mem_free=5G GMAP_filt_trinity_run.sh
```
