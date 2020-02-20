# Now, let's consider only nuclear contigs

### Excluding plastid contigs

```
cd /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280

#Copy plastid fasta files

scp -r augustold@143.107.54.134:/projects/yutaka/Analysis/Genome/AlignmentSugarcaneGenomevsChloroplastMitochondria/Mitochondrial.contigs.LC107874.1.fa .

scp -r augustold@143.107.54.134:/projects/yutaka/Analysis/Genome/AlignmentSugarcaneGenomevsChloroplastMitochondria/Mitochondrial.contigs.LC107875.1.fa .

scp -r augustold@143.107.54.134:/projects/yutaka/Analysis/Genome/AlignmentSugarcaneGenomevsChloroplastMitochondria/Chloroplast.contigs.NC_005878.2.fa .

#Get contig IDs
grep -e ">" Mitochondrial.contigs.LC107874.1.fa > Mitochondrial.contigs.LC107874.1.txt
grep -e ">" Mitochondrial.contigs.LC107875.1.fa > Mitochondrial.contigs.LC107875.1.txt
grep -e ">" Chloroplast.contigs.NC_005878.2.fa > Chloroplast.contigs.NC_005878.2.txt

cat Mitochondrial.contigs.LC107874.1.txt Mitochondrial.contigs.LC107875.1.txt Chloroplast.contigs.NC_005878.2.txt > Plastid.contigs.v1.txt

sed 's/>//g' Plastid.contigs.v1.txt > Plastid.contigs.v2.txt
sed 's/ 1//g' Plastid.contigs.v2.txt > Plastid.contigs.txt

# Sort contig IDs manualy :(

grep -e ">" sc.mlc.cns.sgl.utg.scga7.importdb.fa > scga7.contigs.v1.txt

sed 's/>//g' scga7.contigs.v1.txt > scga7.contigs.txt

### Use R!
R

scga = read.table("scga7.contigs.txt")
plastid = read.table("Plastid.contigs.txt")

df = as.data.frame(rbind(scga, plastid))
df2 = as.data.frame(table(df))
df3 = subset(df2, Freq == 1)

write.table(df3$df, file = "Nuclear.scga7.contigs.v1.txt", row.names = F, col.names = F)

quit()

### Back to terminal!

sed 's/"//g' Nuclear.scga7.contigs.v1.txt > Nuclear.scga7.contigs.txt

/sonas-hs/ware/hpc_norepl/data/programs/maker_v3_updated/maker/bin/fasta_tool --select Nuclear.scga7.contigs.txt sc.mlc.cns.sgl.utg.scga7.importdb.fa > Nuclear.scga7.fa
```

## Map full-length unique isoforms back to the NEW genome reference

### Build GMPA index
script: GMAP_build.sh
```
#!/bin/bash
cd /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/Nuclear.scga7

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load GMAP-GSNAP/2019-03-15
module load SAMtools/1.9

date

gmap_build -d gmap_index -D . -k 13 ../Nuclear.scga7.fa

date
```
```
qsub -cwd -pe threads 1 -l m_mem_free=16G GMAP_build.sh
```

### Considering unique high-quality isoforms across all samples 

script: GMAP_pacbio_hq_transcripts_run.sh

```
#!/bin/bash
cd /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/Nuclear.scga7

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load GMAP-GSNAP/2019-03-15
module load SAMtools/1.9

date

gmap \
-D /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280/Nuclear.scga7 \
-d gmap_index \
-f gff3_gene \
-t 16 \
-n 1 \
--min-trimmed-coverage=0.70 \
--min-identity=0.95 \
../pacbio_hq_transcripts.fasta > pacbio_hq_transcripts.gff3

date
```
```
qsub -cwd -pe threads 16 -l m_mem_free=2G GMAP_pacbio_hq_transcripts_run.sh
```

## Map full-length unique isoforms back to the genome reference - Scaffolds from pacbio data

