## Map full-length unique isoforms back to the genome reference

script: GMAP_run.sh

```
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load GMAP-GSNAP/2019-03-15
module load SAMtools/1.9

date

gmap \
-D /sonas-hs/ware/hpc_norepl/data/diniz/analysis/SP80-3280 \
-d gmap_index \
-f gff3_gene \
-t 16 \
-n 1 \
--min-trimmed-coverage=0.70 \
--min-identity=0.95 \
/sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/
pacbio_hq_transcripts.fasta > pacbio_hq_transcripts.gff3

date
```
```
qsub -cwd -pe threads 16 -l m_mem_free=2G GMAP_run.sh
```
