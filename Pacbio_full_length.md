## Map full-length unique isoforms back to the genome reference

### Considering unique high-quality isoforms across all samples 

script: GMAP_pacbio_hq_transcripts_run.sh

```
cd /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280

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
pacbio_hq_transcripts.fasta > pacbio_hq_transcripts.gff3

date
```
```
qsub -cwd -pe threads 16 -l m_mem_free=2G GMAP_pacbio_hq_transcripts_run.sh
```

# Running Iso-Seq 3 from scratch

Instructions: https://github.com/PacificBiosciences/IsoSeq_SA3nUP/wiki/Tutorial:-Installing-and-Running-Iso-Seq-3-using-Conda

### Installing

```
cd /sonas-hs/ware/hpc/home/diniz/FL_Pacbio

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load IntelPython/2.7.14
module load Anaconda2/5.3.0

conda create -n anaCogent5.2 python=2.7 anaconda
source activate anaCogent5.2

conda install -n anaCogent5.2 biopython
conda install -n anaCogent5.2 -c bioconda bx-python
conda install -n anaCogent5.2 -c bioconda pbccs=4.0
conda install -n anaCogent5.2 -c bioconda isoseq3=3.2

#Optional
conda install -n anaCogent5.2 -c bioconda pbcoretools # for manipulating PacBio datasets
conda install -n anaCogent5.2 -c bioconda bamtools    # for converting BAM to fasta
conda install -n anaCogent5.2 -c bioconda pysam       # for making CSV reports

#Check version
isoseq3 --version # isoseq3 3.2.x (commit v3.2.x)
ccs --version # ccs 4.0.x (commit v4.0.x)
lima --version # lima 1.9.0 (commit v1.9.0)
```

### Raw files in:
```
/isilon/seq/smrt3/data_root/runs/54333U/r54333U_20200123_202438/1_A01/
```

### Running Iso-Seq 3
```
```
