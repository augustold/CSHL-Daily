# PASA tool for annotation update.

## Install PASA

```
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6
module load Anaconda2/5.3.0

cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/PASA_run
cp /mnt/grid/ware/hpc_norepl/data/data/kapeel/NAM/NAM_Canu1.8_new_runs/PASA_runs/environment.yml .
conda env create -f environment.yml
source activate mypasa
conda install gmap

wget https://github.com/PASApipeline/PASApipeline/releases/download/pasa-v2.3.3/PASApipeline-v2.3.3.tar.gz
tar -zxvf PASApipeline-v2.3.3.tar.gz
rm -rf PASApipeline-v2.3.3.tar.gz
cd PASApipeline-v2.3.3; make
```

```
# Prepare working directory

cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis
mkdir PASA_run; cd PASA_run; mkdir scripts

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6
module load Anaconda2/5.3.0
module load git/2.19.1

conda create -n pasa
source activate pasa
conda install -c bioconda pasa
conda install -c bioconda/label/cf201901 gmap
conda install -c bioconda/label/cf201901 ucsc-blat 

wget https://github.com/PASApipeline/PASApipeline/releases/download/pasa-v2.4.1/PASApipeline.v2.4.1.FULL.tar.gz

tar -zxvf PASApipeline.v2.4.1.FULL.tar.gz
rm -rf PASApipeline.v2.4.1.FULL.tar.gz
cd PASApipeline.v2.4.1; make
```

## Input
- Repeatmasked genome
- Trasncripts
	- ESTs
	- Full-Lengths
	- Mikado CDS

## Step 1: combine transcripts

```
cd /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280

cat \
ESTs_SP80-3280.fasta \
sugarcane.fulllength.analysis.all.fasta \
/sonas-hs/ware/hpc_norepl/data/diniz/analysis/mikado_3rd_test/mikado.loci.cds.fasta \
> /sonas-hs/ware/hpc_norepl/data/diniz/analysis/PASA_run/SP80.est.flc.mikado.combined.fasta
```

## Step 2: cleaning the transcript sequences

script: seqclean.sh
```
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/PASA_run 
 
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6
module load Anaconda2/5.3.0
source activate pasa

PASApipeline.v2.4.1/bin/seqclean SP80.est.flc.mikado.combined.fasta -c 16

mkdir cleaning
mv cleaning_* /cleaning
```
```
qsub -cwd -pe threads 16 -l m_mem_free=1G seqclean.sh 
```
  
Output summary:

```
less seqcl_SP80.est.flc.mikado.combined.fasta.log

**************************************************
Sequences analyzed:    577771
-----------------------------------
                   valid:    577694  (5503 trimmed)
                 trashed:        77
**************************************************
----= Trashing summary =------
            by 'low_qual':        4
               by 'short':       33
              by 'shortq':       24
                by 'dust':       16
------------------------------
Output file containing only valid and trimmed sequences: SP80.est.flc.mikado.combined.fasta.clean
For trimming and trashing details see cleaning report  : SP80.est.flc.mikado.combined.fasta.cln
--------------------------------------------------
```

## Step 3: Transcript alignments followed by alignment assembly

Before running the transcript alignment step, copy and edit the copy files annotCompare.config and alignAssembly.config from PASA installation folder.

```
#cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/PASA_run

cp /sonas-hs/ware/hpc_norepl/data/kapeel/NAM/NAM_Canu1.8_new_runs/PASA_runs/B97/alignAssembly.config .
cp /sonas-hs/ware/hpc_norepl/data/kapeel/NAM/NAM_Canu1.8_new_runs/PASA_runs/B97/annotCompare.config .
```

Edit the database path :

```
# database settings
DATABASE=/sonas-hs/ware/hpc_norepl/data/diniz/analysis/PASA_run/sqlite/SP80.sqlite
```

Create a sqlite folder and set the permission

```
mkdir sqlite; chmod -R 777 sqlite
```

Run the PASA alignment assembly pipeline

script: alignAssembly.sh
```
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/PASA_run 
 
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6
module load Anaconda2/5.3.0
module load GMAP-GSNAP/2019-03-15

source activate pasa

PASApipeline.v2.4.1/Launch_PASA_pipeline.pl \
-c alignAssembly.config \
-C -R -T \
-g /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg.scga7.importdb.masked.fa \
-t SP80.est.flc.mikado.combined.fasta.clean \
-u SP80.est.flc.mikado.combined.fasta \
--ALIGNERS gmap \
--CPU 16
```
```
qsub -cwd -pe threads 16 -l m_mem_free=2G alignAssembly.sh 
```

## Step 4: Annotation Comparisons and Annotation Updates

Loading your preexisting protein-coding gene annotations and performing an annotation comparison and generating an updated gene set

script: annotLoadandCompare.sh
```
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/PASA_run 
 
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6
module load Anaconda2/5.3.0
source activate pasa

PASApipeline.v2.4.1/scripts/Load_Current_Gene_Annotations.dbi \
-c alignAssembly.config \
-g /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg.scga7.importdb.masked.fa \
-P /sonas-hs/ware/hpc_norepl/data/diniz/analysis/mikado_3rd_test/mikado.loci.gff3

PASApipeline.v2.4.1/Launch_PASA_pipeline.pl \
-c annotCompare.config \
-A \
-g /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg.scga7.importdb.masked.fa \
-t SP80.est.flc.mikado.combined.fasta.clean
```
```
qsub -cwd -pe threads 1 -l m_mem_free=16G annotLoadandCompare.sh 
```
