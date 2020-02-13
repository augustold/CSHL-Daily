# PASA tool for annotation update.

## Install PASA

```
# Prepare working directory

cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/
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
conda install -c bioconda/label/cf201901 blat
conda install cdbtools
conda install openssl=1.0
conda install -c bioconda/label/cf201901 perl-dbi
conda install -c bioconda perl-dbd-sqlite
conda install -c bioconda perl-uri
conda install -c bioconda perl-db-file

wget https://github.com/PASApipeline/PASApipeline/releases/download/pasa-v2.4.1/PASApipeline.v2.4.1.FULL.tar.gz

tar -zxvf PASApipeline.v2.4.1.FULL.tar.gz
rm -rf PASApipeline.v2.4.1.FULL.tar.gz
cd PASApipeline.v2.4.1; make
```

## Input
- Repeatmasked genome
- Trasncripts
	- ESTs
	- Full-Lengths (Nishiyama 2014)
	- Pacbio Full-Lengths
	- Mikado CDS

## Step 1: combine transcripts

```
cd /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280

cat \
ESTs_SP80-3280.fasta \
sugarcane.fulllength.analysis.all.fasta \
pacbio_hq_transcripts.v2.fasta \
/sonas-hs/ware/hpc_norepl/data/diniz/analysis/mikado_3rd_test/mikado.loci.TErmv.cds.fasta \
> /sonas-hs/ware/hpc_norepl/data/diniz/analysis/PASA_run/SP80.est.flc.mikado.combined.fasta
```

Also get the full-lengths IDs in a different file
```
grep -e ">" pacbio_hq_transcripts.v2.fasta | sed 's/>//' > /sonas-hs/ware/hpc_norepl/data/diniz/analysis/PASA_run/FL_accs.txt
```

## Step 2: cleaning the transcript sequences

script: seqclean.sh
```
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/PASA_run

date

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6
module load Anaconda2/5.3.0
source activate pasa

PASApipeline.v2.4.1/bin/seqclean SP80.est.flc.mikado.combined.fasta -c 16

mkdir cleaning
mv cleaning_* cleaning/

date
```
```
qsub -cwd -pe threads 16 -l m_mem_free=1G seqclean.sh 
```
  
Output summary:

```
less seqcl_SP80.est.flc.mikado.combined.fasta.log

**************************************************
Sequences analyzed:    838049
-----------------------------------
                   valid:    837874  (6470 trimmed)
                 trashed:       175
**************************************************
----= Trashing summary =------
            by 'low_qual':        4
                by 'dust':       16
               by 'short':      131
              by 'shortq':       24
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

date

PASApipeline.v2.4.1/Launch_PASA_pipeline.pl \
-c alignAssembly.config \
-C -R \
-g /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg.scga7.importdb.masked.fa \
-t SP80.est.flc.mikado.combined.fasta.clean -T \
-u SP80.est.flc.mikado.combined.fasta \
-f FL_accs.txt \
--ALIGNERS gmap \
--CPU 30

date
```
```
qsub -cwd -pe threads 30 -l m_mem_free=0.5G alignAssembly.sh 
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

date

PASApipeline.v2.4.1/scripts/Load_Current_Gene_Annotations.dbi \
-c alignAssembly.config \
-g /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg.scga7.importdb.masked.fa \
-P /sonas-hs/ware/hpc_norepl/data/diniz/analysis/mikado_3rd_test/mikado.loci.TErmv.gff3

PASApipeline.v2.4.1/Launch_PASA_pipeline.pl \
-c annotCompare.config \
-A \
-g /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg.scga7.importdb.masked.fa \
-t SP80.est.flc.mikado.combined.fasta.clean

date
```
```
qsub -cwd -pe threads 1 -l m_mem_free=50G annotLoadandCompare.sh 
```

## Step 5: consider run this
```
https://github.com/PASApipeline/PASApipeline/wiki/PASA_abinitio_training_sets
```

script: abinitiotraining.sh
```
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/PASA_run
 
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6
module load Anaconda2/5.3.0
source activate pasa

date

PASApipeline.v2.4.1/scripts/pasa_asmbls_to_training_set.dbi \
--pasa_transcripts_fasta SP80.sqlite.assemblies.fasta \
--pasa_transcripts_gff3 SP80.sqlite.pasa_assemblies.gff3
```
```
qsub -cwd -pe threads 1 -l m_mem_free=16G abinitiotraining.sh

date
```
