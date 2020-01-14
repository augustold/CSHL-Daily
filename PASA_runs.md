# PASA tool for annotation update.

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

/sonas-hs/ware/hpc/home/diniz/software/PASApipeline.v2.4.1/bin/seqclean SP80.est.flc.mikado.combined.fasta -c 32
```
```
qsub -cwd -pe threads 32 -l m_mem_free=1G seqclean.sh 
```
  
Output summary:

```
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
source activate pasa

/sonas-hs/ware/hpc/home/diniz/software/PASApipeline.v2.4.1/Launch_PASA_pipeline.pl \
-c alignAssembly.config \
-C -R -T \
-g /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg.scga7.importdb.masked.fa \
-t SP80.est.flc.mikado.combined.fasta.clean \
-u SP80.est.flc.mikado.combined.fasta \
--ALIGNERS blat,gmap \
--CPU 32
```
```
qsub -cwd -pe threads 32 -l m_mem_free=1G alignAssembly.sh 
```

## Step 4: Annotation Comparisons and Annotation Updates

Loading your preexisting protein-coding gene annotations

```
/sonas-hs/ware/hpc/home/diniz/software/PASApipeline.v2.4.1/scripts/Load_Current_Gene_Annotations.dbi \
-c alignAssembly.config \
-g M162W.maker.repeatmasked.fasta \
-P M162W.maker.gene_only.gff3

```

Performing an annotation comparison and generating an updated gene set

```
/sonas-hs/ware/hpc/home/diniz/software/PASApipeline.v2.4.1/Launch_PASA_pipeline.pl \
-c annotCompare.config \
-A \
-g M162W.maker.repeatmasked.fasta \
-t maize.flc.iso.est.combined.fasta

```
