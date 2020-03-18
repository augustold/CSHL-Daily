# PASA tool for annotation update.

## Split genome

### Get subset of contigs from mikado GFF

```
cd /mnt/grid/ware/hpc_norepl/data/data/diniz/analysis/SP80-3280/braker_masked_RNA

grep -v "^#" augustus.hints.gff3 | cut -f 1 | uniq > contigs_subset.txt

split -l 20000 contigs_subset.txt

#xaa
#xab
#xac
#xad
#xae
#xaf
#xag
#xah
#xai
#xaj
#xak
#xal
#xam
#xan
#xao
#xap
#xaq
#xar

grep -Fwf xaa augustus.hints.gff3 | awk '$3 !~ /superlocus/' > xaa.gff3
grep -Fwf xab augustus.hints.gff3 | awk '$3 !~ /superlocus/' > xab.gff3
grep -Fwf xac augustus.hints.gff3 | awk '$3 !~ /superlocus/' > xac.gff3
grep -Fwf xad augustus.hints.gff3 | awk '$3 !~ /superlocus/' > xad.gff3
grep -Fwf xae augustus.hints.gff3 | awk '$3 !~ /superlocus/' > xae.gff3
grep -Fwf xaf augustus.hints.gff3 | awk '$3 !~ /superlocus/' > xaf.gff3
grep -Fwf xag augustus.hints.gff3 | awk '$3 !~ /superlocus/' > xag.gff3
grep -Fwf xah augustus.hints.gff3 | awk '$3 !~ /superlocus/' > xah.gff3
grep -Fwf xai augustus.hints.gff3 | awk '$3 !~ /superlocus/' > xai.gff3
grep -Fwf xaj augustus.hints.gff3 | awk '$3 !~ /superlocus/' > xaj.gff3
grep -Fwf xak augustus.hints.gff3 | awk '$3 !~ /superlocus/' > xak.gff3
grep -Fwf xal augustus.hints.gff3 | awk '$3 !~ /superlocus/' > xal.gff3
grep -Fwf xam augustus.hints.gff3 | awk '$3 !~ /superlocus/' > xam.gff3
grep -Fwf xan augustus.hints.gff3 | awk '$3 !~ /superlocus/' > xan.gff3
grep -Fwf xao augustus.hints.gff3 | awk '$3 !~ /superlocus/' > xao.gff3
grep -Fwf xap augustus.hints.gff3 | awk '$3 !~ /superlocus/' > xap.gff3
grep -Fwf xaq augustus.hints.gff3 | awk '$3 !~ /superlocus/' > xaq.gff3
grep -Fwf xar augustus.hints.gff3 | awk '$3 !~ /superlocus/' > xar.gff3
mv x* PASA_run/
```

### Get subset of contigs fasta from original genome

```
cd /mnt/grid/ware/hpc_norepl/data/data/diniz/analysis/SP80-3280/braker_masked_RNA/PASA_run

for i in $(ls x* | sort -u)
do
/sonas-hs/ware/hpc_norepl/data/programs/maker_v3_updated/maker/bin/fasta_tool \
--select ${i} \
/sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg.scga7.importdb.masked.fa \
> genome/${i}.fasta
done

```

### Check fasta files

```
cd /mnt/grid/ware/hpc_norepl/data/data/diniz/analysis/SP80-3280/braker_masked_RNA/PASA_run

for i in $(ls x* | sort -u)
do
grep -c "^>" genome/${i}.fasta
done
```

## Step 1: Get transcripts files

```
cd /mnt/grid/ware/hpc_norepl/data/data/diniz/analysis/SP80-3280/braker_masked_RNA/PASA_run

cp /sonas-hs/ware/hpc_norepl/data/diniz/analysis/PASA_run/SP80.est.flc.mikado.combined.fasta .
cp /sonas-hs/ware/hpc_norepl/data/diniz/analysis/PASA_run/FL_accs.txt .
```

## Step 2: cleaning the transcript sequences

script: seqclean.sh
```
#!/bin/bash
cd /mnt/grid/ware/hpc_norepl/data/data/diniz/analysis/SP80-3280/braker_masked_RNA/PASA_run

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6
module load Anaconda2/5.3.0
source activate pasa

date

/sonas-hs/ware/hpc/home/diniz/.conda/envs/pasa/opt/pasa-2.4.1/bin/seqclean SP80.est.flc.mikado.combined.fasta -c 16

mkdir cleaning
mv cleaning_* cleaning/

date
```
```
qsub -cwd -pe threads 16 -l m_mem_free=1G seqclean.sh 
```

## Step 3: Transcript alignments followed by alignment assembly

Make separete directories for running the aligment for each subset of contigs
```
cd /mnt/grid/ware/hpc_norepl/data/data/diniz/analysis/SP80-3280/braker_masked_RNA/PASA_run
for i in $(ls x* | sort -u)
do
mkdir PASA_${i}
done
```

Before running the transcript alignment step, copy and edit the copy files annotCompare.config and alignAssembly.config from PASA installation folder.

```
cd /mnt/grid/ware/hpc_norepl/data/data/diniz/analysis/SP80-3280/braker_masked_RNA/PASA_run

cp /sonas-hs/ware/hpc_norepl/data/kapeel/NAM/NAM_Canu1.8_new_runs/PASA_runs/B97/alignAssembly.config .
cp /sonas-hs/ware/hpc_norepl/data/kapeel/NAM/NAM_Canu1.8_new_runs/PASA_runs/B97/annotCompare.config .
```

Edit the database path :

```
# database settings
DATABASE=/mnt/grid/ware/hpc_norepl/data/data/diniz/analysis/SP80-3280/braker_masked_RNA/PASA_run/PASA_xaa/sqlite/SP80.sqlite
```

Copy files and create a sqlite folder and set the permission
```
cd /mnt/grid/ware/hpc_norepl/data/data/diniz/analysis/SP80-3280/braker_masked_RNA/PASA_run
for i in $(ls x* | sort -u)
do
cd PASA_${i}
cp ../alignAssembly.config .
cp ../annotCompare.config .
cp ../alignAssembly.sh .
cp ../annotLoadandCompare.sh .
mkdir sqlite; chmod -R 777 sqlite
cd ..
done
```

Run the PASA alignment assembly pipeline

script: alignAssembly.sh
```
#!/bin/bash
cd /mnt/grid/ware/hpc_norepl/data/data/diniz/analysis/SP80-3280/braker_masked_RNA/PASA_run/PASA_xaa
 
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6
module load Anaconda2/5.3.0
module load GMAP-GSNAP/2019-03-15

#source activate pasa

date

/sonas-hs/ware/hpc/home/diniz/.conda/envs/pasa/opt/pasa-2.4.1/Launch_PASA_pipeline.pl \
-c alignAssembly.config \
-C -R \
-g ../genome/xaa.fasta \
-t ../SP80.est.flc.mikado.combined.fasta.clean -T \
-u ../SP80.est.flc.mikado.combined.fasta \
-f ../FL_accs.txt \
--ALIGNERS gmap \
--CPU 16

date
```
```
qsub -cwd -pe threads 16 -l m_mem_free=1G alignAssembly.sh
```

## Step 4: Annotation Comparisons and Annotation Updates

Loading your preexisting protein-coding gene annotations and performing an annotation comparison and generating an updated gene set

script: annotLoadandCompare.sh
```
#!/bin/bash
cd /mnt/grid/ware/hpc_norepl/data/data/diniz/analysis/SP80-3280/braker_masked_RNA/PASA_run/PASA_xaa

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6
module load Anaconda2/5.3.0
source activate pasa

date

/sonas-hs/ware/hpc/home/diniz/.conda/envs/pasa/opt/pasa-2.4.1/scripts/Load_Current_Gene_Annotations.dbi \
-c alignAssembly.config \
-g ../genome/xaa.fasta \
-P ../gff3/xaa.gff3

/sonas-hs/ware/hpc/home/diniz/.conda/envs/pasa/opt/pasa-2.4.1/Launch_PASA_pipeline.pl \
-c annotCompare.config \
-A \
-g ../genome/xaa.fasta \
--CPU 32 -t ../SP80.est.flc.mikado.combined.fasta.clean

date
```
```
qsub -cwd -pe threads 32 -l m_mem_free=1G annotLoadandCompare.sh 
```

## Replace path for each subset folder
```
cd /mnt/grid/ware/hpc_norepl/data/data/diniz/analysis/SP80-3280/braker_masked_RNA/PASA_run
for i in $(ls x* | sort -u)
do
cd PASA_${i}
sed -i "s/xaa/${i}/g" alignAssembly.config
sed -i "s/xaa/${i}/g" alignAssembly.sh
sed -i "s/xaa/${i}/g" annotCompare.config
sed -i "s/xaa/${i}/g" annotLoadandCompare.sh
cd ..
done
```

```
#1st batch
cd PASA_xaa/ ; qsub -cwd -pe threads 16 -l m_mem_free=1G alignAssembly.sh
cd ../PASA_xab/ ; qsub -cwd -pe threads 16 -l m_mem_free=1G alignAssembly.sh
cd ../PASA_xac/ ; qsub -cwd -pe threads 16 -l m_mem_free=1G alignAssembly.sh
cd ../PASA_xad/ ; qsub -cwd -pe threads 16 -l m_mem_free=1G alignAssembly.sh
cd ../PASA_xae/ ; qsub -cwd -pe threads 16 -l m_mem_free=1G alignAssembly.sh
cd ../PASA_xaf/ ; qsub -cwd -pe threads 16 -l m_mem_free=1G alignAssembly.sh
cd ../PASA_xag/ ; qsub -cwd -pe threads 16 -l m_mem_free=1G alignAssembly.sh
cd ../PASA_xah/ ; qsub -cwd -pe threads 16 -l m_mem_free=1G alignAssembly.sh
cd ../PASA_xai/ ; qsub -cwd -pe threads 16 -l m_mem_free=1G alignAssembly.sh
cd ../PASA_xaj/ ; qsub -cwd -pe threads 16 -l m_mem_free=1G alignAssembly.sh

#2nd batch
cd PASA_xak/ ; qsub -cwd -pe threads 16 -l m_mem_free=1G alignAssembly.sh
cd ../PASA_xal/ ; qsub -cwd -pe threads 16 -l m_mem_free=1G alignAssembly.sh
cd ../PASA_xam/ ; qsub -cwd -pe threads 16 -l m_mem_free=1G alignAssembly.sh
cd ../PASA_xan/ ; qsub -cwd -pe threads 16 -l m_mem_free=1G alignAssembly.sh
cd ../PASA_xao/ ; qsub -cwd -pe threads 16 -l m_mem_free=1G alignAssembly.sh
cd ../PASA_xap/ ; qsub -cwd -pe threads 16 -l m_mem_free=1G alignAssembly.sh
cd ../PASA_xaq/ ; qsub -cwd -pe threads 16 -l m_mem_free=1G alignAssembly.sh
cd ../PASA_xar/ ; qsub -cwd -pe threads 16 -l m_mem_free=1G alignAssembly.sh
```

```
cd /mnt/grid/ware/hpc_norepl/data/data/diniz/analysis/SP80-3280/braker_masked_RNA/PASA_run
for i in $(ls x* | sort -u)
do
cd PASA_${i}
sed -i "s/#source/source/g" alignAssembly.sh
cd ..
done
```

```
cd PASA_xaa/ ; qsub -cwd -pe threads 32 -l m_mem_free=1G annotLoadandCompare.sh
cd ../PASA_xab/ ; qsub -cwd -pe threads 32 -l m_mem_free=1G annotLoadandCompare.sh
cd ../PASA_xac/ ; qsub -cwd -pe threads 32 -l m_mem_free=1G annotLoadandCompare.sh
cd ../PASA_xad/ ; qsub -cwd -pe threads 32 -l m_mem_free=1G annotLoadandCompare.sh
cd ../PASA_xae/ ; qsub -cwd -pe threads 32 -l m_mem_free=1G annotLoadandCompare.sh
cd ../PASA_xaf/ ; qsub -cwd -pe threads 32 -l m_mem_free=1G annotLoadandCompare.sh
cd ../PASA_xag/ ; qsub -cwd -pe threads 32 -l m_mem_free=1G annotLoadandCompare.sh
cd ../PASA_xah/ ; qsub -cwd -pe threads 32 -l m_mem_free=1G annotLoadandCompare.sh
cd ../PASA_xai/ ; qsub -cwd -pe threads 32 -l m_mem_free=1G annotLoadandCompare.sh
cd ../PASA_xaj/ ; qsub -cwd -pe threads 32 -l m_mem_free=1G annotLoadandCompare.sh
cd ../PASA_xak/ ; qsub -cwd -pe threads 32 -l m_mem_free=1G annotLoadandCompare.sh
cd ../PASA_xal/ ; qsub -cwd -pe threads 32 -l m_mem_free=1G annotLoadandCompare.sh
cd ../PASA_xam/ ; qsub -cwd -pe threads 32 -l m_mem_free=1G annotLoadandCompare.sh
cd ../PASA_xan/ ; qsub -cwd -pe threads 32 -l m_mem_free=1G annotLoadandCompare.sh
cd ../PASA_xao/ ; qsub -cwd -pe threads 32 -l m_mem_free=1G annotLoadandCompare.sh
cd ../PASA_xap/ ; qsub -cwd -pe threads 32 -l m_mem_free=1G annotLoadandCompare.sh
cd ../PASA_xaq/ ; qsub -cwd -pe threads 32 -l m_mem_free=1G annotLoadandCompare.sh
cd ../PASA_xar/ ; qsub -cwd -pe threads 32 -l m_mem_free=1G annotLoadandCompare.sh
```

## Combine 'gene_structures_post_PASA_updates' files
```

cd /mnt/grid/ware/hpc_norepl/data/data/diniz/analysis/SP80-3280/braker_masked_RNA/PASA_run

# mkdir gene_structures_post_PASA_updates

for i in $(ls x* | sort -u)
do
mv PASA_${i}/SP80.sqlite.gene_structures_post_PASA_updates.*.gff3 PASA_${i}/SP80.sqlite.gene_structures_post_PASA_updates.${i}.gff3
cp PASA_${i}/SP80.sqlite.gene_structures_post_PASA_updates.${i}.gff3 gene_structures_post_PASA_updates/
done

cat \
SP80.sqlite.gene_structures_post_PASA_updates.xaa.gff3 \
SP80.sqlite.gene_structures_post_PASA_updates.xab.gff3 \
SP80.sqlite.gene_structures_post_PASA_updates.xac.gff3 \
SP80.sqlite.gene_structures_post_PASA_updates.xad.gff3 \
SP80.sqlite.gene_structures_post_PASA_updates.xae.gff3 \
SP80.sqlite.gene_structures_post_PASA_updates.xaf.gff3 \
SP80.sqlite.gene_structures_post_PASA_updates.xag.gff3 \
SP80.sqlite.gene_structures_post_PASA_updates.xah.gff3 \
SP80.sqlite.gene_structures_post_PASA_updates.xai.gff3 \
SP80.sqlite.gene_structures_post_PASA_updates.xaj.gff3 \
SP80.sqlite.gene_structures_post_PASA_updates.xak.gff3 \
SP80.sqlite.gene_structures_post_PASA_updates.xal.gff3 \
SP80.sqlite.gene_structures_post_PASA_updates.xam.gff3 \
SP80.sqlite.gene_structures_post_PASA_updates.xan.gff3 \
SP80.sqlite.gene_structures_post_PASA_updates.xao.gff3 \
SP80.sqlite.gene_structures_post_PASA_updates.xap.gff3 \
SP80.sqlite.gene_structures_post_PASA_updates.xaq.gff3 \
SP80.sqlite.gene_structures_post_PASA_updates.xar.gff3 \
> BRAKER_PASA.gff3

# Get fasta
gffread BRAKER_PASA.gff3 -g /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg.scga7.importdb.fa -x BRAKER_PASA.gff3.cds.fasta

# Summary stats
/sonas-hs/ware/hpc_norepl/data/kapeel/scripts/gff2gene_stats_withUTRs_canonical.pl BRAKER_PASA.gff3
```

## TE-filtering

```
#Copy fasta to Helix

cd /projects/augustold/CSHL/TEsorter/SP80_braker
conda activate py2
python ../TEsorter.py BRAKER_PASA.gff3.cds.fasta -eval 1e-6 -p 30
```

## Filtering MIKADO gff3
```
cd /mnt/grid/ware/hpc_norepl/data/data/diniz/analysis/SP80-3280/braker_masked_RNA/PASA_run/gene_structures_post_PASA_updates
# get TEsorter output from HELIX
scp -r augustold@143.107.54.134:/projects/augustold/CSHL/TEsorter/SP80_braker/BRAKER_PASA.gff3.cds.fasta.rexdb.cls.tsv .

grep -v "^#" BRAKER_PASA.gff3.cds.fasta.rexdb.cls.tsv | cut -f 1 | sort | uniq | cut -d "." -f 1-2 | sort |uniq > TE-genes.txt
grep -Fwvf TE-genes.txt BRAKER_PASA.gff3 | awk '$3 !~ /superlocus/' > temp.gff3

```
