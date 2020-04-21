# Extensive de novo TE Annotator (EDTA) on sugarcane genome

## Installation
```
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/EDTA
```

Script: EDTA_install.sh
```
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/EDTA
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6
module load Anaconda2/5.3.0
module load git/2.19.1

conda create -n EDTA
conda activate EDTA
conda config --env --add channels anaconda --add channels conda-forge --add channels bioconda
conda install -n EDTA -y cd-hit repeatmodeler muscle mdust blast java-jdk perl perl-text-soundex multiprocess regex tensorflow=1.14.0 keras=2.2.4 scikit-learn=0.19.0 biopython pandas glob2 python=3.6 tesorter genericrepeatfinder genometools-genometools ltr_retriever
git clone https://github.com/oushujun/EDTA
./EDTA/EDTA.pl
```

```
qsub -cwd -pe threads 1 -l m_mem_free=8G EDTA_install.sh 
```

## 1. Get raw libraries from a genome (specify -type ltr|tir|helitron in different runs)

```
mkdir ltr tir helitron
```

### LTR
Script: EDTA_ltr.sh
```
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/EDTA/ltr
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6
module load Anaconda2/5.3.0
source activate EDTA

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sonas-hs/ware/hpc/home/diniz/.conda/envs/EDTA/lib

perl /sonas-hs/ware/hpc_norepl/data/diniz/analysis/EDTA/EDTA/EDTA_raw.pl \
--genome	/sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg.scga7.importdb.fa \
--type ltr \
--threads 30
```
```
qsub -cwd -pe threads 30 -l m_mem_free=0.5G EDTA_ltr.sh
```

### TIR
Script: EDTA_tir.sh
```
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/EDTA/tir
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6
module load Anaconda2/5.3.0
source activate EDTA

perl /sonas-hs/ware/hpc_norepl/data/diniz/analysis/EDTA/EDTA/EDTA_raw.pl \
-genome	/sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg.scga7.importdb.fa \
-type tir \
-threads 16
```
```
qsub -cwd -pe threads 16 -l m_mem_free=1G EDTA_tir.sh
```

### Helitron
Script: EDTA_helitron.sh
```
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/EDTA/tir
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6
module load Anaconda2/5.3.0
source activate EDTA

perl /sonas-hs/ware/hpc_norepl/data/diniz/analysis/EDTA/EDTA/EDTA_raw.pl \
-genome	/sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg.scga7.importdb.fa \
-type helitron \
-threads 16
```
```
qsub -cwd -pe threads 16 -l m_mem_free=1G EDTA_helitron.sh
```

## Finish the rest of the EDTA analysis

```
perl EDTA.pl -overwrite 0 [options]
```

## On Helix

```
cd /projects/augustold/CSHL/EDTA

conda activate EDTA

perl EDTA/EDTA_raw.pl \
--genome scga7.fa \
--species others \
--type tir \
--threads 50

helitron


#############

awk -v size=10000 -v pre=subset -v pad=3 '
   /^>/ { n++; if (n % size == 1) { close(fname); fname = sprintf("%s.%0" pad "d", pre, n) } }
   { print >> fname }
' scga7.fa

ls subset* | sort -u > list

for i in $(cat list)
do
mv ${i} ${i}.fa
done

for i in $(cat list)
do
mkdir ${i}
done

for i in $(cat list)
do
mv ${i}.fa ${i}
done

cd /projects/augustold/CSHL/EDTA/subset.001
perl ../EDTA/EDTA_raw.pl \
--genome subset.001.fa \
--species others \
--type tir \
--threads 20

for i in $(ls -d subset.1*)
do
cd /projects/augustold/CSHL/EDTA/${i}
perl ../EDTA/EDTA_raw.pl \
--genome ${i}.fa \
--species others \
--type tir \
--threads 20
done

for i in $(ls -d subset.2*)
do
cd /projects/augustold/CSHL/EDTA/${i}
perl ../EDTA/EDTA_raw.pl \
--genome ${i}.fa \
--species others \
--type tir \
--threads 20
done

for i in $(ls -d subset.3*)
do
cd /projects/augustold/CSHL/EDTA/${i}
perl ../EDTA/EDTA_raw.pl \
--genome ${i}.fa \
--species others \
--type tir \
--threads 20
done

for i in $(ls -d subset.4*)
do
cd /projects/augustold/CSHL/EDTA/${i}
perl ../EDTA/EDTA_raw.pl \
--genome ${i}.fa \
--species others \
--type tir \
--threads 20
done
```
