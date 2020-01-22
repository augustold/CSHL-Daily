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
source activate EDTA
conda config --env --add channels anaconda --add channels conda-forge --add channels biocore --add channels bioconda --add channels cyclus
conda install -n EDTA -y cd-hit repeatmodeler muscle mdust blast-legacy java-jdk perl perl-text-soundex multiprocess regex tensorflow=1.14.0 keras=2.2.4 scikit-learn=0.19.0 biopython pandas glob2 python=3.6 tesorter
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

perl /sonas-hs/ware/hpc_norepl/data/diniz/analysis/EDTA/EDTA/EDTA_raw.pl \
-genome	/sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg.scga7.importdb.fa \
-type ltr \
-threads 16
```
```
qsub -cwd -pe threads 16 -l m_mem_free=1G EDTA_ltr.sh
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
