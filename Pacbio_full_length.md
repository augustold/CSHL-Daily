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

## Step 2 SUPPA/astalavista


This step is to identify the intron-retention isoforms using SUPPA/astalavista.

### Convert GFF3 to GTF

Note: An annotation file in GTF format is required.

```
gffread -E pacbio_hq_transcripts.gff3 -T -o pacbio_hq_transcripts.gtf
```

### Run suppa.py to identify the intron-retention isoforms

script: suppa.sh
```
#!/bin/bash
cd /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6

date

python /sonas-hs/ware/hpc/home/bwang/software/SUPPA-BUILD/suppa.py generateEvents -i pacbio_hq_transcripts.gtf -o cdna_RI_p -e RI -p 

date
```
```
qsub -cwd -pe threads 1 -l m_mem_free=16G suppa.sh
```

## Map full-length unique isoforms back to the genome reference - Scaffolds from pacbio data

### Build GMPA index
script: GMAP_build.sh
```
#!/bin/bash
cd /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/scaffolds

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load GMAP-GSNAP/2019-03-15
module load SAMtools/1.9

date

gmap_build -d gmap_index -D . -k 13 scaffolds.fasta

date
```
```
qsub -cwd -pe threads 1 -l m_mem_free=16G GMAP_build.sh
```
### Map to scaffolds
script: GMAP_pacbio_hq_transcripts_run.sh
```
#!/bin/bash
cd /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/scaffolds

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load GMAP-GSNAP/2019-03-15
module load SAMtools/1.9

date

gmap \
-D /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/scaffolds \
-d gmap_index \
-f gff3_gene \
-t 16 \
-n 1 \
--min-trimmed-coverage=0.80 \
--min-identity=0.90 \
../pacbio_hq_transcripts.fasta > pacbio_hq_transcripts_to_scaffolds.gff3

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
conda install -c bioconda pbccs
conda install -n anaCogent5.2 -c bioconda isoseq3=3.2

#Optional
conda install -n anaCogent5.2 -c bioconda pbcoretools # for manipulating PacBio datasets
conda install -n anaCogent5.2 -c bioconda bamtools    # for converting BAM to fasta
conda install -n anaCogent5.2 -c bioconda pysam       # for making CSV reports

#Check version
isoseq3 --version # isoseq3 3.2.2 (commit v3.2.2)
ccs --version # ccs 4.0.x (commit v4.0.x)
lima --version # lima 1.9.0 (commit v1.9.0)
```

### Raw files in:
```
/isilon/seq/smrt3/data_root/runs/54333U/r54333U_20200123_202438/1_A01/

# Sara's results are in:
/isilon/seq/smrt4/smrtlink/userdata/jobs_root/0000/0000001/0000001589
```

### Get flnc for each barcode

```
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load SAMtools/1.9

cd /sonas-hs/ware/hpc/home/diniz/FL_Pacbio

# Subset reads per barcode using R

R

data = read.csv("/seq/smrt4/smrtlink/userdata/jobs_root/0000/0000001/0000001589/outputs/flnc.report.csv", header = T)
data2 = as.data.frame(table(data$primer))
write.csv(data2, file = "barcode_count.txt")

bc1001 = subset(data, primer == "bc1001_5p--bc1001_3p")
bc1002 = subset(data, primer == "bc1002_5p--bc1002_3p")
bc1003 = subset(data, primer == "bc1003_5p--bc1003_3p")
bc1004 = subset(data, primer == "bc1004_5p--bc1004_3p")

write.table(bc1001$id, file = "bc1001_temp.txt", row.names = F, col.names = F)
write.table(bc1002$id, file = "bc1002_temp.txt", row.names = F, col.names = F)
write.table(bc1003$id, file = "bc1003_temp.txt", row.names = F, col.names = F)
write.table(bc1004$id, file = "bc1004_temp.txt", row.names = F, col.names = F)

quit()

# Back to terminal

sed 's/"//g' bc1001_temp.txt > bc1001.txt
sed 's/"//g' bc1002_temp.txt > bc1002.txt
sed 's/"//g' bc1003_temp.txt > bc1003.txt
sed 's/"//g' bc1004_temp.txt > bc1004.txt

rm -rf *_temp.txt
```

Use this python script: idfilter.py
```
#!/usr/bin/env python3

import sys

with open(sys.argv[1], 'r') as indexfile:
    ids = set(l.rstrip('\r\n') for l in indexfile)

for line in sys.stdin:
    qname, _ = line.split('\t', 1)
    if qname in ids:
        sys.stdout.write(line)
```

Get flnc for each barcode
```
samtools view /seq/smrt4/smrtlink/userdata/jobs_root/0000/0000001/0000001589/outputs/flnc.bam | ./idfilter.py bc1001.txt > bc1001.flnc.sam
samtools view -S -b bc1001.flnc.sam > bc1001.flnc.bam
rm -rf bc1001.flnc.sam

samtools view /seq/smrt4/smrtlink/userdata/jobs_root/0000/0000001/0000001589/outputs/flnc.bam | ./idfilter.py bc1002.txt > bc1002.flnc.sam
samtools view -S -b bc1002.flnc.sam > bc1002.flnc.bam
rm -rf bc1002.flnc.sam

samtools view /seq/smrt4/smrtlink/userdata/jobs_root/0000/0000001/0000001589/outputs/flnc.bam | ./idfilter.py bc1003.txt > bc1003.flnc.sam
samtools view -S -b bc1003.flnc.sam > bc1003.flnc.bam
rm -rf bc1003.flnc.sam

samtools view /seq/smrt4/smrtlink/userdata/jobs_root/0000/0000001/0000001589/outputs/flnc.bam | ./idfilter.py bc1004.txt > bc1004.flnc.sam
samtools view -S -b bc1004.flnc.sam > bc1004.flnc.bam
rm -rf bc1004.flnc.sam
```

### Running Iso-Seq 3
script: isoseq3_cluster_bc1001.sh
```
#!/bin/bash
cd /sonas-hs/ware/hpc/home/diniz/FL_Pacbio

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load IntelPython/2.7.14
module load Anaconda2/5.3.0

source activate anaCogent5.2

date

isoseq3 cluster bc1001.flnc.bam bc1001.polished.bam --verbose --use-qvs -j 16

date
```
```
qsub -cwd -pe threads 16 -l m_mem_free=2G isoseq3_cluster_bc1001.sh
```

script: isoseq3_cluster_bc1002.sh
```
#!/bin/bash
cd /sonas-hs/ware/hpc/home/diniz/FL_Pacbio

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load IntelPython/2.7.14
module load Anaconda2/5.3.0

source activate anaCogent5.2

date

isoseq3 cluster bc1002.flnc.bam bc1002.polished.bam --verbose --use-qvs -j 16

date
```
```
qsub -cwd -pe threads 16 -l m_mem_free=2G isoseq3_cluster_bc1002.sh
```

script: isoseq3_cluster_bc1003.sh
```
#!/bin/bash
cd /sonas-hs/ware/hpc/home/diniz/FL_Pacbio

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load IntelPython/2.7.14
module load Anaconda2/5.3.0

source activate anaCogent5.2

date

isoseq3 cluster bc1003.flnc.bam bc1003.polished.bam --verbose --use-qvs -j 16

date
```
```
qsub -cwd -pe threads 16 -l m_mem_free=2G isoseq3_cluster_bc1003.sh
```

script: isoseq3_cluster_bc1004.sh
```
#!/bin/bash
cd /sonas-hs/ware/hpc/home/diniz/FL_Pacbio

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load IntelPython/2.7.14
module load Anaconda2/5.3.0

source activate anaCogent5.2

date

isoseq3 cluster bc1004.flnc.bam bc1004.polished.bam --verbose --use-qvs -j 16

date
```
```
qsub -cwd -pe threads 16 -l m_mem_free=2G isoseq3_cluster_bc1004.sh
```




#### 0. Generate CCS
script: generate_ccs.sh
```
#!/bin/bash
cd /sonas-hs/ware/hpc/home/diniz/FL_Pacbio

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load IntelPython/2.7.14
module load Anaconda2/5.3.0

source activate anaCogent5.2

date

/seq/smrt4/smrtlink/smrtcmds/bin/ccs \
/isilon/seq/smrt3/data_root/runs/54333U/r54333U_20200123_202438/1_A01/m54333U_200123_203332.subreads.bam \
m54333U_200123_203332.ccs.bam --min-rq 0.9

date
```
```
qsub -cwd -pe threads 16 -l m_mem_free=2G generate_ccs.sh
```

#### 1. Classify full-length reads:
Script: lima.sh
```
#!/bin/bash
cd /sonas-hs/ware/hpc/home/diniz/FL_Pacbio

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load IntelPython/2.7.14
module load Anaconda2/5.3.0

source activate anaCogent5.2

date

/seq/smrt4/smrtlink/smrtcmds/bin/lima \
--isoseq --dump-clips --no-pbi --peek-guess -j 16 \
m54333U_200123_203332.ccs.bam \
IsoSeq_Primers_12_Barcodes_v1.fasta \
m54333U_200123_203332.demux.bam       

date
```
```
qsub -cwd -pe threads 16 -l m_mem_free=2G generate_ccs.sh
```
