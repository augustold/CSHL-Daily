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

## Step 3 EMBOSS

Identify ORF for each isoform, select the longest ORF, and define the 5'UTR, CDS, 3'UTR of.

### Find ORF for each isoform
```
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load EMBOSS/6.6.0

getorf -sequence pacbio_hq_transcripts.fasta -outseq pacbio_hq_transcripts_orf_nc_3_NOreverse.fa -find 3 -reverse N

```
### Identify candidate NMDs

script: find_longest_ORF_and_NMDs_seqId_singleExonGeneNA.py
```
#!/usr/bin/env python 
# Created by: Xiaofei Wang
# Date: 06/30/2020
# Descript: find the longest orf for each isoform and find the distance to the last exon junction from stop codon
# Inputs: fasta and annotation file (e.g. gff3)
# Out: a fasta file for longest ORF, a TXT file for longest ORF information 
# (id, length, lenght of exons exclude last exon), a TXT for potential NMDs 
# (id, distance between stop codon and last exon junciton)

# Usage: python ******.py input.fa input.gff3 out.fa out_info.txt NMD.txt
# Example: python find_longest_ORF_and_NMDs_seqId_singleExonGeneNA.py test.fa t2.gff3 t2_out.fa t2_longestORFs.txt t2NMDs.txt
#----------------------------------------------------------------------------------------
#===========================================================================================================
#Imports:

import sys
import os
import pyfasta
from operator import itemgetter, attrgetter
import gffutils
#===========================================================================================================
# Functions:


# posInput = open(sys.argv[1], 'rU')
# scaffoldFa = 'test.fa'
# scaffoldFa = 'CS_IsoSeq.CISIs_orf_nc_3_NOreverse.fa'
scaffoldFa = str(sys.argv[1])
anno = str(sys.argv[2])

print(scaffoldFa)
print(anno)


outFa = open(sys.argv[3], 'w')
outLongestORFs = open(sys.argv[4], 'w')
outNMDs = open(sys.argv[5], 'w')

# start = int(sys.argv[2])
# end = int(sys.argv[3])

orfDict={}
longestORFkeys=[]
longestORFlengthDict={}
NMDs=[]
# 
# attrList = []
# attrList.append('KOS_ref_%s_%s' %(start, end))

# print attrList

# seqList = []

# fastaSeq = pyfasta.Fasta(scaffoldFa, key_fn=lambda key: key.split('|')[3])
fastaSeq = pyfasta.Fasta(scaffoldFa)

print(len(fastaSeq.keys()))
print("The length of the sequence is", len(str(fastaSeq[sorted(fastaSeq.keys())[0]])))





for item in sorted(fastaSeq.keys()):
	# print item
	orfDictKey=item.split()[0].rsplit('_',1)[0]
	# print orfDictKey
	if orfDictKey not in orfDict:
		orfDict[orfDictKey]=[]
	orfDict[orfDictKey].append([len(str(fastaSeq[item])),item])

# print orfDict
# print len(orfDict)
print(len(orfDict))


for item in orfDict:
	# print item
	# print max([x[0] for x in orfDict[item]]) 
	longestORF=max(orfDict[item])
	longestORFkeys.append(longestORF[1])
	longestORFlengthDict[item]=[longestORF[0]]

# print longestORFkeys
print(len(longestORFkeys))
# print(longestORFlengthDict)

for i in longestORFkeys:
	# print i
	# print fastaSeq[i]
	outFa.write('>'+i+'\n')
	seq=fastaSeq[i]
	outFa.write(seq[:]+'\n')


# db3 = gffutils.create_db("CS_IsoSeq.CISIs.gff3", dbfn='test3.db', force=True, keep_order=True,merge_strategy='merge', sort_attribute_values=True)
# db3 = gffutils.create_db("test.gff3", dbfn='test3.db', force=True, keep_order=True,merge_strategy='merge', sort_attribute_values=True)
db3 = gffutils.create_db(anno, dbfn='test3.db', force=True, keep_order=True,merge_strategy='merge', sort_attribute_values=True)

# sorted(exons, key= lambda exon: int(exon.start))
for i in db3.features_of_type('gene'):
	exons = list(db3.children(i, featuretype='exon'))
	# print exons
	if i.strand == '+':
		sortedExons=sorted(exons, key=attrgetter('start'))
	elif i.strand == '-':
		sortedExons=sorted(exons, key=attrgetter('start'), reverse=True) # it should be end coordiates conceptly, but start is also fine
	# print sortedExons
	lengths=[len(j) for j in sortedExons]
	# print lengths
	lastExonJunc=sum(lengths[:-1])
	# print i
	geneName=i.attributes['Name']
	# print(geneName)
	# print(geneName[0])
	if geneName[0] in longestORFlengthDict and len(lengths)>1:
		longestORFlengthDict[geneName[0]].append(str(i.seqid))
		longestORFlengthDict[geneName[0]].append(lastExonJunc)
		stopToLastExonDis=lastExonJunc-longestORFlengthDict[geneName[0]][0]
		# print(stopToLastExonDis)
		# print(longestORFlengthDict[geneName[0]])
		longestORFlengthDict[geneName[0]].append(stopToLastExonDis)
		if stopToLastExonDis > 55:
			NMDs.append([geneName[0], str(stopToLastExonDis),str(i.seqid)])
	elif geneName[0] in longestORFlengthDict and len(lengths)==1:
		longestORFlengthDict[geneName[0]].append(str(i.seqid))
		longestORFlengthDict[geneName[0]].append("NA")
		longestORFlengthDict[geneName[0]].append("NA")

	else:
		print("%s not in longest ORF of input fasta" %(geneName[0],))

# print(longestORFlengthDict)
print("The length of the NMDs is", len(NMDs))
# print(NMDs)

for item in NMDs:
	# print(item)
	outNMDs.write("\t".join(item)+'\n')
# 
for item in longestORFlengthDict:
	outLongestORFs.write(item+"\t"+"\t".join(str(i) for i in (longestORFlengthDict[item]))+"\n") 
	# print("%s"+"\n") %(item,)



print("Done!!!")
```

script: run_find_longest_ORF_and_NMDs.sh
```
#!/bin/bash
cd /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load EMBOSS/6.6.0
module load Anaconda2/5.3.0
#conda create -n gffutils
source activate gffutils
#conda install -c bioconda gffutils
#conda install -c bioconda pyfasta

date
python find_longest_ORF_and_NMDs.py longestORF_FL_SP80.fa LongestORFinfo_FL_SP80.txt NMD_FL_SP80.txt
date
```
```
qsub -cwd -pe threads 1 -l m_mem_free=16G run_find_longest_ORF_and_NMDs.sh
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
