#Pfam Scan

### Install

```
conda create -n pfam
source activate pfam

conda install -c bioconda pfam_scan

mkdir Pfam; cd Pfam
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm

```

### Run
```
cd /projects/augustold/CSHL/mikado_SP80
source activate pfam
hmmsearch --tblout out.txt -E 1e-5 --cpu 20 \
/home/augustold/Pfam/Pfam-A.hmm \
uniprot_sprot_plants.fasta

pfam_scan.pl -fasta uniprot_sprot_plants.fasta -dir /projects/ftencaten/database/pfam
```
