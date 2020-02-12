#Pfam Scan

### Install

```
conda create -n pfam
source activate pfam

conda install -c bioconda pfam_scan

cd /home/augustold/Pfam
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
```

### Run
```
cd /projects/augustold/CSHL/mikado_SP80
source activate pfam
hmmsearch --tblout mikado.protein.pfam.out.txt -E 1e-5 --cpu 20 \
/home/augustold/Pfam/Pfam-A.hmm \
mikado.loci.protein_newheader.v3.fasta
```

### Get list of IDs
```
cd /projects/augustold/CSHL/mikado_SP80
grep -v "^#" mikado.protein.pfam.out.txt | awk -F" " '{print $1}' | sort | uniq | wc -l
```
