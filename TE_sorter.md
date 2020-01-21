# Filtering MIKADO transcripts for TE related genes

## Installing TEsorter

This step was run on 'Helix'.
```
conda create --name py2 python=2.7

conda activate py2

pip install biopython
pip install pp
conda install -c biocore hmmer
#conda install -c asmeurer glibc
#conda install -c dan_blanchard glibc
conda install -c bioconda tesorter

cd /projects/augustold/CSH
git clone https://github.com/zhangrengang/TEsorter
cd TEsorter
sh build_database.sh 
```

## Testing TEsorter

```
/projects/augustold/CSH/TEsorter/test
# run
python ../TEsorter.py rice6.9.5.liban
```

## Running TEsorter on SP80-3280 MIKADO transcripts

```
cd /projects/augustold/CSH/TEsorter
mkdir SP80; cd SP80
#copy mikado.loci.cds.fasta from BNB

python ../TEsorter.py mikado.loci.cds.fasta -eval 1e-6 -p 30
```

## Filtering MIKADO gff3

```
cd /projects/augustold/CSH/TEsorter/SP80
#copy mikado.loci.gff3 from BNB

grep -v "^#" mikado.loci.cds.fasta.rexdb.cls.tsv | cut -f 1 | sort | uniq | cut -d "." -f 1-2 | sort |uniq > TE-genes.txt
grep -Fwvf TE-genes.txt mikado.loci.gff3 | awk '$3 !~ /superlocus/' > temp.gff3
grep -Fwf TE-genes.txt mikado.loci.gff3 | awk '$3 !~ /superlocus/' > temp.other.gff3

gt gff3 -sort -tidy -retainids temp.gff3  > mikado.loci.TErmv.gff3
gt gff3 -sort -tidy -retainids temp.other.gff3  > mikado.loci.maker.withTE.gff3
```

## Get summary stats from filtered gff3

```
cd /projects/augustold/CSH/TEsorter/SP80

/projects/augustold/CSHL/gff2gene_stats_canonical.pl mikado.loci.TErmv.gff3
```
