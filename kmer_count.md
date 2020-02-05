# Install Jellyfish

```
wget https://github.com/gmarcais/Jellyfish/releases/download/v2.3.0/jellyfish-2.3.0.tar.gz
tar -zxvf jellyfish-2.3.0.tar.gz ; rm -rf jellyfish-2.3.0.tar.gz
cd jellyfish-2.3.0

./configure --prefix=$HOME
make -j 4
make install
```

https://github.com/gmarcais/Jellyfish/tree/master/doc
### count all k-mers

```
cd /projects/augustold/CSHL/jellyfish

jellyfish count -m 21 -s 100M -t 10 -C \
/projects/augustold/CSHL/Saccharum_genome_refs/SP803280/sc.mlc.cns.sgl.utg.scga7.importdb.fa
```

### compute the histogram of the k-mer occurrences
```
jellyfish histo mer_counts.jf
```

### query the counts of a particular k-mer
```
jellyfish query mer_counts.jf AACGTTG
```

### output all the counts for all the k-mers in the file
```
jellyfish dump mer_counts.jf > mer_counts_dumps.fa
```
