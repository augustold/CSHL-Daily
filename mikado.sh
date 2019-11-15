cd /sonas-hs/ware/hpc/home/diniz/mikado_2nd_test

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6  
module load TransDecoder/5.5.0-Perl-5.28.0

#Creating the configuration file for Mikado
mikado configure \
--list list.txt \
--reference /sonas-hs/ware/hpc/home/diniz/Saccharum_genome_refs/SP803280/sc.mlc.cns.sgl.utg.scga7.importdb.fa \
-t 20 \
--mode permissive \
--scoring plant.yaml  \
--copy-scoring plant.yaml \
--junctions portcullis_all.junctions.bed \
-bt uniprot_sprot_plants.fasta \
configuration.yaml

#Mikado prepare
mikado prepare --json-conf configuration.yaml

#BLAST of the candidate transcripts
makeblastdb -in uniprot_sprot_plants.fasta -dbtype prot -parse_seqids > blast_prepare.log

blastx -max_target_seqs 5 -num_threads 20 -query mikado_prepared.fasta -outfmt 5 -db uniprot_sprot_plants.fasta -evalue 0.000001 2> blast.log | sed '/^$/d' | gzip -c - > mikado.blast.xml.gz

TransDecoder.LongOrfs -t mikado_prepared.fasta
TransDecoder.Predict -t mikado_prepared.fasta

#Mikado serialise
mikado serialise --json-conf configuration.yaml --xml mikado.blast.xml.gz --orfs mikado_prepared.fasta.transdecoder.bed --blast_targets uniprot_sprot_plants.fasta --transcripts mikado_prepared.fasta

#Mikado pick
mikado pick --json-conf configuration.yaml --subloci-out mikado.subloci.gff3 --procs 20

