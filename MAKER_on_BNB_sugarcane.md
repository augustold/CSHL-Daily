# MAKER run on BNB

The below protocol describes the steps to run MAKER-P pipeline on CSHL BNB cluster. The protocol uses maize genome for annotation using MAKER-P

### Prerequisites:
* MAKER-P installed along with dependencies. The current installation is at `/sonas-hs/ware/hpc_norepl/data/programs/maker_v3_updated/maker`
* To make this simpler source the bash profile file `source /sonas-hs/ware/hpc/home/kchougul/.bash_profile`. This will set all the paths and environment variables for tools needed to run MAKER-P
* Genome fasta to annotate. Make sure the fasta headers are proper i.e there is no missing newline character. You can check if there are any issues using `grep "[ATGC]>[ATGC]*" testing.fasta`. Also if there are non ATGC characters present in the sequence please inspect or remove it using `sed -e '/^[^>]/s/[^ATGCatgc]/N/g' testing.fasta`
* Evidence data for aligning to genome (EST/mRNA/cDNA/proteins). A baseline evidence data used for b73 annotation can be found at `/sonas-hs/ware/hpc_norepl/data/kapeel/b73_NAM/standarized_evidence_set`


### Step 1. Create a working directory for genome and evidence files
Before running maker please check the avaible disk space. For a standard MAKER-P run using maize genome requires ~300Gb for a complete run. Check the disk space by loggin into filezone1 server and running the below command.
```
$ cd /sonas-hs/ware/hpc_norepl/data
$ griddf

                         Block Limits                                    |     File Limits
Filesystem type         blocks      quota      limit   in_doubt    grace |    files   quota    limit in_doubt    grace  Remarks
grid       FILESET      29.74T        32T        32T      42.5G     none | 42515845       0        0   137287     none 

```
You may choose to run MAKER-P on `/sonas-hs/ware/hpc/data` or `/sonas-hs/ware/hpc_norepl/data`. You should have a directory named for both the locations. If not contact Peter to create one. 

I will choose to create a directory under `/sonas-hs/ware/hpc_norepl/data`  

```
$ mkdir /sonas-hs/ware/hpc_norepl/data/diniz/MAKER_SP80
$ cd /sonas-hs/ware/hpc_norepl/data/diniz/MAKER_SP80
```

Geting evidence files

```
$ mkdir sugarcane_evidence_set
$ cd sugarcane_evidence_set

# ESTs from SUCEST project
cp /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/ESTs_SP80-3280.fasta .

# Full-lengths 
cp /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sugarcane.fulllength.analysis.all.fasta .

# Protien files from GRAMENE database

#S.spontaneum
wget http://www.life.illinois.edu/ming/downloads/Spontaneum_genome/Sspon.v20190103.protein.fasta.gz
#Sorghum
wget ftp://ftp.gramene.org/pub/gramene/release-62/fasta/sorghum_bicolor/pep/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.pep.all.fa.gz 
#Maize
wget ftp://ftp.gramene.org/pub/gramene/release-62/fasta/zea_mays/pep/Zea_mays.B73_RefGen_v4.pep.all.fa.gz
#Rice
wget ftp://ftp.gramene.org/pub/gramene/release-62/fasta/oryza_sativa/pep/Oryza_sativa.IRGSP-1.0.pep.all.fa.gz
#Brachypodium
wget ftp://ftp.gramene.org/pub/gramene/release-62/fasta/brachypodium_distachyon/pep/Brachypodium_distachyon.Brachypodium_distachyon_v3.0.pep.all.fa.gz
#Arabidopsis
wget ftp://ftp.gramene.org/pub/gramene/release-62/fasta/arabidopsis_thaliana/pep/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz
```

### Step 2. Splitting genome fasta file
A general observation if you have a assembly with scaffold/contig > 1000 no need of splitting continue to step ??? but if you have < 100 scaffold/contig split the genome fasta as below 
```
$ mkdir contig_fasta
$ cp <YOUR_GENOME_FILE> contig_fasta
$ contig_fasta
$ fasta_tool split cml333consensus_v2.fasta
```
This particular assembly has 53 contigs named Contig[0..52].fasta

```diff
! Sugarcane current assembly has > 400 K contigs. No need to run this step.
```

Instead, copy the assembly fasta here
```
$ mkdir contig_fasta
$ cd contig_fasta

cp /sonas-hs/ware/hpc_norepl/data/diniz/Saccharum_genome_refs/SP80-3280/sc.mlc.cns.sgl.utg.scga7.importdb.masked.fa .

```

### Step 3. Creating maker control files

`$ cd ..; maker -CTL`

this will create the following files 
* maker_opts.ctl - gives location of input files (genome and evidence) and sets options that affect MAKER behavior
* maker_exe.ctl - gives path information for the underlying executables.
* maker_bopt.ctl - sets parameters for filtering BLAST and Exonerate alignment results
* maker_evm.ctl - sets parameters adding weights to evidence modeller tool

Provide the file and tools paths within `maker_opts.ctl` and `maker_exe.ctl

Here is full control file:

maker_exe.ctl
```
#-----Location of Executables Used by MAKER/EVALUATOR
makeblastdb=/mnt/grid/ware/hpc_norepl/data/data/programs/maker_v3_updated/maker/bin/../exe/blast/bin/makeblastdb #location of NCBI+ makeblastdb executable
blastn=/mnt/grid/ware/hpc_norepl/data/data/programs/maker_v3_updated/maker/bin/../exe/blast/bin/blastn #location of NCBI+ blastn executable
blastx=/mnt/grid/ware/hpc_norepl/data/data/programs/maker_v3_updated/maker/bin/../exe/blast/bin/blastx #location of NCBI+ blastx executable
tblastx=/mnt/grid/ware/hpc_norepl/data/data/programs/maker_v3_updated/maker/bin/../exe/blast/bin/tblastx #location of NCBI+ tblastx executable
formatdb= #location of NCBI formatdb executable
blastall= #location of NCBI blastall executable
xdformat= #location of WUBLAST xdformat executable
blasta= #location of WUBLAST blasta executable
prerapsearch= #location of prerapsearch executable
rapsearch= #location of rapsearch executable
RepeatMasker=/mnt/grid/ware/hpc_norepl/data/data/programs/maker_v3_updated/maker/bin/../exe/RepeatMasker/RepeatMasker #location of RepeatMasker executable
exonerate=/mnt/grid/ware/hpc_norepl/data/data/programs/maker_v3_updated/maker/bin/../exe/exonerate/bin/exonerate #location of exonerate executable

#-----Ab-initio Gene Prediction Algorithms
snap=/mnt/grid/ware/hpc_norepl/data/data/programs/maker_v3_updated/maker/bin/../exe/snap/snap #location of snap executable
gmhmme3= #location of eukaryotic genemark executable
gmhmmp= #location of prokaryotic genemark executable
augustus=/sonas-hs/ware/hpc/home/mcampbel/applications/augustus-3.1/bin/augustus #location of augustus executable
fgenesh=/sonas-hs/ware/hpc_norepl/data/programs/fgenesh_v8/fgenesh_suite_v8.0.0a/fgenesh #location of fgenesh executable
evm= #location of EvidenceModeler executable
tRNAscan-SE=/sonas-hs/ware/hpc/home/mcampbel/applications/tRNAscan-SE-1.23/bin/tRNAscan-SE #location of trnascan executable
snoscan= #location of snoscan executable

#-----Other Algorithms
probuild= #location of probuild executable (required for genemark)
```
Note: most of the paths will be prefilled by MAKER when you run `maker -CTL` command. The path that needs to be added is for fgenesh tool


maker_opts.ctl
```
#-----Genome (these are always required)
genome=/sonas-hs/ware/hpc_norepl/data/diniz/MAKER_SP80/contig_fasta/sc.mlc.cns.sgl.utg.scga7.importdb.masked.fa #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est=/sonas-hs/ware/hpc_norepl/data/diniz/MAKER_SP80/sugarcane_evidence_set/ESTs_SP80-3280.fasta,/sonas-hs/ware/hpc_norepl/data/diniz/MAKER_SP80/sugarcane_evidence_set/sugarcane.fulllength.analysis.all.fasta #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=/sonas-hs/ware/hpc_norepl/data/diniz/MAKER_SP80/sugarcane_evidence_set/Sspon.v20190103.protein.fasta.gz,/sonas-hs/ware/hpc_norepl/data/diniz/MAKER_SP80/sugarcane_evidence_set/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.pep.all.fa.gz,/sonas-hs/ware/hpc_norepl/data/diniz/MAKER_SP80/sugarcane_evidence_set/Zea_mays.B73_RefGen_v4.pep.all.fa.gz,/sonas-hs/ware/hpc_norepl/data/diniz/MAKER_SP80/sugarcane_evidence_set/Oryza_sativa.IRGSP-1.0.pep.all.fa.gz,/sonas-hs/ware/hpc_norepl/data/diniz/MAKER_SP80/Brachypodium_distachyon.Brachypodium_distachyon_v3.0.pep.all.fa.gz,/sonas-hs/ware/hpc_norepl/data/diniz/MAKER_SP80/sugarcane_evidence_set/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz #protein sequence file in fasta format (i.e. from mutiple organisms)
protein_gff=  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org= #select a model organism for RepBase masking in RepeatMasker
rmlib=/sonas-hs/ware/hpc_norepl/data/kapeel/b73_NAM/standarized_evidence_set/Wessler-Bennetzen_2.fasta #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein=/sonas-hs/ware/hpc_norepl/data/programs/maker_v3_updated/maker/data/te_proteins.fasta #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff= #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm= #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species=sugarcane_sp80 #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
run_evm=0 #run EvidenceModeler, 1 = yes, 0 = no
est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=0 #infer predictions from protein homology, 1 = yes, 0 = no
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
snoscan_meth= #-O-methylation site fileto have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=1 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=1 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=1 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=1 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
min_intron=20 #minimum intron length (used for alignment polishing)
single_exon=1 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=1 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP=/sonas-hs/ware/hpc_norepl/data/diniz/temp #specify a directory other than the system default temporary directory for temporary files
```
Note: 
* Provide full paths somehoe MAKER-P on BNB does not recognize relative paths
* multiple proteins/EST files can be added separated by `,` and tagged using `:<tag>` at the end of file name
* the TMP location can be empty which will default to /tmp on the systems else choose a desired location with enough space
* CPU=1, this setting is desired if one is running MAKER-P in MPI mode

### Step 4. preparing qsub script 

To submit maker jobs we need to prepare a qsub script to schedule jobs. More details on how to submit a job on BNB see this documentation: http://intranet.cshl.edu/administration/information-technology/hpcc/uge.

For our maize run we will use MPI and submit parallel array job submission. Documentation on array jobs: https://github.com/warelab/annotation_workflows/wiki/maker_run_BNB_cluster

cat submit_maker_here.sh
```
#!/bin/bash
# specify BASH shell
#$ -S /bin/bash
# run job in the current working directory where qsub is executed from
#$ -cwd
# Tell the SGE that this is an array job, with "tasks" to be numbered 
#$ -t 1-53
#  specify that the job requires 3GB of memory for each process
#$ -pe mpi 16 -l m_mem_free=2G
# make sure that there is at least 20G of tmp space on the node
#$ -l tmp_free=200G

module load openmpi-x86_64
source /sonas-hs/ware/hpc/home/kchougul/.bash_profile
export LD_PRELOAD=/sonas-hs/ware/hpc/home/mcampbel/applications/openmpi-1.8.8/build/lib/libmpi.so

let i=$SGE_TASK_ID-1
/sonas-hs/ware/hpc/home/mcampbel/applications/openmpi-1.8.8/build/bin/mpiexec -mca btl ^openib -n 16 /sonas-hs/ware/hpc/home/mcampbel/applications/maker/bin/maker -g contig_fasta/Contig$SGE_TASK_ID.fasta -fix_nucleotides -TMP $TMPDIR
``` 
array job specification is represernted by -t option. Here we use 1-53 instead of 0-52. 

Note the -t parameter can be modified to submit 10 array jobs at a time to avoid overload on the system : 
`-t 1-53:10`

### Step 5. submitting qsub job and monitoring status of the job.

Submit qsub script
` $ qsub submit_maker_here.sh `

To check status of array jobs:
`qstat -u kchougul`

To check if MAKER finishes for each contig/array job there will 2 log file submit_maker_here.sh.e[0-9]* submit_maker_here.sh.oe[0-9]*. Once MAKER-P finishes check the status of submit_maker_here.sh.e[0-9]*  file. You should see the following message

```
$ tail -f submit_maker_here.sh.e42180609
adding statistics to annotations
Calculating annotation quality statistics
choosing best annotation set
Choosing best annotations
processing chunk output
processing contig output
Maker is now finished!!!
Start_time: 1509292679
End_time:   1509303084
Elapsed:    10405
```
You can also choose to see the maker output directory if the annotation gff and fasta files are been created.
` $ ls Contig*/*out*/*data*/*/*/*/*gff | wc -l `
This should match the number of contigs in the assembly. If for some reason a contig fails or gets stuck for more than 2 days kill the job using `qdel <Job-ID>` and resubmit by editing the qsub script and specifying the contig number under -t option.

### Step 6. Merge gff and fasta files generated from MAKER-P run

A typical set of outputs for a contig looks like this:
```
kchougul@bnbdev2:/sonas-hs/ware/hpc_norepl/data/diniz/MAKER_SP80$ ls Contig*/*out*/*data*/*/*/Contig0
Contig0.gff
Contig0.maker.augustus_masked.proteins.fasta
Contig0.maker.augustus_masked.transcripts.fasta
Contig0.maker.fgenesh_masked.proteins.fasta
Contig0.maker.fgenesh_masked.transcripts.fasta
Contig0.maker.non_overlapping_ab_initio.proteins.fasta
Contig0.maker.non_overlapping_ab_initio.transcripts.fasta
Contig0.maker.proteins.fasta
Contig0.maker.transcripts.fasta
run.log
theVoid.Contig0
```

To merge gff's and extract only gene features and remove nucleotide seq

`$ gff3_merge -g -n -o CML33.maker.gene_only.gff Contig*/*out*/*data*/*/*/*/*gff`

To merge gff's and extract all features and remove nucleotide seq

`$ gff3_merge -n -o CML33.maker.all.gff Contig*/*out*/*data*/*/*/*/*gff`

To merge fasta
```
$ cat Contig*/*out*/*data*/*/*/*/*maker.augustus_masked.proteins.fasta > CML333.maker.augustus_masked.proteins.fasta
$ cat Contig*/*out*/*data*/*/*/*/*maker.augustus_masked.transcripts.fasta > CML333.maker.augustus_masked.transcripts.fasta
$ cat Contig*/*out*/*data*/*/*/*/*maker.fgenesh_masked.proteins.fasta > CML333.maker.fgenesh_masked.proteins.fasta
$ cat Contig*/*out*/*data*/*/*/*/*maker.fgenesh_masked.transcripts.fasta > CML333.maker.fgenesh_masked.transcripts.fasta
$ cat Contig*/*out*/*data*/*/*/*/*maker.non_overlapping_ab_initio.proteins.fasta > CML333.maker.non_overlapping_ab_initio.proteins.fasta
$ cat Contig*/*out*/*data*/*/*/*/*maker.non_overlapping_ab_initio.transcripts.fasta > CML333.maker.non_overlapping_ab_initio.transcripts.fasta
$ cat Contig*/*out*/*data*/*/*/*/*maker.proteins.fasta > CML333.maker.proteins.fasta
$ cat Contig*/*out*/*data*/*/*/*/*maker.transcripts.fasta > CML333.maker.transcripts.fasta
```

To extract repeatmasked genome seq and gff
```
$ cat Contig*/*out*/*data*/*/*/*/*/*query.masked.fasta > CML333.maker.repeatmasked.fasta
$ cat Contig*/*out*/*data*/*/*/*/*/*query.masked.gff  > CML333.maker.repeatmasked.gff
```

### Step 7. move all the output results to final output directory

```
$ mkdir maker_final_annotations
$ mv *gff *fasta maker_final_annotations
```
