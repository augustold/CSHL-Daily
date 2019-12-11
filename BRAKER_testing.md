### Load modules

module load GCC/7.3.0-2.30
module load GCCcore/7.3.0
module load OpenMPI/3.1.1
module load BRAKER/2.1.2
module load AUGUSTUS/3.3
module load GeneMark-ET/4.46
module load BamTools/2.5.1
module load SAMtools/1.9
module load GenomeThreader/1.7.1-Linux_x86_64-64bit
module load BLAST+/2.9.0


Missing modules
Spaln 2.3.1 
Exonerate 2.2.0
DIAMOND
cdbfasta
cdbyank


########################
## BRAKER
########################

cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker

module load Anaconda2/5.3.0

conda create -n braker #only once
source activate braker

conda install -c anaconda perl
conda install -c bioconda perl-app-cpanminus
conda install -c bioconda perl-hash-merge
conda install -c bioconda perl-parallel-forkmanager
conda install -c bioconda perl-scalar-util-numeric
conda install -c bioconda perl-yaml
conda install -c bioconda perl-class-data-inheritable
conda install -c bioconda perl-exception-class
conda install -c bioconda perl-test-pod
conda install -c anaconda biopython
#conda install -c bioconda perl-file-homedir
conda install -c bioconda perl-file-which # skip if you are not comparing to reference annotation
cpanm Logger::Simple

module load GCCcore/7.3.0
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load git/2.19.1
module load BamTools/2.5.1
module load SAMtools/1.9
module load CMake/3.13.4
module load Python/3.6.6

git clone https://github.com/Gaius-Augustus/BRAKER.git

#Install GeneMark
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker
wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_Zid0r/gm_et_linux_64.tar.gz
tar xvzf gm_et_linux_64.tar.gz
cd gm_et_linux_64
wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_Zid0r/gm_key_64.gz
gunzip gm_key_64.gz
cp gm_key_64 ~/.gm_key

for f in bet_to_gff.pl bp_seq_select.pl build_mod.pl calc_introns_from_gtf.pl \
change_path_in_perl_scripts.pl gc_distr.pl get_sequence_from_GTF.pl \
gmes_petap.pl histogram.pl hmm_to_gtf.pl make_nt_freq_mat.pl \
parse_by_introns.pl parse_ET.pl parse_gibbs.pl parse_set.pl predict_genes.pl \
reformat_fasta.pl reformat_gff.pl rescale_gff.pl rnaseq_introns_to_gff.pl \
run_es.pl run_hmm_pbs.pl scan_for_bp.pl star_to_gff.pl verify_evidence_gmhmm.pl;
do
   cat $f | perl -pe 's/\/usr\/bin\/perl/\/usr\/bin\/env perl/' > $f.tmp
   mv $f.tmp $f
   chmod u+x $f
done

#Install Augustus
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker
git clone https://github.com/Gaius-Augustus/Augustus.git
cd Augustus/
make

#Install diamond
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker
conda install -c bioconda diamond
conda update diamond

#Install cdbfasta
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker
git clone https://github.com/gpertea/cdbfasta.git


#Example
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker/BRAKER/example
wget http://bioinf.uni-greifswald.de/bioinf/braker/RNAseq.bam


export GENEMARK_PATH=/sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker/gm_et_linux_64
export CDBTOOLS_PATH=/sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker/cdbfasta

braker.pl --genome=genome.fa --bam=RNAseq.bam --softmasking








PATH=/sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker/BRAKER/scripts:$PATH
export PATH

export AUGUSTUS_CONFIG_PATH=/sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker/Augustus/config
export AUGUSTUS_BIN_PATH=/sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker/Augustus/bin
export AUGUSTUS_SCRIPTS_PATH=/sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker/Augustus/scripts
export BAMTOOLS_BIN_PATH=/sonas-hs/it/hpc/home/easybuild/install_prod/software/MPI/GCC/7.3.0-2.30/OpenMPI/3.1.1/BamTools/2.5.1/bin
export SAMTOOLS_PATH=/sonas-hs/it/hpc/home/easybuild/install_prod/software/MPI/GCC/7.3.0-2.30/OpenMPI/3.1.1/SAMtools/1.9/bin
export DIAMOND_PATH=/sonas-hs/ware/hpc/home/bwang/software/diamond
export PYTHON3_PATH=/sonas-hs/it/hpc/home/easybuild/install_prod/software/MPI/GCC/7.3.0-2.30/OpenMPI/3.1.1/Python/3.6.6/bin









########################
## BRAKER 2
########################

cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker2

module load Anaconda2/5.3.0

#conda create -n braker2-2.1.2 -c conda-forge -c bioconda -c defaults braker2=2.1.2=2
conda create -n braker2 -c conda-forge -c bioconda -c defaults braker2

# To activate this environment, use:
source activate braker2

which perl
which braker.pl
ls ~/.conda/envs/braker2/bin/helpMod.pm

# To deactivate an active environment, use:
# conda deactivate

module load GCCcore/7.3.0
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load git/2.19.1

#Install Braker
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker2
git clone https://github.com/Gaius-Augustus/BRAKER.git

#Install GeneMark
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker2
wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_Zid0r/gm_et_linux_64.tar.gz
tar xvzf gm_et_linux_64.tar.gz
cd gm_et_linux_64
wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_Zid0r/gm_key_64.gz
gunzip gm_key_64.gz
cp gm_key_64 .gm_key

for f in bet_to_gff.pl bp_seq_select.pl build_mod.pl calc_introns_from_gtf.pl \
change_path_in_perl_scripts.pl gc_distr.pl get_sequence_from_GTF.pl \
gmes_petap.pl histogram.pl hmm_to_gtf.pl make_nt_freq_mat.pl \
parse_by_introns.pl parse_ET.pl parse_gibbs.pl parse_set.pl predict_genes.pl \
reformat_fasta.pl reformat_gff.pl rescale_gff.pl rnaseq_introns_to_gff.pl \
run_es.pl run_hmm_pbs.pl scan_for_bp.pl star_to_gff.pl verify_evidence_gmhmm.pl;
do
   cat $f | perl -pe 's/\/usr\/bin\/perl/\/usr\/bin\/env perl/' > $f.tmp
   mv $f.tmp $f
   chmod u+x $f
done

#Install Augustus
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker2
git clone https://github.com/Gaius-Augustus/Augustus.git
cd Augustus/
make


#Install diamond
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker2
conda install -c bioconda diamond
conda update diamond

#Install cdbfasta
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker2
git clone https://github.com/gpertea/cdbfasta.git
cd cdbfasta
make all

########################
# Wed Dec 11 18:10:17 2019: ERROR: in file /sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker2/BRAKER/scripts/braker.pl at line 3153
Perl module 'Logger::Simple' is required but not installed yet.
########################

#trying
conda install -c bioconda perl-logger-simple

########################
==> WARNING: A newer version of conda exists. <==
  current version: 4.5.11
  latest version: 4.8.0

Please update conda by running

    $ conda update -n base -c defaults conda
########################

#trying
conda update -n base -c defaults conda

#Peter's hint
export PERL5LIB=/sonas-hs/ware/hpc/home/diniz/perl5/lib/perl5:$PERL5LIB
perl -MLogger::Simple -e 1
#no output means it is ok

#Example
module load Anaconda2/5.3.0

source activate braker2

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load Python/3.6.6

cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker2/BRAKER/example

wget http://bioinf.uni-greifswald.de/bioinf/braker/RNAseq.bam

source activate braker2
cd /sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker2/BRAKER/example

export CDBTOOLS_PATH=/sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker2/cdbfasta
export GENEMARK_PATH=/sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker2/gm_et_linux_64

/sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker2/BRAKER/scripts/braker.pl \
--genome=genome.fa \
--bam=RNAseq.bam \
--softmasking \
--AUGUSTUS_CONFIG_PATH=/sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker2/Augustus/config \
--AUGUSTUS_BIN_PATH=/sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker2/Augustus/bin \
--AUGUSTUS_SCRIPTS_PATH=/sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker2/Augustus/scripts \
--CDBTOOLS_PATH=/sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker2/cdbfasta \
--GENEMARK_PATH=/sonas-hs/ware/hpc_norepl/data/diniz/analysis/braker2/gm_et_linux_64






########################
## BRAKER 2 on Mac
########################

cd /Users/diniz/braker2
#conda create -n braker2 -c conda-forge -c bioconda -c defaults braker2

# To activate this environment, use:
source activate braker2

which perl
which braker.pl
ls ~/miniconda3/envs/braker2/bin/helpMod.pm
#testing braker
braker.pl

#Install GeneMark
cd /Users/diniz/braker2
curl http://topaz.gatech.edu/GeneMark/tmp/GMtool_aNOK4/gm_et_macosx.tar.gz --output gm_et_macosx.tar.gz
tar xvzf gm_et_macosx.tar.gz
rm -rf gm_et_macosx.tar.gz
cd gm_et_macosx/gmes_petap
curl http://topaz.gatech.edu/GeneMark/tmp/GMtool_aNOK4/gm_key_64.gz --output gm_key_64.gz
gunzip gm_key_64.gz
cp gm_key_64 ~/.gm_key

#Install diamond
conda install -c bioconda diamond
conda update diamond

#Install cdbfasta
git clone https://github.com/gpertea/cdbfasta.git

# git clone braker
git clone https://github.com/Gaius-Augustus/BRAKER.git
cd cdbfasta/
make

# example
cd /Users/diniz/braker2/BRAKER/example
curl http://bioinf.uni-greifswald.de/bioinf/braker/RNAseq.bam --output RNAseq.bam

export GENEMARK_PATH=/Users/diniz/braker2/gm_et_macosx/gmes_petap
export CDBTOOLS_PATH=/Users/diniz/braker2/cdbfasta

braker.pl \
--genome=genome.fa \
--bam=RNAseq.bam \
--softmasking


braker.pl --genome=genome.fa --softmasking --esmode















export DIAMOND_PATH=/sonas-hs/ware/hpc/home/bwang/software/diamond













########################
## BRAKER 2 on helix
########################

cd /projects/augustold/CSHL/braker2

#conda create -n braker2 -c conda-forge -c bioconda -c defaults braker2
conda activate braker2

which perl
which braker.pl
ls ~/./miniconda3/envs/braker2/bin/helpMod.pm

# To deactivate an active environment, use:
# conda deactivate

#Install Braker
cd /projects/augustold/CSHL/braker2
git clone https://github.com/Gaius-Augustus/BRAKER.git

#Install GeneMark
cd /projects/augustold/CSHL/braker2
wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_Zid0r/gm_et_linux_64.tar.gz
tar xvzf gm_et_linux_64.tar.gz
rm -rf gm_et_linux_64.tar.gz 
cd gm_et_linux_64
wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_Zid0r/gm_key_64.gz
gunzip gm_key_64.gz
cp gm_key_64 ~/.gm_key

for f in bet_to_gff.pl bp_seq_select.pl build_mod.pl calc_introns_from_gtf.pl \
change_path_in_perl_scripts.pl gc_distr.pl get_sequence_from_GTF.pl \
gmes_petap.pl histogram.pl hmm_to_gtf.pl make_nt_freq_mat.pl \
parse_by_introns.pl parse_ET.pl parse_gibbs.pl parse_set.pl predict_genes.pl \
reformat_fasta.pl reformat_gff.pl rescale_gff.pl rnaseq_introns_to_gff.pl \
run_es.pl run_hmm_pbs.pl scan_for_bp.pl star_to_gff.pl verify_evidence_gmhmm.pl;
do
   cat $f | perl -pe 's/\/usr\/bin\/perl/\/usr\/bin\/env perl/' > $f.tmp
   mv $f.tmp $f
   chmod u+x $f
done

#Install diamond
conda install -c bioconda diamond
conda update diamond

#Install cdbfasta
git clone https://github.com/gpertea/cdbfasta.git
cd cdbfasta
make all


#Example

source activate braker2

cd /projects/augustold/CSHL/braker2/BRAKER/example
#wget http://bioinf.uni-greifswald.de/bioinf/braker/RNAseq.bam

export GENEMARK_PATH=/projects/augustold/CSHL/braker2/gm_et_linux_64
export CDBTOOLS_PATH=/projects/augustold/CSHL/braker2/cdbfasta


braker.pl \
--genome=genome.fa \
--bam=RNAseq.bam \
--softmasking

