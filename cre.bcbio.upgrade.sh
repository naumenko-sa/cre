#!/bin/bash
#PBS -l walltime=23:00:00,nodes=1:ppn=10
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

# nohups are dying on qlogin nodes, data nodes are better for long data installation runs
#######################################################################
# fresh install of the second bcbio instance:
# mv ~/.conda/environments.txt ~/.conda/environments.default.txt - move back
# export PATH=/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/bin
# export PYTHONPATH=
# which python
# wget https://raw.github.com/bcbio/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py

# echo "Installing to " $1
# install data later
# module load python/2.7.15
# which bash
# python bcbio_nextgen_install.py $1 --tooldir $1 --genomes GRCh37 --aligners bwa --isolate --nodata
########################################################################
# fresh installation for Sam with human and mouse genome

# to check what enviroments were picked up during the installation
# conda info --envs --json
# check file ~/.conda/environments.txt - if it has environments from all installations they could interfere
# wget https://raw.github.com/bcbio/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
# export PYTHONPATH=/hpf/largeprojects/lauryl/bcbio110/anaconda/lib/python2.7

# PATH=/hpf/largeprojects/lauryl/bcbio110/anaconda/bin
# PATH=${PATH}:/usr/local/bin:/opt/moab/bin:/home/naumenko/cre:/home/naumenko/crt:/home/naumenko/crg:/home/naumenko/tools/mc-4.8.16/bin:/home/naumenko/jkent_tools
# PATH=${PATH}:/home/naumenko/bioscripts:.:/home/naumenko/.aspera/connect/bin:/usr/local/bin:/usr/lib64/qt-3.3/bin:/opt/moab/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin
# PATH=${PATH}:/opt/ibutils/bin:/sbin:/usr/sbin:/sbin:/usr/sbin
# export PATH

# echo $PATH
# echo $PYTHONPATH
# which python

# python bcbio_nextgen_install.py /hpf/largeprojects/lauryl/bcbio110 --tooldir=/hpf/largeprojects/lauryl/bcbio110 --genomes mm10 --aligners bwa --isolate
# bcbio_nextgen.py upgrade -u skip --tools --tooldir /hpf/largeprojects/lauryl/bcbio110
# bcbio_nextgen.py upgrade -u skip --data --genomes mm10 --datatarget variation --datatarget vep
#########################################################################
# upgrade code to the latest stable version
# bcbio_nextgen.py upgrade -u stable
# upgrade code to development
# bcbio_nextgen.py upgrade -u development
# upgrade tools
# bcbio_nextgen.py upgrade -u skip --tools
#########################################################################
. /hpf/largeprojects/ccmbio/naumenko/tools/bcbio_testing/.test_profile
which bcbio_nextgen.py
bcbio_nextgen.py upgrade -u skip --genomes GRCh37 --aligners star --cores 10
#--aligners hisat2 --aligners rtg
#########################################################################
# upgrades gemini, cadd, rnaseq if they were installed before, for all references
# bcbio_nextgen.py upgrade --data

# VEP is upgraded quite often ~2-3 months - when upgrading tools it looks for new VEP cache
# bcbio_nextgen.py upgrade -u skip --genomes GRCh37 --datatarget vep

# gemini
# bcbio_nextgen.py upgrade -u skip --genomes GRCh37 --datatarget gemini
# --genomes hg38

# cadd
# bcbio_nextgen.py upgrade -u skip --genomes GRCh37 --datatarget cadd
# --genomes hg38

# gnomad
# bcbio_nextgen.py upgrade -u skip --genomes GRCh37 --datatarget gnomad
# bcbio_nextgen.py upgrade -u skip --genomes hg38 --datatarget gnomad

# rnaseq
# bcbio_nextgen.py upgrade -u skip --genomes GRCh37 --datatarget rnaseq

# mouse
# bcbio_nextgen.py upgrade -u skip --genomes mm10 --datatarget rnaseq --cores 5
