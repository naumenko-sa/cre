#!/bin/bash

# compute nodes do not have internet access, 
# qnodes do have internet access, but it is better to run on a data node, because nohups are dying on qlogins

# fresh install
# wget https://raw.github.com/bcbio/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
# python bcbio_nextgen_install.py /usr/local/share/bcbio --tooldir=/usr/local \
#  --genomes GRCh37 --aligners bwa --aligners bowtie2

#upgrade code to stable version
#bcbio_nextgen.py upgrade -u stable
#upgrade code to development
#bcbio_nextgen.py upgrade -u development

#upgrade tools
#bcbio_nextgen.py upgrade -u skip --tools

#check tools
#bcbio_conda list | grep vep

#bcbio_nextgen.py upgrade  --data
#upgrade data - it upgrades gemini, cadd, rnaseq if they were installed before

# VEP is upgraded quite often ~2-3 months - when upgrading tools it looks for new cache
bcbio_nextgen.py upgrade -u skip --genomes GRCh37 --datatarget vep

#gemini
#bcbio_nextgen.py upgrade -u skip --genomes GRCh37 --datatarget gemini

#cadd
#bcbio_nextgen.py upgrade -u skip --genomes GRCh37 --datatarget cadd
#--genomes hg38

#gnomad
#bcbio_nextgen.py upgrade -u skip --genomes GRCh37 --datatarget gnomad

#rnaseq
#bcbio_nextgen.py upgrade -u skip --genomes GRCh37 --datatarget rnaseq
