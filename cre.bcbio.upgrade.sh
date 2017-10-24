#!/bin/bash

#compute nodes don't have internet access

#upgrade code to stable version
#bcbio_nextgen.py upgrade -u stable

#upgrade tools
bcbio_nextgen.py upgrade -u skip --tools

#check tools
#bcbio_conda list | grep vep

#upgrade data
#bcbio_nextgen.py upgrade -u skip --genomes GRCh37 --datatarget vep

#--tools 
#--data \
#--genomes GRCh37 \
#--datatarget gemini \
#--datatarget cadd
#--genomes hg38
#--datatarget rnaseq --datatarget gemini --datatarget cadd --genomes hg38
