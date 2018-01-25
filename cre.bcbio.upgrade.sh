#!/bin/bash

#compute nodes do not have internet access, 
#qnodes do have internet access, but it is better to run on a data node

#upgrade code to stable version
#bcbio_nextgen.py upgrade -u stable

#upgrade tools
#bcbio_nextgen.py upgrade -u skip --tools

#check tools
#bcbio_conda list | grep vep

#upgrade data - it upgrades gemini, cadd, rnaseq if they were installed before
#bcbio_nextgen.py upgrade -u skip --genomes GRCh37 --datatarget vep

#gemini
#bcbio_nextgen.py upgrade -u skip --genomes GRCh37 --datatarget gemini

#cadd
#bcbio_nextgen.py upgrade -u skip --genomes GRCh37 --datatarget cadd
#--genomes hg38

#rnaseq
#bcbio_nextgen.py upgrade -u skip --genomes GRCh37 --datatarget rnaseq
