#!/bin/bash
#PBS -l walltime=240:00:00,nodes=1:ppn=10
#PBS -joe .
#PBS -d .
#PBS -l vmem=31g,mem=31g

. /hpf/largeprojects/ccmbio/naumenko/tools/bcbio_testing/.test_profile
which bcbio_nextgen.py

hostname

# sometimes cannot upgrade STAR on data nodes - memory is low, or cannot take much CPUs
# can take more than a day
export PATH=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio_testing/bin:$PATH && STAR \
--genomeDir /hpf/largeprojects/ccmbio/naumenko/tools/bcbio_testing/genomes/Hsapiens/GRCh37/star \
--genomeFastaFiles /hpf/largeprojects/ccmbio/naumenko/tools/bcbio_testing/genomes/Hsapiens/GRCh37/seq/GRCh37.fa \
--runThreadN 10 --limitGenomeGenerateRAM 30000000000 --genomeChrBinNbits 14 --runMode genomeGenerate --genomeSAindexNbases 14
