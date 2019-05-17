#!/bin/bash
#PBS -l walltime=240:00:00,nodes=1:ppn=30
#PBS -joe .
#PBS -d .
#PBS -l vmem=100g,mem=100g

. /hpf/largeprojects/ccmbio/naumenko/tools/bcbio_1.1.5/.test_profile
which bcbio_nextgen.py

hostname

bcbio_path=/hpf/largeprojects/ccmbio/naumenko/tools/bcbio_1.1.5

# sometimes cannot upgrade STAR on data nodes - memory is low, or cannot take much CPUs
# can take more than a day
export PATH=${bcbio_path}/bin:$PATH && STAR \
--genomeDir ${bcbio_path}/genomes/Hsapiens/GRCh37/star \
--genomeFastaFiles ${bcbio_path}/genomes/Hsapiens/GRCh37/seq/GRCh37.fa \
--runThreadN 30 --limitGenomeGenerateRAM 30000000000 --genomeChrBinNbits 14 --runMode genomeGenerate --genomeSAindexNbases 14
