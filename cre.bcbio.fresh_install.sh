#!/bin/bash
#fresh install with human and mouse genome

# to check what enviroments were picked up during the installation
#conda info --envs --json
# check file ~/.conda/environments.txt - it has environment from all installations - they can interfere

#wget https://raw.github.com/bcbio/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
export PYTHONPATH=/hpf/largeprojects/lauryl/bcbio110/anaconda/lib/python2.7

PATH=/hpf/largeprojects/lauryl/bcbio110/anaconda/bin
PATH=${PATH}:/usr/local/bin:/opt/moab/bin:/home/naumenko/cre:/home/naumenko/crt:/home/naumenko/crg:/home/naumenko/tools/mc-4.8.16/bin:/home/naumenko/jkent_tools
PATH=${PATH}:/home/naumenko/bioscripts:.:/home/naumenko/.aspera/connect/bin:/usr/local/bin:/usr/lib64/qt-3.3/bin:/opt/moab/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin
PATH=${PATH}:/opt/ibutils/bin:/sbin:/usr/sbin:/sbin:/usr/sbin

export PATH

echo $PATH
echo $PYTHONPATH
which python

#python bcbio_nextgen_install.py /hpf/largeprojects/lauryl/bcbio110 --tooldir=/hpf/largeprojects/lauryl/bcbio110 --genomes mm10 --aligners bwa --isolate
				  
#bcbio_nextgen.py upgrade -u skip --tools --tooldir /hpf/largeprojects/lauryl/bcbio110
#bcbio_nextgen.py upgrade -u skip --data --genomes mm10 --datatarget variation --datatarget vep
