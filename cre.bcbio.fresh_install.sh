#!/bin/bash

# fresh install with human and mouse genome

#wget https://raw.github.com/bcbio/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
export PYTHONPATH=
export PATH=/hpf/largeprojects/lauryl/bcbio110/anaconda/bin:/usr/local/bin:/opt/moab/bin:/home/naumenko/cre:/home/naumenko/crt:/home/naumenko/crg:/home/naumenko/tools/mc-4.8.16/bin:/home/naumenko/jkent_tools:/home/naumenko/bioscripts:.:/home/naumenko/.aspera/connect/bin:/usr/local/bin:/usr/lib64/qt-3.3/bin:/opt/moab/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/ibutils/bin:/sbin:/usr/sbin:/sbin:/usr/sbin

python bcbio_nextgen_install.py /hpf/largeprojects/lauryl/bcbio110 \
				  -u stable \
				  --genomes mm10 --aligners bwa --cores 8 --tooldir /hpf/largeprojects/lauryl/bcbio110 --isolate
				  
				  
#bcbio_nextgen.py upgrade -u skip --tools --tooldir /hpf/largeprojects/lauryl/bcbio110

