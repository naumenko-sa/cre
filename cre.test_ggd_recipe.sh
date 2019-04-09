#!/bin/bash
date
export PYTHONPATH=/hpf/largeprojects/ccmbio/naumenko/tools/cloudbiolinux:$PYTHONPATH
python -c 'from cloudbio.biodata.ggd import install_recipe; install_recipe("/hpf/largeprojects/ccmbio/naumenko/test", "/hpf/largeprojects/ccmbio/naumenko/tools/bcbio", "gnomad_exome.yaml", "hg38")'
date