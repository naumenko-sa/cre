#!/bin/bash

module load python/3.5.2
python3 ~/cre/pc.upload_vcf.py $1 $2  ~/work/project_cheo/pc.txt phenomecentral.org
