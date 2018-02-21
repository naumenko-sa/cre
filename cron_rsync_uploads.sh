#!/bin/bash

# push files from the uploader to hpf
# ssh cron1

echo "$(date)"

magic_path=/hpf/tools/centos6/python/2.7.11

export PATH=$magic_path/bin:$PATH
export LD_LIBRARY_PATH=$magic_path/lib/:$LD_LIBRARY_PATH
export PYTHONPATH=$magic_path:$PYTHONPATH
export LD_LIBRARY_PATH=$magic_path/usr/lib64/atlas:$LD_LIBRARY_PATH

python /hpf/largeprojects/ccm_dccforge/dccforge/cron/cron_rsync_uploads.py
