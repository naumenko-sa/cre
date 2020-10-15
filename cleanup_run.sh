#!/bin/bash

#PBS -l walltime=01:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=5g,mem=5g

# script to move the regenerate folder files up one level after completed
# deletes empty family subfolder

# parameters: family = familyid

# input: familyid 

if [ -z "$family" ]
then
	family=$1
fi

mv "$family"/* .
rm -rf "$family"
