#!/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/anaconda/bin/python

"""
Looks for a specific sample
"""

import re
import sys
import os
import os.path

sample = sys.argv[1]

family,sample_only = sample.split("_")

match = re.match('\d*',family)

if match:
    prefix=str(int(match.group(0))/100)
    report_path = prefix+'x/'+family
    
    report=0
    bam=0
    
    errors = []
    
    if os.path.isfile(report_path+'/'+family+'.csv'):
	#print("Report exists")
	report=1
    else:
	errors.append('Error: no report')
	
    if os.path.isfile(report_path+'/'+sample+'.bam'):
	#print("Bam exists")
	bam=1
    else:
	errors.append(' ERROR: no bam')
	
    if (bam==1 and report==1):
        print(sample+'\t'+os.getcwd()+"/"+report_path+"\t"+os.getcwd()+"/"+report_path+'/'+sample+'.bam')
    else:
	print(sample+'\t'+' '.join(errors))
else:
    print("Family ID is not starting with digital")
