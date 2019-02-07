#!/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/anaconda/bin/python

import csv
import sys
import re

family = sys.argv[1]
report = sys.argv[2]

with open('samples.txt', 'rb') as f_samples:
    samples = f_samples.readlines()
samples = [x.strip() for x in samples]
    
n_samples = len(samples)

for sample in samples:
    with open(report, 'rb') as f_csv:
	reader = csv.DictReader(f_csv)
	with open(family + '.' + sample + '.c4r', 'w') as f_sample:
	    zygosity_field = 'Zygosity.' + sample
	    #because R when building report substitutes - with _ in column names
	    zygosity_field = zygosity_field.replace("-", "_")
	    fieldnames = ['Position', 'Ref', 'Alt', 'Variation', zygosity_field, 'Refseq_change',
			  'Gene', 'Conserved_in_20_mammals', 'Sift_score', 'Polyphen_score', 'Cadd_score', 'Gnomad_af']
	    for row in reader:
		l = []
		for key in fieldnames:
		    l.append(row[key])
		zygosity = row[zygosity_field]
		if (zygosity != '-' and not re.search('Insufficient', zygosity)):
		    f_sample.write(','.join(l)+','+sample)
		    f_sample.write('\n')

f_sample.close()
f_samples.close()
f_csv.close()
