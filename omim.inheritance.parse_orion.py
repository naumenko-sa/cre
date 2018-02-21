#!/usr/bin/python

# extracts information about omim inheritance modes from https://www.cs.toronto.edu/~buske/cheo/

import re
from os.path import expanduser

home = expanduser('~')

inheritance = {}
inheritance['Autosomal recessive']='AR'
inheritance['Autosomal dominant'] = 'AD'
inheritance['X-linked recessive'] = 'XLR'
inheritance['X-linked dominant'] = 'XLD'
inheritance['Isolated cases'] = 'IC'
inheritance['Mitochondrial'] = 'Mi'
inheritance['X-linked'] = 'XL'
inheritance['Y-linked'] = 'YL'
inheritance['Digenic recessive'] = 'DR'
inheritance['Digenic dominant'] = 'DD'
inheritance['Multifactorial'] = 'Mu'
inheritance['Somatic mosaicism'] = 'Smo'
inheritance['Somatic mutation'] = 'Smu'

genes = {}

f1 = open(home+'/cre/ensembl_w_description.txt','r')
for line in f1:
    ar = line.split('\t')
    genes[ar[0]] = ar[1]
f1.close()

f = open('genemap2.txt','r')

print 'Ensembl_gene_id Gene_name2 Omim_inheritance'

for line in f:
    #print line,
    match = re.search(r'ENSG[0-9]{11}',line)
    if match:
	modes = []
	#print 'found', match.group()
	for key in inheritance.keys():
	    match1 = re.search(key,line)
	    if match1:
		#print 'found inheritance', match1.group()
		modes.append(inheritance[key])
	if len(modes)>0:
	    print match.group(), genes.get(match.group()), ','.join(modes)

f.close()
