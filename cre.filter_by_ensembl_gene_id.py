#!/bin/env python

# reports variants only in a gene panel using csv report on rare variants as input

import csv
import sys
import re

original_report = sys.argv[1]
panel_report = original_report.replace("csv","panel.csv")

gene_panel_file_name = sys.argv[2]
with open(gene_panel_file_name,'rb') as f_gene_panel:
    genes = f_gene_panel.readlines()
f_gene_panel.close()

genes = [x.strip() for x in genes]


with open(original_report,'rb') as f_original_report:
    header_line = f_original_report.readline().strip().replace('"','')
    header = header_line.split(',')
f_original_report.close()

with open(original_report,'rb') as f_original_report:
	reader = csv.DictReader(f_original_report)
	with open(panel_report,'w') as f_panel_report:

	    f_panel_report.write('"'+'","'.join(header)+'"')
	    f_panel_report.write('\n')
	    for row in reader:
		if (row['Ensembl_gene_id'] in genes):
		    row['UCSC_Link'] = 'UCSC_Link' #a bug with quotes
		    values = []
		    for column in header:
			values.append(row[column])
		    f_panel_report.write('"'+'","'.join(values)+'"')
		    f_panel_report.write('\n')

f_original_report.close()
f_panel_report.close()
