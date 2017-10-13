#!/usr/bin/python

#annotates cre report with genes list from PhenomeCentral

import csv
import sys

pc_number_of_occurrences = {}
pc_features = {}
pc_hpo_ids = {}

with open(sys.argv[2],'r') as f_genes:
    reader = csv.DictReader(f_genes, delimiter='\t')
    for row in reader:
	gene_id=row['Gene ID']
	pc_number_of_occurrences[gene_id] = row['Number of occurrences']
	pc_features[gene_id] = row['Features']
	pc_hpo_ids[gene_id] = row['HPO IDs']

f_genes.close()

with open(sys.argv[1],'rb') as f_csv:
    with open(sys.argv[3],'wb') as f_out:
	writer = csv.writer(f_out,quotechar='"', quoting = csv.QUOTE_ALL, lineterminator='\n')
        reader = csv.DictReader(f_csv)
	fieldnames = reader.fieldnames
	new_fields = ['PC_number_of_occurrences','PC_features','PC_HPO_IDs']
	for field in new_fields:
    	    fieldnames.append(field)

	writer.writerow(fieldnames)
	
	for row in reader:
	    ensembl_gene = row['Ensembl_gene_id']
	    row['PC_number_of_occurrences'] = pc_number_of_occurrences.get(ensembl_gene,'NA')
	    row['PC_features'] = pc_features.get(ensembl_gene,'NA')
	    row['PC_HPO_IDs'] = pc_hpo_ids.get(ensembl_gene,'NA')
	    ll = []
	    for field in fieldnames:
		ll.append(row[field])
	    writer.writerow(ll)

f_csv.close()
f_out.close()
