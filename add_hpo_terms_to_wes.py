import pandas as pd
import sys
from os.path import basename, splitext

#Usage: python3 add_hpo_terms_to_wes.py <phenotips_hpo_terms> <wes.report.csv>

HPO_DF = pd.read_csv(sys.argv[1], comment='#', skip_blank_lines=True,\
	sep="\t",  encoding="ISO-8859-1", engine='python').set_index("Gene ID")\
	.drop(columns=["Gene Symbol", "HPO IDs"])
WES_REPORT = pd.read_csv(sys.argv[2], encoding="ISO-8859-1").set_index("Position")
OUT = "%s.w_hpoterms.tsv" % splitext(basename(sys.argv[2]))[0]

OUT_DF = WES_REPORT.join(HPO_DF, on="Ensembl_gene_id")\
	.rename(columns={"Number of occurrences": "HPO_count", "Features": "HPO_terms"}).fillna("NA")
OUT_DF.to_csv(OUT, sep="\t")