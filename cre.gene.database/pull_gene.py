#!/usr/bin/env python3

import pandas as pd
import sys
from datetime import datetime
from os import mkdir
from os.path import exists, join

GENE = sys.argv[1]
DB_DATE = sys.argv[2]
DAY = datetime.now().strftime("%Y-%m-%d")

SAMPLE_WISE_DB = "%s.c4r.sample_wise.csv" % DB_DATE
VARIANT_WISE_DB = "%s.c4r.variant_wise.csv" % DB_DATE

SAMPLE_WISE_GENE_REPORT = "%s.%s.sample_wise.csv" % (DAY, GENE)
VARIANT_WISE_GENE_REPORT = "%s.%s.variant_wise.csv" % (DAY, GENE)

sample_df = pd.read_csv(SAMPLE_WISE_DB)
variant_df = pd.read_csv(VARIANT_WISE_DB)

if not exists(DAY):
	mkdir(DAY)

sample_df[sample_df['Gene'] == GENE].to_csv(join(DAY, SAMPLE_WISE_GENE_REPORT), index=False)
variant_df[variant_df['Gene'] == GENE].to_csv(join(DAY, VARIANT_WISE_GENE_REPORT), index=False)