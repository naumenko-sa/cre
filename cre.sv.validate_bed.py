""" Calculates sensitivity and precision for 2 bed files, modified from bcbio.structural.validate
    $1 = bed file to evaluate
    $2 = truth set
    $3 = caller
    outpus: caller.csv
"""

import csv
import sys

from bcbio.structural import validate

data='ploidy'
svtype = 'DEL'
vcaller = sys.argv[3]

with open(vcaller+".csv","w") as out_handle:
    with open(vcaller+".df.csv", "w") as df_out_handle:
        writer = csv.writer(out_handle)
        dfwriter = csv.writer(df_out_handle)
        writer.writerow(["svtype", "size", "caller", "sensitivity", "precision"])
        dfwriter.writerow(["svtype", "size", "caller", "metric", "value", "label"])
        for size in validate.EVENT_SIZES:
            str_size = "%s-%s" % size
            evalout = validate._evaluate_one(vcaller,svtype,size,sys.argv[1],sys.argv[2],data)
            writer.writerow([svtype, str_size, vcaller,
                             evalout["sensitivity"]["label"], evalout["precision"]["label"]])
            for metric in ["sensitivity", "precision"]:
                dfwriter.writerow([svtype, str_size, vcaller, metric,
                                   evalout[metric]["val"], evalout[metric]["label"]])