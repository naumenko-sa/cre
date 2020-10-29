import argparse
import csv
from datetime import date
import pandas as pd
import sys 

def db_output_to_dict(db_output):
    db1 = {}
    db2 = {}
    with open(db_output) as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            alt_depths = [row[col] for col in row if 'gt_depths' in col]
            if row['db'] == args.prefix1:
                variant = row['chrom'] + ':' +  row['start'] + ':' + row['end'] + ':' + row['ref'] + ':' + row['alt']
                db1[variant] = {'impact_severity': row['impact_severity'], 'clinvar_sig': row['clinvar_sig'], 'clinvar_pathogenic': row['clinvar_pathogenic'], 'alt_depths': alt_depths}
            elif row['db'] == args.prefix2:
                variant = row['chrom'] + ':' +  row['start'] + ':' + row['end'] + ':' + row['ref'] + ':' + row['alt']
                db2[variant] = {'impact_severity': row['impact_severity'], 'clinvar_sig': row['clinvar_sig'], 'clinvar_pathogenic': row['clinvar_pathogenic'], 'alt_depths': alt_depths}
    return db1, db2

def get_explanations(report1_var, report2_var):
    explanation = {}
    for variant in report1_var:
        report1_var[variant]['alt_depths'] = [int(ad.strip(' ')) for ad in report1_var[variant]['alt_depths']]
        if variant not in report2_var:
            explanation[variant] = 'Variant not present in comparison database'
        elif report1_var[variant]['clinvar_sig'] != 'None' and report2_var[variant]['clinvar_sig'] == 'None':
            explanation[variant] = 'Change in clinvar_sig from %s to None'%report1_var[variant]['clinvar_sig']
        elif report1_var[variant]['clinvar_pathogenic'] != 'None' and report2_var[variant]['clinvar_pathogenic'] == 'None':
            explanation[variant] = 'Change in clinvar_pathogenic from %s to None'%report1_var[variant]['clinvar_pathogenic']
        elif report1_var[variant]['impact_severity'] != 'LOW' and report2_var[variant]['impact_severity'] == 'LOW':
            explanation[variant] = 'Change in impact_severity from %s to LOW'%report1_var[variant]['impact_severity']
        elif 3 <= max(report1_var[variant]['alt_depths']) < 10:
            explanation[variant] = 'Alt depth less than 10 but greater than 3'

        else:
            explanation[variant] = 'Cannot explain'
    return explanation


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Determine cause for inclusion/exclusion of variants between two reports')
    parser.add_argument('-db_output1', help='output in format <prefix>.uniq.db.txt from gemini_compare.sh for first report', required=True)
    parser.add_argument('-db_output2', help='output in format <prefix>.uniq.db.txt from gemini_compare.sh for second report', required=True)
    parser.add_argument('-prefix1', help='prefix of first report, e.g. 428.wes.regular.2020-10-14', required=True)
    parser.add_argument('-prefix2', help='prefix of first report, e.g. 428.wes.regular.2020-10-19', required=True)
    args = parser.parse_args()

    # For variants unique to report 1, determine reason they were not included in report 2
    db1_unique = db_output_to_dict(args.db_output1)
    report1_var = db1_unique[0]
    report2_var = db1_unique[1]
    explanation_1 = get_explanations(report1_var,report2_var)
    
    today = date.today()
    today = today.strftime("%Y-%m-%d")
    explanation_1_df = pd.DataFrame.from_dict(explanation_1, orient='index').reset_index()
    explanation_1_df.to_csv('validation_summary_unique_in_%s_%s.csv'%(args.prefix1,today), header=['Variant', 'Explanation'], index=False)

    # For variants unique to report 2, determine reason they were not included in report 1
    db2_unique = db_output_to_dict(args.db_output2)
    report1_var = db2_unique[0]
    report2_var = db2_unique[1]
    explanation_2 = get_explanations(report2_var, report1_var)
    
    explanation_2_df = pd.DataFrame.from_dict(explanation_2, orient='index').reset_index()
    explanation_2_df.to_csv('validation_summary_unique_in_%s_%s.csv'%(args.prefix2,today), header=['Variant', 'Explanation'], index=False)






