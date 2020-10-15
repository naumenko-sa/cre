import argparse
import pandas as pd
import sys

def import_variants(report):
    report = pd.read_csv(report,  encoding='ISO-8859-1')
    report = report.set_index(['Position', 'Ref', 'Alt'])
    return report


def get_missing_variants(report1, report2):
    report1_joined_report2 = report1.join(report2, how='left', lsuffix='_1', rsuffix='_2')
    not_in_report1 = report2[~report2.index.isin(report1_joined_report2.index)]
    return not_in_report1

def select_columns(report, columns):
    report = report[columns]
    return report

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Outputs variants that differ between reports')
    parser.add_argument('-old', help='old variant report', required=True)
    parser.add_argument('-new', help='new variant report', required=True)
    parser.add_argument('-columns', nargs= '+', help='list of columns of interest', required=True)
    args = parser.parse_args()

    old = import_variants(args.old)
    new = import_variants(args.new)

    if old.columns.tolist() == new.columns.tolist():
        print('Columns in old report are the same as columns in the new report')
    else:
        new_columns = list(set(new.columns.tolist()) - set(old.columns.tolist()))
        missing_columns = list(set(old.columns.tolist()) - set(new.columns.tolist()))
        print('Columns in old and new report are different')
        print('New columns:')
        print(new_columns)
        print('Missing columns:')
        print(missing_columns)

    not_in_old = get_missing_variants(old, new)
    not_in_new = get_missing_variants(new, old)

    columns = args.columns
    not_in_old = select_columns(not_in_old, columns)
    not_in_new = select_columns(not_in_new, columns)

    not_in_old.to_csv('variants_not_in_old.csv',  encoding='ISO-8859-1')
    not_in_new.to_csv('variants_not_in_new.csv',  encoding='ISO-8859-1')



