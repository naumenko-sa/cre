import argparse
import pandas as pd
import sys

def import_variants(report):
    report = pd.read_csv(report,  encoding='ISO-8859-1')
    report = report.set_index(['Position', 'Ref', 'Alt'])
    report.columns = [col.lower() for col in report.columns.tolist()]
    return report


def get_missing_variants(report1, report2):
    not_in_report1 = report2[~report2.index.isin(report1.index)]
    not_in_report1 = not_in_report1.reset_index()[['Position', 'Ref', 'Alt']]
    not_in_report1['Chrom'] = [pos.split(':')[0] for pos in not_in_report1['Position'].tolist()]
    not_in_report1['Start'] = [int(pos.split(':')[1])-1 for pos in not_in_report1['Position'].tolist()]
    not_in_report1['End'] = not_in_report1['Start']
    not_in_report1 = not_in_report1[['Chrom', 'Start', 'End', 'Ref', 'Alt']]
    return not_in_report1


def select_columns(report, columns):
    report = report[columns]
    return report


def get_diff_columns(columns1, columns2):
    not_in_1 = list(set(columns1)- set(columns2))
    not_in_2 = list(set(columns2)- set(columns1))
    return not_in_1, not_in_2


def get_diff_annotations(report1, report2, diff_columns, left_suffix, right_suffix):
    # identifies non-equivalent annotations for variants present in both reports
    # returns dataframe with divergent annotations for variants present in both reports
    # note that comparisons are by column (annotation),not by variant
    report1_common = report1[report1.index.isin(report2.index)].drop_duplicates()
    report2_common = report2[report2.index.isin(report1.index)].drop_duplicates()
    columns = list(set(report1_common.columns.tolist() + report2_common.columns.tolist()) - set(diff_columns))
    diff_annotations = []
    for col in columns: 
        report1_col = report1_common[col].copy().sort_index().tolist()
        report2_col = report2_common[col].copy().sort_index().tolist()
        if report1_col != report2_col:
            diff_annotations.append(col)
    merged_report = report1_common[diff_annotations].join(report2_common[diff_annotations], how='inner', lsuffix=left_suffix, rsuffix=right_suffix)
    return merged_report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Assess differences in variants and annotations between two variant reports')
    parser.add_argument('-old', help='old variant report', required=True)
    parser.add_argument('-new', help='new variant report', required=True)
    args = parser.parse_args()

    old_name = args.old.split('/')[-1]
    new_name = args.new.split('/')[-1]
    family = args.new.split('/')[-1].split('.')[0]
    
    old = import_variants(args.old)
    old_cols = old.columns.tolist()
    new = import_variants(args.new)
    new_cols = new.columns.tolist()

    # determine if column names in new report are equivalent to column names in old report
    if old_cols == new_cols:
        new_columns, missing_columns = [], []
    else:
        new_columns, missing_columns = get_diff_columns(new_cols, old_cols)

    # get variants unique to new and old reports
    not_in_old = get_missing_variants(old, new)
    not_in_old.to_csv('%s.uniq.pos'%new_name.strip('.csv'), index=False, header=False, sep='\t')
    not_in_new = get_missing_variants(new, old)
    not_in_new.to_csv('%s.uniq.pos'%old_name.strip('.csv'), index=False, header=False, sep='\t')

    # for variants common to new and old reports, test for annotation equivalence
    # include new and old annotations in csv with common variants if annotations are not equivalent
    different_columns = new_columns + missing_columns
    diff_annotations = get_diff_annotations(new, old, different_columns, '_new', '_old')
    diff_annotations = diff_annotations[sorted(diff_annotations)]
    diff_annotations.to_csv('%s.different_annotations.csv'%family, encoding='ISO-8859-1')

    # write summary file
    summary = open('summary.txt', 'w')
    summary.write(
        ('Comparing reports: {}, {}\nColumns unique to {}: {}\nColumns unique to {}: {}\nNumber of variants unique to {}: {}\nNumber of variants unique to {}: {}\nNumber of variants shared: {}').format(
            old_name, new_name, old_name, missing_columns, new_name, new_columns, old_name, len(not_in_new), new_name, len(not_in_old), len(diff_annotations))
        )