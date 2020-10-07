import pandas as pd
import sys

#reports generated from cre report_changes branch
# Synonymous: no filter applied on the severity of variants
# new report: alt_depth at least 3 in any of three samples, alt_depth of 1 and in clinvar
# variants in old report should all be included in the new report, and any additional variants will have alt_depth of at least 3 in any samples or 1 if in clinvar

#old reports
wes = sys.argv[1]
wes = pd.read_csv(wes, encoding='latin-1')

#reports from report_changes branch (new reports)
wes_regular = sys.argv[2]
wes_regular = pd.read_csv(wes_regular,  encoding='latin-1')


print('Comparing wes.regular to old wes report')
#test that all variants in wes (old report) are also in wes_regular (new report)
wes_regular_variants = wes_regular['Position'].tolist()
not_in_wes_regular = []
for variant in wes['Position'].tolist():
    if variant not in wes_regular_variants:
        print(wes[wes['Position'] == variant]['Depth'].tolist()[0])
        if int(wes[wes['Position'] == variant]['Depth'].tolist()[0]) >= 10:
            pass
        else:
            print(variant)
            print('Variant not present in new report that was present in old report')
    else:
        pass

#get variants in wes_regular (new report) that are not in wes (old report)
wes_variants = wes['Position'].tolist()
not_in_wes = []
for variant in wes_regular['Position'].tolist():
    if variant not in wes_variants:
        not_in_wes.append(variant)
    else:
        pass

#make dataframe of variants that were not present in old report but are present in new
not_in_wes_list = [wes_regular[wes_regular['Position'] == position] for position in not_in_wes]

#If no additional variants in new report, nothing to test
if len(not_in_wes_list) == 0:
    print('No additional variants with new filters in wes report')
#If there are additional variants with the new filters, test that they comply with stated new filters
else:
    not_in_wes_df = pd.concat(not_in_wes_list)

    col_names = not_in_wes_df.columns.tolist()
    alt_depth = [col_name for col_name in col_names if 'Alt_depths' in col_name]

    #check that variants present in new report that are not present in the old reports
    #have alt_depth > 3 in at least one sample or are present in ClinVar
    for index,row in not_in_wes_df.iterrows():
        alt_depth_all = [int(row[ad]) for ad in alt_depth]
        if min(alt_depth_all) >= 3:
            pass
        elif row['Clinvar'] != 'None':
            pass
        else:
             print('Variants present in new report do not comply with new filters')
