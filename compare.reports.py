import os
import sys
import csv
import re

skip_columns = ['UCSC_Link','GNOMAD_Link']
compare_columns = []

def read_report(report_file):
    
    report_lines = []
    header_pos_dict = {}
    with open(report_file, 'r+',encoding="ISO-8859-1") as report_fh:
        
        csv_reader = csv.reader(report_fh,delimiter=',')
        index = 0
        for line in csv_reader:
            if index == 0:
                line[0] = line[0].lstrip('\ufeff')
                line[0] = re.sub(r'\W+','',line[0])
                print(line[0])
                header_pos_dict = {index:col for index,col in enumerate(line)}
            else:
                report_lines.append({header_pos_dict[index]:col for index,col in enumerate(line)})
            index = index+1
    return report_lines

def compare_reports(report_arr, first_arg, second_arg, comparison_order):
   
    
    if int(comparison_order) == 1:
        first_file, second_file, first_file_name, second_file_name = report_arr[0], report_arr[1], first_arg, second_arg
    else:
        first_file, second_file, first_file_name, second_file_name = report_arr[1], report_arr[0], second_arg, first_arg

    print("Comparing {} to {} and reporting entires that are different or absent in {}".format(first_file_name, second_file_name, second_file_name))
    print("Not comparing columns "+ ", ".join(skip_columns)+"\n\n\n")
    print("-"*100)
    for report_one_line in first_file:
        var_name = "_".join([report_one_line['Position'],report_one_line['Ref'],report_one_line['Alt']])
        diff_columns = {}
        for report_two_line in second_file:
            if report_one_line['Position'] == report_two_line['Position'] and report_one_line['Ref'] == report_two_line['Ref'] and report_one_line['Alt'] == report_two_line['Alt']:
                for column in report_one_line:
                    column_prefix = column[0:column.index('.')] if '.' in column else column
                    if len(compare_columns) >0:
                        if column_prefix not in compare_columns:
                            continue
                    #if column not in skip_columns and column in report_two_line and 'Burden' not in column: #not comparing burden for now 
                    if column not in skip_columns and column in report_two_line: #not comparing burden for now 
                        if not (str(report_one_line[column]).lower() in ('','0','na') and str(report_two_line[column]).lower() in ('','0','na')):
                            if str(report_one_line[column]).lower() != str(report_two_line[column]).lower():
                                diff_columns[column] = {'report_one': report_one_line[column], 'report_two': report_two_line[column]}
                if len(diff_columns.keys()) > 0:
                    print("Variant "+ var_name+ " differs in reports for following columns\n")
                    for column in diff_columns:
                        print("Field: {}\n\nValue in report {}: {}\nValue in report {}: {}\n".format(column,first_file_name,diff_columns[column]['report_one'].encode('utf-8'),second_file_name,diff_columns[column]['report_two'].encode('utf-8')))
                    print("-"*100)
                break
        else: #this else corresponds to for loop. this will get executed if and only if the for loop reaches the end of the loop
            print("variant "+var_name+" is not present in report {}".format(second_file_name))

if __name__ == "__main__" :

    if len(sys.argv) !=4:
        sys.exit("Need these command line arguments 1) path to report one 2) path to report two 3) 1 to compare report one to two or 2 to compare report two to one ")
    
    report_arr = []
    for report_file in sys.argv[1:3]:
        if not os.path.exists(report_file):
            sys.exit(report_file + "doesn't exist")
        report_arr.append(read_report(report_file))
    if int(sys.argv[3]) <1 or int(sys.argv[3]) >2:
        sys.exit("Please enter either 1 or 2 for argument number three")
    compare_reports(report_arr, os.path.basename(sys.argv[1]), os.path.basename(sys.argv[2]), sys.argv[3])
    
    
