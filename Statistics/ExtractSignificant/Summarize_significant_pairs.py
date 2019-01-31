# -*- coding: utf-8 -*-
import argparse
import xlsxwriter
import codecs
import re
import numpy as np
import pandas as pd
from pandas import Series, DataFrame

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--pair", help="significant pair file")
parser.add_argument("-o", "--out", help="output file")
parser.add_argument("-g", "--group", help="group")
args = parser.parse_args()


FILE_in = args.pair
FILE_out = args.out
GROUP = [0]*11 + list(map(int, args.group.split(",")))

#==========================================
# load input
#==========================================
DATA = pd.read_csv(FILE_in, sep='\t', header=0)

workbook = xlsxwriter.Workbook(FILE_out)
worksheet = workbook.add_worksheet()

format_sig = workbook.add_format({'bg_color': '#D1E68F', 'font_color': '#000000'})
format_string = workbook.add_format({'align': 'right'})
format_header = workbook.add_format({'bold': True, 'align': 'center', 'text_wrap': True})
format_header2 = workbook.add_format({'bold': True, 'align': 'center', 'text_wrap': True, 'font_color': '#0070C0'})
format_header3 = workbook.add_format({'bold': True, 'align': 'center', 'text_wrap': True, 'font_color': '#FF0000'})
format_number = workbook.add_format({'num_format': '#,##0'})
format_small = workbook.add_format({'num_format': '0.000'})
format_pval = workbook.add_format({'num_format': '0.000E+00'})

#==========================================
# pathway
#==========================================
worksheet.set_column('B:C', 12)
worksheet.set_column('E:G', 12)

string_col = set([0,3])
commmaNum_col = set([1,2,4,5,6])
smallNum_col = set([7,8,10])
pval_col = set([9])
worksheet.write_row(0, 0, DATA.columns, format_header)

# header
for col in range(0,len(DATA.columns)):
    if GROUP[col] == 0:
        worksheet.write(0, col, unicode(DATA.columns[col], 'utf-8'), format_header)
    elif GROUP[col] == 1:
        worksheet.write(0, col, unicode(DATA.columns[col], 'utf-8'), format_header2)
    elif GROUP[col] == 2:
        worksheet.write(0, col, unicode(DATA.columns[col], 'utf-8'), format_header3)


for i in range(0,len(DATA)):
    row=i+1
    for col in range(0,len(DATA.columns)):
        if col in string_col:
            worksheet.write(row, col, unicode(DATA.iloc[i, col], 'utf-8'), format_string)
        elif col in commmaNum_col:
            worksheet.write(row, col, float(DATA.iloc[i, col]), format_number)
        elif col in smallNum_col:
            worksheet.write(row, col, float(DATA.iloc[i, col]), format_small)
        elif col in pval_col:
            worksheet.write(row, col, float(DATA.iloc[i, col]), format_pval)
        else:
            worksheet.write_number(row, col, float(DATA.iloc[i, col]))




Range = "H2:H" + str(len(DATA)+1)
worksheet.conditional_format(Range, {'type': '3_color_scale',
                                         'min_type': 'num',
                                         'max_type': 'num',
                                         'mid_type': 'num',
                                         'min_value': -10,
                                         'max_value': 10,
                                         'mid_value': 0,
                                         'min_color': '#0070C0',
                                         'mid_color': '#FFFFFF',
                                         'max_color': '#FF0000'
})

Range = "J2:J" + str(len(DATA)+1)
worksheet.conditional_format(Range, {'type': 'cell',
                                         'criteria': '<',
                                         'value': 0.05,
                                         'format': format_sig})

workbook.close()
