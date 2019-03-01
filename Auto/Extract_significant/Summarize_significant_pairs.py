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
parser.add_argument("-t", "--title", help="title of excel tab")
parser.add_argument("--image", help="path to image directory")
args = parser.parse_args()


FILE_in = args.pair
FILE_out = args.out
DIR_image = args.image

#==========================================
# load input
#==========================================
DATA = pd.read_csv(FILE_in, sep='\t', header=0)

workbook = xlsxwriter.Workbook(FILE_out)
worksheet = workbook.add_worksheet(args.title)

format_sig = workbook.add_format({'bg_color': '#D1E68F', 'font_color': '#000000', 'valign': 'top'})
format_string = workbook.add_format({'align': 'right', 'valign': 'top'})
format_header = workbook.add_format({'bold': True, 'align': 'center', 'text_wrap': True, 'bg_color': '#DBDBDB'})
format_number = workbook.add_format({'num_format': '#,##0', 'valign': 'top'})
format_small = workbook.add_format({'num_format': '0.000', 'valign': 'top'})
format_pval = workbook.add_format({'num_format': '0.000E+00', 'valign': 'top'})
format_url = workbook.add_format({'underline':  1, 'font_color': 'blue','valign': 'top'})

#==========================================
# pathway
#==========================================
worksheet.set_column('A:B', 6.5)
worksheet.set_column('C:D', 12)
worksheet.set_column('E:E', 6.5)
worksheet.set_column('F:H', 12)
worksheet.set_column('I:I', 13.5)
worksheet.set_column('J:L', 7.5)
worksheet.set_column('M:N', 11)


TITLE = ["rank", "chr1", "start1", "end1", "chr2", "start2", "end2", "distance", "map", "score", "control", "fc", "pval", "qval"]
worksheet.write_row(0, 0, TITLE, format_header)

for i in range(0,len(DATA)):
	row=i+1
	worksheet.set_row(row, 70)
	worksheet.write(row, 0, float(DATA.iloc[i, 0]), format_number)
	worksheet.write(row, 1, unicode(DATA.iloc[i, 1], 'utf-8'), format_string)
	worksheet.write(row, 2, float(DATA.iloc[i, 2]), format_number)
	worksheet.write(row, 3, float(DATA.iloc[i, 3]), format_number)
	worksheet.write(row, 4, unicode(DATA.iloc[i, 4], 'utf-8'), format_string)
	worksheet.write(row, 5, float(DATA.iloc[i, 5]), format_number)
	worksheet.write(row, 6, float(DATA.iloc[i, 6]), format_number)
	worksheet.write(row, 7, float(DATA.iloc[i, 7]), format_number)
 
	FILE_image = DIR_image + "/rank_" + str(DATA.iloc[i, 0]) + ".png"
	worksheet.insert_image(row,8, FILE_image, {'x_offset': 2, 'y_offset': 2, 'x_scale': .17, 'y_scale': .17})

	worksheet.write(row, 9, float(DATA.iloc[i, 8]), format_small)
	worksheet.write(row, 10, float(DATA.iloc[i, 9]), format_small)
	worksheet.write(row, 11, float(DATA.iloc[i, 10]), format_small)
	worksheet.write(row, 12, float(DATA.iloc[i, 11]), format_pval)
	worksheet.write(row, 13, float(DATA.iloc[i, 12]), format_pval)

workbook.close()
