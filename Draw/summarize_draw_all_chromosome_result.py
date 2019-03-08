# -*- coding: utf-8 -*-
import argparse
import xlsxwriter
import codecs
import re
import os
import numpy as np
import pandas as pd
from pandas import Series, DataFrame

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--out", help="output excel file name")
parser.add_argument("-n", "--name", help="sample namaes")
parser.add_argument("-c", "--chromosome", help="chromosome namaes")
parser.add_argument("--image", help="path to image directory")
args = parser.parse_args()


FILE_out = args.out
DIR_image = args.image
NAMEs=args.name.split(",")
CHROMOSOMEs=args.chromosome.split(",")

workbook = xlsxwriter.Workbook(FILE_out)
worksheet = workbook.add_worksheet()

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
worksheet.set_column('A:A', 13)
worksheet.set_column(1,len(NAMEs)+1, 42.5)

TITLE = ["chromosome"] + NAMEs
worksheet.write_row(0, 0, TITLE, format_header)

for i in range(0,len(CHROMOSOMEs)):
	row=i+1
	worksheet.set_row(row, 228)
	worksheet.write(row, 0, unicode(CHROMOSOMEs[i], 'utf-8'), format_string)

	for j in range(0, len(NAMEs)):
		FILE_image = DIR_image + "/" + NAMEs[j] + "_" +  CHROMOSOMEs[i] + ".png"
		if os.path.isfile(FILE_image):
			worksheet.insert_image(row,j+1, FILE_image, {'x_offset': 2, 'y_offset': 2, 'x_scale': .58, 'y_scale': .58})

workbook.close()
