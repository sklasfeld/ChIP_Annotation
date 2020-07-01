#!/usr/bin/env python

import os
import argparse
import sys
import pandas as pd
import numpy as np


parser = argparse.ArgumentParser(description="reorder columns")
parser.add_argument('in_table', help='file that contains table that needs to be reorganized')
parser.add_argument('list_file', help='file that contains new-line delimited list in the \
	order that columns should be output. By default, this list should be organized using \
	the column name. Set `order_by_number` if this list is written using the column number.')
parser.add_argument('out', help='output file')
parser.add_argument('-n', '--noheader', action='store_true', default=False, \
					help='Set if `in_table` has no header. If this is set, \
					user can set label the columns with `colnames` parameter')
parser.add_argument('-c', '--colnames', nargs='+', help='`Set or Rename the columns.')
parser.add_argument('--sep', '-s', default="\t", help='table delimiter of \
					`in_table`. By default, the table is expected to be tab-delimited')
parser.add_argument('--out_sep', '-so', default="\t", \
					help='table delimiter of `out_table`. By default, \
					the out table will be tab-delimited')
parser.add_argument('--indexCol', '-i', \
					help='Column(s) to use as the row labels of the \
					`in_table`, either given as string name or column index.')
parser.add_argument('-obn', '--order_by_number', action='store_true', default=False, \
					help='Set if the values of `list_file` are the column numbers. \
					Note that the first column is columns 0, the second is column 1, etc.')
args = parser.parse_args()

# check file paths
if not os.path.isfile(args.in_table):
	sys.exit(("ERROR: cannot locate the `in_table` file: %s") % (args.in_table))
if not os.path.isfile(args.list_file):
	sys.exit(("ERROR: cannot locate the `list_file` file: %s") % (args.list_file))

# get the order of the columns
with open(args.list_file) as f:
    order_of_columns = f.readlines()
# remove whitespace characters like `\n` at the end of each line
order_of_columns = [x.strip() for x in order_of_columns]


# 1. Read input table
read_table_param={}
read_table_param["sep"]=args.sep
if args.noheader:
	read_table_param["header"]=None
if args.indexCol:
	read_table_param["index_col"]=args.indexCol


in_df = pd.read_csv(args.in_table, dtype=str, **read_table_param)

if args.colnames:
	if len(in_df.columns) != len(args.colnames):
		sys.exit(("ValueError: Length mismatch: Expected axis " +
			"has %i elements, new values have %i elements") %
		(len(in_df.columns), len(args.colnames)))
	in_df.columns = args.colnames

# check to make sure each column is available
if not args.order_by_number:
	for col in order_of_columns:
		if col not in list(in_df.columns):
			sys.exit("ERROR: Column name `%s` is not in the table" % col)
	out_df = in_df.loc[:,order_of_columns].copy()
else:
	ints_of_columns=[int(x) for x in order_of_columns]
	out_df = in_df.iloc[:,order_of_columns].copy()

out_param={}
out_param["sep"]=args.out_sep
if not args.indexCol:
	out_param["index"]=False

	out_df.to_csv(args.out, **out_param)