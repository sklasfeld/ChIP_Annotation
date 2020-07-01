#!/usr/bin/env python3

# -*- coding: iso-8859-15 -*-
# 2017, Samantha Klasfeld, the Wagner Lab
# the Perelman School of Medicine, the University of Pennsylvania
# Samantha Klasfeld, 12-21-2017


import argparse
import sys
import pandas as pd
import numpy as np


parser = argparse.ArgumentParser(description="this script takes \
								in a 2 tables and performs a \
								joins them to create a merged table")
parser.add_argument('left_table', help='left table file name')
parser.add_argument('right_table', help='right table file name')
parser.add_argument('out_table', help='output table file name')
parser.add_argument('-w','--how', help='Type of merge to be performed: \
	`left`,`right`,`outer`,`inner`, `antileft`. Default:`inner`', 
	choices=['left', 'right', 'outer', 'inner', 'antileft'], default='inner')
parser.add_argument('-j','--on', help='Column or index level names \
	to join on. These must be found in both DataFrames. If on is None \
	and not merging on indexes then this defaults to the intersection \
	of the columns in both DataFrames.', nargs='+')
parser.add_argument('-lo','--left_on', help='Column or index level names \
	to join on in the left DataFrame. Can also be an array or list of arrays \
	of the length of the left DataFrame. These arrays are treated as if \
	they are columns.', nargs='+')
parser.add_argument('-ro','--right_on', help='Column or index level names \
	to join on in the right DataFrame. Can also be an array or list of arrays \
	of the length of the left DataFrame. These arrays are treated as if \
	they are columns.', nargs='+')
parser.add_argument('-ml','--merge_left_index', help='Use the index from the left \
	DataFrame as the join key(s). If it is a MultiIndex, the number of keys \
	in the other DataFrame (either the index or a number of columns) must \
	match the number of levels.', action='store_true', default=False)
parser.add_argument('-mr','--merge_right_index', help='Use the index from the right \
	DataFrame as the join key(s). If it is a MultiIndex, the number of keys \
	in the other DataFrame (either the index or a number of columns) must \
	match the number of levels.', action='store_true', default=False)
parser.add_argument('-or','--order', help='Order the join keys \
	lexicographically in the result DataFrame. If False, the \
	order of the join keys depends on the join type (how keyword).', \
	action='store_true', default=False)
parser.add_argument('-su','--suffixes', help='Tuple of (str,str). Each str is a \
	Suffix to apply to overlapping column names in the left and right side, \
	respectively. To raise an exception on overlapping columns \
	use (False, False). Default:(`_x`,`_y`)', nargs=2)
parser.add_argument('-nl', '--noheader_l', action='store_true', default=False, \
					help='Set if `left_table` has no header. If this is set, \
					user must also set `colnames_l`')
parser.add_argument('-nr', '--noheader_r', action='store_true', default=False, \
					help='Set if `right_table` has no header. If this is set, \
					user must also set `colnames_r`')
parser.add_argument('-cl', '--colnames_l', nargs='+', \
					help='`If `noheader_l` is set, add column names \
					to `left_table`. Otherwise, rename the columns.')
parser.add_argument('-cr', '--colnames_r', nargs='+', \
					help='`If `noheader_r` is set, add column names \
					to `right_table`. Otherwise, rename the columns.')
parser.add_argument('--left_sep', '-sl', default="\t", \
					help='table delimiter of `left_table`. By default, \
					the table is expected to be tab-delimited')
parser.add_argument('--right_sep', '-sr', default="\t", \
					help='table delimiter of `right_table`. By default, \
					the table is expected to be tab-delimited')
parser.add_argument('--out_sep', '-so', default="\t", \
					help='table delimiter of `out_table`. By default, \
					the out table will be tab-delimited')
parser.add_argument('--left_indexCol', '-il', \
					help='Column(s) to use as the row labels of the \
					`left_table`, either given as string name or column index.')
parser.add_argument('--right_indexCol', '-ir', \
					help='Column(s) to use as the row labels of the \
					`right_table`, either given as string name or column index.')
parser.add_argument('-clc','--change_left_cols', nargs='+', 
	help='list of specific column names you want to change in left table. \
	For example, if you want to change columns `oldColName1` and \
	`oldColName2` to `newColName1` \
	and `newColName2`, respectively, then set this to \
	`oldColName2,newColName1 oldColName2,newColName2`')
parser.add_argument('-crc','--change_right_cols', nargs='+', 
	help='list of specific column names you want to change in right table. \
	For example, if you want to change columns `oldColName1` and \
	`oldColName2` to `newColName1` \
	and `newColName2`, respectively, then set this to \
	`oldColName2,newColName1 oldColName2,newColName2`')

#parser.add_argument('--header','-H', action='store_true', default=False, \
#					help='true if header in table')
args = parser.parse_args()

if args.noheader_l and not args.colnames_l:
	sys.exit("Error: If `noheader_l` is set, user must also set `colnames_l`\n")
if args.noheader_r and not args.colnames_r:
	sys.exit("Error: If `noheader_r` is set, user must also set `colnames_r`\n")

if args.change_left_cols and args.colnames_l:
	sys.exit("Error: Can only set one of these parameters:\n" +
		"\t* change_left_cols\n"+
		"\t* colnames_l\n")
if args.change_right_cols and args.colnames_r:
	sys.exit("Error: Can only set one of these parameters:\n" +
		"\t* change_right_cols\n"+
		"\t* colnames_r\n")
if not args.on:
	if not args.left_on and not args.right_on:
		sys.exit("Error: must set columns to join on.")

# 1. Read input files
read_ltable_param={}
read_rtable_param={}
read_ltable_param["sep"]=args.left_sep
read_rtable_param["sep"]=args.right_sep
if args.noheader_l:
	read_ltable_param["header"]=None
if args.noheader_r:
	read_rtable_param["header"]=None
if args.left_indexCol:
	read_ltable_param["index_col"]=args.left_indexCol
if args.right_indexCol:
	read_rtable_param["index_col"]=args.right_indexCol

left_df = pd.read_csv(args.left_table, **read_ltable_param)
right_df = pd.read_csv(args.right_table, **read_rtable_param)

# 2. Change/Update column names of the input tables
if args.colnames_l:
	if len(left_df.columns) != len(args.colnames_l):
		sys.exit(("ValueError: Length mismatch: Expected axis " +
			"has %i elements, new values have %i elements") %
		(len(left_df.columns), len(args.colnames_l)))
	left_df.columns = args.colnames_l
if args.colnames_r:
	if len(right_df.columns) != len(args.colnames_r):
		sys.exit(("ValueError: Length mismatch: Expected axis " +
			"has %i elements, new values have %i elements") %
		(len(right_df.columns), len(args.colnames_r)))
	right_df.columns = args.colnames_r

if args.change_left_cols:
	for left_changeCol_param in args.change_left_cols:
		if len(left_changeCol_param.split(",")) != 2:
			sys.exit("ERROR: values set to `change_left_cols` must " +
				"be in the format [old_col_name],[new_column_name]")
	rename_left_cols = dict(x.split(",") for x in args.change_left_cols)
	left_df = left_df.rename(columns=rename_left_cols)
if args.change_right_cols:
	for right_changeCol_param in args.change_right_cols:
		if len(right_changeCol_param.split(",")) != 2:
			sys.exit("ERROR: values set to `change_right_cols` must " +
				"be in the format [old_col_name],[new_column_name]")
	rename_right_cols = dict(x.split(",") for x in args.change_right_cols)
	right_df = right_df.rename(columns=rename_right_cols)






# 3. Set merge parameters
merge_param={}
if args.how == "antileft":
	merge_param['how']="left"
else:
	merge_param['how']=args.how
if args.on:
	merge_param['on']=args.on
if args.left_on:
	merge_param['left_on']=args.left_on
if args.right_on:
	merge_param['right_on']=args.right_on
if args.merge_left_index:
	merge_param['left_index']=args.merge_left_index
if args.merge_right_index:
	merge_param['right_index']=args.merge_right_index
if args.order:
	merge_param['sort']=args.order
if args.suffixes:
	merge_param['suffixes']=args.suffixes


# 4. Perform Merge
merge_df = left_df.merge(
	right_df, **merge_param) 

# 4B. There is an extra step for a left anti-join
# 5.  Export merged table
out_param={}
out_param["sep"]=args.out_sep
if not args.left_indexCol:
	out_param["index"]=False

if args.how == "antileft":
	antimerge_df = left_df.loc[merge_df.index,:].copy()
	antimerge_df.to_csv(args.out_table, **out_param)
else:
	merge_df.to_csv(args.out_table, **out_param)