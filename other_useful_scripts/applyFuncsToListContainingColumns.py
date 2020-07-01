#!/usr/bin/env python

import os
import argparse
import sys
import pandas as pd
import numpy as np

# function to find the maximum/minimum
# in a list that is currently in 
#string format
def max_min(str_list, str_del,func, abval=False): 
	if len(str_list) == 0:
		return np.NaN
	elif str_list == "False":
		return np.NaN
	elif str_list == "None":
		return np.NaN
	else:
		split_list = str_list.split(str_del)
		if 'None' in split_list:
			split_list.remove('None')
		split_list = [np.float64(x) for x in split_list]
		if abval:
				split_list = np.abs(split_list)
		if func == "max":
			return np.max(split_list)
		else:
			return np.min(split_list)

# function to find the most
# frequent value in a list
# that is currently in string
# format
def most_frequent(str_list, str_del): 
    split_list = str_list.split(str_del)
    counter = 0
    most_freq_obj = ""
      
    for obj in split_list: 
        curr_frequency = split_list.count(obj) 
        if curr_frequency > counter: 
            counter = curr_frequency 
            most_freq_obj = obj
        elif curr_frequency== counter:
        	if len(most_freq_obj) > 0:
	        	most_freq_obj = most_freq_obj + str_del + obj
	        else:
	        	most_freq_obj = obj
  
    return obj


parser = argparse.ArgumentParser(description="Given a table that contains lists as values, \
	add columns that contain the maximum, minimum, mode, or anyTrue of the user-specified columns")
parser.add_argument('in_table', help='file that contains table that needs to be reorganized. \
	The table MUST CONTAIN A HEADER with the column names')
parser.add_argument('cmd_table', help='file that contains table where the first column \
	contains a column name from `in_table`. The second column contains one of the \
	values: `max` (the maximum),`min` (the minimum), `absmax` (the maximum of the \
	absolute values), `absmin` (the minimum of the absolute values), `mode`(most frequent value(s)), \
	`anyTrue` (True if list contains True),. The third column contains the new column \
	name. Note that this table should NOT contain a header.')
parser.add_argument('out', help='output file')
parser.add_argument('--tsep1', '-t1', default="\t", help='table delimiter of \
					`in_table`. By default, the table is expected to be tab-delimited')
parser.add_argument('--tsep2', '-t2', default="\t", help='table delimiter of \
					`cmd_table`. By default, the table is expected to be tab-delimited')
parser.add_argument('--lsep', '-l', help='list delimiter of \
					values in `in_table`. By default, the lists are expected to be \
					delimited by semicolons', default=";")
parser.add_argument('--out_sep', '-so', default="\t", \
					help='table delimiter of `out_table`. By default, \
					the out table will be tab-delimited')
parser.add_argument('--indexCol', '-i', \
					help='Column(s) to use as the row labels of the \
					`in_table`, either given as string name or column index.')
args = parser.parse_args()

# check file paths
if not os.path.isfile(args.in_table):
	sys.exit(("ERROR: cannot locate the `in_table` file: %s") % (args.in_table))
if not os.path.isfile(args.cmd_table):
	sys.exit(("ERROR: cannot locate the `cmd_table` file: %s") % (args.cmd_table))

# 1. Read input table
read_table_param={}
read_table_param["sep"]=args.tsep1
if args.indexCol:
	read_table_param["index_col"]=args.indexCol


in_df = pd.read_csv(args.in_table, **read_table_param, dtype=str)

# read command table

read_cmdtable_param={}
read_cmdtable_param["sep"]=args.tsep2


cmd_df = pd.read_csv(args.cmd_table, sep=args.tsep2, header=None)
if len(cmd_df.columns) != 3:
	sys.exit(("ERROR: The `cmd_table` has %i columns,"+ 
		" but should have EXACTLY THREE columns.") % (len(cmd_df.columns)))
cmd_df.columns=['old_colname','cmd','new_colname']

for index, row in cmd_df.iterrows():
	old_cname = row["old_colname"]
	new_cname = row["new_colname"]
	if row['cmd'] == "max":
		in_df[new_cname] = in_df.apply(
			lambda x: max_min(str(x[old_cname]), args.lsep, "max"), 
			axis=1)
	elif row['cmd'] == "min":
		in_df[new_cname] = in_df.apply(
			lambda x: max_min(str(x[old_cname]), args.lsep, "min"), 
			axis=1)
	elif row['cmd'] == "absmax":
		in_df[new_cname] = in_df.apply(
			lambda x: max_min(str(x[old_cname]), args.lsep, "max", True), 
			axis=1)
	elif row['cmd'] == "absmin":
		in_df[new_cname] = in_df.apply(
			lambda x: max_min(str(x[old_cname]), args.lsep, "min", True), 
			axis=1)
	elif row['cmd'] == "mode":
		in_df[new_cname] = in_df.apply(
			lambda x: most_frequent(str(x[old_cname]),args.lsep), 
			axis=1)
	elif row['cmd'] == "anyTrue":
		in_df[new_cname] = in_df.apply(
			lambda x: (str(x[old_cname]).split(args.lsep)).count('True') > 0, 
			axis=1)
	else:
		sys.exit(("ERROR: Row %i in `cmd_table` (%s) does not " + 
			"contain any one of the valid arguments:\n" + 
			"\t1. `max` :the maximum\n" +
			"\t2. `min` :the minimum\n" +
			"\t3. `mode`:most frequent value(s)\n" +
			"\t4. `anyTrue`: reports True if any value " +
			"in list is True\n\n" +
			"The argument %s is not available at this time.") % 
			(index+1, args.cmd_table, row['cmd']))

out_param={}
out_param["sep"]=args.out_sep
if not args.indexCol:
	out_param["index"]=False

in_df.to_csv(args.out, **out_param)

