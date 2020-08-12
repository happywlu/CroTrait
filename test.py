# -*- coding: utf-8 -*-
# @Author: Wang Lu
# @Date:   2020-08-11 23:20:52
# @Last Modified by:   Lenovo1
# @Last Modified time: 2020-08-12 14:38:56


import os
import datetime
import argparse

today=datetime.date.today().strftime('%Y%m%d')

parser=argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='input the original file')
parser.add_argument('-o', '--out', help='the out file')
args=parser.parse_args()

def find_galF(i):
	pt=os.getcwd()

	if args.out:
		wf=open(os.path.join(pt,args.out),'w')
	else:
		wf=open(os.path.join(pt,today),'w')

	with open(os.path.join(pt,args.input)) as f:
		for line in f:
			if "galF" in line:
				wf.write(line)
				wf.write(next(f))
	wf.close()


#if __name__=='__main__':
#	find_galF(args.input)


