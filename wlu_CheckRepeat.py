# -*- coding:utf-8 -*-

"""
__author__ = "Wanglu"
__license__ = "NanKai University"
__version__ = "1.0"
__maintainer__ = "Wanglu"
__email__ = "wlubio@sina.com" or "wluhappy@gmail.com"

This script is write to check the repeat in the table.

"""

import os
import datetime
import argparse

today=datetime.date.today().strftime('%Y%m%d')

parser=argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='input the original file')
parser.add_argument('-p', '--prefix', help='add the prefix in the final result')
args=parser.parse_args()

def find_repeat (i):
	pt=os.getcwd()+'/'

	if args.prefix:
		b=open(pt+today+args.prefix+'.txt','w')
	else:
		b=open(pt+today+'.txt','w')

	list1=[]
	list2=[]
	for line in open(pt+i,'r'):
		list1.append(line.strip('\n'))
	for element in list1 :
		if element not in list2 :
			list2.append(element)
		else:
			b.write(element+'\n')
	b.close()


if __name__=='__main__':
	find_repeat(args.input)
