# -*- coding: utf-8 -*-
# @Author: Lenovo1
# @Date:   2020-08-12 14:54:15
# @Last Modified by:   Lenovo1
# @Last Modified time: 2020-08-12 15:35:51


import os
import argparse

parser = argparse.ArgumentParser(description='wanwglu')

parser.add_argument('--cluster', '-c', type=str, required=True)
parser.add_argument('--extract', '-o', type=str, required=True)
args = parser.parse_args()

pt = os.getcwd()
wf = open(os.path.join(pt, args.extract), "w")

with open(os.path.join(pt, args.cluster)) as f:
	for line in f:
		if "galF" in line:
			wf.write(line)
			wf.write(next(f))

wf.close()
