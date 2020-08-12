# -*- coding: utf-8 -*-
# @Author: Wang Lu
# @Date:   2020-08-11 23:20:52
# @Last Modified by:   Lenovo1
# @Last Modified time: 2020-08-12 13:56:15


from os import getcwd

from os import listdir

wf = open("test.fasta", "w")
with open("Cronobacter_OACs.fasta") as f:
	for line in f:
		if ">" in line and "galF" in line:
			m = line.strip("\n").split("_")[1]
			wf.write(">"+m+"\n")
			wf.write(next(f))

wf.close()