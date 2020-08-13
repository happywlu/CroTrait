# -*- coding: utf-8 -*-
# @Author: wanglu
# @Date:   2020-08-12 14:54:15
# @Last Modified by:   Lenovo1
# @Last Modified time: 2020-08-14 00:54:38
# 
"""
discover currently known O-antigen biosynthesis gene clusters or
predict new O-antigen biosynthesis gene clusters.

Written by Wanglu
contacted at wlubio@sina.com
"""

from __future__ import print_function

import os
# import argparse

# parser = argparse.ArgumentParser(description='wanwglu')

# parser.add_argument('--cluster', '-c', type=str, required=True)
# parser.add_argument('--extract', '-o', type=str, required=True)
# args = parser.parse_args()

pt = os.getcwd()

speies_O_antigen = {"condimenti":["CO1"],
					"dublinensis":["DO1a", "DO1b", "DO2"],
					"malonaticus":["MaO1", "MaO2", "MaO3", "MaO4"],
					"muytjensii":["MuO1", "MuO2"],
					"sakazakii":["SO1", "SO2", "SO3", "SO4", "SO6", "SO7"],
					"turicensis":["TO1", "TO3", "TO4"],
					"universalis":["UO1"]}

species = ["condimenti","dublinensis","malonaticus","muytjensii",
		   "sakazakii","turicensis","universalis"]

O_types = ["CO1","DO1a","DO1b","DO2","MaO1","MaO2","MaO3","MaO4",
		   "MuO1","MuO2","SO1","SO2","SO3","SO4","SO6","SO7",
		   "TO1","TO3","TO4","UO1"]


galF_gnd = ["galF", "gnd"]

def seperate_sequence(OAC):
	contig_sequence = dict()
	contig_name = list()

	contig_length = dict()
	with open(OAC, "r") as f:
		for line in f:
			if ">" in line:
				m = line.strip("\n")[1:]
				seq = next(f).strip("\n")
				contig_sequence[m] = seq
				contig_name.append(m)
				contig_length[m] = len(seq)
	return(contig_sequence, contig_name, contig_length)



def generate_blast(OACs, gene, genome):
	wf = open(gene+".fasta", "w")
	m1 = list(OACs)[0]
	m2 = list(OACs)[1]
	for i in m2:
		if gene in i:
			wf.write(">"+i+"\n")
			wf.write(m1.get(i)+"\n")
	wf.close()
	os.system("makeblastdb -in "+gene+".fasta "+"-dbtype nucl "+"-out "+gene)
	os.system("blastn -query "+genome+" -db "+gene+" -outfmt 6 -out "+gene+".txt")
	result = open(gene+".txt").readline()
	os.remove(gene+".txt")
	os.remove(gene+".fasta")
	os.remove(gene+".nhr")
	os.remove(gene+".nin")
	os.remove(gene+".nsq")

	if len(result) == 0:
		return("bad")
	else:
		return(result)

def tran_genome(genome):
	g_name = genome[::-1]
	g_name = g_name[g_name.find(".")+1:][::-1]
	wf = open(g_name+"1.fasta", "w")
	trait = "wlu"
	with open(genome) as f:
		for line in f:
			if ">" in line:
				contig_name = line.strip("\n").split(" ")[0][1:]
				if trait == "wlu":
					wf.write(">"+contig_name+"\n")
					trait = "wlei"
				else:
					wf.write("\n>"+contig_name+"\n")
			else:
				wf.write(line.strip("\n"))
	wf.write("\n")
	wf.close()



# tran_genome("MaO2_LMG23826.fna")




wl = seperate_sequence(OAC="Cronobacter_OACs.fasta")
wl1 = generate_blast(OACs = wl, genome = "MaO2_LMG23826.fna", gene = "galF")

print(wl1)





#def generate_blast(sequence, genmoe, gene):


#	wf = open("galF.fasta", "w")
#	with open(OAGCs) as f:
#		for line in f:
#			if "galF" in line:
#				wf.write(line)
#				wf.write(next(f))
#	os.system("makeblastdb -in galF -dbtype nucl -out galF")