# -*- coding: utf-8 -*-
# @Author: Lenovo1
# @Date:   2020-08-14 10:20:38
# @Last Modified by:   Lenovo1
# @Last Modified time: 2020-08-14 10:40:35


# -*- coding: utf-8 -*-
# @Author: wanglu
# @Date:   2020-08-12 14:54:15
# @Last Modified by:   Lenovo1
# @Last Modified time: 2020-08-14 09:49:18
# 
"""
discover currently known O-antigen biosynthesis gene clusters or
predict new O-antigen biosynthesis gene clusters.

Written by Wanglu
contacted at wlubio@sina.com
"""

from __future__ import print_function

import os
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
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


def tran_sequence(sequence):
	g_name = sequence[::-1]
	g_name = g_name[g_name.find(".")+1:][::-1]
	wf = open(g_name+"1.fasta", "w")
	trait = "wlu"
	with open(sequence) as f:
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


# for genome, run tran_sequence before running seperate_sequence
def seperate_sequence(sequence):
	contig_sequence = dict()
	contig_name = list()
	contig_length = dict()
	with open(sequence, "r") as f:
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



def O_antigen_cluster(OACs, genome, outfmt = 0):
	O_AGCs = seperate_sequence(OACs)

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

	new_genome = seperate_sequence(g_name+"1.fasta")

	galF_result = generate_blast(OACs=O_AGCs, genome=g_name+"1.fasta", gene="galF")
	gnd_result = generate_blast(OACs=O_AGCs, genome=g_name+"1.fasta", gene="gnd")

	if galF_result=="bad" or gnd_result=="bad":
		return("bad1")
	else:
		galF_split = galF_result.split("\t")
		gnd_split  = gnd_result.split("\t")
		if galF_split[0] == gnd_split[0]:
			galF_identity = float(galF_split[2])
			gnd_identity = float(gnd_split[2])
			galF_coverage = float(galF_split[3])/891
			gnd_coverage = float(gnd_split[3])/1407
			if galF_coverage>0.6 and gnd_coverage>0.6 and galF_identity>80 and gnd_identity>80:
				aim_sequence = list(new_genome)[0].get(galF_split[1])
				if float(galF_split[6]) < float(gnd_split[6]) and float(galF_split[8]) < float(galF_split[9]:
					my_sequence = aim_sequence[int(galF_split[6]):int(gnd_split[7])]
					print(my_sequence)
		#	return(galF_identity, galF_coverage, gnd_identity, gnd_coverage)
		else:
			return("bad2")


wl = O_antigen_cluster(OACs="Cronobacter_OACs.fasta", genome = "MaO2_LMG23826.fna", outfmt = 0)
print(wl)