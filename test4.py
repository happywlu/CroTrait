# -*- coding: utf-8 -*-
# @Author: Wanglu
# @Date:   2020-08-12 14:54:15
# @Last Modified by:   Lenovo1
# @Last Modified time: 2020-08-15 00:39:50
# 
"""
discover currently known O-antigen biosynthesis gene clusters or
predict new O-antigen biosynthesis gene clusters in Cronobacter genus.

Written by Wanglu
contacted at wlubio@sina.com
"""

from __future__ import print_function
import os
import sys
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import subprocess
from shutil import copy
# import argparse

# parser = argparse.ArgumentParser(description='wanwglu')

# parser.add_argument('--cluster', '-c', type=str, required=True)
# parser.add_argument('--extract', '-o', type=str, required=True)
# args = parser.parse_args()

species_Oserotype = {"condimenti":["CO1"],
					"dublinensis":["DO1a", "DO1b", "DO2"],
					"malonaticus":["MaO1", "MaO2", "MaO3", "MaO4"],
					"muytjensii":["MuO1", "MuO2"],
					"sakazakii":["SO1", "SO2", "SO3", "SO4", "SO6", "SO7"],
					"turicensis":["TO1", "TO3", "TO4"],
					"universalis":["UO1"]}

def tran_sequence(sequence):
	g_name = sequence[::-1]
	g_name = g_name[g_name.find(".")+1:][::-1]
	wf = open(g_name+"_1.fasta", "w")
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
	g_name = os.path.split(genome)[1]

	g_name = g_name[::-1]
	g_name = g_name[g_name.find(".")+1:][::-1]
	wf = open(g_name+"_1.fasta", "w")
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

	new_genome = seperate_sequence(g_name+"_1.fasta")

	galF_result = generate_blast(OACs=O_AGCs, genome=g_name+"_1.fasta", gene="galF")
	gnd_result = generate_blast(OACs=O_AGCs, genome=g_name+"_1.fasta", gene="gnd")
	os.remove(g_name+"_1.fasta")

	if galF_result=="bad" or gnd_result=="bad":
		return("bad")
	else:
		galF_split = galF_result.split("\t")
		gnd_split  = gnd_result.split("\t")
		if galF_split[0] == gnd_split[0]:
			galF_identity = float(galF_split[2])
			gnd_identity = float(gnd_split[2])
			galF_coverage = float(galF_split[3])/891
			gnd_coverage = float(gnd_split[3])/1407
			if galF_coverage > 0.6 and gnd_coverage > 0.6 and galF_identity > 80 and gnd_identity > 80:
				aim_sequence = list(new_genome)[0].get(galF_split[0])
				lst = [int(galF_split[6]), int(galF_split[7]), int(gnd_split[6]), int(gnd_split[7])]
				my_sequence = aim_sequence[min(lst)-1:max(lst)]
				if float(galF_split[6]) < float(gnd_split[6]) and float(galF_split[8]) < float(galF_split[9]):
					if outfmt == 0:
						return("yes")
					elif outfmt == 1:
						with open(g_name+"_OAGCs.fasta", "w") as wf:
							wf.write(">"+g_name+"\n")
							wf.write(my_sequence+"\n")
				elif float(galF_split[6]) > float(gnd_split[6]) and float(galF_split[8]) > float(galF_split[9]):
					my_sequence = Seq(my_sequence, IUPAC.unambiguous_dna)
					my_sequence = my_sequence.reverse_complement()
					my_sequence = str(my_sequence)
					if outfmt == 0:
						return("yes")
					elif outfmt == 1:
						with open(g_name+"_OAGCs.fasta", "w") as wf:
							wf.write(">"+g_name+"\n")
							wf.write(my_sequence+"\n")
				else:
					return("bad")
		else:
			return("bad")

# os.path.split()
def O_serotype(OACs_sequence, genome, species):
	g_name = genome[::-1]
	g_name = g_name[g_name.find(".")+1:][::-1]

	serotype = species_Oserotype.get(species)
	OACs = list(seperate_sequence(OACs_sequence))
	contig_name = OACs[1]
	contig_sequence = OACs[0]
	wf = open("specific_gene.fasta", "w")
	for i in contig_name:
		for sero in serotype:
			if sero in i:
				if "wzx" in i or "wzm" in i:
					wf.write(">"+i+"\n")
					wf.write(contig_sequence.get(i)+"\n")
	wf.close()

	wzx_wzm = seperate_sequence("specific_gene.fasta")

	os.system("makeblastdb -in specific_gene.fasta -dbtype nucl -out wzx_or_wzm")
	os.system("blastn -query "+genome+" -db wzx_or_wzm -outfmt 6 -out wzx_or_wzm.txt")
	result = open("wzx_or_wzm.txt").readline()
	os.remove("wzx_or_wzm.txt")
	os.remove("wzx_or_wzm.nhr")
	os.remove("wzx_or_wzm.nin")
	os.remove("wzx_or_wzm.nsq")
	os.remove("specific_gene.fasta")

	if len(result) == 0:
		new_type = O_antigen_cluster(OACs = OACs_sequence, genome = genome, outfmt = 0)
		if new_type == "yes":
			return([species,g_name,"new"])
		else:
			return([species,g_name,"bad"])
	else:
		result_split = result.split("\t")
		result_identity = float(result_split[2])
		coverage = float(result_split[3])/float(list(wzx_wzm)[2].get(result_split[1]))
		if result_identity > 75:
			if coverage > 0.4:
				idetified_serotype = result_split[1].split("_")[0]
				return([species,g_name,idetified_serotype])
			else:
				return([species,g_name,"bad"])
		else:
			new_type = O_antigen_cluster(OACs = OACs_sequence, genome = genome, outfmt = 0)
			if new_type == "yes":
				return([species,g_name,"new"])
			else:
				return([species,g_name,"bad"])


def main(OACs_sequence, genome_directory, species, blast, out_file, new_O_file):

	star_dir = os.getcwd()
	dir_path = os.path.abspath("%s" % genome_directory)
	if blast == "blastn":
		# windows subprocess.call(["where", "blastn"]); linux subprocess.call(["which", blastn"])
		ab = subprocess.call(["where", "blastn"])
		if ab == 0:
			print("welcome to O antigen cluster finder, please citate: XXXXX")
		else:
			print("blastn don't exist in path, please install or add it !")
			sys.exit()
	else:
		print("please check the right spelling of blastn !")
		sys.exit()
	wf1 = open(os.path.join(star_dir, out_file), "w")
	wf2 = open(os.path.join(star_dir, new_O_file), "w")
	for i in os.listdir(dir_path):
		copy(os.path.join(dir_path,i), star_dir)
		m = O_serotype(OACs_sequence=OACs_sequence, genome=i, species=species)
		wf1.write("\t".join(m)+"\n")
		if m == "new":
			gn = i[::-1]
			gn = gn[gn.find(".")+1:][::-1]
			new_fasta = O_antigen_cluster(OACs_sequence=OACs_sequence, genome=i, outfmt=1)
			print(new_fasta)
			wf2.write(">"+gn+"\n")
			wf2.write(new_fasta+"\n")
		os.remove(os.path.join(star_dir,i))
	wf1.close()
	wf2.close()
	
if __name__ == "__main__":
	main(OACs_sequence="Cronobacter_OACs.fasta", genome_directory="example_sequence",
		species="sakazakii", out_file="wanglu.txt", blast="blastn", new_O_file="new_O_serotype.fasta")