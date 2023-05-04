#!/usr/bin/env python3.6
def GCstrand(strand):
	CG = 0
	total = 0
	percent = []
	
	for a in strand:
		total +=1
		if a == 'C' or a == 'G':
			CG+=1
	global GCcontents
	GCcontents = CG/total*100
	return(GCcontents)

def GCDB():
	import sys
	global seq_dic
	seq_dic = {}
	#part1
	#path = sys.argv[1]
	
	database = "./database.fa"
	dic_aln = {}
	seq = ""
	seqname = ""
	sn_list = []
	seq_n = 0
	#part2
	
	for line in open(database):
		if line[0] == ">":
			if len(seq)>0:
				dic_aln[seqname] = seq
			seqname = line.rstrip()[1:]
			seq = ""
			seq_n += 1
		elif seq_n == 1:
			seq += line.rstrip().upper()
			seq_n = 0
			#append just 1 linse - DNA strand
	
	#last species
	dic_aln[seqname] = seq


	if len(seq)>0:
		seq_dic[seqname] = seq
	seq_dic = dic_aln
	
	CG = 0
	total = 0
	percent = []

	for a in seq_dic.keys():
		for b in seq_dic[a]:
			total += 1
			if b == 'C' or b == 'G':
				CG += 1
		seq_dic[a] = CG/total * 100
		percent.append(CG/total * 100)
		CG = 0
		total = 0
	return(seq_dic)

def GC():
	import sys
	DNAstrand = sys.argv[1]
	strand_percent = GCstrand(DNAstrand)
	DNA_dic = GCDB()
	b = 100
	similar_n = ""
	for a in DNA_dic.keys():
		if abs(DNA_dic[a] - strand_percent) < b:
			b = abs(DNA_dic[a] - strand_percent)
			similar_n = a
	print("Species of input strand is similar with", similar_n)
	print("input percent:",strand_percent, "%/", "DB percent:", DNA_dic[similar_n],"%\n")
	


if __name__ == '__main__':
	GC()


