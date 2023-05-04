#!/usr/bin/env python3.6
import pandas as pd
import sys
import random
def DP():
#seqY는 DB내 전부!
	def DP_class(seqX,seqY): #from DB
	# Create a DP matrix
		global max_score
		dp_matrix = []
		arrow_matrix = []
		max_score = 0
		for i in range(0,len(seqX)+1):
			scores = []
			directions = []
			for j in range(0,len(seqY)+1):
				if i == 0 and j == 0:
					scores.append(0)
					directions.append('e')
				elif i == 0:
					scores.append(j*g_score)
					directions.append('l')
				elif j == 0:
					scores.append(i*g_score)
					directions.append('u')
				else:
					scores.append(0)
					directions.append('')
			dp_matrix.append(scores)
			arrow_matrix.append(directions)
		# Fill in the DP matrix
		for i in range(1, len(seqX)+1):
			for j in range(1, len(seqY)+1):
				d_score, u_score, l_score = 0,0,0
				# Calculate F(i-1,j-1) + s(xi,yi)
				if seqX[i-1] == seqY[j-1]:
					d_score = dp_matrix[i-1][j-1] + m_score 
				else:
					d_score = dp_matrix[i-1][j-1] + mm_score
				# Calculate F(i-1,j) + gap_score
				u_score = dp_matrix[i-1][j] + g_score
			# Calculate F(i,j-1) + gap_score
				l_score = dp_matrix[i][j-1] + g_score
				
				max_score = max([d_score, u_score, l_score])
				dp_matrix[i][j] = max_score

				max_dir = ''
				if d_score == max_score: max_dir += 'd'
				if u_score == max_score: max_dir += 'u'
				if l_score == max_score: max_dir += 'l'
				arrow_matrix[i][j] = max_dir
	#Tracback
		i = len(seqX)	
		j = len(seqY)
		aln_seqX = ''
		aln_seqY = ''
		while i != 0 or j != 0:
			max_dirs = list(arrow_matrix[i][j])
			tb_dir = random.choice(max_dirs)
			if tb_dir == 'd':
				aln_seqX = seqX[i-1] + aln_seqX
				aln_seqY = seqY[j-1] + aln_seqY
				i -= 1
				j -= 1
			elif tb_dir == 'u':
				aln_seqX = seqX[i-1] + aln_seqX
				aln_seqY = '-' + aln_seqY
				i -= 1
			else:
				aln_seqX = '-' + aln_seqX
				aln_seqY = seqY[j-1] + aln_seqY
				j -= 1	
		'''
		print('Optimal alignment score:',max_score)
		print(aln_seqX)
		print(aln_seqY)
		# Print the pairwise alignment and the optimal alignment score
		'''
		return([max_score, aln_seqX, aln_seqY])


	#seq_f = sys.argv[1]
	#m_score = int(sys.argv[2])
	#mm_score = int(sys.argv[3])
	#g_score = int(sys.argv[4])
	database = "./database.fa"
	m_score = int(input("write match score : "))
	mm_score = int(input("write mismatch score : "))
	g_score = int(input("write gap score : "))
	


	# Read a sequence file
	
	#main code
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
			sn_list.append(seqname)
			seq = ""
			seq_n += 1
		elif seq_n == 1:
			seq += line.rstrip().upper()
			seq_n = 0
			#append just 1 linse - DNA strand
	
	if len(seq)>0:
		dic_aln[seqname] = seq
	
	X = sys.argv[1]
	Y = ""
	#seqY는 DB내 전부!
	max_list = []
	for a in sn_list:
		Y = dic_aln[a]
		max_list.append(DP_class(X,Y)[0])
	for a in sn_list:
		Y = dic_aln[a]
		if max(max_list) == DP_class(X,Y)[0]:
			K = DP_class(X,Y)
			print('\noptimal alignment score :', K[0])
			print('\n> strand that you entered')
			print(K[1])
			print('\n>', a)
			print(K[2])
			print("\nThis strand is similar with", a, "!\n")

if __name__ == '__main__':
	DP()

