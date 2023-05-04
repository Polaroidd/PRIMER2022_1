#!/usr/bin/env python3

def PSI_BLAST():
	# PSI-BLAST scoring Table function
	def apply_protein_Table(sequence_length):
		F = [[0 for i in range(sequence_length)] for j in range(20)]
		#row: A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y
		return F


	#####


	#making Protein position probility * score matrix
	def scoring_protein(List, alignment_number, query, iteration_cnt, loc, score):
		A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y = 0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05
		for i in range(0, alignment_number):
			if List[i] == 'A': A += 1
			elif List[i] == 'C': C += 1
			elif List[i] == 'D': D += 1
			elif List[i] == 'E': E += 1
			elif List[i] == 'F': F += 1
			elif List[i] == 'G': G += 1
			elif List[i] == 'H': H += 1
			elif List[i] == 'I': I += 1
			elif List[i] == 'K': K += 1
			elif List[i] == 'L': L += 1
			elif List[i] == 'M': M += 1
			elif List[i] == 'N': N += 1
			elif List[i] == 'P': P += 1
			elif List[i] == 'Q': Q += 1
			elif List[i] == 'R': R += 1
			elif List[i] == 'S': S += 1
			elif List[i] == 'T': T += 1
			elif List[i] == 'V': V += 1
			elif List[i] == 'W': W += 1
			else: Y += 1    #List[i] == 'Y'
	
	
		if query == 'A': idx = 0
		elif query == 'C': idx = 1
		elif query == 'D': idx = 2
		elif query == 'E': idx = 3
		elif query == 'F': idx = 4
		elif query == 'G': idx = 5
		elif query == 'H': idx = 6
		elif query == 'I': idx = 7
		elif query == 'K': idx = 8
		elif query == 'L': idx = 9
		elif query == 'M': idx = 10
		elif query == 'N': idx = 11
		elif query == 'P': idx = 12
		elif query == 'Q': idx = 13
		elif query == 'R': idx = 14
		elif query == 'S': idx = 15
		elif query == 'T': idx = 16
		elif query == 'V': idx = 17
		elif query == 'W': idx = 18
		else: idx = 19    #List[i] == 'Y'
	

		#To prevent zero score, add 0.05 each. TOTAL 1
		i += 1

		if iteration_cnt == 1:
			A = A/i * score[0][idx]
			C = C/i * score[1][idx]
			D = D/i * score[2][idx]
			E = E/i * score[3][idx]
			F = F/i * score[4][idx]
			G = G/i * score[5][idx]
			H = H/i * score[6][idx]
			I = I/i * score[7][idx]
			K = K/i * score[8][idx]
			L = L/i * score[9][idx]
			M = M/i * score[10][idx]
			N = N/i * score[11][idx]
			P = P/i * score[12][idx]
			Q = Q/i * score[13][idx]
			R = R/i * score[14][idx]
			S = S/i * score[15][idx]
			T = T/i * score[16][idx]
			V = V/i * score[17][idx]
			W = W/i * score[18][idx]
			Y = Y/i * score[19][idx]
	
		else:
			A = A/i * score[0][loc]
			C = C/i * score[1][loc]
			D = D/i * score[2][loc]
			E = E/i * score[3][loc]
			F = F/i * score[4][loc]
			G = G/i * score[5][loc]
			H = H/i * score[6][loc]
			I = I/i * score[7][loc]
			K = K/i * score[8][loc]
			L = L/i * score[9][loc]
			M = M/i * score[10][loc]
			N = N/i * score[11][loc]
			P = P/i * score[12][loc]
			Q = Q/i * score[13][loc]
			R = R/i * score[14][loc]
			S = S/i * score[15][loc]
			T = T/i * score[16][loc]
			V = V/i * score[17][loc]
			W = W/i * score[18][loc]
			Y = Y/i * score[19][loc]
	
		return A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y


	#####


	def finding_consensus(F, length):
		best_score = 0
		temp_nt = ''
		consensus_motif = ''
		other_seq_num = 1
	
		for i in range(length):
			best_score = max(F[0][i],F[1][i],F[2][i],F[3][i],F[4][i],F[5][i],F[6][i],F[7][i],F[8][i],F[9][i],F[10][i],F[11][i],F[12][i],F[13][i],F[14][i],F[15][i],F[16][i],F[17][i],F[18][i],F[19][i])
		
			if best_score == F[0][i]:
				temp_nt += 'A'
			elif best_score == F[1][i]:
				temp_nt += 'C'
			elif best_score == F[2][i]:
				temp_nt += 'D'
			elif best_score == F[3][i]:
				temp_nt += 'E'
			elif best_score == F[4][i]:
				temp_nt += 'F'
			elif best_score == F[5][i]:
				temp_nt += 'G'
			elif best_score == F[6][i]:
				temp_nt += 'H'
			elif best_score == F[7][i]:
				temp_nt += 'I'
			elif best_score == F[8][i]:
				temp_nt += 'K'
			elif best_score == F[9][i]:
				temp_nt += 'L'
			elif best_score == F[10][i]:
				temp_nt += 'M'
			elif best_score == F[11][i]:
				temp_nt += 'N'
			elif best_score == F[12][i]:
				temp_nt += 'P'
			elif best_score == F[13][i]:
				temp_nt += 'Q'
			elif best_score == F[14][i]:
				temp_nt += 'R'
			elif best_score == F[15][i]:
				temp_nt += 'S'
			elif best_score == F[16][i]:
				temp_nt += 'T'				        
			elif best_score == F[17][i]:
				temp_nt += 'V'
			elif best_score == F[18][i]:
				temp_nt += 'W'
			elif best_score == F[19][i]:
				temp_nt += 'Y'
	
			other_seq_num *= len(temp_nt)
			consensus_motif += temp_nt
			temp_nt = ''
		
		return consensus_motif, other_seq_num


	#####

	def update_score_matrix(alignment, seq_length, seq_num, query_seq, cnt, score): #F is score matrix
		Temp = []
		for i in range(0, seq_length):
			for j in range(0, seq_num):
				Temp.append(alignment[j][i])
	
			protein = query_seq[i]
			F[0][i],F[1][i],F[2][i],F[3][i],F[4][i],F[5][i],F[6][i],F[7][i],F[8][i],F[9][i],F[10][i],F[11][i],F[12][i],F[13][i],F[14][i],F[15][i],F[16][i],F[17][i],F[18][i],F[19][i] = scoring_protein(Temp, seq_num, protein, cnt, i, score)
			Temp = []

		
		PSSM = {'A':F[0],'C':F[1],'D':F[2],'E':F[3],
			'F':F[4],'G':F[5],'H':F[6],'I':F[7],
			'K':F[8],'L':F[9],'M':F[10],'N':F[11],
			'P':F[12],'Q':F[13],'R':F[14],'S':F[15],
			'T':F[16],'V':F[17],'W':F[18],'Y':F[19]}
		origin_col = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
	
		return F, PSSM, origin_col
	


	#####	


	def calculate_score(alignment, query_seq, seq_length, seq_num, match_score, mismatch_score):
		score = 0
		for j in range(0, seq_length):
			if alignment[j] == query_seq[j]: score += match_score
			else: score += mismatch_score #alignment[i][j] != query_seq[j]
		
		return score
	

	#####

	#main code
	import pandas as pd
	import sys
	import random

	#database = sys.argv[1]
	#score_matrix = sys.argv[2]

	database = "./database.fa"
	score_matrix = "./BLOSUM62.fa"

	###load given database
	DB = {}
	species = []
	values = []
	key = ''
	value = []
	for line in open(database):
		line = line.rstrip()
		if line[0] == ">":
			if key != '':
				DB[key] = value
				species.append(key)
				values.append(value)
				value = []
			key = line[1:]	
		else:
			value.append(line)
	
	DB[key] = value
	species.append(key)
	values.append(value)

	
	# loading only the amino acid sequences
	data = []
	seq_temp = ''
	seq_num = len(species)
	
	for i in range(seq_num):
		data.append(values[i][1])


	# loading the score matrix
	score = [[0 for i in range(20)] for j in range(20)]
	line_num = 0


	for line in open(score_matrix):
		temp = line.rstrip().split() #split by ' '

		for i in range(20):
			score[line_num][i] = float(temp[i])

		line_num += 1



	#Printing Sequence
	print("<Sequences>")
	for i in range(0, seq_num):
		print("%3d sequence : %s" %(i+1, data[i]))
	print()



	#Choosing a Query Sequence(motif)
	#query = int(input("Please choose a Query Sequence(motif) : "))
	#query_seq = data[query-1]
	query_seq = sys.argv[1]
	
	match_score = int(input("Score of match(int) : "))
	mismatch_score = int(input("Score of mismatch(int) : "))
	threshold = float(input("Threshold score(float) : "))

	length = len(query_seq)
	F = apply_protein_Table(length)
	cnt = 1

	while (1):
		#Protein PSSM
		F, PSSM, origin_col = update_score_matrix(data, length, seq_num, query_seq, cnt, score)

		#Printing PSSM
		sj = pd.DataFrame(PSSM, index=[i for i in range(1, length+1)], columns=origin_col)
		sj.index.name = 'seq loc'
		print("\n<New Score Matrix>")
		print(sj)


		#printing the Best score and Consensus motif
		consensus_motif, other_seq_num = finding_consensus(F, length)
		con_score = calculate_score(consensus_motif, query_seq, length, 1, match_score, mismatch_score)
		Raw_score = []

		print("\n<Score>\n")
		for i in range(0, seq_num):
			Raw_score.append(calculate_score(data[i], query_seq, length, 1, match_score, mismatch_score))
			if Raw_score[i] >= con_score*(1-threshold):
				print("%3d sequence - %s : %d" %(i+1, data[i],Raw_score[i]))

		print("\n<RESULT %d>\n" %cnt)
		print("Query sequence  :", query_seq)
		print("Consensus motif :", consensus_motif)
		print("Consensus score :", con_score,"\n")

		if query_seq == consensus_motif:
			for query in range(seq_num):
				if query_seq == values[query][1]:
					break
			print("This is same with %s's sequence" %species[query])

		else:
			maximum = con_score
			index = 0
			for i in range(seq_num):
				if maximum <= Raw_score[i]:
					maximum = Raw_score[i]
					index = i
			print("This is highly similar with %s's sequence" %species[index])

		iteration = input("If you want run one more time, write 'yes': ")
		if iteration != 'yes': break
	
		threshold = float(input("Threshold score(float) : "))
		cnt += 1
		score = F

if __name__ == '__main__':
	PSI_BLAST()
