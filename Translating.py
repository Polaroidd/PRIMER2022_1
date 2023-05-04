#!/usr/bin/env python3

def Translation():
	#function that return the index
	def apply_table(string):
		if string == 'A':
			temp = 0
		elif string == 'U':
			temp = 1
		elif string == 'G':
			temp = 2
		elif string == 'C':
			temp = 3

		return temp

	import pandas as pd
	import sys

	#need to programming 'enter key' recognization
	data = sys.argv[1]
	
	#3 dimensional list to make Table which is ordered in 'ATGC' 3 times
	F = [[['K', 'N', 'K', 'N'], ['I', 'I', 'M', 'I'], ['R', 'S', 'R', 'S'], ['T', 'T', 'T', 'T']],
	     [['/', 'Y', '/', 'Y'], ['L', 'F', 'L', 'F'], ['/', 'C', 'W', 'C'], ['S', 'S', 'S', 'S']],
		 [['E', 'D', 'E', 'D'], ['V', 'V', 'V', 'V'], ['G', 'G', 'G', 'G'], ['A', 'A', 'A', 'A']],
		 [['Q', 'H', 'Q', 'H'], ['L', 'L', 'L', 'L'], ['R', 'R', 'R', 'R'], ['P', 'P', 'P', 'P']]]
		 

	#printing the table
	print("\n<Table of mRNA codon to Protein>")

	table = pd.DataFrame(data={
			'AA':F[0][0],'AU':F[0][1],'AG':F[0][2],'AC':F[0][3],
			'UA':F[1][0],'UU':F[1][1],'UG':F[1][2],'UC':F[1][3],
			'GA':F[2][0],'GU':F[2][1],'GG':F[2][2],'GC':F[2][3],
			'CA':F[3][0],'CU':F[3][1],'CG':F[3][2],'CC':F[3][3]
			},index=['A','U','G','C'], columns=['AA','AU','AG','AC','UA','UU','UG','UC','GA','GU','GG','GC','CA','CU','CG','CC'])
	table.index.name='1st base of codon'
	table.columns.name='2+3 base of codon'
		
	print(table)

	print("\n")
	idx = 0


	#finding 'ATG' codon(starting point)
	while idx<len(data):
		if data[idx] == 'A':
			idx += 1
			if data[idx] == 'U':
				idx += 1
				if data[idx] == 'G':
					break
		idx += 1
		
	#after found 'ATG', translating codon to protein
	for times in range(0, len(data)):
		l = times*3+idx-2	#first amino acid
		m = times*3+idx-1	#second amino acid
		n = times*3+idx		#third amino acid
		
		#index out of the range(the length of data)
		if n >= len(data) or m >= len(data) or l >= len(data):
			break
		
		#to find in table, we need to know the index of table
		i = apply_table(data[l])
		j = apply_table(data[m])
		k = apply_table(data[n])
		
		#stop codon
		if F[i][j][k] == '/':
			print()
			break
		#other protein codon
		else:
			print(F[i][j][k], end='')
	
	

if __name__ == '__main__':
	Translation()
