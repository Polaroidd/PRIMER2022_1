#!/usr/bin/env python3

def Transcription():
	import sys

	data = list(sys.argv[1])
	seq = ''

	for i in range(0,len(data)):
		if data[i] == 'T':
			data[i] = 'U'
		seq += data[i]
	
	print(seq)

if __name__ == '__main__':
	Transcription()
