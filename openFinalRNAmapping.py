#!/usr/bin/env python3

def printFinal():
	import sys

	final = "./{}_stemloop_final_result.txt".format(sys.argv[1][0:-4])
	for line in open(final):
		line = line.rstrip()
		if line[0] == ">":
			print("Name : %s\n<Final Mapping>\n" %line[1:])
		else:
			print(line)

if __name__ == '__main__':
	printFinal()
