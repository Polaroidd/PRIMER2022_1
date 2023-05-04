#!/usr/bin/env python3
def write_txt():
	import sys
	
	RNA = sys.argv[1]
	each_line_100 = len(RNA)//100

	name = input("Write the name of species : ")

	newfile = open("./Gene.txt", 'w')
	
	newfile.write(">%s\n" %name)

	letter = -1
	for i in range(each_line_100):
		for j in range(100):
			letter += 1
			newfile.write(RNA[letter])
		newfile.write("\n")
	
	if len(RNA)%100 != 0:
		for m in range(len(RNA)%100):
			letter += 1
			newfile.write(RNA[letter])
	
	newfile.close()

if __name__ == '__main__':
	write_txt()
