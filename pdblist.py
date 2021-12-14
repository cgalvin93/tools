#writes file names of all pdb files in working directory to a text file
#that can be fed to rosetta for scoring multiple structures
#USAGE:
#ipython ~/Desktop/tools/pdblist.py pdblist.txt
import sys
import os

ofilename=sys.argv[1]
ofile=open(ofilename,'w')

pdbfiles=[file for file in os.listdir(os.getcwd()) if file[-3:]=='pdb']

for strc in pdbfiles:
	s = str(strc).strip()
	ofile.write(s+'\n')

ofile.close()
