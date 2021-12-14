#writes a resfile for an input pdb file, requires text file specifying which
#residues to allow packing or design, fixes current sidechain of other residues

#ipython ~/desktop/tools/resfile.py protein.pdb options.txt out.resfile

#options file format is desiredcode,num1chain,num2chain,etc. with diff code on newline
#ex.
#ALLAA,1A,2B
#NATAA,3A,4A,1X
'''
appears to be an issue with the script where the earliest residue specified
in options file just gets natro instead of what its supposed to
'''
#time ipython ~/desktop/tools/resfile.py protein.pdb options.txt out.resfile
#time ipython ~/desktop/tools/resfile.py clean_3s6f.pdb 3s6f_res_opt.txt 3s6f.resfile
#time ipython ~/desktop/tools/resfile.py clean_1zk4.pdb 1zk4resopts.txt 1zk4.resfile
#time ipython ~/desktop/tools/resfile.py clean_1F4P_wt.pdb opts.txt 1f4p.resfile

'''
ALLAA ................ allow all 20 amino acids INCLUDING cysteine
POLAR ................ allow only canonical polar amino acids (DEHKNQRST)
APOLAR ............... allow only canonical non polar amino acids
PROPERTY <property> .. disallow any residue type that lacks the given property
NOTAA <list of AAs> .. disallow only the specified amino acids ( Use one letter codes, undelimited like ACFYRT. For NCAAs, use X[<full name>]. )
PIKAA <list of AAs> .. allow only the specified amino acids ( Use one letter codes, undelimited like ACFYRT. For NCAAs, use X[<full name>].)
NATAA ................ allow only the native amino acid (NATive Amino Acid) - repack without design
NATRO ................ preserve the input rotamer ( do not pack at all) (NATive ROtamer)
'''
import sys

ptnfile = sys.argv[1]
pdbfile = open(ptnfile, 'r')

optfile = sys.argv[2]
options = open(optfile, 'r')

outfilename = sys.argv[3]
ofile = open(outfilename,'w')

#parse pdb file to get list of str w resnum, icode, chain, for each res
rawres=[]
for line in pdbfile.readlines():
	if line[0:4]=='ATOM':
		resnum = line[22:27]
		chain = line[21:22]
		rawres.append((resnum,chain))
residues=[]
for i,x in rawres:
	if (i,x) not in residues:
		residues.append((i,x))

#parse options file to get lists of res w certain options
allaa=[]
nataa=[]
for line in options.readlines():
	if line[0:5]=='ALLAA':
		l = line.split(',')
		l2=[i.strip() for i in l[1:]]
		for i in l2:
			allaa.append(i)
	elif line[0:5]=='NATAA':
		l = line.split(',')
		l2=[i.strip() for i in l[2:]]
		for i in l2:
			nataa.append(i)

#write header
ofile.write('USE_INPUT_SC')
ofile.write('\nstart'+'\n')

#write line for each residue + it's corresponding code
for i,x in residues:
	res_id = str(i.strip())+str(x)
	if res_id in allaa:
		s = str(i) + ' ' + str(x) + ' ALLAA'	#resnum of diff length puts char in diff place on line!!!!
		ofile.write(s+'\n')
	elif res_id in nataa:
		s = str(i) + ' ' + str(x) + ' NATAA'
		ofile.write(s+'\n')
	else:
		s = str(i) + ' ' + str(x) + ' NATRO'
		ofile.write(s+'\n')

ofile.close()

'''
ALLAA           # allow all amino acids
EX 1 EX 2       # allow extra chi rotamers at chi-id 1 and 2 (note: multiple commands can be on the same line.)
USE_INPUT_SC    # allow the use of the input side chain conformation   ( see below for more detailed description of commands)
start
<PDBNUM>[<ICODE>] <CHAIN>  <COMMANDS> 	#<PDBNUM>[<ICODE>] corresponds to columns 22-26
40A Q ALLAA   # Residue 40, insertion code A, on chain Q, use any residue type
'''
