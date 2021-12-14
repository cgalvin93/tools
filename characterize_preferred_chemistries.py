#time ipython characterize_preferred_chemistries.py fasta.fasta pdbID

import sys
import Bio
from Bio import SeqUtils
from Bio.SeqUtils import seq3

ptn_pdb_id=sys.argv[2]
#define binding site residues
if ptn_pdb_id == '2xbn':
    binding_site_res=[73, 133, 134, 135, 138, 159,161, 202, 231, 233,
                      234, 262,264, 265]
elif ptn_pdb_id == '1f4p':
    binding_site_res=[10, 11, 12, 13, 14, 15, 58, 59,
                      60, 61, 62, 68, 93, 94, 95, 98,
                      100, 101, 102]
elif ptn_pdb_id == '1zk4':
    binding_site_res=[13, 14, 15, 16, 17, 18, 19, 36,
                      37, 38, 62, 63, 89, 90, 91, 92,
                      112, 140, 142, 155, 159]
elif ptn_pdb_id == '3dk9':
    binding_site_res=[26, 27, 28, 29, 30, 31, 49, 50,
                      51, 52, 56, 57, 62,
                      66, 129, 130, 155, 156, 157,
                      177, 181, 197, 198, 201, 202,
                      291, 294, 298, 330, 331]
elif ptn_pdb_id == '3dlc':
    binding_site_res=[48, 49, 50,
                      51, 52, 53, 55, 72, 73, 74, 77,
                      100, 101, 102, 117, 118, 119,
                      122, 123]
elif ptn_pdb_id == '3r2q':
    binding_site_res=[9, 10, 11, 33, 34, 48, 49, 50,
                      62, 63, 64]
elif ptn_pdb_id == '3s6f':
    binding_site_res=[78, 79, 80, 85, 86, 87, 88, 90, 91, 114, 115, 118, 119, 121]
elif ptn_pdb_id == '1zk4':
    binding_site_res=[13, 14, 15, 16, 17, 18, 19, 36, 37, 38, 62, 63, 89, 90,
                      91, 92, 112, 140, 142, 155, 159, 187, 189, 190, 192, 194,
                      195, 205]
elif ptn_pdb_id == '2ij2':
    binding_site_res=[69,75,86,87,96,100,107,153,260,261,264,265,268,269,272,
                      322,327,328,331,357,392,393,394,398,399,400,401,402,405,406]
else:
    print('the binding site residues of the query protein cannot be identified')

#get sequences from fasta file
f=open(sys.argv[1],'r')
lines=[line for line in f.readlines()]
f.close()
seqs=[]
for line in lines:
    if line[0]!='>':
        seqs.append(line)

#convert 1 letter to 3 letter and return as list rather than single string
def one_to_three(s):
    tls=Bio.SeqUtils.seq3(s.upper().strip('\n'))
    tls_list=[]
    for i in range(0,len(tls),3):
        res = tls[i:i+3]
        tls_list.append(res)
    return tls_list

#convert to 3letter code
all_3l_seqs=[]
for seq in seqs:
    all_3l_seqs.append(one_to_three(seq))

#make a dictionary position=[res1,res2....] for all positions in seq profile
n_positions=len(all_3l_seqs[0])
positional_seqs={}
for i in range(n_positions):
    position_seqs=[]
    for x in all_3l_seqs:
        position_seqs.append(x[i])
    positional_seqs[binding_site_res[i]]=position_seqs

#define aa chemistries
polar=['SER','THR','ASN','GLN']
hydrophobic=['ALA','ILE','LEU','VAL','MET','CYS']
aromatic=['PHE','TYR','TRP']
charged=['GLU','ASP','ARG','LYS','HIS']
gly=['GLY']
pro=['PRO']

chemslist=[]
#go through each position, find preferred chemistry, strength of preference
for key in positional_seqs.keys():
    npolar=1;nnonpolar=1;naromatic=1;ncharged=1;ngly=1;npro=1 #pseudocounts
    for res in positional_seqs[key]:
        if res.upper() in polar:
            npolar+=1
        elif res.upper() in hydrophobic:
            nnonpolar+=1
        elif res.upper() in aromatic:
            naromatic+=1
        elif res.upper() in charged:
            ncharged+=1
        elif res.upper() in gly:
            ngly+=1
        elif res.upper() in pro:
            npro+=1
    chems=[npolar,nnonpolar,naromatic,ncharged,ngly,npro]
    if chems.index(max(chems))==0:
        preferred_chem='polar'
    elif chems.index(max(chems))==1:
        preferred_chem='hydrophobic'
    elif chems.index(max(chems))==2:
        preferred_chem='aromatic'
    elif chems.index(max(chems))==3:
        preferred_chem='charged'
    elif chems.index(max(chems))==4:
        preferred_chem='glycine'
    elif chems.index(max(chems))==5:
        preferred_chem='proline'
    sorted_chems=sorted(chems,reverse=True)
    n_second_preference=sorted_chems[1]
    if chems.index(n_second_preference)==0:
        second_preference='polar'
    elif chems.index(n_second_preference)==1:
        second_preference='hydrophobic'
    elif chems.index(n_second_preference)==2:
        second_preference='aromatic'
    elif chems.index(n_second_preference)==3:
        second_preference='charged'
    elif chems.index(n_second_preference)==4:
        second_preference='glycine'
    elif chems.index(n_second_preference)==5:
        second_preference='proline'
    preference_strength=float(max(chems)/n_second_preference)
    print(preferred_chem);print(second_preference)
    print(preference_strength)
    chemslist.append((preferred_chem,preference_strength))

#write results
ofile=open('preferred_chems.txt','a')
ofile.write('\n')
ofile.write(ptn_pdb_id)
ofile.write('\n')
ofile.write(str(chemslist))
ofile.close()
