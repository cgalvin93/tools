#time ipython find_redundant_structures.py paramsfile
#time ipython ~/desktop/tools/find_redundant_structures.py

from pyrosetta import *
init()
import os
import sys
from Levenshtein import distance as lev_dist
import Bio
from Bio import SeqIO as seqio


allseqs={}
pdbfiles=[file for file in os.listdir(os.getcwd()) if file[-3:]=='pdb']
#make fasta files for all pdb files in working directory, put them in fragment_picker dir
for file in pdbfiles:
    seqs=[]
    for record in seqio.parse(file,'pdb-atom'):
        s=record.seq
        seqs.append((file[:-4],s))
    count=1
    for name, sequence in seqs:
        allseqs[name+str(count)]=sequence
        count+=1



#store names of all structures, remove homologous as they are IDed
allchains=[]
for key in allseqs.keys():
    if key[:-1] not in allchains:
        allchains.append(key[:-1])

#remove homologs from allchains
homologs=[]
for chain1 in allseqs.keys():
    for chain2 in allseqs.keys():
        if chain1[:-1]!=chain2[:-1]:
            ldistance=lev_dist(str(allseqs[chain1]),str(allseqs[chain2]))
            alength=max([len(str(allseqs[chain1])),len(str(allseqs[chain2]))])
            if ldistance==0:
                if (chain2[:-1],chain1[:-1]) not in homologs:
                    try:
                        allchains.remove(chain2[:-1])
                        homologs.append((chain1[:-1],chain2[:-1]))
                    except:
                        pass
            else:
                lratio=ldistance/alength
                if lratio<0.2:
                    if (chain2[:-1],chain1[:-1]) not in homologs:
                        try:
                            allchains.remove(chain2[:-1])
                            homologs.append((chain1[:-1],chain2[:-1]))
                        except:
                            pass

print(homologs)
print(len(allchains))
os.mkdir('nonhomologous')
for strc in allchains:
    os.rename(strc+'.pdb','nonhomologous/'+strc+'.pdb')



'''
                chainstoremove=[]
                for key in allseqs.keys():
                    if key[:-1]==chain2[:-1]:
                        chainstoremove.append(key)
                for chain in chainstoremove:
                    allseqs.pop(chain)
'''
