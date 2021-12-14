#this script parses all pdb files in working directory
#creates a new subdirectory with new versions of these pdb files,
#containing only the ligand of interest (specified in first arg on command line)
#as well as the residues with atoms within 4.1 angstroms of the ligand
#the number of lines in a pdbfile for the ligand must also be specified
#as the second argument on command line. this is used to help grab only one
#copy of the ligand from the input pdb files, in case they contain multiple
#copies bound to duplicate chains
#time ipython ~/desktop/tools/contact_pdb.py pmp 16

#how about I just make a contact fuzzball for the lig from native structures?
'''
1.identify ONE copy of the ligand
2.find residues with atoms within cutoff of that copy
3.export the lig with contact sc, just do it as a bunch of files
'''
import os
import math
import sys
#define function to get indices of residues in list of pdbfile lines
def return_residue_indices(fuzzball_lines):
    fuzzball_residue_indices=[];fuzzball_seen_resnums=[]
    starts=[];lasts=[]
    for index, line in enumerate(fuzzball_lines): #collecting start/end indices for unique res
            resnum=int(fuzzball_lines[index][22:26])
            resname=line[17:20]
            try:
                lastresnum=int(fuzzball_lines[index-1][22:26])
                lastresname=fuzzball_lines[index-1][17:20]
                if resnum!=lastresnum or resname!=lastresname:
                    start=index
                    starts.append(start)
            except:
                start=index
                starts.append(start)
            try:
                nextresname=fuzzball_lines[index+1][17:20]
                next_resnum=int(fuzzball_lines[index+1][22:26])
                if resnum!=next_resnum or resname!=nextresname:
                    last=index+1
                    lasts.append(last)
            except:
                last=len(fuzzball_lines)
                lasts.append(last)
    for index,start in enumerate(starts): #put the indices together for each res
        fuzzball_residue_indices.append((start,lasts[index]))
    fuzzball_residue_indices=list(set(fuzzball_residue_indices)) #ensure no redundant indices
    fuzzball_residue_indices=sorted(fuzzball_residue_indices, key=lambda first: first[0])
    return fuzzball_residue_indices

#set of functions to get distance between two points
def displace(p1,p2):
	x = p1[0] - p2[0]
	y = p1[1] - p2[1]
	z = p1[2] - p2[2]
	return (x,y,z)
def norm(x):
    return math.sqrt(sum(i**2 for i in x))
def dist(p1,p2):
	v = displace(p1,p2)
	return norm(v)

os.mkdir('contact_pdbs')
pdbfiles=[file for file in os.listdir(os.getcwd()) if file[-3:]=='pdb']
for file in pdbfiles:
    lig_atoms=[]
    ptn_atoms=[]
    f=open(file,'r')
    for line in f.readlines():
        if line[0:6]=='HETATM':
            if line[17:20]==sys.argv[1] or line[17:20]==sys.argv[1].upper():
                lig_atoms.append(line)
        elif line[0:4]=='ATOM':
            ptn_atoms.append(line)
    f.close()
    lig_atoms=lig_atoms[0:int(sys.argv[2])]
    output_pdb_lines=[]
    res_indices=return_residue_indices(ptn_atoms)
    for (a,b) in res_indices:
        current_residue_lines=[i for i in ptn_atoms[a:b]]
        for line in current_residue_lines:
            ptn_atom_coords=( float(line[31:39]),float(line[39:47]),float(line[47:54]) )
            for atom in lig_atoms:
                lig_atom_coords=( float(atom[31:39]),float(atom[39:47]),float(atom[47:54]) )
                d=dist(ptn_atom_coords, lig_atom_coords)
                print(d)
                if d<4.1:
                    output_pdb_lines.append(current_residue_lines)
    clean_output_pdb_lines=[]
    for l in output_pdb_lines:
        if l not in clean_output_pdb_lines:
            clean_output_pdb_lines.append(l)
    clean_output_pdb_lines.append(lig_atoms)
    of=open('contacts_'+file,'w')
    for l in clean_output_pdb_lines:
        for line in l:
            of.write(line)
    of.close()
    os.rename('contacts_'+file,'contact_pdbs/'+'contacts_'+file)





'''
#imports
import pandas as pd
import os
import numpy as np
import sys
import Bio
from Bio import SeqUtils
from Bio.SeqUtils import IUPACData
#list pdbfiles in working directory
pdbfiles=[file for file in os.listdir(os.getcwd()) if file[-3:]=='pdb']
#pdbfile to dataframe
def pdb_to_df(pdbfile):
    l=[]
    f=open(pdbfile,'r')
    for line in f.readlines():
        if line[0:4]=='ATOM':
            d={}
            d['recordname']=line[0:6]
            d['atomnumber']=line[6:11]
            d['atomname']=line[11:16]
            d['altloc']=line[16:17]
            d['resname']=line[17:20]
            d['chain']=line[20:22]
            d['resnum']=line[22:29]
            d['achar']=line[29:31]
            d['x']=line[31:39]
            d['y']=line[39:47]
            d['z']=line[47:54]
            d['occupancy']=line[54:60]
            d['temp_factor']=line[60:66]
            d['seg']=line[66:76]
            d['element']=line[76:78]
            d['q']=line[78:80]
            l.append(d)
        elif line[0:6]=='HETATM':
            if line[17:20]=='pmp' or line[17:20]=='pmp'.upper():
                d={}
                d['recordname']=line[0:6]
                d['atomnumber']=line[6:11]
                d['atomname']=line[11:16]
                d['altloc']=line[16:17]
                d['resname']=line[17:20]
                d['chain']=line[20:22]
                d['resnum']=line[22:29]
                d['achar']=line[29:31]
                d['x']=line[31:39]
                d['y']=line[39:47]
                d['z']=line[47:54]
                d['occupancy']=line[54:60]
                d['temp_factor']=line[60:66]
                d['seg']=line[66:76]
                d['element']=line[76:78]
                d['q']=line[78:80]
                l.append(d)
    df=pd.DataFrame(l)
    return df
#convert dataframe to pdb file
def df_to_pdb(dataframe,ofile):
    of=open(ofile,'w')
    for i in range(dataframe.shape[0]):
        line=''.join(dataframe.iloc[i])
        of.write(line)
    of.close()
#do stuff to pdbfiles
for file in pdbfiles[0:1]:
    df=pdb_to_df(file)
    chains=list(set(list(df['chain'])))
    for chain in chains:
        seq1=''
        chain_indices=np.where(df['chain']== chain)
        seq3=list(df.iloc[chain_indices]['resname'])



    for i in range(df.shape[0]):
        chainid=df.iloc[i][5]

######################################################################
import os
pdbfiles=[file for file in os.listdir(os.getcwd()) if file[-3:]=='pdb']

fuzzball_lines=[]
for file in pdbfiles[0:1]: #
    f=open(file,'r')
    for line in f.readlines():
        if line[0:4]=='ATOM':
            fuzzball_lines.append(line)
        elif line[0:6]=='HETATM':
            if line[17:20]!='HOH':
                fuzzball_lines.append(line)
#get indices of unique residues in fuzzball_lines
fuzzball_residue_indices=[];fuzzball_seen_resnums=[]
starts=[];lasts=[]
for index, line in enumerate(fuzzball_lines): #collecting start/end indices for unique res
        resnum=int(fuzzball_lines[index][22:26])
        resname=line[17:20]
        try:
            lastresnum=int(fuzzball_lines[index-1][22:26])
            lastresname=fuzzball_lines[index-1][17:20]
            if resnum!=lastresnum or resname!=lastresname:
                start=index
                starts.append(start)
        except:
            start=index
            starts.append(start)
        try:
            nextresname=fuzzball_lines[index+1][17:20]
            next_resnum=int(fuzzball_lines[index+1][22:26])
            if resnum!=next_resnum or resname!=nextresname:
                last=index+1
                lasts.append(last)
        except:
            last=len(fuzzball_lines)
            lasts.append(last)
for index,start in enumerate(starts): #put the indices together for each res
    fuzzball_residue_indices.append((start,lasts[index]))
fuzzball_residue_indices=list(set(fuzzball_residue_indices)) #ensure no redundant indices
fuzzball_residue_indices=sorted(fuzzball_residue_indices, key=lambda first: first[0])
'''

#####okay just separate the residue indices by chain
#####then iterate resnames for each chain and store as strings to compare
######may have to handle ligands independently since no standard naming convention
######ehhh idk this is so much time for something that should seemingly be simple

'''
CANT FIGURE OUT HOW TO COPY DESIRED CHAINS AND DUMP TO PDB

from pyrosetta import *
init()
from Levenshtein import distance as lev_dist

pdbfiles=[file for file in os.listdir(os.getcwd()) if file[-3:]=='pdb']

for file in pdbfiles:
    p=pose_from_pdb(file)

#get the chain sequences
chain_seqs=[]
for x in pyrosetta.rosetta.core.pose.get_chains(p):
    l=list(pyrosetta.rosetta.core.pose.get_resnums_for_chain_id(p,x))
    first_res=l[0]
    last_res=l[len(l)-1]
    chain_sequence=p.sequence(first_res,last_res)
    chain_seqs.append(chain_sequence)

#get unique chains only, if levenstein distance less than 10 call chains duplicates
unique_chains=[]
for i in chain_seqs:
    unique_chains.append(i)

for ind,i in enumerate(chain_seqs):
    try:
        for k in chain_seqs[ind+1:]:
            if i==k or lev_dist(i,k)<10:
                print(i,k)
                try:
                    unique_chains.remove(k)
                except:
                    print('?')
    except:
        print('all done')



ligand_pose=p.residue(p.total_residue()).clone()
res_pose=p.residue(resnum).clone()
contact_pose.append_residue_by_jump(res_pose, 1)
contact_pose.append_residue_by_jump(ligand_pose, 1)

p.residue(1).name()
p.dump_pdb(desired name)
'''
