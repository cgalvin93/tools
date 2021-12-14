import os
import Bio
from Bio import SeqUtils
from Bio.SeqUtils import IUPACData

def threeto1(s):
    olseq=''
    for i in range(0,len(s),3):
        res = s[i:i+3]
        threelc=res[0]+res[1:3].lower()
        onelc = Bio.SeqUtils.IUPACData.protein_letters_3to1[threelc]
        olseq+=onelc
    return olseq

def return_bs_res(ptn_pdb_id):
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
    return binding_site_res

def return_bs_seq(pdbfiles):
    seqs=[]
    for file in pdbfiles:
        fullseq=[]
        f=open(file,'r');lines=[line for line in f.readlines()]
        for line in lines:
            if line[0:4]=='ATOM':
                resname=line[17:21].strip()
                resnum = line[22:27].strip()
                fullseq.append((resnum, resname))
        fullseq=list(set(fullseq))
        s=''
        for i in binding_site_res:
            for num, name in fullseq:
                if i==int(num):
                    s+=name
        s2=threeto1(s)
        seqs.append((file,s2))
    return seqs



names=['1f4p','2xbn','3dk9','3s6f','1zk4']
for dir in names:
  os.chdir(dir)
  os.chdir('cm2')
  binding_site_res=return_bs_res(dir)
  f=open('pdblist.txt','r')
  pdbnames=[line for line in f.readlines()]
  f.close()
  cleannames=[]
  for strc in pdbnames:
      cleannames.append(strc.replace('\n',''))
  seqs=return_bs_seq(cleannames)
  ofile=open(dir+'_bsseqs.txt','w')
  for name, sequence in seqs:
    ofile.write('>'+name + '\n')
    length=len(sequence)
    remainder=length%60; n_seg=length/60
    indices=[]
    if length-remainder!=0:
        for i in range(0,60,length-remainder):
            indices.append(i)
    else:
        for i in range(0,60,length):
            indices.append(i)
    for i,x in enumerate(indices[:-1]):
        start=x
        end=indices[i+1]
        s=sequence[start:end]+'\n'
        ofile.write(s)
    last_lim=indices[-1]
    last_s=sequence[last_lim:length]
    ofile.write(str(last_s) +'\n')
  ofile.close()
  os.system('cp '+dir+'_bsseqs.txt ~/cm2seqs/'+dir+'.txt')
  os.chdir('..')
  os.chdir('..')
