#time ipython seq_entropy.py pdbid

import os
import sys
import numpy as np
import scipy
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#must put binding site res numbers here manually
#store pdb numbering of binding site residues for query protein
ptn_pdb_id=sys.argv[1]
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

from pyrosetta import *
# init()
init('-ignore_unrecognized_res -load_PDB_components False')
pdbs=[i for i in os.listdir() if i[-3:]=='pdb' and i[-13:]!='0001_0001.pdb']
all_bss=[]
bss_length=len(binding_site_res)


for design in pdbs:##########
    p=pose_from_pdb(design)
    s=p.sequence()
    bss=''
    for res in binding_site_res:
        if ptn_pdb_id=='1f4p':
            bss+=s[int(res)-2]#starts at res 2
        elif ptn_pdb_id=='2xbn':
            bss+=s[int(res)-22]#starts at res 22
        elif ptn_pdb_id=='3dk9':
            bss+=s[int(res)-17]#starts at res 17
        else:
            bss+=s[int(res)-1]
    all_bss.append(bss)

def p_dist(seqs):
    n_positions=len(seqs[0])
    n_seqs=len(seqs)
    a=np.zeros(shape=(20,n_positions))
    ala=[];asp=[];asn =[];arg =[];cys=[];phe =[];gly=[];glu =[];gln =[];his =[];
    ile =[];leu=[];lys =[];met =[];pro =[];ser =[];trp =[];tyr =[];thr=[];val=[]
    for i in range(0,n_positions):
        nala=0;nasp=0;nasn=0;narg=0;ncys=0;nphe=0;ngly=0;nglu=0;ngln=0;nhis=0;
        nile=0;nleu=0;nlys=0;nmet=0;npro=0;nser=0;ntrp=0;ntyr=0;nthr=0;nval=0;
        ngap=0
        for seq in seqs:
            aa=seq[i]
            if aa=='A':
                nala+=1
            elif aa=='C':
                ncys+=1
            elif aa=='D':
                nasp+=1
            elif aa=='E':
                nglu+=1
            elif aa=='F':
                nphe+=1
            elif aa=='G':
                ngly+=1
            elif aa=='H':
                nhis+=1
            elif aa=='I':
                nile+=1
            elif aa=='K':
                nlys+=1
            elif aa=='L':
                nleu+=1
            elif aa=='M':
                nmet+=1
            elif aa=='N':
                nasn+=1
            elif aa=='P':
                npro+=1
            elif aa=='Q':
                ngln+=1
            elif aa=='R':
                narg+=1
            elif aa=='S':
                nser+=1
            elif aa=='T':
                nthr+=1
            elif aa=='V':
                nval+=1
            elif aa=='W':
                ntrp+=1
            elif aa=='Y':
                ntyr+=1
            else:
                ngap+=1
        ntot=float((n_seqs-ngap));fala=nala/ntot
        fasp=nasp/ntot;fasn=nasn/ntot;farg=narg/ntot;
        fcys=ncys/ntot;fphe=nphe/ntot;fgly=ngly/ntot;
        fglu =nglu/ntot;fgln=ngln/ntot;fhis=nhis/ntot;
        ffile=nile/ntot;fleu=nleu/ntot;flys=nlys/ntot;
        fmet=nmet/ntot;fpro=npro/ntot;fser=nser/ntot;
        ftrp=ntrp/ntot;ftyr=ntyr/ntot;fthr=nthr/ntot;
        fval=nval/ntot
        ala.append(fala);asp.append(fasp);asn.append(fasn);arg.append(farg);
        cys.append(fcys);phe.append(fphe);gly.append(fgly);glu.append(fglu);
        gln.append(fgln);his.append(fhis);ile.append(ffile);leu.append(fleu);
        lys.append(flys);met.append(fmet);pro.append(fpro);ser.append(fser);
        trp.append(ftrp);tyr.append(ftyr);thr.append(fthr);val.append(fval)
    a[0]=ala;a[1]=asp;a[2]=asn;a[3]=arg;a[4]=cys;a[5]=phe;a[6]=gly
    a[7]=glu;a[8]=gln;a[9]=his;a[10]=ile;a[11]=leu;a[12]=lys;a[13]=met
    a[14]=pro;a[15]=ser;a[16]=trp;a[17]=tyr;a[18]=thr;a[19]=val
    return a



p=p_dist(all_bss)


all_pos_freqs=[]#dist of amino acid frequencies at each position
for i in range(p.shape[1]):
    positional_freqs=[]
    for x in range(p.shape[0]):
        positional_freqs.append(p[x][i])
    all_pos_freqs.append(positional_freqs)

entropies=[]
for i in all_pos_freqs:
    entropies.append(scipy.stats.entropy(i))

ofile=open('ogentropies.txt','a')
ofile.write('\n')
ofile.write(str(ptn_pdb_id))
ofile.write('\n')
ofile.write(str(entropies))
ofile.close()
