#time ipython pps.py scorefile native_seqs outfile pdbID
#must be run in directory with design pdb files so it can get sequences
#to calculate pps with
#produces pps plots for all metrics

#time ipython ~/tools/pps.py 1f4p_all.sc ~/filt/nat/native_seqs/1f4p_binding_site.txt pps1f4pfdgenpot.pdf 1f4p

#time ipython ~/desktop/tools/pps.py cm1f4p.sc ~/desktop/prj/filter_comparison/native_seqs/1f4p_binding_site.txt pps1f4pcm.pdf 1f4p
#time ipython ~/tools/pps.py 2xbnfdgenpot.sc ~/filt/nat/native_seqs/2xbn_binding_site.txt pps2xbnfdgenpot.pdf 2xbn
#time ipython ~/tools/pps.py 1zk4fdgenpot.sc ~/filt/nat/native_seqs/1zk4_binding_site.txt pps1zk4fdgenpot.pdf 1zk4


import os
import sys
import Bio
from Bio import SeqUtils
from Bio.SeqUtils import IUPACData
import numpy as np
import scipy
from scipy import spatial
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#must put binding site res numbers here manually
#store pdb numbering of binding site residues for query protein
ptn_pdb_id=sys.argv[4]
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

#need this later
p_nat_letters_dict={'A':0,
                    'D':1,
                    'N':2,
                    'R':3,
                    'C':4,
                    'F':5,
                    'G':6,
                    'E':7,
                    'Q':8,
                    'H':9,
                    'I':10,
                    'L':11,
                    'K':12,
                    'M':13,
                    'P':14,
                    'S':15,
                    'W':16,
                    'Y':17,
                    'T':18,
                    'V':19}

#store lines of scorefile
scorefile=open(sys.argv[1],'r')
lines=[line for line in scorefile.readlines()[1:]]
scorefile.close()
#acquire the names of strc in scorefile
strcs=[]
scores=[]
for line in lines[1:]:
    strc=line.split()[-1]
    score=line.split()[1]
    strcs.append(strc+'.pdb')
    scores.append(score)

#convert 3 letter aa code to 1 letter
def threeto1(s):
    olseq=''
    for i in range(0,len(s),3):
        res = s[i:i+3]
        threelc=res[0]+res[1:3].lower()
        onelc = Bio.SeqUtils.IUPACData.protein_letters_3to1[threelc]
        olseq+=onelc
    return olseq


#from a list of pdb files, return binding site sequences from those files
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
        seqs.append(s2)
    return seqs


#return matrix of probability distribution for each amino acid at each position
#uses pseudo counts to avoid any P=0
#rows = amino acids
#columns = frequencies at each position
def p_dist(seqs):
    n_positions=len(seqs[0])
    n_seqs=len(seqs)
    a=np.zeros(shape=(20,n_positions))
    ala=[];asp=[];asn =[];arg =[];cys=[];phe =[];gly=[];glu =[];gln =[];his =[];
    ile =[];leu=[];lys =[];met =[];pro =[];ser =[];trp =[];tyr =[];thr=[];val=[]
    for i in range(0,n_positions):
        nala=1;nasp=1;nasn=1;narg=1;ncys=1;nphe=1;ngly=1;nglu=1;ngln=1;nhis=1;
        nile=1;nleu=1;nlys=1;nmet=1;npro=1;nser=1;ntrp=1;ntyr=1;nthr=1;nval=1;
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
        ntot=float((n_seqs-ngap)+20);fala=nala/ntot
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

#get all the sequences from a fasta file
def capture_seqs(file):
    indices=[]
    f=open(file,'r')
    lines=[line for line in f.readlines()]
    f.close()
    for line in lines:
        if line[0]=='>':
            indices.append(lines.index(line))
    strs=[]
    for i in range(0,len(indices)-1):
        start=indices[i]
        end=indices[i+1]
        s=''
        for line in lines[start+1:end]:
            s+=str(line.strip('\n'))
        strs.append(s.upper())
    return list(set(strs))


#calculate the pps score for each position in 2 distributions
def calc_pps(p,q):
    n_pos = p.shape[1]
    pps_scores=[]
    for i in range(0,n_pos):
        jsd = (scipy.spatial.distance.jensenshannon(p[:,i],q[:,i]))**2
        pps=1-jsd
        pps_scores.append(pps)
    return pps_scores


#open the output pdf file to put all the graphs on
outfilename=sys.argv[3]
pdf = PdfPages(outfilename)

#do pps stuff
native_seqs=capture_seqs(sys.argv[2])
p_nat=p_dist(native_seqs)
seqs=return_bs_seq(strcs)
p_mut=p_dist(seqs)
pps_scores=calc_pps(p_nat,p_mut)

#do bs probs
all_pos_prob_avgs=[]
for seq in seqs:
    positional_probs=[]
    for index,aaletter in enumerate(seq):
        p_nat_index_1=p_nat_letters_dict[aaletter]
        probability_aa_pos=p_nat[p_nat_index_1][index]
        positional_probs.append(probability_aa_pos)
    sum_p=sum(positional_probs)
    avg_pos_prob=sum_p/float(len(seq))
    all_pos_prob_avgs.append(avg_pos_prob)
'''
needs path to original structures
r'/wynton/home/kortemme/cgalvin/filt/fd/pdb/'
/Users/Student/Desktop/amandas/
'''
amandas_path=r'/wynton/home/kortemme/cgalvin/filt/fd/pdb/'+ptn_pdb_id
amandas_strc_1=[i for i in os.listdir(amandas_path) if i[-4:]=='.pdb']
amandas_strc_2=[amandas_path+'/'+str(i) for i in amandas_strc_1]
amandas_seqs=return_bs_seq(amandas_strc_2)
p_amanda=p_dist(amandas_seqs)
pps_amanda=calc_pps(p_nat,p_amanda)

#
term_scores=[]
term_scores.append(pps_scores)
term_scores.append(pps_amanda)


#plot pps
fig,ax=plt.subplots()
bp_dict = ax.boxplot(term_scores)
ax.set_title(str(ptn_pdb_id))
ax.set_ylabel('PPS')
xlocs=[x+1 for x in range(len(term_scores))]
xlabs=['Gen_pot Designs','OG Designs']
plt.xticks(xlocs, xlabs)
for line in bp_dict['medians']:
    x,y = line.get_xydata()[1]
    plt.text(x,y, '%.3f' % y,verticalalignment='center')
pdf.savefig()
plt.clf()

#plot bs probs
y=[float(i) for i in scores];x=[float(k) for k in all_pos_prob_avgs]

current_corr=scipy.stats.pearsonr(x,y)
rval=current_corr[0];pval=current_corr[1]
plt.scatter(x,y)
plt.xlabel('Average Probability Across Designed Positions')
plt.ylabel('Total_Score')
plt.title('R='+'{:.4}'.format(current_corr[0]))
pdf.savefig()
plt.clf()

pdf.close()
