#input scorefile path 
sf_path='~/desktop/genpottest2.sc'


#parse the scorefile to create a list of tuples where each element is (structure_name, scoreterm_name,scoreterm_value)
sf=open(sf_path,'r') 
lines=[line for line in sf.readlines()]
sf.close()
terms=[]
for i in lines[1].split()[1:]:
    terms.append(i)
strc_names=[]
datas=[]
for i,line in enumerate(lines[2:]):
    properties=line.split()[1:]
    try:
        name=properties[-1]
        if name not in strc_names:
            strc_names.append(name)
            for index,term in enumerate(terms):
                val=properties[index]
                datas.append((name,term,val))
    except:
        print('problem with entry in scorefile, index '+str(i))

print(terms) #print the names of the terms in the scorefile
#an explanation of many terms evaluated by the interfaceanalyzermover can be found in the 'expected outputs' section of:
#https://new.rosettacommons.org/docs/latest/application_documentation/analysis/interface-analyzer

print(len(strc_names)) #should be the number of designs in your scorefile
print(len(properties))#should be the number of metrics in your scorefile
print(len(datas))#should be the number of designs times the number of properties

#for a given term, print all values sorted in ascending order and
#plot their distribution
import matplotlib.pyplot as plt
import numpy as np

def plot_dist(datas, desired_term):
    scores=[]
    for name,term,val in datas:
        if term==desired_term:
            scores.append(float(val))
    print(sorted(scores))
    fig,ax=plt.subplots()
    ax.violinplot(scores)
    plt.title(desired_term)
    plt.ylabel('Score (REU)')
    plt.xlabel('designs')
    mean=np.mean(scores)
    text_s='MEAN: '+str(mean)
    plt.text(1,1,text_s,verticalalignment='top',transform=ax.transAxes)


plot_dist(datas, 'total_score') #input desired score term


#identify designs which meet some specified criteria
def filter(names,datas,desired_term,condition,threshold):
    passing_strc=[]
    for strc_name in names:
        for name,term,val in datas:
            if strc_name==name:
                if desired_term==term:
                    if condition=='>=':
                        if float(val)>=threshold:
                            passing_strc.append(strc_name)
                    elif condition=='<=':
                        if float(val)<=threshold:
                            passing_strc.append(strc_name)
    return passing_strc

#example:return designs with negative total score
negative_total_score=filter(strc_names,datas,'total_score','<=',0.0)
print(len(negative_total_score))
print(negative_total_score)

#pass the results of one filter through the function again to identify designs meeting multiple criteria
#eg return designs with negative total_score and good packstat score (>0.6):
negative_total_score_good_packstat=filter(negative_total_score,datas,'packstat','>=',0.6)
print(len(negative_total_score_good_packstat))
print(negative_total_score_good_packstat)

#a promising design might have a good total score (say, lower than the mean total_score),
#good binding energy (dG_separated)
#a high packstat score (>=0.6),
#no buried unstatified hbonds (buns_bb_heavy<=0),
#high shape complementarity between ligand and binding site (sc>=0.55),
#high buried nonpolar surface area (bsa>=mean), etc.
#good fragment quality


#make a pdf showing all score distributions
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

outfilename='scorefile_2-chains_Wynton_dists.pdf'
pdf = PdfPages(outfilename)
#for each term, plot the distribution of values
for term in terms:
    allscores=[]#score
    if term != 'description':
        for entry in datas:
            if entry[1]==term:
                allscores.append(float(entry[2]))
    fig,ax=plt.subplots()
    if len(allscores)!=0:
        ax.hist(allscores)#bins=int(len(allscores)/20)
        ax.set_title(term)
        ax.set_ylabel('frequency')
        ax.set_xlabel('score')
        pdf.savefig()
        plt.clf()
pdf.close()


#returning top nmetrics
sf=open('simon_ib2_designs_iam.sc','r')
lines=[line for line in sf.readlines()]
sf.close()
terms=[]
for i in lines[1].split()[1:-1]:
    terms.append(i)
strc_names=[]
datas=[]
for line in lines[2:]:
    properties=line.split()[1:]
    name=properties[-1]
    strc_names.append(name)
    for index,term in enumerate(terms):
        val=properties[index]
        datas.append((name,term,float(val)))


#terms in the score file that are more favorable when value higher must be
#specified:
higher_better=['packstat','sc','norm_bsa','bsa','dSASA_hphobic','dSASA_int','dSASA_polar',
               'hbonds_int','nres_int','acc','sc','design_fragqual_good/total',
               'design_fragqual_coverage','wt_diff_rmsd']


#find the lowest and highest values for each term,
#store in a list (term, low, high) that I can then feed to
#'scan feature' function that will do five intervals between low and hi
def scan(datas,terms):
    to_scan=[]
    for term in terms:
        if term != 'description':
            vals=[]
            for entry in datas:
                if entry[1]==term:
                    vals.append(float(entry[2]))
            high=max(vals); low=min(vals)
            to_scan.append((term,low,high))
    return to_scan

to_scan=scan(datas,terms)

#find four equally spaced values between high and low value
def foursteps(min, max):
    d = max - min
    int_size = d/5
    p1=min+int_size;p2=min+2*int_size;p3=min+3*int_size
    p4=min+4*int_size
    threshvals=[p1,p2,p3,p4]
    return threshvals

#
ntop={}
nametop={}
conditions={}
for theterm,lowval,hival in to_scan:
    vals=foursteps(lowval,hival)
    if theterm in higher_better:
        conditions[theterm]=vals[len(vals)-1]
    else:
        conditions[theterm]=vals[0]

for strc_name in strc_names:
    ntopmetrics=0
    namestopmetrics=[]
    for condition in conditions.keys():
        for a,b,c in datas:
            if a==strc_name:
                if b==condition:
                    if condition in higher_better:
                        if c>conditions[condition]:
                            ntopmetrics+=1
                            namestopmetrics.append(b)
                    else:
                        if c<conditions[condition]:
                            ntopmetrics+=1
                            namestopmetrics.append(b)
    ntop[strc_name]=ntopmetrics
    nametop[strc_name]=namestopmetrics

#
sd=sorted(ntop.items(), key=lambda x: x[1], reverse=True)
print(sd[:50])
print(nametop[:50])
#
max=sd[0][1]
perc90val=max*0.9
ninety_perc=[i for i in sd if i[1]>=perc90val]
for i in sd:
    if i[1]>=22.0:
len(ninety_perc)
select_terms=['packstat','sc','buns_bb_heavy','dG_separated','rama_prepro']
npterms={}
for i in ninety_perc:
    name=i[0]
    for key in nametop.keys():
        if name==key:
            count=0
            for z in nametop[key]:
                if z in select_terms:
                    count+=1
            npterms[name]=count
nptsorted=sorted(npterms.items(), key=lambda x: x[1], reverse=True)
nptsorted
topstrc=[i for i in nptsorted if i[1]>=3.0]#change number of select metrics to filter on






####pdb to fasta file
from pyrosetta import *
init('-ignore_unrecognized_res')
def fasta(pdb):
    p=pose_from_pdb(pdb)
    sequence=str(p.sequence())[:-1]
    ofile=open(pdb[:-4]+'.fasta','w')
    ofile.write('>'+pdb+ '\n')
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

for a,b in topstrc:
    name=a[:-5]+'.pdb'
    try:
        fasta(name)
    except:
        print('failure')



'''
mv *.fasta ~/ib2fastas

import os
for a,b in topstrc:
    try:
        name=a[:-5]+'.pdb'
        os.system('cp '+name+' ~/ib2fastas')
    except:
        pass

In [36]: len(ninety_perc)
    ...:
Out[36]: 93
In [38]: len(topstrc)
Out[38]: 16

In [39]: topstrc
Out[39]:
[('UM_1_F31W35W57_1_model_87920_IBP_BSFF_v2_new_1_0062_0001', 5),
 ('UM_1_F31W35W57_1_model_87920_IBP_BSFF_v2_new_1_0591_0001', 4),
 ('UM_1_F31W35W57_1_model_87920_IBP_BSFF_v2_new_1_0891_0001', 4),
 ('UM_1_F31W35W53_1_model_403374_IBP_BSFF_v2_new_1_0002_0001', 3),
 ('UM_1_F31W35W53_1_model_403374_IBP_BSFF_v2_new_1_0240_0001', 3),
 ('UM_1_F31W35W53_1_model_403374_IBP_BSFF_v2_new_1_0389_0001', 3),
 ('UM_1_F31W35W53_1_model_403374_IBP_BSFF_v2_new_1_0489_0001', 3),
 ('UM_1_F41W81W73_1_model_136340_IBP_BSFF_v2_new_1_0376_0001', 3),
 ('UM_1_F31W35W57_1_model_87920_IBP_BSFF_v2_new_1_0232_0001', 3),
 ('UM_1_F31W35W53_1_model_403374_IBP_BSFF_v2_new_1_0024_0001', 3),
 ('UM_1_F31W35W53_1_model_403374_IBP_BSFF_v2_new_1_0083_0001', 3),
 ('UM_1_F31W35W53_1_model_403374_IBP_BSFF_v2_new_1_0138_0001', 3),
 ('UM_1_F31W35W53_1_model_403374_IBP_BSFF_v2_new_1_0216_0001', 3),
 ('UM_1_F31W35W53_1_model_403374_IBP_BSFF_v2_new_1_0486_0001', 3),
 ('UM_1_F31W35W53_1_model_403374_IBP_BSFF_v2_new_1_0421_0001', 3),
 ('UM_1_F41W81W73_1_model_136340_IBP_BSFF_v2_new_1_0572_0001', 3)]
'''











#pareto filtering
#time ipython ~/desktop/bsd_filters/pareto.py

#okay do it for top terms within each method for sure
#maybe top 2-5 terms, pps compared to ts along with nstrc
#extend to 10 terms if time maybe
#compare with using top terms from all on all strc


import numpy as np
import pandas as pd
import scipy
from scipy import spatial


rel_terms=['norm_bsa', 'design_avg_rmsd'] ##########################

higher_better=['norm_bsa','ps']
#a function to find the pareto front of a set of points
#returns indices of pareto points in input list
#LOWER VAL = BETTER VAL
def pareto(points):
    size=points.shape[0]
    ids=np.arange(size)
    pareto_front=np.ones(size,dtype=bool)
    for i in range(size):
        for j in range(size):
            if all(points[j]<=points[i])and any(points[j]<points[i]):
                pareto_front[i]=0
                break
    return ids[pareto_front]


def return_scores(sc):
    scorefile=open(sc,'r')
    lines=[line for line in scorefile.readlines()[1:]]
    scorefile.close()
    terms=[]
    for i in lines[0].split()[1:]:
        terms.append(i)
    del lines[0]
    strc_names=[]
    datas=[]
    for line in lines:
        properties=line.split()[1:]
        name=properties[-1]
        strc_names.append(name)
        for index,term in enumerate(terms[:-1]):
            if term!='description':
                if term in higher_better:
                    ogval=float(properties[index])
                    val=0-ogval
                    datas.append((name,term,val))
                else:
                    val=float(properties[index])
                    datas.append((name,term,val))
    df=pd.DataFrame()
    rel_terms.sort()
    df['terms']=pd.Series(rel_terms)
    names=[a for a,b,c in datas]
    names=list(set(names))
    for name in names:
        vals = [c for a,b,c in datas if a==name and b in rel_terms]
        df[name]=pd.Series(vals)
    scores=[]
    for name in names:
        x,y = df[name]      #####################
        scores.append((x,y))        ######################
    scores=np.array(scores)
    return datas,scores, names

def return_pareto_pps(scores,names):
    pareto_indices=pareto(scores)
    pareto_opt=[i[:-5]+'.pdb' for i in np.array(names)[pareto_indices]]
    n_strc=len(pareto_opt)
    return pareto_opt,n_strc



sf='master_sf.sc'
datas,scores,names=return_scores(sf)
par_opt,nstrc=return_pareto_pps(scores,names)

#plotting pareto front for 2 terms
import matplotlib.pyplot as plt
allvals={}
par_vals={}
for term in rel_terms:
    vals=[]
    parvals=[]
    for a,b,c in datas:
        if b==term:
            if a[:-5]+'.pdb' in par_opt:
                parvals.append(c)
            else:
                vals.append(c)
    allvals[term]=vals
    par_vals[term]=parvals
fig,ax=plt.subplots()
allx=[i for i in allvals['design_avg_rmsd']];ally=[i for i in allvals['norm_bsa']]
parx=[i for i in par_vals['design_avg_rmsd']];pary=[i for i in par_vals['norm_bsa']]
ax.scatter(allx,ally)
ax.scatter(parx,pary,c='r')
ax.set_title('Pareto Front for Fragment RMSD and Buried NPSA')
ax.set_xlabel('Average Fragment RMSD(A)')
ax.set_ylabel('Buried NPSA per Residue (A^2)')
