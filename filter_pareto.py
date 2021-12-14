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


#########################
#########################
#########################
#########################
#########################
sf='master_sf.sc'
datas,scores,names=return_scores(sf)
par_opt,nstrc=return_pareto_pps(scores,names)
#########################
#########################
#########################
#########################
#########################

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
