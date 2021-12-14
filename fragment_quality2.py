# time ipython ~/desktop/tools/fragment_quality2.py design.pdb
# time ipython ~/desktop/tools/fragment_quality2.py 1F4P_0366.pdb
# these files must be in working directory:
# design.pdb design9mers
# the 9mers file should have the following format:
# 'PDBID_designnumber_frags.200.9mers'
# so if the design file is 2xbn_0346.pdb, the 9mers should be:
# 2xbn_0346_frags.200.9mers
# if the 9mer files you generate have a different suffix, you can change it
# in the 'design_path' variable around line 27
# designed residue positions must be manually input on line 34
# or just put designed_positions=[1] if you don't want the design
# position specific metrics

import sys
import numpy as np
from pyrosetta import *
init()
import os
import pandas as pd

ptn_pdb_id=sys.argv[1][:4]
cwd=os.getcwd()


#load ninemers
fragset9 = pyrosetta.rosetta.core.fragment.ConstantLengthFragSet(9)
design_path=cwd+'/'+sys.argv[1][:-4].upper()+'_frags.200.9mers'
fragset9.read_fragment_file(design_path)

#load design pdb
designs=[sys.argv[1:]]

#manually input designed positions
designed_positions=[1]

#open text file to write results
outfilename=cwd+'/'+ptn_pdb_id+'_frag_filters.txt'
of=open(outfilename,'a')
#######
resultslist=[]
resultsdict={}
resultsdict['SCORE:']='SCORE:'

#do stuff
for design in designs:
    of.write('\n'+'\n'+'DESIGN'+'\n')
    of.write('\n'+str(design[:-4].upper())+'\n')
    p=pose_from_pdb(design)
    #run fragqual filter on the whole ptn for n good frags + fraction of ptn covered by good frags
    calc9=pyrosetta.rosetta.protocols.pose_metric_calculators.FragQualCalculator(fragset9)
    calc9.begin(1)
    calc9.end(p.total_residue()-1)
    calc9.rmsd_cutoff(1.0)
    fragscore9=calc9.get('num_goodfrag',p)
    total_frags=fragset9.size() #total number of fragments
    of.write('The number of good fragments is: '+str(fragscore9)+' out of: '+str(total_frags)+' total fragments')
    resultsdict['design_fragqual_good/total']=float(fragscore9)/total_frags
    ###############coverage
    calc9r=pyrosetta.rosetta.protocols.pose_metric_calculators.FragQualCalculator(fragset9)
    calc9r.begin(1)
    calc9r.end(p.total_residue()-1)
    calc9r.ratio_cutoff(0.3)
    fragscore9r=calc9r.get('coverage',p)
    of.write('\nThe fraction of positions covered by good fragments is : '+str(fragscore9r)+'\n')
    resultsdict['design_fragqual_coverage']=float(fragscore9r)
    ############### design position specific metrics for 9mers
    of.write('\nDesign position specific metrics: '+'\n')
    fr=pyrosetta.rosetta.core.fragment.FragmentRmsd(fragset9)
    designed_position_fragments=[]
    for i in designed_positions:
        for k in range(i-8,i+1):
            if k>0 and k<=p.total_residue()-9:
                designed_position_fragments.append(k)
    designed_position_fragments=list(set(designed_position_fragments))
    #first calculate rmsd for all unique fragments containg des pos, to avoid redundant calculation
    all_rmsd_data={}
    print('calculating fragment rmsd values for designed positions')
    for fragment_number in designed_position_fragments:
        design_fraglib_rmsd=[]
        for x in range(1,201):
            curr_rmsd=fr.rmsd(fragment_number,x,p)
            design_fraglib_rmsd.append(curr_rmsd)
            print(str(curr_rmsd))
        all_rmsd_data[fragment_number]=design_fraglib_rmsd
    #now go thru each design position and get the rmsd values for frags its a part of
    for i in designed_positions:
        frags=[]
        pos_all_rmsd=[]
        for k in range(i-8,i+1):
            if k>0 and k<=p.total_residue()-9:
                frags.append(k)
        for k in frags:
            for x in all_rmsd_data[k]:
                pos_all_rmsd.append(x)
        position_avg_rmsd=np.mean(pos_all_rmsd)
        position_max_rmsd=np.max(pos_all_rmsd)
        position_min_rmsd=np.min(pos_all_rmsd)
        position_med_rmsd=np.median(pos_all_rmsd)
        of.write('\nFor Design Position: '+str(i))
        of.write('\nThe Average RMSD is: '+str(position_avg_rmsd))
        of.write('\nThe Median RMSD is: '+str(position_med_rmsd))
        of.write('\nThe Min RMSD is: '+str(position_min_rmsd))
        of.write('\nThe Max RMSD is: '+str(position_max_rmsd)+'\n')
    #now get the average rmsd for all fragments that contain a design position
    all_design_position_means=[]
    for key in all_rmsd_data.keys():
        for value in all_rmsd_data[key]:
            all_design_position_means.append(np.mean(value))
    of.write('\nThe Average RMSD for all fragments containing a designed position is: '+str(np.mean(all_design_position_means)))
    resultsdict['design_avg_rmsd_des_pos']=float(np.mean(all_design_position_means))

of.close()


resultsdict['description']=sys.argv[1][:-4]+'_0001'#always this way in scorefiles
resultslist.append(resultsdict)
df=pd.DataFrame(resultslist)
dfname=sys.argv[1][:4]+'_fragment_filters.sc'
dffnew=open(dfname,'a')
dffnew.write('SEQUENCE:\n')
df.to_csv(dffnew, index=False, sep='\t', mode='a')
dffnew.close()
