#time ipython ~/desktop/tools/fragment_quality.py des9mers des.pdb params wt9mers wt.pdb pdbID
#script assumes 200 fragments in frag libraries
#time ipython ~/desktop/tools/fragment_quality.py frags2009mers 1F4P_0395.pdb FMN.params 1f4pwt_9mers 1f4p_wt.pdb test

#time ipython ~/desktop/tools/fragment_quality.py frags2009mers 1F4P_0395.pdb FMN.params 1f4pwt_9mers 1f4p_wt.pdb 1f4p
#time ipython ~/desktop/tools/fragment_quality.py 1F4P_0366_frags.200.9mers 1F4P_0366.pdb FMN.params 1f4pwt_9mers 1f4p_wt.pdb 1f4p

import sys
# import os
import numpy as np
import pandas as pd
from pyrosetta import *
init()
#from  pyrosetta.rosetta.core.scoring import *

#manually input designed positions
#store pdb numbering of binding site residues for query protein
ptn_pdb_id=sys.argv[6]
if ptn_pdb_id == '2xbn':
    designed_positions=[73, 133, 134, 135, 138, 159,161, 202, 231, 233,
                        234, 262,264, 265]
elif ptn_pdb_id == '1f4p':
    designed_positions=[10, 11, 12, 13, 14, 15, 58, 59,
                        60, 61, 62, 68, 93, 94, 95, 98,
                        100, 101, 102]
elif ptn_pdb_id == '1zk4':
    designed_positions=[13, 14, 15, 16, 17, 18, 19, 36,
                        37, 38, 62, 63, 89, 90, 91, 92,
                        112, 140, 142, 155, 159]
elif ptn_pdb_id == '3dk9':
    designed_positions=[26, 27, 28, 29, 30, 31, 49, 50,
                        51, 52, 56, 57, 62,
                        66, 129, 130, 155, 156, 157,
                        177, 181, 197, 198, 201, 202,
                        291, 294, 298, 330, 331]
elif ptn_pdb_id == '3dlc':
    designed_positions=[48, 49, 50,
                        51, 52, 53, 55, 72, 73, 74, 77,
                        100, 101, 102, 117, 118, 119,
                        122, 123]
elif ptn_pdb_id == '3r2q':
    designed_positions=[9, 10, 11, 33, 34, 48, 49, 50,
                        62, 63, 64]
elif ptn_pdb_id == '3s6f':
    designed_positions=[78, 79, 80, 85, 86, 87, 88, 90, 91, 114, 115, 118, 119, 121]
elif ptn_pdb_id == '1zk4':
    designed_positions=[13, 14, 15, 16, 17, 18, 19, 36, 37, 38, 62, 63, 89, 90,
                        91, 92, 112, 140, 142, 155, 159, 187, 189, 190, 192, 194,
                        195, 205]
elif ptn_pdb_id == '2ij2':
    designed_positions=[69,75,86,87,96,100,107,153,260,261,264,265,268,269,272,
                        322,327,328,331,357,392,393,394,398,399,400,401,402,405,406]
elif ptn_pdb_id == 'test':
    designed_positions=[1]
else:
    print('the binding site residues of the query protein cannot be identified')


#load ninemers
fragset9 = pyrosetta.rosetta.core.fragment.ConstantLengthFragSet(9)
fragset9.read_fragment_file(sys.argv[1])
#load design pdb and params file
designs=[sys.argv[2]]
paramsfile=[sys.argv[3]]
#load wildtype scaffold pdb and 9mers
wt_ninemers=pyrosetta.rosetta.core.fragment.ConstantLengthFragSet(9)
wt_ninemers.read_fragment_file(sys.argv[4])
wt_scaffolds=[sys.argv[5]]  #this doesnt need to be a list im just lazy

#open text file to write results
outfilename=sys.argv[2][:-4]+'_frag_filters.txt'
of=open(outfilename,'a')



#to store stuff for pandas df
resultslist=[]
resultsdict={}
resultsdict['SCORE:']='SCORE:'
#get the rmsd of design position containing fragments in wildtype structure so
#that they may be compared with designs
#mostly just copying code from below (i wrote the block for design evaluation first)
#and changing names to wt 9mers/pdb
wt_values_for_comparison=[] #(n_good_frags/total_frags,coverage,(pos,avg_rmsd_pos),avg_rmsd_all_pos)
for wt_scaffold in wt_scaffolds:
    of.write('\n'+'WILD TYPE SCAFFOLD'+'\n')
    of.write('\n'+str(wt_scaffold[:-4].upper())+'\n')
    p=Pose()
    generate_nonstandard_residue_set(p,paramsfile)
    pose_from_file(p, wt_scaffold)
    #run fragqual filter on the whole ptn for n good frags + fraction of ptn covered by good frags
    calc9=pyrosetta.rosetta.protocols.pose_metric_calculators.FragQualCalculator(wt_ninemers)
    calc9.begin(1)
    calc9.end(p.total_residue()-1)
    calc9.rmsd_cutoff(1.0)
    fragscore9=calc9.get('num_goodfrag',p)
    total_frags=wt_ninemers.size() #total number of fragments
    of.write('The number of good fragments is: '+str(fragscore9)+' out of: '+str(total_frags)+' total fragments')
    wt_values_for_comparison.append(float(fragscore9)/total_frags)
    ###############coverage
    calc9r=pyrosetta.rosetta.protocols.pose_metric_calculators.FragQualCalculator(wt_ninemers)
    calc9r.begin(1)
    calc9r.end(p.total_residue()-1)
    calc9r.ratio_cutoff(0.3)
    fragscore9r=calc9r.get('coverage',p)
    of.write('\nThe fraction of positions covered by good fragments is : '+str(fragscore9r)+'\n')
    wt_values_for_comparison.append(float(fragscore9r))
    ############### design position specific metrics for 9mers
    of.write('\nDesign position specific metrics: '+'\n')
    fr=pyrosetta.rosetta.core.fragment.FragmentRmsd(wt_ninemers)
    designed_position_fragments=[]
    for i in designed_positions:
        for k in range(i-8,i+1):
            if k>0 and k<=p.total_residue()-9:
                designed_position_fragments.append(k)
    designed_position_fragments=list(set(designed_position_fragments))
    #first calculate rmsd for all unique fragments containg des pos, to avoid redundant calculation
    all_rmsd_data={}
    print('\n'+'\n'+'\n'+'\n'+'calculating fragment rmsd values for designed positions for wildtype'.upper())
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
        wt_values_for_comparison.append((i,position_avg_rmsd))
    #now get the average rmsd for all fragments that contain a design position
    all_design_position_means=[]
    for key in all_rmsd_data.keys():
        for value in all_rmsd_data[key]:
            all_design_position_means.append(np.mean(value))
    of.write('\nThe Average RMSD for all fragments containing a designed position is: '+str(np.mean(all_design_position_means)))
    wt_values_for_comparison.append(float(np.mean(all_design_position_means)))


#evaluate designs and compare with wt
for design in designs:
    des_values_for_comparison=[]
    of.write('\n'+'\n'+'DESIGNS'+'\n')
    of.write('\n'+str(design[:-4].upper())+'\n')
    p=Pose()
    generate_nonstandard_residue_set(p,paramsfile)
    pose_from_file(p, design)
    #run fragqual filter on the whole ptn for n good frags + fraction of ptn covered by good frags
    calc9=pyrosetta.rosetta.protocols.pose_metric_calculators.FragQualCalculator(fragset9)
    calc9.begin(1)
    calc9.end(p.total_residue()-1)
    calc9.rmsd_cutoff(1.0)
    fragscore9=calc9.get('num_goodfrag',p)
    total_frags=fragset9.size() #total number of fragments
    of.write('The number of good fragments is: '+str(fragscore9)+' out of: '+str(total_frags)+' total fragments')
    des_values_for_comparison.append(float(fragscore9)/total_frags)
    resultsdict['design_fragqual_good/total']=float(fragscore9)/total_frags
    ###############coverage
    calc9r=pyrosetta.rosetta.protocols.pose_metric_calculators.FragQualCalculator(fragset9)
    calc9r.begin(1)
    calc9r.end(p.total_residue()-1)
    calc9r.ratio_cutoff(0.3)
    fragscore9r=calc9r.get('coverage',p)
    of.write('\nThe fraction of positions covered by good fragments is : '+str(fragscore9r)+'\n')
    des_values_for_comparison.append(float(fragscore9r))
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
    print('\n'+'\n'+'\n'+'\n'+'calculating fragment rmsd values for designed positions'.upper())
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
        des_values_for_comparison.append((i,position_avg_rmsd))
    #now get the average rmsd for all fragments that contain a design position
    all_design_position_means=[]
    for key in all_rmsd_data.keys():
        for value in all_rmsd_data[key]:
            all_design_position_means.append(np.mean(value))
    of.write('\nThe Average RMSD for all fragments containing a designed position is: '+str(np.mean(all_design_position_means)))
    des_values_for_comparison.append(float(np.mean(all_design_position_means)))
    resultsdict['design_avg_rmsd_des_pos']=float(np.mean(all_design_position_means))
    #now comparing the design with the wildtype scaffold
    of.write('\n'+'\n'+'COMPARISON OF DESIGN '+str(design[:-4].upper())+' WITH WILDTYPE (wt_value - design_value)'+'\n')
    of.write('\n'+'n_good_frags/total_frags, coverage, (pos,pos,avg_rmsd_pos), avg_rmsd_all_pos')
    comparisonlist=[]
    for i,e in enumerate(wt_values_for_comparison):
        if type(e)==float:
            diff=e-des_values_for_comparison[i]
            of.write('\n'+str(diff))
            comparisonlist.append(diff)
        elif type(e)==tuple:
            diff=e[1]-des_values_for_comparison[i][1]
            of.write('\n'+str(e[0])+','+str(des_values_for_comparison[i][0])+': '+str(diff))
    resultsdict['wt_diff_fragqual']=comparisonlist[0]
    resultsdict['wt_diff_coverage']=comparisonlist[1]
    resultsdict['wt_diff_rmsd']=comparisonlist[2]

of.close()

resultsdict['description']=sys.argv[2][:-4]+'_0001'#always this way in scorefiles
resultslist.append(resultsdict)
df=pd.DataFrame(resultslist)
dfname=ptn_pdb_id+'fragment_filters.sc'
if os.path.exists(dfname):
    df.to_csv(dfname, header=None,index=False, sep='\t', mode='a')
else:
    dffnew=open(dfname,'w')
    dffnew.write('SEQUENCE:\n')
    dffnew.close()
    df.to_csv(dfname, index=False, sep='\t', mode='a')



'''
DEVELOPMENT


calc9r.ratio_cutoff(0.3)
#my guess at what this is:
#???in order for a position to be considered covered by 'good' fragments, the fraction of
#fragments it is contained in (typically 200 for each position) that have less than 1 angstrom rmsd
#with fragments in fragment library must be greater than or equal to THIS VALUE.
#so, when this value is 0, a position must have 0 good fragments out of the N fragments it
#is a part of to be considered covered, thus all positions will pass, and it will return 1
#if it is 1, 100% of fragments containing a position must have <1A rmsd with the corresponding fragments in the design
#   so (likely) no positions pass and the filter returns 0
#The higher you make it, the more stringent the criteria for a position 'passing'
#   as being considered 'covered' by good fragments
#0.3 is the default value, so I will just stick with that


!!!!!!!!!!
rmsd(self: pyrosetta.rosetta.core.fragment.FragmentRmsd, position: int, k: int, reference: pyrosetta.rosetta.core.pose.Pose) → float¶
Returns the RMSD of the kth fragment at the specified position
in the fragment library and pose.
!!!!!!!!!!
!!!!!!!!!SAYS ON https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/fragment-file
!!!!!!!!!THAT THE FRAGMENT NUMBER CORRESPONDS WITH THE POSE NUMBERING
!!!!!!!!!BUT THEN I DONT UNDERSTAND WHY THERE ARE LESS FRAGMENTS
!!!!!!!!!!(139 FRAGS VS 147 RESIDUES FOR 1F4P_0395)
''Where "position" is the pose number of the starting point of the fragment,
and "neighbors" is the number of fragments (ignored on reading).''
    OKAY I THINK WHATS GOING ON IS THAT THE FRAGMENT NUMBERS ARE FOR EACH OF THE N-8
    UNIQUE FRAGMENT SEQUENCES OF LENGTH 9 THAT YOU CAN MAKE OUT OF A PTN SEQUENCE OF LENGTH N
        SO RESIDUE 1 IS ONLY INCLUDED IN FRAGMENT 1 (RESIDUES 123456789)
        RESIDUE 2 IS INCLUDED IN RESIDUE 1 AND 2
        RES 3 = 1,2,3
        ...
        RES9 =1,2,3,4,5,6,7,8,9
        RES10=2,3,4,5,6,7,8,9,10
        ...
        RES n = {n-8 > 0: n<=N-8}
        ....
    AND THERE CAN ONLY BE N-8 FRAGMENTS OF LENGTH 9 BECAUSE AT RESIDUE N-7 THE REMAINING SEGMENT
    IS ONLY 8 RESIDUES LONG, SHORTER FOR REMAINING N-6 RESIDUES
    SO LAST RESIDUE IS ALSO ONLY INCLUDED IN 1, SECOND LAST IN 2, ETC....

fragset9 = pyrosetta.rosetta.core.fragment.ConstantLengthFragSet(9)
fragset9.read_fragment_file('1f4p395_9mers')

fr=pyrosetta.rosetta.core.fragment.FragmentRmsd(fragset9)
print(fr.fragment(15,150))
print(fr.rmsd(15,150,p))

calc9=pyrosetta.rosetta.protocols.pose_metric_calculators.FragQualCalculator(fragset9)
# calc9.begin(1)
# calc9.end(p.total_residue())
calc9.begin(10)
calc9.end(15)
calc9.rmsd_cutoff(1.0)
fragscore9=calc9.get('num_goodfrag',p)
print(fragscore9)

calc9r=pyrosetta.rosetta.protocols.pose_metric_calculators.FragQualCalculator(fragset9)
calc9r.begin(1)
calc9r.end(p.total_residue())
# calc9r.begin(10)
# calc9r.end(15)
calc9r.ratio_cutoff(0.3)
fragscore9r=calc9r.get('coverage',p)
print(fragscore9r)


'''

'''
FOR DOING IT WITHOUT PYROSETTA --- UNFINISHED
I THINK ILL NEED THIS FOR TINBERG FILTER
BECAUSE IT DOESNT SEEM THERE IS ANY WAY TO ACCESS RMSD OF FRAGMENTS THEMSELVES IN
PYROSETTA
        JUST KIDDING I DONT NEED THE RMSD AMONGST FRAGMENTS THEMSELVES
        BUT IM GONNA LEAVE THIS HERE IN CASE I EVER DO NEED IT FOR SOME REASON
# threemers=open(sys.argv[1],'r')
# ninemers=open(sys.argv[2],'r')
threemers=open('aat000_03_05.200_v1_3','r') #
ninemers=open('aat000_09_05.200_v1_3','r') #
threemerlines=[line for line in threemers.readlines()[1:] if line!='\n']
ninemerlines=[line for line in ninemers.readlines()[1:] if line!='\n']
threemers.close(); ninemers.close()

#load relevant fragment data (residueID,phi,psi,omega)
threemer_data=[]
ninemer_data=[]
#threemers
for i in range(0,len(threemerlines),3):
    current_fragment=(threemerlines[i:i+3])
    fragment_data=[]
    for k in current_fragment:
        fragment_data.append((k[14],k[17:26],k[27:35],k[36:44]))#x,y,z=k[46:55],k[55:64],k[65:73]
    threemer_data.append(fragment_data)
#ninemers
for i in range(0,len(ninemerlines),9):
    current_fragment=(ninemerlines[i:i+9])
    fragment_data=[]
    for k in current_fragment:
        fragment_data.append((k[14],k[17:26],k[27:35],k[36:44]))#x,y,z=k[46:55],k[55:64],k[65:73]
    ninemer_data.append(fragment_data)

def load_fragment(frag_res_list):
    sequence=''
    for i in frag_res_list:
        sequence+=i[0]
    print(sequence)
    fragment_pose=pose_from_sequence(sequence)
    for i,e in enumerate(frag_res_list):
        fragment_pose.set_phi(i+1,float(e[1]))
        fragment_pose.set_psi(i+1,float(e[2]))
        fragment_pose.set_omega(i+1,float(e[3]))
    return fragment_pose

'''


'''
ROCKLIN

Fragment quality analysis:
Fragments were chosen for each designed protein using the standard Rosetta fragment generation protocol (17), which uses the designed sequence and PSIPRED-predicted secondary structure (75) as input. These metrics quantify the geometric agreement between the selected 9-mer fragments from natural proteins and the corresponding 9-mer segments of the designs (200 9-mer fragments are chosen per designed segment).
avg_all_frags: the average RMSD of all selected fragments to their corresponding segments of the designs, in Å. (200 × (n - 8) fragments in total for protein length n)
avg_best_frags: the average RMSD of the lowest-RMSD fragment for each designed segment, in Å. (n - 8 fragments in total)
sum_best_frags: the sum of the RMSDs of the lowest-RMSD fragment for each designed segment. (n - 8 fragments in total)
worstfrag: among the set of fragments that are the lowest-RMSD fragments for their positions, the highest RMSD found
worst6frags: among the set of fragments that are the lowest-RMSD fragments for their positions, the sum of the RMSDs of the six highest RMSD fragments
'''


'''
TINBERG

"we developed a metric to estimate the impact of design on local backbone
structure and used this metric to discard designs that were predicted to
lead to backbone structure changes. Using the structure prediction modules
of Rosetta19, we generated a set of 9-mer fragment structures for each
designed and wild type scaffold sequence and compared the average RMSD
of these fragments to those of the scaffold backbone structures. If the
average RMSD of conformations predicted in these fragments (200 9-mers)
near any designed position was greater (> 0.8 Å) for the designed sequence
than the wild type scaffold sequence, we flagged that region of the designed
protein as unlikely to adopt the local backbone conformation of the scaffold
protein and rejected that designed protein
'''
