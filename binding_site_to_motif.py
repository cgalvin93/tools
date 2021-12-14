#USAGE:
#time ipython binding_site_to_motif.py binding_site.pdb ligand.params motifname.pdb
#-ensure ligand is last residue in pdb file
import sys
import collections
from pyrosetta import *
init()

#define the score function
sf = ScoreFunction()
from pyrosetta.rosetta.core.scoring import fa_atr,fa_rep,fa_sol,hbond_sc,fa_elec,hbond_bb_sc
sf.set_weight(fa_atr, 1)
sf.set_weight(fa_rep, .55)
sf.set_weight(fa_sol, 1)
sf.set_weight(hbond_sc, 1)
sf.set_weight(fa_elec, 1)
sf.set_weight(hbond_bb_sc,1)

##load the binding site
frag = [sys.argv[2]] #ligand params
p = Pose()
generate_nonstandard_residue_set(p,frag)
pose_from_file(p, sys.argv[1])
ligand_resnum=p.total_residue()

# #load the binding site in terminal
# frag = ['DEX.params'] #ligand params
# p = Pose()
# generate_nonstandard_residue_set(p,frag)
# pose_from_file(p, '3mnp_binding_site_dexH.pdb')
# ligand_resnum=p.total_residue()

#go through bs res and get 2 body energies w/ lig
residue_energies=collections.defaultdict(list)
for i in range(p.total_residue()-1):
    resnum=i+1
    contact_pose=Pose()  #make new pose for ligand and each res, cus we dont want rest of bs to influence energies
    ligand_pose=p.residue(p.total_residue()).clone()
    res_pose=p.residue(resnum).clone()
    contact_pose.append_residue_by_jump(res_pose, 1)
    contact_pose.append_residue_by_jump(ligand_pose, 1)
    rosetta.core.pack.optimizeH(contact_pose, sf)
    '''
    OptH -
    This function will optimize the placement of all movable
    hydrogen atoms. This includes the hydroxyl hydrogens as well as
    the HIS protonation state. If the -flip_HNQ flag is on the command
    line, then it will also consider the flip states of histadine, asparagine
    and glutamine, (nearly) as described by Word et al. 1999.
    '''
    sf(contact_pose)
    w=contact_pose.energies().energy_graph().find_energy_edge(1,2)
    w.fill_energy_map()
    hbsc=w[rosetta.core.scoring.hbond_sc];residue_energies[resnum].append(hbsc)
    hbbbsc=w[rosetta.core.scoring.hbond_bb_sc];residue_energies[resnum].append(hbbbsc)
    faatr=w[rosetta.core.scoring.fa_atr];residue_energies[resnum].append(faatr)
    farep=w[rosetta.core.scoring.fa_rep];residue_energies[resnum].append(farep)
    fasol=w[rosetta.core.scoring.fa_sol];residue_energies[resnum].append(fasol)
    faelec=w[rosetta.core.scoring.fa_elec];residue_energies[resnum].append(faelec)
    total_score=faelec+fasol+farep+faatr+hbbbsc+hbsc;residue_energies[resnum].append(total_score)


'''
okay, now to identify 3 best residues to define motif
i want to prioritize directional interaxns, hbonds and estat
beyond that total score is as good a discriminator as any
SO:
filter to res with favorable total score
of these, ID 3 res with best sum hbsc, hbbbsc, fa elec
'''
total_E_d={}
favorable_residues=[]
for key in residue_energies.keys():
    total_E=residue_energies[key][-1]
    if total_E<0.:
        favorable_residues.append(key)
    total_E_d[key]=total_E
best_total_E=sorted(total_E_d.items(), key=lambda x: x[1])

directional_E={}
for res in favorable_residues:
    dir_E=residue_energies[res][0]+residue_energies[res][1]+residue_energies[res][5]
    directional_E[res]=dir_E
sd=sorted(directional_E.items(), key=lambda x: x[1])

#For top 3 directional, if they are above median total score, add to motif
motif_residues=[]
median_E=len(best_total_E)/2
median_E=int(median_E)
for a,b in sd[:3]:
    for index,i in enumerate(best_total_E):
        if a==i[0]:
            if index<median_E:
                if a not in motif_residues:
                    motif_residues.append(a)
#otherwise, just add the top scoring residues until we have 3
for i in range(len(best_total_E)):
    while len(motif_residues)<3:
        if best_total_E[i][0] not in motif_residues:
            motif_residues.append(best_total_E[i][0])
    else:
        break

#finally, output a pdb file of the motif residues with the ligand
motif_pose=Pose()
for i in motif_residues:
    res_pose=p.residue(i).clone()
    motif_pose.append_residue_by_jump(res_pose, 1)
ligand_pose=p.residue(p.total_residue()).clone()
motif_pose.append_residue_by_jump(ligand_pose, 1)
motif_pose.dump_pdb(sys.argv[3])
