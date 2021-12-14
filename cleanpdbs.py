#ipython ~/desktop/tools/cleanpdbs.py *.pdb
#CLEAN PDB FILES TO WORK BETTER WITH ROSETTA
#erase non atom/ter/hetnam/hetatm lines
#change hetatm to atom records
#erase waters
#change num + chain of a specified ligand to 1 X
#writes ligands to their own pdb file


import sys
for filename in sys.argv[1:]:
    f=open(filename,'r')
    lines=[line for line in f.readlines() if line[0:4]=='ATOM' or line[0:3]=='TER' or line[0:6]=='HETNAM' or line[0:6]=='HETATM']
    ligands=[]
    for line in lines:
        if line[0:6]=='HETATM':
            ligname=line[17:20]
            if ligname!='HOH':
                if ligname not in ligands:
                    ligands.append(ligname)
    for ligname in ligands:
        ofile=open(ligname+'.pdb','w')
        for line in lines:
            if line[17:20]==ligname:
                ofile.write(line)
        ofile.close()
    ofile=open('clean_'+filename,'w')
    for line in lines:
        if line[0:4]=='ATOM':
            ofile.write(line)
        if line[0:6]=='HETATM':
            if line[17:20]=='HOH':
                pass
            # else:
            #     ofile.write(line)
            elif line[17:20]=='DEX':  #gotta put name of ligand you wanna change to 1x for rosetta
                # newline='ATOM  '+line[6:20]+' X'+' 1     '+line[29:]
                newline='ATOM  '+line[6:17]+'DEX'+' X'+' 1     '+line[29:] #TO ALSOP CHANGE LIG NAME
                ofile.write(newline)
            else:
                newline='ATOM  '+line[6:]
                ofile.write(newline)

        else:
            ofile.write(line)
    f.close()
    ofile.close()
