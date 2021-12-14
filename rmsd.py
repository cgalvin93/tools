####pymol script to calculate the RMSD between a thermophile with all of its mesophilic homologues

import __main__
__main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
from time import sleep
import pymol
pymol.finish_launching()
from pymol import cmd
from glob import glob

thermofile = 'clean1pcz_A_100.pdb'
cmd.load(thermofile,"thermo")

mesofiles = 'mesos/*.pdb'
for file in glob(mesofiles):
	print file
	cmd.load(file,"meso")
	rms = cmd.align("thermo","meso")[0]
	print "%s rmsd: %s" % (file,rms)
	cmd.delete("meso")
pymol.cmd.quit()
