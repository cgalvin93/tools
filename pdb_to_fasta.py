#time ipython ~/desktop/tools/pdb_to_fasta.py outfile.fasta *.pdb
#python ~/desktop/tools/pdb_to_fasta.py outfile.fasta onefile.pdb
#will output a fasta file with all sequences for pdb files in working directory

import Bio
from Bio import SeqIO as seqio
import sys


seqs=[]
for file in sys.argv[2:]:
    for record in seqio.parse(file,'pdb-atom'):
        s=record.seq
        seqs.append((file,s))

ofile=open(sys.argv[1],'w')
for name, sequence in seqs:
    ofile.write('>'+name + '\n')
    length=len(sequence)
    remainder=length%60; n_seg=length/60
    indices=[]
    for i in range(0,60,length-remainder):
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
