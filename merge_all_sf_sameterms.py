#!/usr/bin/env python

##time ipython ~/desktop/tools/merge_all_sf_sameterms.py outfilename
##time ipython ~/desktop/tools/merge_all_sf_sameterms.py simon_ib2_designs_iam.sc


import sys
import os

scfiles=[i for i in os.listdir() if i[-3:]=='.sc']

sc1=open(scfiles[0],'r')
sc1lines=[line for line in sc1.readlines()[1:]] #sequence: printed at top
sc1.close()

#acquire the names of the metrics calculated in the scorefile from first line
sc1terms=[]
for i in sc1lines[0].split()[0:]: #score: always first thing on line
    sc1terms.append(i)
ofile=open(sys.argv[1],'w')
ofile.write('SEQUENCE:\n')
term_line='  '.join(sc1terms)
ofile.write(term_line+'\n')


for i in sc1lines[1:]:
    props=[]
    for x in i.split():
        props.append(x)
    propline='  '.join(props)
    ofile.write(propline+'\n')


for i in scfiles[1:]:
    sc=open(i,'r')
    lines=[line for line in sc.readlines()[1:]] #sequence: printed at top
    sc.close()
    terms=[]
    for x in lines[0].split()[0:]: #score: always first thing on line
        terms.append(x)
    matchedterms=[]
    for k,e in enumerate(terms):
        if e==sc1terms[k]:
             matchedterms.append(e)
    if len(matchedterms)==len(terms):
        for z in lines[1:]:
            props=[]
            for y in z.split():
                props.append(y)
            propline='  '.join(props)
            ofile.write(propline+'\n')
    else:
        print('METRICS DO NOT MATCH BETWEEN SCORE FILES')

ofile.close()
