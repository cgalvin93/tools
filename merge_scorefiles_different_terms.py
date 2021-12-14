#time ipython ~/desktop/tools/merge_scorefiles_different_terms.py sc1 sc2 outfilename
#time ipython ~/desktop/tools/merge_scorefiles_different_terms.py 1f4p_score.sc 1f4p_rb_scores.sc 1f4pall.sc
#time ipython ~/desktop/tools/merge_scorefiles_different_terms.py normalized_scores.sc fixedfragqual.sc master_sf.sc

import sys

sc1=open(sys.argv[1],'r')
sc2=open(sys.argv[2],'r')

sc1lines=[line for line in sc1.readlines()[1:]] #sequence: printed at top
sc2lines=[line for line in sc2.readlines()[1:]]
sc1.close()
sc2.close()
#acquire the names of the metrics calculated in the scorefile from first line
sc1terms=[]
for i in sc1lines[0].split()[1:]: #score: always first thing on line
    sc1terms.append(i)
sc2terms=[]
for i in sc2lines[0].split()[1:]:
    sc2terms.append(i)
ogsc2terms=[i for i in sc2terms]
#make master list of terms,description last, no duplicate terms
#dict description:[term vals in order of terms list]
#line 1 is terms.join
#for every dict entry, line = value.join + key at end

#make a master list of terms, excluding those in second which
#are already in first scorefile
all_terms=[]
all_terms.append('SCORE:')
for a in sc1terms:
    if a!='description':
        all_terms.append(a)
for i in range(10): #idk why but i have to do this multiple times to erase all instances of duplicates
    for k in sc2terms:
        if k in sc1terms:
            sc2terms.remove(k)
for i in sc2terms:
    if i!='description':
        all_terms.append(i)
all_terms.append('description')# making sure to put this one last as is standard in .sc from rosetta

#parse the score file to create a list where each element is (structure name,property name, property value)
#for all structures and properties in the score file
def return_data(lines,terms):
    strc_names=[]
    datas=[]
    for line in lines[1:]:
        try:
            properties=line.split()[1:]
            name=properties[-1]
            if name not in strc_names:
                strc_names.append(name)
                for index,term in enumerate(terms[:-1]):
                    val=properties[index]
                    datas.append((name,term,val))
        except:
            continue
    return datas, list(set(strc_names))

sc1data,sc1names=return_data(sc1lines,sc1terms)
sc2data,sc2names=return_data(sc2lines,ogsc2terms)
sc1data=list(set(sc1data))
sc2data=list(set(sc2data))

allnames=[i for i in sc1names if i in sc2names]
if len(allnames)<len(sc1names) or len(allnames)<len(sc2names):
    print('\n\n\n\nALL STRUCTURES NOT PRESENT IN BOTH SCOREFILES')
    print(len(sc2names))
dict={}
for strc in allnames:#assuming both files have same structures
    l=[]
    l.append('SCORE:')
    for term in sc1terms:
        for a,b,c in sc1data:
            if strc==a:
                if term==b:
                    l.append(c)
    for term2 in sc2terms:
        for d,e,f in sc2data:
            if strc==d:
                if term2==e:
                    l.append(f)
    l.append(strc)#add description as last
    dict[strc]=l


ofile=open(sys.argv[3],'w')
ofile.write('SEQUENCE:\n')
term_line='  '.join(all_terms)
ofile.write(term_line+'\n')
for key in dict.keys():
    newline='  '.join(dict[key])
    ofile.write(newline+'\n')
ofile.close()
