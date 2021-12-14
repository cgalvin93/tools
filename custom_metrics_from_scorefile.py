#parse the scorefile to create a list of tuples where each element is (structure_name, scoreterm_name,scoreterm_value)

sf=open('stability_scores.sc','r') #input scorefile name
lines=[line for line in sf.readlines()]
sf.close()
terms=[]
for i in lines[1].split()[1:]:
    terms.append(i)
strc_names=[]
datas=[]
for line in lines[2:]:
    properties=line.split()[1:]
    name=properties[-1]
    strc_names.append(name)
    for index,term in enumerate(terms):
        if term!='description':
            val=properties[index]
            datas.append((name,term,val))

#get the residue counts handy for each structure
rescount_dict={}
for strc_name in strc_names:
    for name,term,val in datas:
        if name==strc_name:
            if term=='rescount':
                rescount_dict[strc_name]=val
#get the total sasa handy for each structure
tsasa_dict={}
for strc_name in strc_names:
    for name,term,val in datas:
        if name==strc_name:
            if term=='total_sasa':
                tsasa_dict[strc_name]=val

#total_score/rescount
#bsa/rescount
#exphyd/rescount
#hphobe_sasa/total_sasa
newterms=[]
newterms.append('SCORE:')
for i in terms[:-1]:
    newterms.append(i)
newterms.append('norm_total_score')
newterms.append('norm_bsa')
newterms.append('norm_exphyd')
newterms.append('ratio_sasa')
newterms.append('description')
newdatas=[i for i in datas]
for strc_name in rescount_dict.keys():
    for name,term,val in datas:
        if name==strc_name:
            if term=='total_score':
                norm_score=float(val)/float(rescount_dict[strc_name])
                newdatas.append((strc_name,'norm_total_score',norm_score))
            elif term=='bsa':
                norm_score=float(val)/float(rescount_dict[strc_name])
                newdatas.append((strc_name,'norm_bsa',norm_score))
            elif term=='exphyd':
                norm_score=float(val)/float(rescount_dict[strc_name])
                newdatas.append((strc_name,'norm_exphyd',norm_score))
for strc_name in tsasa_dict.keys():
    for name,term,val in datas:
        if name==strc_name:
            if term=='hphobe_sasa':
                norm_score=float(val)/float(tsasa_dict[strc_name])
                newdatas.append((strc_name,'ratio_sasa',norm_score))

#write to new scorefile
dict={}
for strc in strc_names:
    l=[]
    l.append('SCORE:')
    for a,b,c in newdatas:
        if strc==a:
            l.append(str(c))
    l.append(strc)#add description as last
    dict[strc]=l
#
ofile=open('normalized_scores.sc','w')
ofile.write('SEQUENCE:\n')
term_line='  '.join(newterms)
ofile.write(term_line+'\n')
for key in dict.keys():
    newline='  '.join(dict[key])
    ofile.write(newline+'\n')
ofile.close()


'''
#fixing names on frag sc so they match norm scores sc
sf=open('fragment_filters.sc','r')
lines=[line for line in sf.readlines()]
sf.close()
newlines=[]
newlines.append(lines[0])
newlines.append(lines[1])
for line in lines[2:]:
    properties=line.split()
    name=properties[-1][46:]+'\n'
    newprops=[i for i in properties[:-1]]
    newprops.append(name)
    newline='   '.join(newprops)
    newlines.append(newline)

ofile=open('fixedfragqual.sc','w')
for line in newlines:
    ofile.write(line)
ofile.close()

'''
