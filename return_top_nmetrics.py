#time ipython ~/desktop/tools/return_top_nmetrics.py scorefile.sc

# sf=open(sys.argv[1],'r') #input scorefile name
sf=open('master_sf.sc','r') ###################################
lines=[line for line in sf.readlines()]
sf.close()
terms=[]
for i in lines[1].split()[1:-1]:
    terms.append(i)
strc_names=[]
datas=[]
for line in lines[2:]:
    properties=line.split()[1:]
    name=properties[-1]
    strc_names.append(name)
    for index,term in enumerate(terms):
        val=properties[index]
        datas.append((name,term,float(val)))


#terms in the score file that are more favorable when value higher must be
#specified:
higher_better=['ps','norm_bsa','bsa','dSASA_hphobic','dSASA_int','dSASA_polar',
               'hbonds_int','nres_int','acc','sc','design_fragqual_good/total',
               'design_fragqual_coverage','wt_diff_rmsd']


#find the lowest and highest values for each term,
#store in a list (term, low, high) that I can then feed to
#'scan feature' function that will do five intervals between low and hi
def scan(datas,terms):
    to_scan=[]
    for term in terms:
        if term != 'description':
            vals=[]
            for entry in datas:
                if entry[1]==term:
                    vals.append(float(entry[2]))
            high=max(vals); low=min(vals)
            to_scan.append((term,low,high))
    return to_scan

to_scan=scan(datas,terms)

#find four equally spaced values between high and low value
def foursteps(min, max):
    d = max - min
    int_size = d/5
    p1=min+int_size;p2=min+2*int_size;p3=min+3*int_size
    p4=min+4*int_size
    threshvals=[p1,p2,p3,p4]
    return threshvals

ntop={}
nametop={}
conditions={}
for theterm,lowval,hival in to_scan:
    vals=foursteps(lowval,hival)
    if theterm in higher_better:
        conditions[theterm]=vals[len(vals)-1]
    else:
        conditions[theterm]=vals[0]

for strc_name in strc_names:
    ntopmetrics=0
    namestopmetrics=[]
    for condition in conditions.keys():
        for a,b,c in datas:
            if a==strc_name:
                if b==condition:
                    if condition in higher_better:
                        if c>conditions[condition]:
                            ntopmetrics+=1
                            namestopmetrics.append(b)
                    else:
                        if c<conditions[condition]:
                            ntopmetrics+=1
                            namestopmetrics.append(b)
    ntop[strc_name]=ntopmetrics
    nametop[strc_name]=namestopmetrics

sd=sorted(ntop.items(), key=lambda x: x[1], reverse=True)
