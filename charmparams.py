f=open('samH_bcc.params','r')
lines=[line for line in f.readlines()]
atoms=[]
for line in lines:
    if line[0:4]=='ATOM':
        atoms.append(line)

def rosettatocharmm(r):
    if r=='Nam':
        newr='NH3  '
    elif r=='CSp':
        newr='CT1  '
    elif r=='CDp':
        newr='C    '
    elif r=='Oat':
        newr='OC   '
    elif r=='CS2':
        newr='CT2  '
    elif r=='Oad':
        newr='O    '
    elif r=='Nad':
        newr='NH1  '
    elif r=='CS1':
        newr='CT1  '
    elif r=='Sth':
        newr='S    '
    elif r=='HN':
        newr='HR3  '
    elif r=='HC':
        newr='HA   '
    elif r=='CS2':
        newr='CT2  '
    elif r=='HS':
        newr='HS   '
    elif r=='SG5':
        newr='S    '
    elif r=='CS3':
        newr='CT3  '
    elif r=='Oet':
        newr='OS   '
    elif r=='Ohx':
        newr='OH1  '
    elif r=='Nad3':
        newr='NH2  '
    elif r=='Nim':
        newr='NR2  '
    elif r=='CR':
        newr='CA   '
    elif r=='CRp':
        newr='CA   '
    elif r=='NG22':
        newr='NH2  '
    elif r=='HO':
        newr='H    '
    elif r=='HR':
        newr='HP   '
    elif r=='CD1':
        newr='CC   '
    else:
        print('\n\nUNDEFINED ATOM TRANSLATIONS! '+str(r)+'\n\n')
        newr=r
    return newr
charmmatoms=[]
for atom in atoms:
    rtype=atom[10:14].strip()
    charmmtype=rosettatocharmm(rtype)
    newline=atom[:15]+charmmtype+atom[16:]
    charmmatoms.append(newline)


of=open('sam_fixed.params','w')
for line in lines:
    if line[0:4]!='ATOM':
        print(lines.index(line))

#
for line in lines[:4]:
    of.write(line)

for newline in charmmatoms:
    of.write(newline)

for line in lines[40:]:
    of.write(line)

of.close()
