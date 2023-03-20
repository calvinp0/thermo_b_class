#molstats.py
##Dictionary of all elements with molecular weight


MM_of_Elements = {'H': 1.00794, 'He': 4.002602, 'Li': 6.941, 'Be': 9.012182, 'B': 10.811, 'C': 12.0107, 'N': 14.0067,
              'O': 15.9994, 'F': 18.9984032, 'Ne': 20.1797, 'Na': 22.98976928, 'Mg': 24.305, 'Al': 26.9815386,
              'Si': 28.0855, 'P': 30.973762, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.0983, 'Ca': 40.078,
              'Sc': 44.955912, 'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938045,
              'Fe': 55.845, 'Co': 58.933195, 'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.409, 'Ga': 69.723, 'Ge': 72.64,
              'As': 74.9216, 'Se': 78.96, 'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90585,
              'Zr': 91.224, 'Nb': 92.90638, 'Mo': 95.94, 'Tc': 98.9063, 'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.42,
              'Ag': 107.8682, 'Cd': 112.411, 'In': 114.818, 'Sn': 118.71, 'Sb': 121.760, 'Te': 127.6,
              'I': 126.90447, 'Xe': 131.293, 'Cs': 132.9054519, 'Ba': 137.327, 'La': 138.90547, 'Ce': 140.116,
              'Pr': 140.90465, 'Nd': 144.242, 'Pm': 146.9151, 'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25,
              'Tb': 158.92535, 'Dy': 162.5, 'Ho': 164.93032, 'Er': 167.259, 'Tm': 168.93421, 'Yb': 173.04,
              'Lu': 174.967, 'Hf': 178.49, 'Ta': 180.9479, 'W': 183.84, 'Re': 186.207, 'Os': 190.23, 'Ir': 192.217,
              'Pt': 195.084, 'Au': 196.966569, 'Hg': 200.59, 'Tl': 204.3833, 'Pb': 207.2, 'Bi': 208.9804,
              'Po': 208.9824, 'At': 209.9871, 'Rn': 222.0176, 'Fr': 223.0197, 'Ra': 226.0254, 'Ac': 227.0278,
              'Th': 232.03806, 'Pa': 231.03588, 'U': 238.02891, 'Np': 237.0482, 'Pu': 244.0642, 'Am': 243.0614,
              'Cm': 247.0703, 'Bk': 247.0703, 'Cf': 251.0796, 'Es': 252.0829, 'Fm': 257.0951, 'Md': 258.0951,
              'No': 259.1009, 'Lr': 262, 'Rf': 267, 'Db': 268, 'Sg': 271, 'Bh': 270, 'Hs': 269, 'Mt': 278,
              'Ds': 281, 'Rg': 281, 'Cn': 285, 'Nh': 284, 'Fl': 289, 'Mc': 289, 'Lv': 292, 'Ts': 294, 'Og': 294,
              'ZERO': 0,'D':2.014,'T':3.01604}


#Class Molecule with name as string in standard format. Characters allowed are all combo of chemical Symbols, 2 digit number and '(' and ')'
class Molecule:
    # init method or constructor 
    def __init__(self, name):
        if not type(name)==str:
            raise ValueError
        else:
            self.data = name
    # Sample Method 
    def calculate_num(self):
        s2=self.data
        s2=s2.rstrip().lstrip()
        err=0
        if s2[-1:].isupper():
            s2=s2+'1'
        if s2[0].isupper() and not(s2[1].islower()) and not (s2[1].isnumeric()):
            s2=s2[0]+'1'+s2[1:]
        s1=s2[::-1]
        n1=[]
        n2=[]
        cnt=0
        ns1=''
        temp_mul=1
        total_mul=1
        mul1=[]
        mul2=[]
        if len(s2)==1:
            n1.append(s2)
            n2.append(1)
        elif len(s2)==2:
            if s2[1].isnumeric():
                n1.append(s2[0])
                n2.append(int(s2[1]))
            elif s2[1].isalpha():
                if s2[1].islower():
                    n1.append(s2)
                    n2.append(1)
                elif s2[1].isupper():
                    n1.append(s2[0])
                    n1.append(s2[1])
                    n2.append(1)
                    n2.append(1)
        else:
            for i in range(1,len(s1)):
                if s1[i]==')':
                    if i>1:
                        if s1[i-1].isnumeric() and s1[i-2].isnumeric():
                            temp_mul=int(s1[i-1]+s1[i-2])
                        elif s1[i-1].isnumeric() and not s1[i-2].isnumeric():
                            temp_mul=int(s1[i-1])
                        elif not s1[i-1].isnumeric():
                            temp_mul=temp_mul
                    else:
                        if s1[i-1].isnumeric():
                            temp_mul=int(s1[i-1])
                        elif not s1[i-1].isnumeric():
                            temp_mul=temp_mul
                    total_mul=total_mul*temp_mul
                elif s1[i]=='(':
                    total_mul=int(total_mul/temp_mul)
                elif s1[i].isalpha() and s1[i-1]=='(':
                    if s1[i].isupper():
                        n1.append(s1[i])
                        n2.append(1*total_mul)
                elif s1[i].isalpha() and s1[i-1].isalpha():
                    err=5
                    if s1[i-1].islower() and s1[i].isupper():
                        if i<len(s1):
                            if s1[i-2].isnumeric():
                                n1.append(s1[i]+s1[i-1])
                                n2.append(int(s1[i-2])*total_mul)
                                err=2
                            else:
                                n1.append(s1[i]+s1[i-1])
                                n2.append(1*total_mul)
                                err=1
                        else:
                            n1.append(s1[i]+s1[i-1])
                            n2.append(1*total_mul)
                            err=4
                    elif s1[i].isupper() and s1[i-1].isupper():
                        n1.append(s1[i])
                        n2.append(1*total_mul)
                        err=6
                    else:
                        continue
                        err=7
                elif s1[i].isalpha() and s1[i-1].isnumeric():
                    if s1[i].isupper() and i>1:
                        if s1[i-1].isnumeric() and s1[i-2].isnumeric():
                            n1.append(s1[i])
                            n2.append(int(s1[i-1]+s1[i-2])*total_mul)
                        elif s1[i-1].isnumeric() and not s1[i-2].isnumeric():
                            n1.append(s1[i])
                            n2.append(int(s1[i-1])*total_mul)
                    elif s1[i].islower():
                        continue
                    else:
                        n1.append(s1[i])
                        n2.append(int(s1[i-1])*total_mul)
                mul1.append(temp_mul)
                mul2.append(total_mul)
        return n1,n2
    def getStats(self):
        nn1,nn2=self.calculate_num()
        return nn1,nn2
    def getMolecularWeight(self):
        nn1,nn2=self.calculate_num()
        mol_wt=[]
        for k in range(len(nn1)):
            mol_wt.append(MM_of_Elements[nn1[k]]*nn2[k])
        return sum(mol_wt)
    def getNumElements(self):
        nn1,nn2=self.calculate_num()
        num_el=sum(nn2)
        return num_el