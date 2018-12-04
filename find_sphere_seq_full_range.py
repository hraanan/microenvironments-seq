import random

from Bio import PDB
import math
import numpy
from Bio.PDB import *
import itertools

######initial variabels############# 
rootdir ='f:/microfolds_8_2018/all/' # folder of the pdb original filesdownladed form the pdb
 
AA=['PHE','TRP','TYR','ALA','CYS','ASP','GLU','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL']
chl=['BCB','CLA','CHL','BCL','CL0','PMR','PHO'] #chlorophyll cofactors list
heme=['HEA','HAS','2FH','522','89R','DDH','DHE','HES','HDD','HDE','HDM','HEB','HEC','HEM','HEO','HEV','HP5','MH0','N7H','NTE','OBV','SRM','VER']
pyr_atom_list=['C3A','C3B','C3C','C3D']

aa=['c','d','s','q','k','i','p','t','f','n','g','h','l','r','w','a','v','e','y','m']
aa_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
#####################
organic_cofactors_list=[]
cofactors_dict={}
organic_cofactors_pdb_file=open('manual_cofactor_list_with_quinone.txt','r')
for line in organic_cofactors_pdb_file:
    line=line.split('\t')
    organic_cofactors_list.append(line[1])
    if line[1] not in cofactors_dict:
        if len(line[4])>1:
            cofactors_dict[line[1]]=[line[3],line[4][:-1]]
        else:
            cofactors_dict[line[1]]=[line[3]]

pdbl=PDB.PDBList()




def to_ranges(iterable):
    iterable = sorted(set(iterable))
    for key, group in itertools.groupby(enumerate(iterable),
                                        lambda t: t[1] - t[0]):
        group = list(group)
        yield group[0][1], group[-1][1]



def get_main_chain(microen):
 
    
    protein=microen[0:4]
    microen='_'.join(microen.split('_')[:-2])
    cof=microen.split('_')[0].split('.')[1]
    parser = PDB.PDBParser(PERMISSIVE=1,get_header=1,QUIET=1)
    structure = parser.get_structure(protein,rootdir+cof+'/'+microen+'.pdb')
    model = structure[0]
    chain = model[microen.split('_')[1]]
    res_num_list=[]
    for residue in chain:
        if residue.id[1]==int(microen.split('_')[2]) and residue.get_resname()==cof: continue
        res_num_list.append(residue.id[1])
    return (list(to_ranges(res_num_list)))
                                                 
#print(get_main_chain('11gs.GSH_B_210_transferase_2.5.1.18'))
#y=get_main_chain('11gs.GSH_B_210_transferase_2.5.1.18')
#
##x=to_ranges([1,10)
#x=range(1,10)
#xs=set(x)
#
#for i in y[1]:
#    y=range(i[0],i[1])
#
#    print(y)
#    print(xs.intersection(y))