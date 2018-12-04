# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 15:36:03 2018

@author: hraanan
"""
import Bio
import os
import sys
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import three_to_one
from Bio.PDB.Polypeptide import is_aa
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



from Bio import PDB
from Bio.PDB import PDBIO

import math
import numpy
from collections import Counter
import random 
from Bio.PDB import *
#import gzip







split=False
 
rootdir ='f:/pdb_download_8_2018/pdb/' # folder of the pdb original filesdownladed form the pdb
in_file=open('microenvironment_seq_with_neibor_res_pos_form_merge_metal_organic_chl_ca_10_md_3_rmsd_5_factor_correct_EC_1.4_1.3_1.15_no_ADE.txt','r')
neighbors_res_range_file=open('neighbors_res_range_seq_full.txt','w')
#neighbors_res_range_file=open('neighbors_res_range_seq_splited.txt','w')
#sphere_res_range_file=open('sphere_res_range_seq_full.txt','w')
in_file.readline()

pdbl=PDB.PDBList()
microen_list=[]
for line in in_file:
    line=line.split('\t')
    microen_list.append(line)
in_file=open('microenvironment_seq_with_neibor_res_pos_form_merge_metal_organic_chl_ca_10_md_3_rmsd_5_factor_correct_EC_1.4_1.3_1.15_no_ADE.txt','r')
in_file.readline()
for line in in_file:
    line=line.split('\t')
    #if line[1] != '7': continue
    #print (line[0],line[1])
    main_chain=line[4]
    res_range=line[5].split('-')
    protein=line[0].split('.')[0]
    overlap=False
    if split==True:
        for microen in microen_list:
            prot=microen[0].split('.')[0]
            chain=microen[4]
               
            if microen[0]!=line[0] and prot==protein and main_chain==chain:
                microen_range=microen[5].split('-')
                #print(res_range,microen_range)
                
                if int(res_range[1])>int(microen_range[0]) and int(res_range[0])<int(microen_range[0]) :
                    overlap=True
                    print(res_range,microen_range)
                    res_range[1]=str(int(microen_range[0])-1)
                    
                elif int(res_range[0])<int(microen_range[1]) and int(res_range[1])>int(microen_range[1]):
                    overlap=True
                    print(res_range,microen_range)
                    res_range[0]=str(int(microen_range[1])+1)
                    #print(overlap)
                    #if overlap==True: print('double overlap',line[0])     
    #print ('pdb_code:'+protein)
    protein=protein.lower()
    parser = PDB.PDBParser(PERMISSIVE=1,get_header=1,QUIET=1)
   # print(rootdir+protein[1:3]+'/pdb'+protein+'.ent')
    structure = parser.get_structure(protein,rootdir+protein[1:3]+'/pdb'+protein+'.ent')
     
    i=0
    for model in structure:
        
        if i<1:
            i=i+1 #if model!=0: not working for unknown reason
        #if model!=0:
            #print('yes') 
            for chain in model:
                
                seq = list()
                #print(chain.get_id(),)
                if chain.get_id()==main_chain:
                 
                    for residue in chain:
                        if residue.id[1]<1:continue
                        ## The test below checks if the amino acid
                        ## is one of the 20 standard amino acids
                        ## Some proteins have "UNK" or "XXX", or other symbols
                        ## for missing or unknown residues
                        if is_aa(residue.get_resname(), standard=True):
                            seq.append(three_to_one(residue.get_resname()))
                        else:
                            continue
                            #seq.append("X")
                     
                    ## This line is used to display the sequence from each chain
                     
                    
                    
                    #print (">Chain_" + chain.get_id() + "\n" + str("".join(seq)))
                     
                    ## The next two lines create a sequence object/record
                    ## It might be useful if you want to do something with the sequence later
                     
    #                myProt = Seq(str(''.join(seq)), IUPAC.protein)
    #                seqObj = SeqRecord(myProt, id=chainID, name="", description="")
    #                out_file.write(line[0]+'\t'+line[1]+'\t'+line[2]+'\t'+line[4]+'\t'+res_range[1][0]+'\t'+'-'.join(list(map(str,res_range[2])))+'\t'+'-'.join(list(map(str,res_range[3])))+'\n')
                    
                    neighbors_res_range_file.write('>'+line[0]+':'+line[1]+':'+'-'.join(res_range)+'\n')
                    #neighbors_res_range_file.write(''.join(seq)+'\n')
                    neighbors_res_range_file.write(''.join(seq[int(res_range[0])-1:int(res_range[1])])+'\n')
    #                sphere_res_range_file.write('>'+line[0]+':'+'-'.join(list(map(str,res_range[3])))+'\n')
    #                sphere_res_range_file.write(res_range[5]+'\n')

#out_file.close()
neighbors_res_range_file.close()
print('end')         
## The end
 