# -*- coding: utf-8 -*-
"""
Created on Tue Aug 21 17:01:33 2018

@author: hraanan
"""
import itertools

from find_sphere_seq_full_range import get_main_chain

in_file=open('microenvironment_seq_with_neibor_res_pos_form_merge_metal_organic_chl_ca_10_md_3_rmsd_5_factor_correct_EC_1.4_1.3_1.15_no_ADE.txt','r')
in_file.readline()
out_file=open('sphere_overlap_score.txt','w')
out_file.write('source\ttarget\tsource_group\ttarget_group\tsource_size\ttarget_size\tminimum_size\toverlap\toverlap_residues\n')

#sphere_res_range_file=open('sphere_res_range_seq_full.txt','w')

microen_dict={}
microen_list=[]
for line in in_file:
    line=line.split('\t')
    #if line[1] != '7': continue
    microen_list.append(line[0])
    if line[0] not in microen_dict:
        microen_dict[line[0]]=[line[1],line[6].split('-')[0],line[6].split('-')[1]] #dict of {microen[group,sphere start,sphere end]}

i=0
for microen_i in microen_list:
    try:
        prot_i=microen_i.split('.')[0]
        chain_i=microen_i.split('_')[1]
        i=i+1
        j=0
        #print(i)
        for microen_j in microen_list:
           if microen_dict.get(microen_i)[0]==microen_dict.get(microen_j)[0]:continue # if its from the same group continue
           overlap_list=[]     
           j=j+1
           if j<i:continue
           prot_j=microen_j.split('.')[0]
           chain_j=microen_j.split('_')[1]
           if prot_i==prot_j and chain_i==chain_j: # if two microen on the same chain check for overlap
                
                value_i=list(map(int,microen_dict.get(microen_i)[1:]))
                value_j=list(map(int,microen_dict.get(microen_j)[1:]))
                #print('value_1',value_i)
                x=range(value_i[0],value_i[1])
                y=range(value_j[0],value_j[1])
                #print('x',x)
                xs=set(x)
                #print('xs',xs)
                overlap_list=[]
                if len(xs.intersection(y))>0:
                    range_list_i=get_main_chain(microen_i)
                    range_list_j=get_main_chain(microen_j)
                    #print(range_list_i)
                    ii=0
                    xs=[]
                    ys=[]
                    for x in range_list_i:
                        ii=ii+1
                        x=list(range(x[0],x[1]))
                        xs=xs+x
                        #xs=set.union(x)
                        jj=0
                    xs=set(xs)
                    #print('xs',xs)
                    for y in range_list_j:
                        jj=jj+1
                        #if jj<ii:continue
                        y=list(range(y[0],y[1]))
                        ys=ys+y
                    ys=set(ys)
                    #print('ys',ys)
                    overlap_list=overlap_list+list(xs.intersection(ys))
                    overlap_list=list(set(overlap_list))
                    #print(len(overlap_list),overlap_list)
                   # out_file.write(microen_i+'\t'+microen_j+'\t'+str(len(overlap_list))+'\t'+';'.join(list(map(str,overlap_list)))+'\n')            
                    small_sphere_size=min([len(xs),len(ys)])
                    size='\t'.join([str(len(xs)),str(len(ys)),str(small_sphere_size)])
                    out_file.write(microen_i+'\t'+microen_j+'\t'+str(microen_dict.get(microen_i)[0])+'\t'+str(microen_dict.get(microen_j)[0])+'\t'+size+'\t'+str(len(overlap_list))+'\t'+';'.join(list(map(str,overlap_list)))+'\n')            
    except:
        continue
                       
##out_file.close()
print('end')