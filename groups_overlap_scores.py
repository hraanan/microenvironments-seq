# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 10:40:01 2018

@author: hraanan
"""
import statistics
from collections import defaultdict

groups_file=open('microenvironment_seq_with_neibor_res_pos_form_merge_metal_organic_chl_ca_10_md_3_rmsd_5_factor_correct_EC_1.4_1.3_1.15_no_ADE.txt','r')
groups_file.readline()
group_dict={}
#for line in groups_file:
#    line=line.split('\t')
#    group_dict[line[0]]=int(line[2])





dates_dict = defaultdict(list)

group_dict={}
microens=[]
microen_groups={}
ec_groups={}
prot=[]
for line in groups_file:
    line=line.split('\t')
    microen_groups[line[0]]=line[1]
in_file=open('microenvironment_seq_with_neibor_res_pos_form_merge_metal_organic_chl_ca_10_md_3_rmsd_5_factor_correct_EC_1.4_1.3_1.15_no_ADE.txt','r')
for line in in_file:
    line=line.split('\t')
    if line[1] not in group_dict:
        group_dict.setdefault(line[1],[]).append(line[0])
    if line[1] in group_dict:
        microens=group_dict.get(line[1])
        if line[0] not in microens:
            microens.append(line[0])
            microens=sorted(microens)


edge_file=open('sphere_overlap_score.txt','r')
out_file=open('edges_ave_overlap_score.txt','w')
out_file.write('groups\tave overlap\tstdev overlap\tave small sphere size\tstdev small sphere size\tnumber of overlaps\tsmall group size\tnormalized to small group\n')
edge_file.readline()
edge_dict={}
for line in edge_file:
        line=line.split('\t')
        if line[2]==line[3]:continue # if from the same group 
        
        if int(line[7])<1:continue     # if no overlap
        edge=line[2]+'-'+line[3]
        score_list=[]
        if edge in edge_dict:
            #print(edge_dict)
            score_list=edge_dict.get(edge)
            #print('sco',score_list)
            #print(type(score_list))
            score_list.append([int(line[6]),int(line[7])])
            #print('sco',score_list)
           
            edge_dict[edge]=score_list
        if edge not in edge_dict:
            edge_dict[edge]=[[int(line[6]),int(line[7])]]
        
        #print(len(edge_dict))
    
for key,value in edge_dict.items():
    
    size_list=[]
    overlap_list=[]
   #print(value)
    for i in value:
        size_list.append(i[0])
        overlap_list.append(i[1])
    #print(key,len(size_list),len(overlap_list))
    group_1_size=0
    group_2_size=0
    if len(overlap_list)>1:
        size_ave=sum(size_list) / float(len(size_list))
        size_std=statistics.stdev(size_list)
        
        overlap_ave=sum(overlap_list) / float(len(overlap_list))
        overlap_std=statistics.stdev(overlap_list)
        
        
    else:
        size_ave=size_list[0]
        size_std=0
        
        overlap_ave=overlap_list[0]
        overlap_std=0
    
    group_1_size=len(group_dict.get(key.split('-')[0]))
    group_2_size=len(group_dict.get(key.split('-')[1]))
    print(group_1_size,group_2_size)
    small_group_size=min([group_1_size,group_2_size])    
    #print((key+'\t'+str(overlap_ave)+'\t'+str(overlap_std)+'\t'+str(size_ave)+'\t'+str(size_std)+'\t'+str(len(overlap_list))+'\t'+str(small_group_size)+'\t'+str(len(overlap_list)/small_group_size)+'\n'))
    out_file.write(key+'\t'+str(overlap_ave)+'\t'+str(overlap_std)+'\t'+str(size_ave)+'\t'+str(size_std)+'\t'+str(len(overlap_list))+'\t'+str(small_group_size)+'\t'+str(len(overlap_list)/small_group_size)+'\n')
#    except:
#        continue
out_file.close()
print('end')
