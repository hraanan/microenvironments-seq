# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 16:04:43 2018

@author: hraanan
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 11:38:55 2016

@author: hraanan
"""
import networkx as nx
G=nx.Graph()
giant=nx.Graph()
#mypath='/home/hraanan/pymol/PymolAlign/postaligin/grouping_by_network/factor_1.5/'
in_file=open('HHalign_sphere_101118_woverlap_101918_filter10.txt','r')
out_file=open('seq_scoure_network_distances.txt','w')
#node_file=open('group_48_long_paths_nodes.txt','w')
cofactors={}
groups={}
in_file.readline()
#node_file.write('ID\tLP\n')
#for line in in_file:
#    line=line.split('\t')
#    groups.setdefault(line[1],[]).append(line[0])
##groups_list=['HEM_3','CD_14','SF4_5','FE_2','FE2_0','HEC_1','CU_11','CU_2','CD_22','CD_144','CD_210','FE2_4','FES_0','HEM_8','NI_18','CD_46','HEC_4','HEC_5','F3S_0','SF4_10','CD_32','HEC_2','SF4_2']
#groups_list=['48']
#groups_list_dict={}
#for group,microen_list in groups.items():
#    if group in groups_list:
#        groups_list_dict[group]=microen_list
        
for line in in_file:
    #G.clear() 
    #giant.clear()
    
    #edge_file=open('/home/hraanan/pymol/PymolAlign/postaligin/binding_motifs/align_all_ligands_and_nonredundant_ca_10_md_3_rmsd_5_factor_1.4.txt','r')
    #edge_file.readline()
    #for line in edge_file:
    line=line.split('\t')
    #if line[0] in microen_list:
    G.add_edge(line[0],line[1],weight=float(line[4]))

    #print(G.number_of_nodes())
    #if G.number_of_nodes() < 10:
     #   continue
    #print(group)    
    #print(G.nodes())
       # reps=groups.get(group)
    #print(microen_lis)
    short_paths={}
    long_paths={}     
#    print (G.nodes())
    #giant = max(nx.connected_component_subgraphs(G), key=len)    
    i=0
    #for nod_i in giant.nodes():
    
    for nod_i in G.nodes():
        #min_path=1000000 
        #max_path=0
        i=i+1
        j=0
        for nod_j in G.nodes():
            j=j+1
            if i==j:
                continue
                       
            try:
                short_path=nx.shortest_path(G,source=nod_i,target=nod_j)
            except:# NetworkXNoPath:
                out_file.write(nod_i+'\t'+nod_j+'\tNA\n')
                continue
    #            if len(short_path)<min_path:
    #                short_paths[i]=short_path
    #                min_path=len(short_path)
    #            if len(short_path)>max_path:
    #                long_paths[i]=short_path
    #                max_path=len(short_path)
        
    #    for key,path in long_paths.items():
    #        if len(path)==max_path:
    #            for i in path:
    #                node_file.write(i+'\tLP\n')
            out_file.write(nod_i+'\t'+nod_j+'\t'+str(len(short_path))+'\n')
            
                           
in_file.close()
out_file.close()
##    cofactor= group.split('_')[0]
##    groups_file_name=mypath+'cofactors_groups/'+cofactor+'_groups.txt'
##    edges_file_name=mypath+'cofactors_graphs/'+cofactor+'_weighted_edgelist.edgelist'
##    lines = open(edges_file_name).readlines()
##    open('temp.txt', 'w').writelines(lines[1:])    
##    lines.clear()
#    for node in microen_list:
#        
    #G=nx.read_weighted_edgelist('pry_1_Edges.csv',delimiter='\t')