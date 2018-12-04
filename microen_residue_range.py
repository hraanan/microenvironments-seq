# -*- coding: utf-8 -*-
"""
Created on Tue Aug 21 17:01:33 2018

@author: hraanan
"""

from find_main_chain_of_microen_sphere_range import get_main_chain

in_file=open('merge_metal_organic_chl_ca_10_md_3_rmsd_5_factor_correct_EC_1.4_1.3_1.15_no_ADE.txt','r')
in_file.readline()
#out_file=open('microenvironment_seq_with_neibor_res_pos_form_merge_metal_organic_chl_ca_10_md_3_rmsd_5_factor_correct_EC_1.4_1.3_1.15_no_ADE.txt','w')
#out_file.write('id\tgroup\tcofactor\tcofactor group\tchain\tbindig res range\tsphere res range\n')
dash_seq_file=open('test_dash_seq_file.txt','w')
#sphere_res_range_file=open('sphere_res_range_seq_full.txt','w')
for line in in_file:
    line=line.split('\t')
    if line[1] != '7': continue
    res_range=get_main_chain(line[0])
    #print(res_range)
    if res_range[0]==False and res_range[1]!='err':
        #out_file.write(line[0]+'\t'+line[1]+'\t'+line[2]+'\t'+line[4]+'\t'+res_range[1][0]+'\t'+'-'.join(list(map(str,res_range[2])))+'\t'+'-'.join(list(map(str,res_range[3])))+'\t'+';'.join(res_range[6])+'\n')
#        neighbors_res_range_file.write('>'+line[0]+':'+'-'.join(list(map(str,res_range[2])))+'\n')
#        neighbors_res_range_file.write(res_range[4]+'\n')
#        sphere_res_range_file.write('>'+line[0]+':'+'-'.join(list(map(str,res_range[3])))+'\n')
#        sphere_res_range_file.write(res_range[5]+'\n')
        dash_seq_file.write('>'+line[0]+':'+'-'.join(list(map(str,res_range[2])))+'\n')
        dash_seq_file.write(res_range[7]+'\n')
#out_file.close()
print('end')