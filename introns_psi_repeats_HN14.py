# -*- coding: utf-8 -*-
"""
the program:
1. removes bad introns and clusters from the excised introns count table
2. converts the counts in each repeat to percents
3. gets the significant results
4. plots heatmaps
"""

import pandas as pd
import os, sys
import shutil
sys.path.insert(0, "Z:/Analysis/general_scripts")
from leafcutter_functions_HN2 import read_table
from leafcutter_functions_HN2 import psi_heatmap
from leafcutter_functions_HN2 import add_deltapsi
from leafcutter_functions_HN2 import AS_events
from leafcutter_functions_HN2 import percent
from leafcutter_functions_HN2 import plot_bars
from leafcutter_functions_HN2 import get_limits
import numpy as np
import re 

def main():
    min_expression = 10 #min expression of intron
    delatpsi_cutoff = 0.2 #DS clusters - cluster with this value in at least 1 intron
    min_delta_psi = 0.1 #DS clusters - clusters with at least 2 shared introns with this value or higher
    n = 500 #number of introns to use for heatmap
    version = "_HN14"
    part = 2 # 1 or 2 to run part 1 or part 2
    #human GSE60424
    #name = 'GSE60424'
    #dire = 'Z:/Analysis/ImmGen/bulk_RNA/GSE60424/leafcutter/'
    #dire = 'Z:/Analysis/ImmGen/bulk_RNA/GSE60424/leafcutter_all/'
    #input_expression = dire + 'GSE60424_M10_perind_numers.counts'
    #order_list = ['Bcells', 'CD8', 'CD4', 'NK', 'Monocytes', 'Neutrophils'] #human dataset
    #order_list = ['Bcells', 'CD8', 'CD4', 'Monocytes', 'Neutrophils'] #human dataset
    #disease = ['HealthyControl', 'MSposttreatment', 'Type1Diabetes','MSpretreatment', 'ALS', 'Sepsis']
    #disease = ['MSpretreatment','MSposttreatment']
    #to_replace = ['B-cells', 'Bcells','_', '']
    #folders_list = ['Bcells_CD4','Bcells_CD8','Bcells_Monocytes','Bcells_Neutrophils','Bcells_NK','CD4_CD8','Monocytes_CD4','Monocytes_CD8','Monocytes_Neutrophils','Monocytes_NK','Neutrophils_CD4','Neutrophils_CD8','Neutrophils_NK','NK_CD4','NK_CD8']
    #repeats = 4 #how many repeats there are for each cell
    
    #human GSE64655
    name = 'GSE64655'
    dire = 'Z:/Analysis/ImmGen/bulk_RNA/GSE64655/leafcutter/'
    #dire = 'Z:/Analysis/ImmGen/bulk_RNA/GSE64655/leafcutter/TIV/'
    input_expression = 'Z:/Analysis/ImmGen/bulk_RNA/GSE64655/leafcutter/GSE64655_M10_perind_numers.counts'
    order_list = ['B_0', 'T_0', 'NK_0', 'MO_0', 'NEU_0'] #human dataset with out DC
    #order_list = ['B_0', 'B_1', 'B_3', 'B_7', 'T_0', 'T_1', 'T_3', 'T_7', 'NK_0', 'NK_1', 'NK_3', 'NK_7', 'DC_0', 'DC_1', 'DC_3', 'DC_7', 'MO_0', 'MO_1', 'MO_3', 'MO_7', 'NEU_0', 'NEU_1', 'NEU_3', 'NEU_7'] #
    to_replace = ['.sort.bam', '']
    folders_list = ['B_0_MO_0','B_0_NEU_0','B_0_NK_0','B_0_T_0','MO_0_NEU_0','MO_0_NK_0','MO_0_T_0','NEU_0_NK_0','NEU_0_T_0','NK_0_T_0']
    #folders_list = ['B_0_B_1','B_0_B_3','B_0_B_7','T_0_T_1','T_0_T_3','T_0_T_7', 'NK_0_NK_1','NK_0_NK_3','NK_0_NK_7', 'DC_0_DC_1','DC_0_DC_3','DC_0_DC_7', 'MO_0_MO_1','MO_0_MO_3','MO_0_MO_7', 'NEU_0_NEU_1','NEU_0_NEU_3','NEU_0_NEU_7']
    #repeats = 2
    
    #GSE122597
    #name = 'GSE122597'
    #dire = 'Z:/Analysis/ImmGen/GSE122597/leafcutter/'
    #input_expression = dire + 'GSE122597_M10_perind_numers.counts'
    #order_list = ['B', 'T4', 'Tgd', 'NK', 'MF']
    #to_replace = ['.sort.bam', '', 'T_4_Nve', 'T4']
    #folders_list = ['B_T4','B_Tgd','B_NK','B_MF','T4_Tgd','T4_NK','T4_MF','Tgd_NK','Tgd_MF','MF_NK'] 
    #repeats = 15
    
    #GSE124829 - female male
    #name = 'GSE124829'
    #dire = 'Z:/Analysis/ImmGen/Female_Male/leafcutter/'
    #input_expression = dire + 'Female_Male_M10_perind_numers.counts'
    #order_list = ['B1ab', 'B', 'T4', 'T8', 'Tgd','Treg', 'NKT', 'NK', 'Gn', 'DC','MF']
    #to_replace = ['.sort.bam', '', 'female', 'F', 'male', 'M']
    #folders_list = ['B1ab_DC','B1ab_MF','B1ab_NK','B1ab_NKT','B1ab_T4','B1ab_T8','B1ab_Tgd','B1ab_Treg','B1ab_Gn','B_B1ab','B_DC','B_Gn','B_MF','B_NK','B_NKT','B_T4','B_T8','B_Tgd','B_Treg','DC_Gn','DC_MF','DC_NK','DC_NKT','DC_T4','DC_T8','DC_Tgd','DC_Treg','Gn_MF','Gn_NK','Gn_NKT','MF_NK','MF_NKT','NK_NKT','T4_Gn','T4_MF','T4_NK','T4_NKT','T4_T8','T4_Tgd','T4_Treg','T8_Gn','T8_MF','T8_NK','T8_NKT','T8_Tgd','T8_Treg','Tgd_Gn','Tgd_MF','Tgd_NK','Tgd_NKT','Treg_Gn','Treg_MF','Treg_NK','Treg_NKT','Treg_Tgd']
    #repeats = 2 
    
    leafcutter_output = dire + 'leafcutter_expression_table_' + name + version + '.txt'
    psi_persent_output = dire + 'psi_persent_' + name + version + '.xlsx'
    comparison_output = dire + 'comparisons_deltapsi_' + name + version + '.xlsx'
    AS_event_output = dire + 'AS_events_' + name + version + '.xlsx'
    sashmi_output = dire + 'sashimi_list_' + name + version + '.txt'
    heatmap_output_1 = dire + 'psi_heatmap_' + name + version + '.jpg'
    heatmap_output_2 = dire + 'psi_heatmap_' + name + '_no_col' + version + '.jpg'  
   
    #GSE143943 - edy
    #dire = 'Z:/Analysis/LPS_Brenner_2016/leafcutter/'
    #input_expression = dire + 'LPS_M10_perind_numers.counts'
    #order_list = ['Tgd_','NKT_', 'NK_','MF6Cplus_', 'MF6C__']
    #to_replace = ['.sort.bam', '', 'dotSpNo', '_','dot','']
    #folders_list = ['MF6C__MF6Cplus','MF6C__NK','MF6C__NKT','MF6Cplus_NK','MF6Cplus_NKT','MF6Cplus_Tgd','MF6C__Tgd','NK_NKT','NK_Tgd','NKT_Tgd']
    #order_list = [w.replace('_', '') for w in order_list]
    #order_list[4] = order_list[4].replace('MF6C', 'MF6C_')
    
    if part==1:
        print('load expression table')
        leafcutter_order, order_list_repeats = read_table(input_expression, order_list, to_replace)
        cluster_list = list(leafcutter_order.cluster.unique())  #list of unique clusters
        print('there are ',len(leafcutter_order), 'introns in leafcutter table')
        print('there are ',len(cluster_list), 'clusters in leafcutter table')
    
        leafcutter_order.reset_index(inplace=True)
    #join clumns to get introns names
        leafcutter_order['intron'] = leafcutter_order.iloc[:,0:4].apply(lambda x: ":".join(x.astype(str)), axis=1)
        leafcutter_order.drop(['level_0', 'level_1', 'level_2'], axis=1, inplace=True)
    #remove introns that are not express in any of the samples
        leafcutter_order['max_counts'] = leafcutter_order[order_list_repeats].max(axis=1)
        introns_max = leafcutter_order[['intron','max_counts']] #save max count for each intron       
        leafcutter_order.drop(leafcutter_order.loc[leafcutter_order['max_counts']<min_expression].index, inplace=True)
        leafcutter_order.drop(['max_counts'], axis=1, inplace=True)
        cluster_list = list(leafcutter_order.cluster.unique())  #list of unique clusters
        print('there are ',len(leafcutter_order), 'introns>10 reads in leafcutter table')
        print('there are ',len(cluster_list), 'clusters in leafcutter second table')
    
        print('remove bad introns and cluster and create new leafcutter expression table and psi_percent')
    #psi table for each repeat, based on percent from the expression table
    #remove introns with no shared positions and create new exprssion table
        first = 1
        index = 0
        for one_cluster in cluster_list:
            if index%1000 == 0:
                print(one_cluster)
            index += 1
            nice_cluster_express = leafcutter_order.loc[leafcutter_order['cluster'] == one_cluster]
            nice_cluster_express = get_limits(nice_cluster_express)
            nice_cluster_express = nice_cluster_express[1]       
            if len(nice_cluster_express) > 1: #only clusters with at least 2 introns
                nice_cluster_persent = percent(nice_cluster_express)
                if first:
                    psi_final_percent = nice_cluster_persent
                    leafcutter_final_table = nice_cluster_express
                    first = 0
                else:
                    frames = [psi_final_percent, nice_cluster_persent]
                    psi_final_percent = pd.concat(frames)
                    frames = [leafcutter_final_table, nice_cluster_express]
                    leafcutter_final_table = pd.concat(frames)
        cluster_list = list(psi_final_percent.cluster.unique())  #list of unique clusters
        print('there are ',len(psi_final_percent), 'introns that share location in psi persent table')
        print('there are ',len(cluster_list), 'clusters in psi persent table')
    #make the final leafcutter table look like the original expression table
        leafcutter_final_table.drop(['cluster'], axis=1, inplace=True)
        col_names = list(leafcutter_final_table.columns)
        new_col_names = [col_names[len(col_names)-1]] + col_names[0:len(col_names)-1]
        leafcutter_final_table = pd.DataFrame(leafcutter_final_table, columns = new_col_names)
        leafcutter_final_table.to_csv(leafcutter_output, index = None, sep = ' ') 
        psi_final_percent.to_excel(psi_persent_output) #after part 2 i add to this file more data   
    
        #create group file in each directory
        leafcutter_col = leafcutter_final_table.columns
        for i in range(1,len(disease)):
            for cell in order_list:
                folder = cell + disease[0] + "_" + cell +  disease[i]
                print(folder)
                folder_path = dire + folder
                #print(folder_path)
                #remove exist dir and make new one
                if os.path.exists(folder_path):
                    shutil.rmtree(folder_path)
                os.mkdir(folder_path, 770)
                group_file = folder_path + "/groups_file_HN13.txt"
                fout = open(group_file, 'w') 
                look_for = cell + disease[0]
                r = re.compile(look_for)
                newlist = list(filter(r.match, leafcutter_col))  #Determine if the RE matches at the beginning of the string.
                for groupName in newlist:
                    toWrite = groupName + " " + cell + disease[0] + "\n"
                    fout.write(toWrite)
                look_for = cell + disease[i]
                r = re.compile(look_for)
                newlist = list(filter(r.match, leafcutter_col))  #Determine if the RE matches at the beginning of the string.
                for groupName in newlist:
                    toWrite = groupName + " " + cell + disease[i] + "\n"
                    fout.write(toWrite)
                fout.close()
    #run leafcutter DS analysis with the new expression table first
    if part==2:
        #get folders list if all the folders in the dir are needed
        #files_list = os.listdir(dire)
        #folders_list = []
        #for i in files_list:
        #    if os.path.isdir(dire + i):
        #        folders_list.append(i)   
    #calculate deltapsi for each intron in each comparison
        print('calculate deltapsi')
        psi_final_percent = pd.read_excel(psi_persent_output, index_col=0)   #load the previous table, if i want to run only this part.
        comparison_data = pd.DataFrame(index = psi_final_percent['intron'], columns = folders_list) 
        #make the col names as the names in the directories
        #col_names = psi_final_percent.columns
        #res = list(map(lambda st: str.replace(st, "_", ""), col_names))
        #psi_final_percent.columns = res
        for folder in folders_list:
            DS = folder.split("_")  #if folder is cell_time_cell_time it is not working
            r = re.compile(DS[0])
            cell_col = list(filter(r.search, psi_final_percent.columns)) 
            one_cell_df = psi_final_percent[['intron'] + cell_col]
            one_cell_df = one_cell_df.set_index('intron')
            one_cell_df['mean'] = one_cell_df.mean(axis=1)       
            r = re.compile(DS[2]) #if folder is cell_time_cell_time it is not working when looking for DS[1] ******************************
            cell_col = list(filter(r.search, psi_final_percent.columns)) 
            second_cell_df = psi_final_percent[['intron'] + cell_col]
            second_cell_df = second_cell_df.set_index('intron')
            second_cell_df['mean'] = second_cell_df.mean(axis=1)   
            delta_psi = abs(second_cell_df['mean'] - one_cell_df['mean'])
            comparison_data[folder] = delta_psi
    
        comparison_data['intron'] = comparison_data.index    
        comparison_data['chr_split'] = comparison_data['intron'].str.split(':')
        comparison_data['cluster'] = 'NA'
        for i in range(0,len(comparison_data)):
            cluster_name =  comparison_data['chr_split'][i] 
            cluster_name = cluster_name[3]
            comparison_data['cluster'][i] = cluster_name
        comparison_data.drop(['chr_split'], axis=1, inplace=True)
        comparison_data.to_excel(comparison_output) 
    
        print('add deltapsi')
    #deltapsi table for the introns in each comparison  
        #the name of the input is in the function. i need to remove it from there
        psi_final_percent = add_deltapsi(psi_final_percent, dire, folders_list, psi_persent_output, comparison_data) 
        psi_final_percent_sucess =  psi_final_percent.loc[psi_final_percent['status'] == 'Success']
        cluster_list_sucess = list(psi_final_percent_sucess.cluster.unique())  #li
        print('there are ',len(psi_final_percent_sucess), 'sucess introns')
        print('there are ',len(cluster_list_sucess), 'clusters in sucess table')   
    
        print('AS_table')  #all the introns from AS event intron should remain
    #psi_final_percent = pd.read_excel(psi_persent_output, index_col=0)   #load the previous table, if i want to run only this part
        AS_event_table = AS_events(psi_final_percent, delatpsi_cutoff, min_delta_psi, sashmi_output)
        AS_event_table.drop(['status'], axis=1, inplace=True)
        AS_event_table = pd.merge(AS_event_table, introns_max, on=['intron'])   #intron_max is from part 1. need to find a way to get it
        #order by disease and not by cell
        order_temp = []
        for cell in disease:
            r = re.compile(cell)
            #newlist = list(filter(r.match, table_col))  #Determine if the RE matches at the beginning of the string.
            newlist = list(filter(r.search, AS_event_table.columns)) #Scan through a string, looking for any location where this RE matches.
            order_temp = order_temp + newlist
        new_col = list(AS_event_table.columns)
        new_col[1:-5] = order_temp
        AS_event_table = pd.DataFrame(AS_event_table, columns = new_col)
        AS_event_table.to_excel(AS_event_output, index = False)   
        cluster_list = list(AS_event_table.cluster.unique()) 
        genes_list = list(AS_event_table.genes.unique()) 
        psi_heatmap(AS_event_table, n, heatmap_output_1, heatmap_output_2) 
        AS_high_deltapsi =  AS_event_table.loc[AS_event_table['deltapsi_max'] >= 0.1]
        print('there are ',len(AS_event_table), 'AS events introns')
        print('there are ',len(AS_high_deltapsi), 'introns with deltapsi>=0.1')
        print('there are ',len(cluster_list), 'clusters in AS table')
        print('there are ',len(genes_list), 'genes in AS table')
       
    #plot the cluster
        print('plot data')
        index = 0
        for one_cluster in cluster_list:
            if index%50 == 0:
                print(one_cluster)
            index += 1
            nice_cluster_persent = psi_final_percent.loc[psi_final_percent['cluster'] == one_cluster]
            gene = list(nice_cluster_persent.genes.unique()) 
            introns = nice_cluster_persent['intron']
            nice_cluster_express = pd.merge(leafcutter_order, introns,  on=['intron'])
            nice_cluster_log = nice_cluster_express.replace(0, 1)
            nice_cluster_log = np.log2(nice_cluster_log.loc[:,nice_cluster_log.columns[1:-1]])
            nice_cluster_persent = pd.DataFrame(nice_cluster_persent, columns = nice_cluster_log.columns)
        #bars_output = dire + 'bars/' + one_cluster + '_' + str(gene[0]) + '.eps'
            bars_output = dire + 'bars/' + one_cluster + '_' + str(gene[0]) + '.jpg'
            bars_output = bars_output.replace('NA_', '')
            plot_bars(nice_cluster_persent, nice_cluster_express, nice_cluster_log, one_cluster, gene[0], bars_output, repeats)
    return

if __name__ == "__main__":
    main() 