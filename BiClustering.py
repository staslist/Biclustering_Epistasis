# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 02:54:48 2024

@author: staslist
"""

import re
import csv
import unittest
import numpy as np
import os
#import time
from scipy.stats import hypergeom
from datetime import datetime
from numba import jit

# Bi-Clustering Algorithm Code

# Assume I have a file that lists all indices in the matrix that are 1.

# GLOBAL VARIABLES

# For tested interaction matrix, assume that every interaction has been tested (exhaustive algorithm)

def merge_list(input_list:list):
    output = []
    eliminated_indeces = set()
    outer_index = 0
    while outer_index < len(input_list):
        if(outer_index in eliminated_indeces):
            outer_index += 1 
            continue
        
        #to_merge_indeces = []
        current_set = input_list[outer_index]
        inner_index = 0 
        while inner_index < len(input_list):
            if(inner_index in eliminated_indeces):
                inner_index += 1
                continue
            
            if(input_list[outer_index] == input_list[inner_index]):
                inner_index += 1
                continue
            
            if(len(input_list[outer_index].intersection(input_list[inner_index])) > 0):
                #to_merge_indeces.append(inner_index)
                current_set = current_set.union(input_list[inner_index])
                eliminated_indeces.add(inner_index)
            
            inner_index += 1
            
        output.append(current_set)
        
        outer_index += 1
    
    #print(output)
    if(output == input_list):
        return output
    else:
        return merge_list(output)

def generate_gene_networks(gene_inter_file:str):
    # Convert a collection of gene pairs into a set of gene networks
    # The interacting gene pairs are stored in a text file, one gene pair per line
    
    gene_pairs = []
    gene_networks = []
    with open(gene_inter_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            gene_pairs.append({row[0], row[1]})
            
    gene_networks = merge_list(gene_pairs)       
    
    return gene_networks
                
def check_collection_overlap(collection1, collection2):
    for ele in collection1:
        if(ele in collection2):
            return True
        
    return False

def generate_gene_blocks(fname:str, distance:int):
    # Take in a plink range file and aggregate together genes/regulatory elements 
    # that are located closely together into blocks. 
    # Write out blocks, one per line, into a new file. 
    
    ranges = dict()
    pairs = []
    blocks = []
    
    with open(fname) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        for row in csv_reader:
            ranges[row[3]] = ((row[0]), int(row[1]), int(row[2]))
            
    
    for k,v in ranges.items():
        chrom1 = v[0]
        start1 = v[1]
        end1 = v[2]
        paired = False
        for k2,v2 in ranges.items():
            if(v == v2):
                continue
            
            #print(k, k2)
            
            chrom2 = v2[0]
            start2 = v2[1]
            end2 = v2[2]
            
            # Detect one element within another element
            if(chrom1 == chrom2 and ((start1 < start2 < end1) or (start1 < end2 < end1)
                                       or (start2 < start1 < end2) or (start2 < end1 < end2))):
                pairs.append((k,k2))
                paired = True
            # Detect neighboring elements
            elif(chrom1 == chrom2 and (abs(start2 - end1)<=distance or abs(start1 - end2)<=distance)):
                pairs.append((k,k2))
                paired = True
                
        if(not paired):
            blocks.append({k})
            
    #print("Pairs: ", pairs)
    
    
    used_pairs = []
    # Now conglomerate the pairs into blocks
    for pair in pairs:
        reverse_pair = (pair[1], pair[0])
        if(pair in used_pairs or reverse_pair in used_pairs):
            continue
        
        curr_block = {pair[0], pair[1]}
        for pair2 in pairs:
            if(pair2 in used_pairs):
                pass
            elif(pair == pair2):
                pass
            elif(check_collection_overlap(curr_block, pair2)):
                curr_block.add(pair2[0])
                curr_block.add(pair2[1])
                
                used_pairs.append(pair2)
                
        
        used_pairs.append(pair)
        blocks.append(curr_block)
                
    return blocks
               

def summation_func(n:int):
    # Sum (n-1) + (n-2) + ... + 2 + 1 + 0
    assert(n > 1)
    result = 0
    n = n-1
    while n > 0:
        result += n
        n -= 1
        
    return result

def generate_multicpu_files(anno_files:list, bim_file:str, total_markers:int, num_cpus:int, 
                            out_dir:str, chrom1:str, chrom2:str, pval_cutoff:float):
    interval_length = total_markers//num_cpus
    
    i_start = 0 
    i_end = total_markers
    current_i = 0
    i_intervals = []
    while current_i < total_markers:
        if( (current_i + interval_length) < total_markers):
            i_intervals.append((current_i, current_i + interval_length))
            current_i = current_i + interval_length
        else:
            i_intervals.append((current_i, total_markers))
            current_i = total_markers
            
    counter = 0
    for i_interval in i_intervals:
        i_s = str(i_interval[0])
        i_e = str(i_interval[1])
        tot_marks = str(total_markers)
        fname_out_py = out_dir + 'BiClustering_Launcher_chr' + chrom1 + '_chr' + chrom2 
        fname_out_py += '_' + str(counter) + '.py'
        with open(fname_out_py, 'w') as writer:
            writer.write('from BiClustering import *\n')
            writer.write('import sys\n')
            writer.write('if __name__ == "__main__":\n')
            writer.write('\tassert sys.version_info.major == 3\n')
            writer.write('\tassert sys.version_info.minor >= 7\n')
            writer.write("\tindir = '/gpfs/group/home/slistopad/BiClustering/'\n")
            r = 1
            for anno_file in anno_files:
                writer.write("\tmarker_file" + str(r) + " = indir + '" + str(anno_file) + "'\n")
                r += 1
            r = 1
            writer.write('\tmarker_files = [')
            while r <= len(anno_files):
                writer.write('marker_file' + str(r))
                if(r != len(anno_files)):
                    writer.write(',')
                r += 1
            writer.write(']\n')
            writer.write("\tbim_file = indir + '" + str(bim_file) + "'\n")
            writer.write("\tchrom1,chrom2 = '" + chrom1 + "','" + chrom2 + "'\n")
            writer.write('\tN,n,inter_matrix,up_tri = initialize_matrices2(bim_file,marker_files, chrom1, chrom2, ' + str(pval_cutoff) + ')\n')
            writer.write("\tout_dir = '/gpfs/group/home/slistopad/BiClustering/'\n")
            writer.write('\ti_start, i_end = '+i_s+', '+i_e+'\n')
            writer.write('\tkm_results = compute_k_m_parallel(60, inter_matrix, up_tri, i_start, i_end, N, n)\n')
            writer.write('\tcompute_interval_pval_parallel(km_results, N, n, i_start, i_end, 60, out_dir, chrom1=chrom1, chrom2=chrom2)\n')
            
            if(counter == 0):
                writer.write('\tpval_results = dict()\n')
                writer.write('\ti_intervals = [')
                counter_inner = 0
                for i_inter in i_intervals:
                    counter_inner += 1
                    writer.write('('+str(i_inter[0])+','+str(i_inter[1])+')')
                    if(counter_inner < len(i_intervals)):
                        writer.write(',')
                    else:
                        writer.write(']\n')
                writer.write('\tfor i_pair in i_intervals:\n')
                writer.write("\t\tfilename = out_dir + 'chr" + chrom1 + "_chr" + chrom2 + "_pval_results_' + str(i_pair[0]) + '_' + str(i_pair[1]) + '.csv'\n")
                writer.write('\t\twith open(filename) as csv_file:\n')
                writer.write("\t\t\tcsv_reader = csv.reader(csv_file, delimiter=',')\n")
                writer.write("\t\t\tfor row in csv_reader:\n")
                writer.write("\t\t\t\tpval_results[(int(row[0]),int(row[1]),int(row[2]),int(row[3]))] = float(row[6])\n")
                        
                writer.write("\ts_inter_pairs = sorted(pval_results.items(), key = lambda x: abs(x[1]), reverse = False)\n")            
                writer.write("\ttrimmed_inter_pairs = trim_intervals(s_inter_pairs)\n")
                writer.write("\tfilename = out_dir + 'biclustering_results_chr" + chrom1 + "_chr" + chrom2 + ".txt'\n")
                writer.write("\twith open(filename, 'w') as writer:\n")
                writer.write("\t\tfor inter_pair in trimmed_inter_pairs:\n")
                writer.write("\t\t\twriter.write(str(inter_pair[0][0])+','+str(inter_pair[0][1])+','+str(inter_pair[0][2])+','+str(inter_pair[0][3])+','+str(inter_pair[1]) + '\\n')\n")
        
        counter += 1
        
    fname_out_sh = out_dir + 'BiClustering_Launcher_chr' + chrom1 + '_chr' + chrom2 + '.sh'
    with open(fname_out_sh, 'w') as writer:
        writer.write("#!/bin/sh\n")
        writer.write("#SBATCH --job-name=SL_BiClustering_Launcher_chr" + chrom1 + '_chr' + chrom2 + "\n")
        writer.write("#SBATCH --array=1-" + str(len(i_intervals)-1) + "\n")
        writer.write("#SBATCH --ntasks=1\n")
        writer.write("#SBATCH --mem=48gb\n")
        writer.write("#SBATCH --time=240:00:00\n")
        writer.write("#SBATCH --partition=shared\n\n")
        writer.write("cd $SLURM_SUBMIT_DIR\n")
        writer.write("module load use.own\n")
        writer.write("module load python/3.8.3\n")
        writer.write("export PYTHONPATH=/gpfs/home/slistopad/.local/lib/python3.8/site-packages:$PYTHONPATH\n")
        writer.write("python3 " + 'BiClustering_Launcher_chr' + chrom1 + '_chr' + chrom2 + '_${SLURM_ARRAY_TASK_ID}.py' + '\n')
        writer.write("python3 " + 'BiClustering_Launcher_chr' + chrom1 + '_chr' + chrom2 + '_0.py' + '\n')
          
        
def get_number_of_interacting_snps_in_interval(biclustering_file:str, bim_file:str, gene_location_file:str,
                                               chrom1:str, chrom2:str, interaction_files:list, pval_cutoff:float):
    interacting_pairs_per_interval = dict()
    
    chrom1_snps, chrom2_snps = parse_plink_bim_file(bim_file, chrom1, chrom2)
    
    chrom1_ranges, chrom2_ranges = [], []
    
    with open(biclustering_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            range1 = (int(row[0]), int(row[0]) + int(row[2]))
            chrom1_ranges.append(range1)
            range2 = (int(row[1]), int(row[1]) + int(row[3]))
            chrom2_ranges.append(range2)
            
            
    interactions = parse_remma_interaction_data(interaction_files, pval_cutoff, chrom1, chrom2)
            
    #print('goodbye')
            
    print('chrom1_ranges: ', chrom1_ranges)
    print('chrom2_ranges: ', chrom2_ranges)
    
    gene_reg_locations = []

    with open(gene_location_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        for row in csv_reader:
            gene_reg_locations.append((row[0],row[1],row[2],row[3]))
    
    i = 0
    num_interval_pairs = len(chrom1_ranges)
    while i < num_interval_pairs:
        interval1 = chrom1_snps[chrom1_ranges[i][0]:chrom1_ranges[i][1]]
        interval2 = chrom2_snps[chrom2_ranges[i][0]:chrom2_ranges[i][1]]
        #print(len(interval1))
        #print(len(interval2))
        
        #print(interval1)
        #print(interval2)
        
        #print(interactions)
        
        for snp in interval1:
            for snp2 in interval2:
                if( (snp, snp2) in interactions or (snp2, snp) in interactions ):
                    
                    chrom1_loc1 = snp.split(':')
                    chrom1 = chrom1_loc1[0]
                    loc1 = chrom1_loc1[1]
                    element1 = ''
                    
                    for ele_loc in gene_reg_locations:
                        if(chrom1 == ele_loc[0] and int(loc1) >= int(ele_loc[1]) and int(loc1) <= int(ele_loc[2])):
                            element1 += ele_loc[3] + ' '
                            
                    #print(element1)
                    
                    chrom2_loc2 = snp2.split(':')
                    chrom2 = chrom2_loc2[0]
                    loc2 = chrom2_loc2[1]
                    element2 = ''
                    
                    
                    for ele_loc in gene_reg_locations:
                        if(chrom2 == ele_loc[0] and int(loc2) >= int(ele_loc[1]) and int(loc2) <= int(ele_loc[2])):
                            element2 += ele_loc[3] + ' '
                            
                    #print(element2)
                            
                    if (element1, element2) not in interacting_pairs_per_interval:
                        interacting_pairs_per_interval[(element1, element2)] = 1 
                    else:
                        interacting_pairs_per_interval[(element1, element2)] += 1
        i += 1
        print(interacting_pairs_per_interval)
                    
    return interacting_pairs_per_interval
    
def get_msigdb_enrichment(gene_pairs:set, out_dir:str):
    #print(pairs)
    pairs_map = dict()
    fname = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/Pathways/c2.all.v2023.2.Hs.symbols.gmt'
    with open(fname) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            pathway_name = row[0]
            pathway_genes = row[2:]
            for pair in gene_pairs:
                if(pair[0] in pathway_genes and pair[1] in pathway_genes):
                    try:
                        pairs_map[pair].append(pathway_name)
                    except KeyError:
                        pairs_map[pair] = [pathway_name]
                        
    fname_out = out_dir + 'biclustering_results_msigdb_enrichment.txt'
    with open(fname_out, 'w') as writer:
        for k,v in pairs_map.items():
            writer.write(str(k) + ' : ' + str(v) + '\n')
       
def parse_biclustering_annotation(bicluster_file_annot:str):
    element_pairs_dict = dict()
    element_pairs_annotated = []
    element_pairs_per_interval = []
    unique_interacting_snps = []
    interacting_snp_pairs_all_intervals = []
    regulatory_elements_all = set()
    with open(bicluster_file_annot) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            # If the row is not empty
            if(len(row) == 0):
                continue
            if(len(row[0]) > 0 and 'Interacting snp pairs:' not in row[0] and 
               'Interacting num snp pairs:' not in row[0] and 'Total number of unique' not in row[0]):
                # There must be at least one gene/reg element pair
                # It must be surrounded by (" ")
                interval_pair = row[0]
                interval_pair_split = interval_pair.split("),")
                pval = row[1]
                element_pairs_annotated2 = []
                for element_pair in interval_pair_split:
                    #print(element_pair)
                    element2_begin = element_pair.find('element2:')
                    r = 0
                    record_element = False
                    element1_set = set()
                    element1 = ''
                    while r < (element2_begin-1):
                        #print(element_pair)
                        #print(i)
                        char = element_pair[r]
                        if(char == "'"):
                            if(record_element):
                                if('GH' in element1 and len(element1) == 11):
                                    regulatory_elements_all.add(element1)
                                    element1 = ''
                                else:
                                    element1_set.add(element1)
                                    element1 = ''
                            record_element = not record_element
                        if(record_element and char!= "'"):
                            element1 += char
                        r += 1
                        
                    #print('Element1: ', element1_list)
                        
                    r = element2_begin
                    record_element = False
                    element2_set = set()
                    element2 = ''
                    while r < len(element_pair):
                        #print(interval_pair)
                        #print(r)
                        #print(element_pair_end-2)
                        #print(char)
                        char = element_pair[r]
                        #print(char)
                        if(char == "'"):
                            if(record_element):
                                if('GH' in element2 and len(element2) == 11):
                                    regulatory_elements_all.add(element2)
                                    element2 = ''
                                else:
                                    element2_set.add(element2)
                                    element2 = ''
                            record_element = not record_element
                        if(record_element and char!= "'"):
                            element2 += char
                        r += 1
                        
                    #print('Element2: ', element2_list)
                    element_pair_tuple = (element1_set, element2_set)
                    if(element_pair_tuple not in element_pairs_annotated):
                        element_pairs_annotated.append((element1_set, element2_set))
                    if(element_pair_tuple not in element_pairs_annotated2):
                        element_pairs_annotated2.append((element1_set, element2_set))
                    
                    element_pair_str = str((element1_set, element2_set))
                    element_pair_str_reverse = str((element2_set, element1_set))
                    if(element_pair_str in element_pairs_dict):
                        if(pval < element_pairs_dict[element_pair_str]):
                            element_pairs_dict[element_pair_str] = pval
                    elif(element_pair_str_reverse in element_pairs_dict):
                        if(pval < element_pairs_dict[element_pair_str_reverse]):
                            element_pairs_dict[element_pair_str_reverse] = pval
                    else:
                        element_pairs_dict[element_pair_str] = pval
                element_pairs_per_interval.append(element_pairs_annotated2)
                        
            elif('Total number of unique' in row[0]):
                split_row = row[0].split(';')[1:]
                #print(split_row)
                unique_interacting_snps_in_interval = set()
                for unique_snp in split_row:
                    if(':' in unique_snp):
                        unique_interacting_snps_in_interval.add(unique_snp)
                unique_interacting_snps.append(unique_interacting_snps_in_interval)
                
            elif('Interacting snp pairs:' in row[0]):
                current_row = row[0]
                pattern = "('\d+:\d+', '\d+:\d+')"
                matches = re.findall(pattern, current_row)
                interacting_snp_pairs_all = set()
                for snp_pair_string in matches:
                    snp_pair_string = snp_pair_string.replace("'", "")
                    match_split = snp_pair_string.split(', ')
                    #print(match_split)
                    if((match_split[1], match_split[0]) not in interacting_snp_pairs_all):
                        interacting_snp_pairs_all.add( (match_split[0], match_split[1]) )
                interacting_snp_pairs_all_intervals.append(interacting_snp_pairs_all)
                    
    return element_pairs_dict, element_pairs_annotated, unique_interacting_snps, interacting_snp_pairs_all_intervals, element_pairs_per_interval, regulatory_elements_all

def parse_biclustering_results(biclustering_file:str, bim_file:str, gene_location_file:str,
                               gene_hancer_file:str, chrom1:str, chrom2:str, interaction_files:list,
                               pval_cutoff:float, out_dir:str, single_file_output:bool = False,
                               snp_list_output:bool = False, strict_annotation:bool = True):
    # WARNING, The pval_cutoff should match the p-value cutoff used to conduct the biclustering analysis.
    # Otherwise, the mapping of the interacting intervals to interacting genes will not make sense. 
    # strict annotation = only map interacting snp pairs to genes/reg elements, 
    # loose annotation = map all snp pairs contained within interacting interval to genes/reg elements
    
    chrom1_snps, chrom2_snps = parse_plink_bim_file(bim_file, chrom1, chrom2)
    
    chrom1_ranges, chrom2_ranges = [], []
    p_values = []
    
    with open(biclustering_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            range1 = (int(row[0]), int(row[0]) + int(row[2]))
            chrom1_ranges.append(range1)
            range2 = (int(row[1]), int(row[1]) + int(row[3]))
            chrom2_ranges.append(range2)
            p_values.append(row[4])
    
            
    # modify this section depending on the gene_location_file used
    # herein I assume usage of ucsc-hg19
    # first read in and map the transcripts to the gene names
    transcript_to_gene_name = dict()
    kgXref_file = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/Gene_Annotations/ucsc-hg19-kgXref.txt'
    with open(kgXref_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            transcript_to_gene_name[row[0]] = row[1]
    
    gene_reg_locations = []
    with open(gene_location_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            gene_reg_locations.append((row[1][3:],row[2],row[3],transcript_to_gene_name[row[0]]))
            
    # We also need to account for all the regulatory elements located within the genes 
    # that were used in epistasis analysis
    
    with open(gene_hancer_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        first_line = True
        for row in csv_reader:
            if(first_line):
                first_line = False
                continue
            gene_reg_locations.append((row[0][3:],row[1],row[2],row[3]))
            
    #print("gene_reg_locations: ", gene_reg_locations)
            
    #print('chrom1_ranges: ', chrom1_ranges)
    #print('chrom2_ranges: ', chrom2_ranges)
    
    interactions = parse_remma_interaction_data(interaction_files, pval_cutoff, chrom1, chrom2)
    
    #print("Interactions: ", interactions)
    
    if(not single_file_output):
        if(not snp_list_output):
            
            if(strict_annotation):
                filename = out_dir + 'biclustering_results_chr' + chrom1 + '_chr' + chrom2 + '_annotated.txt'
            else:
                filename = out_dir + 'biclustering_results_chr' + chrom1 + '_chr' + chrom2 + '_annotated_loosely.txt'
            with open(filename, 'w') as writer:
                i = 0
                num_interval_pairs = len(chrom1_ranges)
                while i < num_interval_pairs:
                    
                    gene_pair_to_num_interacting_snp_pairs = dict()
                    gene_pair_to_interacting_snp_pairs = dict()
                    
                    gene_pairs = set()
                    interval1 = chrom1_snps[chrom1_ranges[i][0]:chrom1_ranges[i][1]]
                    interval2 = chrom2_snps[chrom2_ranges[i][0]:chrom2_ranges[i][1]]
                    p_val = p_values[i]
                    
                    #print("Interval1: ", interval1)
                    #print("Interval2: ", interval2)
                    
                    if(strict_annotation):
                        for snp in interval1:
                            for snp2 in interval2:
                                #print(snp, snp2)
                                if( (snp, snp2) in interactions or (snp2, snp) in interactions ):
                                    # Need to identify number of interacting snp pairs belonging to each gene pair
                                    
                                    #print("Found SNP pair in interactions!")
                                    chrom1_loc1 = snp.split(':')
                                    chrom1 = chrom1_loc1[0]
                                    loc1 = chrom1_loc1[1]
                                    element1 = set()
                                    
                                    #print('chrom1: ', chrom1)
                                    #print('loc1: ', loc1)
                                    
                                    element_name1 = 'intergenic'
                                    for ele_loc in gene_reg_locations:
                                        if(chrom1 == ele_loc[0] and int(loc1) >= int(ele_loc[1]) and int(loc1) <= int(ele_loc[2])):
                                            # There are a small number of overlapping genes, these are the 
                                            # poorly defined genes, such as MICB and HLA-C, for whom the exact 
                                            # loci are not yet fully established.
                                            element_name1 = ele_loc[3]
                                            element1.add(element_name1)
                                            
                                    #print('Element1: ', element1)
                                    if(len(element1) == 0):
                                        element1.add(element_name1)
                                    
                                    chrom2_loc2 = snp2.split(':')
                                    chrom2 = chrom2_loc2[0]
                                    loc2 = chrom2_loc2[1]
                                    element2 = set()
                                    
                                    #print('chrom2: ', chrom2)
                                    #print('loc2: ', loc2)
                                    
                                    element_name2 = 'intergenic'
                                    for ele_loc in gene_reg_locations:
                                        if(chrom2 == ele_loc[0] and int(loc2) >= int(ele_loc[1]) and int(loc2) <= int(ele_loc[2])):
                                            element_name2 = ele_loc[3]
                                            element2.add(element_name2)
                                            
                                    #print('Element2: ', element2)
                                    if(len(element2) == 0):
                                        element2.add(element_name2)
                                    
                                    for ele_name1 in element1:
                                        for ele_name2 in element2:
                                            if((ele_name1, ele_name2) not in gene_pair_to_num_interacting_snp_pairs and 
                                               (ele_name2, ele_name1) not in gene_pair_to_num_interacting_snp_pairs):
                                                gene_pair_to_num_interacting_snp_pairs[(ele_name1, ele_name2)] = 1
                                                gene_pair_to_interacting_snp_pairs[(ele_name1, ele_name2)] = [(snp, snp2)]
                                            elif((ele_name1, ele_name2) in gene_pair_to_num_interacting_snp_pairs):
                                                gene_pair_to_num_interacting_snp_pairs[(ele_name1, ele_name2)] += 1
                                                gene_pair_to_interacting_snp_pairs[(ele_name1, ele_name2)].append((snp, snp2))
                                            elif((ele_name2, ele_name1) in gene_pair_to_num_interacting_snp_pairs):
                                                gene_pair_to_num_interacting_snp_pairs[(ele_name2, ele_name1)] += 1
                                                gene_pair_to_interacting_snp_pairs[(ele_name2, ele_name1)].append((snp2, snp))
                                      
                                    #print(gene_pair_to_num_interacting_snp_pairs)
                                    gene_pairs.add( ('element1: '+ str(element1), 'element2: ' + str(element2)) )
                        writer.write(str(gene_pairs) + '\t')
                        writer.write(p_values[i] + '\n')
                        writer.write('Interacting num snp pairs: ' + str(gene_pair_to_num_interacting_snp_pairs) + '\n')
                        writer.write('Interacting snp pairs: ' + str(gene_pair_to_interacting_snp_pairs) + '\n')
                        current_union = set()
                        for k,v in gene_pair_to_interacting_snp_pairs.items():
                            current_union = current_union | set(v)
                        writer.write('Total number of unique interacting pairs: ' + str(len(current_union)) + ';')
                        unique_interacting_snps = set()
                        for unique_snp_pair in current_union:
                            unique_interacting_snps.add(unique_snp_pair[0])
                            unique_interacting_snps.add(unique_snp_pair[1])
                        for unique_snp in unique_interacting_snps:
                            writer.write(str(unique_snp) + ';')
                        writer.write('\n\n')
                    else:
                        elements1 = set()
                        for snp in interval1:
                            chrom1_loc1 = snp.split(':')
                            chrom1 = chrom1_loc1[0]
                            loc1 = chrom1_loc1[1]
                            element1 = set()
                            
                            for ele_loc in gene_reg_locations:
                                if(chrom1 == ele_loc[0] and int(loc1) >= int(ele_loc[1]) and int(loc1) <= int(ele_loc[2])):
                                    element1.add(ele_loc[3])
                                    #print(element1)
                                    
                            elements1.add(str(element1))
                        elements2 = set()
                        for snp2 in interval2:
                            chrom2_loc2 = snp2.split(':')
                            chrom2 = chrom2_loc2[0]
                            loc2 = chrom2_loc2[1]
                            element2 = set()
                            #print(loc2)
                             
                            for ele_loc in gene_reg_locations:
                                if(chrom2 == ele_loc[0] and int(loc2) >= int(ele_loc[1]) and int(loc2) <= int(ele_loc[2])):
                                    element2.add(ele_loc[3])
                            elements2.add(str(element2))
                                 
                        writer.write(str( (elements1,elements2) ) + '\t')
                        writer.write(p_values[i] + '\n')
                    
                    i += 1
        else:
            #print("HELLO!")
            filename = out_dir + 'biclustering_results_chr' + chrom1 + '_chr' + chrom2 + '_snplist_interacting_pairs.txt'
            with open(filename, 'w') as writer:
                num_interval_pairs = len(chrom1_ranges)
                i = 0
                while i < num_interval_pairs:
                    gene_pairs = set()
                    interval1 = chrom1_snps[chrom1_ranges[i][0]:chrom1_ranges[i][1]]
                    interval2 = chrom2_snps[chrom2_ranges[i][0]:chrom2_ranges[i][1]]
                    
                    # format for liftover tool
                    snp1_info = interval1[0].split(':')
                    snp1_chrom,snp1_loci = snp1_info[0],snp1_info[1]
                    snp2_info = interval1[-1].split(':')
                    snp2_chrom,snp2_loci = snp2_info[0],snp2_info[1]
                    
                    #writer.write(snp1_chrom + ' ' + snp1_loci + ' ' + snp2_loci + '\n')
                    
                    snp1_info = interval2[0].split(':')
                    snp1_chrom,snp1_loci = snp1_info[0],snp1_info[1]
                    snp2_info = interval2[-1].split(':')
                    snp2_chrom,snp2_loci = snp2_info[0],snp2_info[1]
                    
                    #writer.write(snp1_chrom + ' ' + snp1_loci + ' ' + snp2_loci + '\n')
                    
                    for snp in interval1:
                        for snp2 in interval2:
                            #print(snp, snp2)
                            if( (snp, snp2) in interactions or (snp2, snp) in interactions ):
                                writer.write(snp + ' : ')
                                writer.write(snp2 + '\n')
                                
                                snp1_info = snp.split(':')
                                snp1_chrom,snp1_loci = snp1_info[0],snp1_info[1]
                                
                                snp2_info = snp2.split(':')
                                snp2_chrom,snp2_loci = snp2_info[0],snp2_info[1]
                                
                                print("plink --bfile /gpfs/group/home/slistopad/REMMA/data/native_american/", end='')
                                print('NA3_Combined_Fixed_Set_score900_db001_and_exp001_1_and_genehancer', end='')
                                print('_score25_maf001_region_qc3 --ld ' + snp + ' ' + snp2 + ' --r2 ', end='')
                                print('--out /gpfs/group/home/slistopad/REMMA/data/native_american/Biclustering_SNPLists/', end='')
                                print('Combined_Fixed2_Set_Biclustering/biclustering_results_ld_' + snp1_chrom, end='')
                                print('_' + snp1_loci + '_' + snp2_chrom + '_' + snp2_loci)
                     
                    i += 1
                    
    else:
        # Note, what is put into the file, are the locations of the genes, which contain 
        # interacting SNPs within the interacting intervals. 
        
        # This file is made for generation of Circo plots. In it we don't want to see the same gene pairs even if they 
        # belong to different interval pairs. 
        filename = out_dir + 'biclustering_results_all.txt'
        with open(filename, 'a') as writer:
            i = 0
            tracker = set()
            num_interval_pairs = len(chrom1_ranges)
            while i < num_interval_pairs:
                gene_pairs = set()
                interval1 = chrom1_snps[chrom1_ranges[i][0]:chrom1_ranges[i][1]]
                interval2 = chrom2_snps[chrom2_ranges[i][0]:chrom2_ranges[i][1]]
                p_val = p_values[i]
                
                #print('interval1: ', interval1)
                #print('interval2: ', interval2)
                
                for snp in interval1:
                    for snp2 in interval2:
                        if( (snp, snp2) in interactions or (snp2, snp) in interactions ):
                
                            chrom1_loc1 = snp.split(':')
                            chrom1 = chrom1_loc1[0]
                            loc1 = chrom1_loc1[1]
                            element1 = ''
                            
                            for ele_loc in gene_reg_locations:
                                if(chrom1 == ele_loc[0] and int(loc1) >= int(ele_loc[1]) and int(loc1) <= int(ele_loc[2])):
                                    element1 = 'chr' + ele_loc[0] + ',' + ele_loc[1] + ',' + ele_loc[2]
                            
                            chrom2_loc2 = snp2.split(':')
                            chrom2 = chrom2_loc2[0]
                            loc2 = chrom2_loc2[1]
                            element2 = ''
                            
                            for ele_loc in gene_reg_locations:
                                if(chrom2 == ele_loc[0] and int(loc2) >= int(ele_loc[1]) and int(loc2) <= int(ele_loc[2])):
                                    element2 = 'chr' + ele_loc[0] + ',' + ele_loc[1] + ',' + ele_loc[2]
                            
                            gene_pairs.add((element1, element2))
                
                #print(gene_pairs)
                #print(p_val)
                for pair in gene_pairs:
                    if(pair not in tracker):
                        writer.write(pair[0] + ',' + pair[1] + ',' + str(p_val) + '\n')
                    tracker.add(pair)
                i += 1
    
    
def parse_remma_interaction_data(interaction_files:list, pval_cutoff:float, chrom1:str, chrom2:str):
    #print(pval_cutoff)
    #print(chrom1)
    #print(chrom2)
    
    interactions = dict()
    for interaction_file in interaction_files:
        with open(interaction_file) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ')
            header = True
            for row in csv_reader:
                if(header):
                    header = False
                    continue
                # Check if entry passes p-value cutoff and matches target chromosomes.
                if(float(row[18]) < pval_cutoff and ((row[1] == chrom1 and row[8] == chrom2) or 
                   (row[8] == chrom1 and row[1] == chrom2))):
                    
                    # If this entry already is recorded in interactions (presumably from another 
                    # interaction file), then overwrite the p-value for the record only if the 
                    # new p-value is smaller.
                    if((row[2], row[9]) in interactions):
                        if(float(row[18]) < interactions[(row[2], row[9])]):
                            interactions[(row[2], row[9])] = float(row[18])
                    elif((row[9], row[2]) in interactions):
                        if(float(row[18]) < interactions[(row[9], row[2])]):
                            interactions[(row[9], row[2])] = float(row[18])
                    else:
                        interactions[(row[2], row[9])] = float(row[18])
    return interactions

def parse_plink_bim_file(bim_file:str, chrom1:str, chrom2:str):
    chrom1_snps = []
    chrom2_snps = []
    
    with open(bim_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            if(chrom1 == row[0]):
                chrom1_snps.append(row[1])
            if(chrom2 == row[0]):
                chrom2_snps.append(row[1])
                
    return chrom1_snps, chrom2_snps
    
  
def initialize_matrices2(bim_file:str, interaction_files:list, chrom1:str, chrom2:str, 
                         pval_cutoff:float):
    # Assume that testing is always done using all relevant SNPs on one chromosome vs 
    # all relevant SNPs on another chromosome. A chromosome can be also tested against itself. 
    # The matrix size is number of SNPs on chromosome 1 vs number of SNPs on chromosome 2.
    # Relevant SNPs means SNPs that were included in the epistasis analysis.
    
    # Bim file informs us the number of relevant SNPs on each chromosome, and which SNPs pair 
    # a given (i,j) value in the matrix represents. 
    # Interaction file informs us which SNP pairs were found to interact.
    # Generally it is expected that the number of interactions is small. 
    # Whenever the chromosomes are different, all SNP pairs are assumed to have been tested 
    # for interaction.
    # If chromsome is the same, the tested interaction matrix is an upper triangular matrix. 
    
    # Important assumption regarding bim file, for each chromosome all SNPs are ordered based 
    # on their location, from beginning of chromosome to end of chromosome. 
    
    # The interaction file is assumed to be a REMMA .anno file. 
    
    chrom1_snps, chrom2_snps = parse_plink_bim_file(bim_file, chrom1, chrom2)
                
    interactions = parse_remma_interaction_data(interaction_files, pval_cutoff, chrom1, chrom2)
        
    #print(interactions)    
    
    upper_triangular = False
    if(chrom1 == chrom2):
        upper_triangular = True
        
    N = 0
    n = 0

    total_markers1 = len(chrom1_snps)
    total_markers2 = len(chrom2_snps)
    
    if(upper_triangular):
        i = 0
        while i < total_markers1:
            N += i
            i += 1
    else:
        N = total_markers1 * total_markers2
    
    inter_array = np.zeros((total_markers1, total_markers2), dtype=int)
    i = 0
    while i < total_markers1:
        j = 0
        while j < total_markers2:
            if((chrom1_snps[i], chrom2_snps[j]) in interactions or
               (chrom2_snps[j], chrom1_snps[i]) in interactions):
                if(upper_triangular):
                    if(i >= j):
                        j += 1
                        continue
                #print(i,j)
                inter_array[i][j] = 1
                n += 1
                
            j += 1
            
        i += 1
        
    print("Interaction matarix initialized.")
    print("N: ", N)
    print("n: ", n)
    print("Chromosome1: ", chrom1)
    print("Chromosome2: ", chrom2)
    print("P-Value Cutoff: ", pval_cutoff)
    
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)
                
    return N, n, inter_array, upper_triangular
    
  
def initialize_matrices(total_markers:int, marker_file:str, upper_triangular:bool = True):
    # Also, we assume that all values in the interaction matrices are 1s and 0s.
    
    # We assume that all interaction and tested-interaction matrices are square.
    
    # We assume that all tested interaction matrices are either upper triangular or full (of 1s).
    # Upper triangular tested interaction matrix represents a set of marker being tested against 
    # itself for epistasis.
    # Full tested interaction matrix represents two completely different sets of markers being
    # tested against each other for epistasis. 
    
    N = 0
    n = 0
    inter_matrix = dict()
    
    if(upper_triangular):
        i = 0
        while i < total_markers:
            N += i
            i += 1
        #print("N = ", N)
    else:
        N = total_markers**2

    i = 0
    j = 0
    while i < total_markers:
        while j < total_markers:
            inter_matrix[(i,j)] = 0
            j += 1
              
        j = 0
        i += 1
    
    with open(marker_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            i,j = int(row[0]), int(row[1])
            inter_matrix[(i,j)] = 1
            n += 1
    
    return N, n, inter_matrix, upper_triangular

def check_overlap(i:int, j:int, a:int, b:int, i2:int, j2:int, a2:int, b2:int):
    # Note if two matrices are identical we return false. 
    if(i == i2 and j == j2 and a==a2 and b==b2):
        return False
    
    assert(i >= 0 and j >= 0 and i2 >= 0 and j2 >= 0)
    assert(a >= 1 and b>= 1 and a2 >= 1 and b2 >= 1)
    
    origin_test = ( ((i+a-1) >= i2 >= i) and ((j+b-1) >= j2 >= j) ) or ( ((i2+a2-1) >= i >= i2) and ((j2+b2-1) >= j >= j2) )
    if(origin_test):
        return True
    bottom_left_test =  ( (i2 <= (i+a-1) <= (i2+a2-1)) and ((j2+b2-1) >= j >= j2) ) or ( (i <= (i2+a2-1) <= (i+a -1)) and ((j+b-1) >= j2 >= j) )
    if(bottom_left_test):
        return True
    bottom_right_test = ( (i2 <= (i+a-1) <= (i2+a2-1)) and ((j2+b2-1) >= (j+b-1) >= j2) ) or (  (i <= (i2+a2-1) <= (i+a-1)) and ((j+b-1) >= (j2+b2-1) >= j)  )
    if(bottom_right_test):
        return True
    upper_right_test = ( ((i+a-1) >= i2 >= i) and ((j2+b2-1) >= (j+b-1) >= j2) ) or ( ((i2+a2-1) >= i >= i2) and ((j+b-1) >= (j2+b2-1) >= j) )
    if(upper_right_test):
        return True
    
    return False

'''
def comp(i:int, j:int, a:int, b:int, val:str, inter_matrix:dict, tested_inter_matrix:dict):
    #print("i =",i,";j = ",j,";a = ", a, ";b = ", b)
    if(i < 0):
        raise ValueError("Invalid i value.")
    elif(j < 0):
        raise ValueError("Invalid j value.")
    elif(a < 0):
        raise ValueError("Invalid a value.")
    elif(b < 0):
        raise ValueError("Invalid b value.")
        
    if(a == 0 or b == 0):
        return 0
    elif(a == 1 and b == 1):
        if(val == 'k'):
            return inter_matrix[(i,j)]
        elif(val == 'm'):
            return tested_inter_matrix[(i,j)]
        else:
            raise ValueError("This parameter must be either 'k' or 'm'.")
    else:
        return comp(i,j,a-1,b-1,val,inter_matrix,tested_inter_matrix) + comp(i+a-1,j,1,b-1,val,inter_matrix,tested_inter_matrix) + comp(i,j+b-1,a-1,1,val,inter_matrix,tested_inter_matrix) + comp(i+a-1,j+b-1,1,1,val,inter_matrix,tested_inter_matrix)
'''

def convert_dict_to_2D_numpy_array(input_dict:dict):
    # This function converts a dictionary that maps (row_index, column_index) tuple to integer values.
    # It does so using the correct ordering.
    s_matrix_values = sorted(input_dict.items(), key = lambda x: x[0], reverse = False)
    TwoDList = []
    curr_inner = []
    prev_row = 0
    for ele in s_matrix_values:
        if(ele[0][0] > prev_row):
            TwoDList.append(curr_inner)
            curr_inner = [ele[1]]
            prev_row = ele[0][0]
        else:
            curr_inner.append(ele[1])
            
    return np.array(TwoDList)

@jit(nopython=True)
def comp_alt(i:int, j:int, a:int, b:int, inter_matrix:'np.array', upper_triangular:bool = True,
             compute_k:bool = True):
    # Assume that tested interaction matrix is either upper triangular or full of 1s.
    # In that case m_val can be computed quickly. Assume that one of the two is always true.
    #print("i =",i,";j = ",j,";a = ", a, ";b = ", b)
    if(i < 0):
        raise ValueError("Invalid i value.")
    elif(j < 0):
        raise ValueError("Invalid j value.")
    elif(a < 0):
        raise ValueError("Invalid a value.")
    elif(b < 0):
        raise ValueError("Invalid b value.")
        
    k_val = 0
    m_val = 0
    
    
    if(compute_k):
        if(not upper_triangular):
            m_val = a * b
        else:
            i2 = i   
            while i2 < (i+a):
                # At which index do 1s start in this row. That is simpy row_index + 1.
                if(j < (i2 + 1)):
                    m_val += max(j + b - (i2 + 1), 0)
                else:
                    m_val += b
                i2 += 1
        
        k_val = np.count_nonzero(inter_matrix[i:(i+a),j:(j+b)])
    
    '''
    i2 = i   
    while i2 < (i+a):
        j2 = j
        while j2 < (j+b):
            k_val += inter_matrix[i2,j2]
            j2 += 1
        i2 += 1
    '''
    
    return k_val, m_val

@jit(nopython=True)
def compute_k_m_parallel(max_inter_length:int, inter_array,
                         upper_triangular:bool, i_start:int, i_end:int, N:int, n:int):
    
    # Our only knowledge about inter_matrix is that it consists of 1s and 0s 
    # and that 1s are typically very sparce.
    
    # We have stronger assumptions about tested_inter_matrix which we deem to be either 
    # upper trinagular (with 1s and 0s) or completely filled with 1s. 
    # Given this we do not need to actually parse tested_inter_matrix to compute m.
    
    # Convert inter_matrix and tested_inter_matrix into numpy arrays.
    dimensions = inter_array.shape
    total_markers1 = dimensions[0]
    total_markers2 = dimensions[1]
    
    i = i_start
    j = 0
    a = max_inter_length
    b = max_inter_length
    
    km_results = dict()
    
    # filename = out_dir + 'km_results_' + str(i_start) + '_' + str(i_end) + '.csv'
    # with open(filename, 'w') as writer:
        # writer.write('#i,j,a,b,k,m')
    
    while i < i_end:
        #skipped_comp_k = 0
        #computed_k = 0
        while j < total_markers2:
            # If we the biggest matrix at (i,j) has k = 0, then all of the matrices within it 
            # will also have k = 0. 
            k_is_zero = False
            # a_thresh and b_thresh represent biggest submatrix that is all zeros
            a_thresh = 0
            b_thresh = 0
            while a >= 1:
                while b >= 1:
                    if((i + a) <= total_markers1 and (j+b) <= total_markers2):
                        if(k_is_zero and a <= a_thresh and b <= b_thresh):
                            #skipped_comp_k += 1
                            k,m = comp_alt(i,j,a,b, inter_array, upper_triangular, compute_k=False)
                        else:
                            k,m = comp_alt(i,j,a,b, inter_array, upper_triangular)
                            #computed_k += 1
                            if(k == 0):
                                k_is_zero = True
                                a_thresh = a
                                b_thresh = b
                        # Do not store the k,m pair if k <= expected value (m*n/N)
                        if(k > (m*n/N)):
                            km_results[(i,j,a,b)] = (k,m)
                        else:
                            pass
                        
                        # writer.write(str(i)+','+str(j)+','+str(a)+','+str(b)+','+str(k)+','+str(m))
                        
                    b -= 1
                b = max_inter_length
                a -= 1
            a = max_inter_length
            b = max_inter_length
            j += 1
            #print('j = ', j)
        #a = max_inter_length
        #b = max_inter_length
        j = 0
        i += 1
        print('i = ', i)
        #print("Number of skipped k computations: ", skipped_comp_k)
        #print("Number of completed k computations: ", computed_k)
    
    print("Length of km_results: ", len(km_results))
    #now = datetime.now()
    #current_time = now.strftime("%H:%M:%S")
    #print("Current Time =", current_time)
    return km_results

#@jit(nopython=True)
def compute_interval_pval_parallel(km_results:dict, N:int, n:int, i_start:int, i_end:int,
                                   max_inter_length:int, out_dir:str, chrom1:str = '-1',
                                   chrom2:str = '-1'):
    # Assume that all k values are bigger than expected k-value of (m*n/N)
    
    computed_pvals = dict()
    # Conservative Bonferonni Correction
    # In practice the number of tests is much smaller
    correction = 0.05 / (N * (max_inter_length**2))
    # this would be correct if we only computed 
    # one chromosome pair, but unfortunately we need to account for all pairs
    
    # hard coded correction:
    #correction = 0.05 / (60055879978 * 3600)
    
    
    if(chrom1 == '-1' or chrom2 == '-1'):
        filename = out_dir + 'pval_results_' + str(i_start) + '_' + str(i_end) + '.csv'
    else:
        filename = out_dir + 'chr' + chrom1 + '_chr' + chrom2 + '_pval_results_'
        filename +=  str(i_start) + '_' + str(i_end) + '.csv'
    with open(filename, 'w') as writer:
    
        pval_results = dict()
        # i = 0
        for key,val in km_results.items():
            k,m = val[0], val[1]
            max_k = min(m, n)
            # min_k = max(0, m-(N-n))
            pval_results[key] = 0
                
            # First check if a p-value for the following k,m pair has already been computed.
            try:
                pval_results[key] = computed_pvals[(val[0],m)]
            except KeyError:
                # If it has not been, then compute the p-value and store it.
                # if(k > (m*n/N)):
                
                pval_results[key] = hypergeom.sf(k-1,N,n,m)
                # while k <= max_k:
                #     pval_results[key] += hypergeom.pmf(k,N,n,m)
                #     k += 1
                
                if(pval_results[key] < correction):
                    #print(pval_results[key])
                    writer.write(str(key[0])+','+str(key[1])+','+str(key[2])+','+str(key[3])+',')
                    writer.write(str(val[0]) + ',' + str(m) + ',' + str(pval_results[key]) + '\n')
                # We do not care about interval pairs that have unusually low 
                # number of interacting pairs
                #else:
                    # while k >= min_k:
                    #     pval_results[key] += hypergeom.pmf(k, N, n, m)
                    #     k -= 1
                    
                computed_pvals[(val[0],m)] = pval_results[key]
                
            # i += 1
            # if(i%1000 == 0):
            #     print(i)
            #     print(len(computed_pvals))
    print("Computed p-values.")
    # now = datetime.now()
    # current_time = now.strftime("%H:%M:%S")
    # print("Current Time =", current_time)
    return pval_results

def trim_intervals(sorted_intervals:list):
    # Now go through interval pairs starting from most significant, and remove all overlapping 
    # intervals for the current interval
    outer_count = 0
    indeces_to_delete = set()
    for tup in sorted_intervals:
        # Skip this interval pair if it is already marked for delition
        if(outer_count in indeces_to_delete):
            outer_count += 1
            #print('SKIPPED!')
            continue
        inner_count = 0
        curr_tup = tup
        i,j,a,b = curr_tup[0][0],curr_tup[0][1],curr_tup[0][2],curr_tup[0][3]
        pval_outer = curr_tup[1]
        for tup in sorted_intervals:
            if(inner_count in indeces_to_delete):
                inner_count += 1
                #print('SKIPPED!')
                continue
            i2,j2,a2,b2 = tup[0][0],tup[0][1],tup[0][2],tup[0][3]
            pval_inner = tup[1]
            #print(i,j,a,b)
            #print(i2,j2,a2,b2)
            if(check_overlap(i,j,a,b,i2,j2,a2,b2)):
                if(pval_inner < pval_outer):
                    indeces_to_delete.add(outer_count)
                else:
                    indeces_to_delete.add(inner_count)
                    #print("bye!")
                    
            inner_count += 1
            
        outer_count += 1
        
    #print(indeces_to_delete)
            
    culled_inter_pairs = []
    i = 0
    while i < len(sorted_intervals):
        if(i not in indeces_to_delete):
            culled_inter_pairs.append(sorted_intervals[i])
        i += 1
        
    return culled_inter_pairs

class TestBiClusterCodeBase(unittest.TestCase):  
    
    def test_get_msigdb_enrichment(self):
        gene_pairs = {('ABAT', 'ACOX1'), ('CSMD1', 'DLGAP1')}
        out_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/Test/'
        get_msigdb_enrichment(gene_pairs, out_dir)
        
        fname = out_dir + 'biclustering_results_msigdb_enrichment.txt'
        with open(fname) as msig_fname:
            line = msig_fname.readline()
        
        exp_line = "('ABAT', 'ACOX1') : ['AFFAR_YY1_TARGETS_UP', "
        exp_line += "'CARRILLOREIXACH_HEPATOBLASTOMA_VS_NORMAL_DN', 'FLECHNER_BIOPSY_KIDNEY_"
        exp_line += "TRANSPLANT_REJECTED_VS_OK_DN', 'LEE_LIVER_CANCER_MYC_TGFA_DN', 'RODRIGUES_"
        exp_line += "DCC_TARGETS_DN', 'RODRIGUES_THYROID_CARCINOMA_ANAPLASTIC_DN', "
        exp_line += "'RODRIGUES_THYROID_CARCINOMA_POORLY_DIFFERENTIATED_DN']\n"
        self.assertEqual(line, exp_line)
    
    def test_parse_remma_interaction_data(self):
        indir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/Test/'
        ad_anno = indir + 'epiAD_ad_NA3_Combined_Strict_Set_Brief.anno'
        dd_anno = indir + 'epiDD_dd_NA3_Combined_Strict_Set_Brief.anno'
        interactions = parse_remma_interaction_data([ad_anno, dd_anno], 1e-6, '1', '20')
        
        self.assertEqual(len(interactions), 7)
        self.assertEqual(interactions[('1:7863293', '20:8515243')], 3.0037315721644386e-07)
        self.assertEqual(interactions[('1:7863293', '20:8515956')], 3.0037315721644386e-07)
        
        self.assertEqual(interactions[('1:7865063', '20:33467717')], 8.801979806253102e-07)
        self.assertEqual(interactions[('1:7865063', '20:33488013')], 7.014568812762179e-07)
        self.assertEqual(interactions[('1:7865063', '20:33514465')], 7.070270435022684e-07)
        self.assertEqual(interactions[('1:7865691', '20:33488013')], 9.85547891823227e-09)
        self.assertEqual(interactions[('1:1725760', '20:9588420')], 7.018367391361476e-07)
        
    
    def test_get_number_of_interacting_snps_in_interval(self):
        indir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/Test/'
        biclustering_file = indir + 'biclustering_results_chr1_chr5.txt'
        bim_file = indir + 'NA3_Combined_Strict_Set_score900_db001_or_exp001_1_and_genehancer_score10_maf001_region_qc3.bim'
        gene_location_file = indir + 'Combined_Strict_Set_score900_db001_or_exp001_1_and_genehancer_score10_ranges.txt'
        marker_files = [indir + 'epiDD_dd_NA3_Combined_Strict_Set.anno']
        inter_pairs_per_interval  = get_number_of_interacting_snps_in_interval(biclustering_file, bim_file, gene_location_file,
                                                                               '1', '5', marker_files, 1e-7)
    
        self.assertEqual(inter_pairs_per_interval[('LHX4 ', 'PPP2R2B ')], 12)
        self.assertEqual(inter_pairs_per_interval[('DISC1 ', 'SLIT3 ')], 5)
        
    def test_parse_biclustering_results(self):
        indir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/REMMA_Results/'
        marker_file1 = 'epiAA_aa1_approx_parallel_merged_NA3_Combined_Strict_Set_score900_db001_or_exp001_1_and_genehancer_score10_maf001_region_qc3.anno'
        marker_file2 = 'epiAA_aa2_approx_parallel_merged_NA3_Combined_Strict_Set_score900_db001_or_exp001_1_and_genehancer_score10_maf001_region_qc3.anno'
        marker_file3 = 'epiAA_aa3_approx_parallel_merged_NA3_Combined_Strict_Set_score900_db001_or_exp001_1_and_genehancer_score10_maf001_region_qc3.anno'
        marker_file4 = 'epiAD_ad_approx_parallel_merged_NA3_Combined_Strict_Set_score900_db001_or_exp001_1_and_genehancer_score10_maf001_region_qc3.anno'
        marker_file5 = 'epiDD_dd_approx_parallel_merged_NA3_Combined_Strict_Set_score900_db001_or_exp001_1_and_genehancer_score10_maf001_region_qc3.anno'
        
        marker_files = [indir + marker_file1, indir + marker_file2, indir + marker_file3, indir + marker_file4, indir + marker_file5]
        
        indir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/Test/'
        biclustering_file = indir + 'biclustering_results_chr1_chr3.txt'
        bim_file = indir + 'NA3_Combined_Strict_Set_score900_db001_or_exp001_1_and_genehancer_score10_maf001_region_qc3.bim'
        gene_location_file = indir + 'Combined_Strict_Set_score900_db001_or_exp001_1_and_genehancer_score10_ranges.txt'
        out_dir = indir

        parse_biclustering_results(biclustering_file, bim_file, gene_location_file, '1', '3', marker_files, 1e-7, out_dir, False) 
        filename = indir + 'biclustering_results_chr1_chr3_annotated.txt'
        rows = []
        with open(filename) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter='\t')
            for row in csv_reader:
                rows.append(row)
                
        self.assertEqual(rows[0], ["{('DTL ', 'SRGAP3 ')}", '9.605131888649506e-62'])
        self.assertEqual(rows[1], ["{('DTL ', 'SRGAP3 ')}", '6.4212833992207265e-27'])
        self.assertEqual(rows[2], ["{('SRGAP2 ', 'ROBO2 ')}", '4.0199636843080785e-19'])
        
        parse_biclustering_results(biclustering_file, bim_file, gene_location_file, '1', '3', marker_files, 1e-7, out_dir, True) 
    
        biclustering_file = indir + 'biclustering_results_chr4_chr10.txt'
        parse_biclustering_results(biclustering_file, bim_file, gene_location_file, '4', '10', marker_files, 1e-7, out_dir, True)
        biclustering_file = indir + 'biclustering_results_chr4_chr5.txt'
        parse_biclustering_results(biclustering_file, bim_file, gene_location_file, '4', '5', marker_files, 1e-7, out_dir, True)
        biclustering_file = indir + 'biclustering_results_chr8_chr17.txt'
        parse_biclustering_results(biclustering_file, bim_file, gene_location_file, '8', '17', marker_files, 1e-7, out_dir, True)
        
        filename = indir + 'biclustering_results_all.txt'
        rows = []
        with open(filename) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter='\t')
            for row in csv_reader:
                rows.append(row)
                
        self.assertEqual(rows[0], ['chr1,212206918,212280187,chr3,9020275,9293369,9.605131888649506e-62'])
        self.assertEqual(rows[1], ['chr1,206514199,206639783,chr3,75984644,77701114,4.0199636843080785e-19'])
        self.assertEqual(rows[2], ['chr4,6320304,6567327,chr10,52748910,54060110,2.993577604860245e-295'])
        self.assertEqual(rows[3], ['chr4,88080213,88143674,chr5,145967066,146463083,5.967210373347171e-208'])
        self.assertEqual(rows[4], ['chr4,20253234,20622788,chr5,15498304,15941900,4.7347537763144387e-20'])
        
        self.assertTrue(['chr8,8957214,8964665,chr17,28519336,28564986,1.1381827143817773e-19'] in rows)
        self.assertTrue(['chr8,8957214,8964665,chr17,28443018,28446019,1.1381827143817773e-19'] in rows)
        
        os.remove(indir + 'biclustering_results_all.txt')
    
    def test_generate_gene_blocks(self):
        fname = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/Biclustering/Test_Input/'
        fname += 'Test_Ranges.txt'
        blocks = generate_gene_blocks(fname, 3000)
        expected_blocks = [{'TEST3'}, {'TEST4'}, {'TEST5'}, {'TEST6'}, {'TEST11'}, {'TEST12'},
                           {'TEST22'}, {'TEST2', 'TEST1'}, {'TEST9', 'TEST10', 'TEST8', 'TEST7'},
                           {'TEST13', 'TEST14', 'TEST16', 'TEST18'},
                           {'TEST15', 'TEST19', 'TEST21', 'TEST20', 'TEST17'}]
        
        self.assertTrue(len(blocks) == len(expected_blocks))
        for block in expected_blocks:
            self.assertTrue(block in blocks)
            
        blocks = generate_gene_blocks(fname, 200)
        expected_blocks = [{'TEST1'}, {'TEST2'}, {'TEST3'}, {'TEST4'}, {'TEST5'},
                           {'TEST6'}, {'TEST8'}, {'TEST9'}, {'TEST11'}, {'TEST12'},
                           {'TEST16'}, {'TEST18'}, {'TEST20'}, {'TEST21'}, {'TEST22'},
                           {'TEST10', 'TEST7'}, {'TEST13', 'TEST14'},
                           {'TEST19', 'TEST15', 'TEST17'}]
        
        self.assertTrue(len(blocks) == len(expected_blocks))
        for block in expected_blocks:
            self.assertTrue(block in blocks)
            
        blocks = generate_gene_blocks(fname, 1000)
        expected_blocks = [{'TEST1'}, {'TEST2'}, {'TEST3'}, {'TEST4'}, {'TEST5'}, {'TEST6'},
                           {'TEST11'}, {'TEST12'}, {'TEST22'}, 
                           {'TEST9', 'TEST10', 'TEST8', 'TEST7'},
                           {'TEST13', 'TEST14', 'TEST18', 'TEST16'},
                           {'TEST15', 'TEST19', 'TEST21', 'TEST20', 'TEST17'}]
        
        self.assertTrue(len(blocks) == len(expected_blocks))
        for block in expected_blocks:
            self.assertTrue(block in blocks)
        
    
    def test_check_overlap(self):
        # Simple Cases
        self.assertTrue(check_overlap(166,621,10,40,167,620,9,42))
        self.assertTrue(check_overlap(0,491,36,11,0,491,35,11))
        self.assertTrue(check_overlap(10,10,5,5,5,12,8,5))
        self.assertTrue(check_overlap(0,0,2,2,1,1,2,2))
        self.assertTrue(check_overlap(2,2,3,6,1,1,3,3))
        self.assertFalse(check_overlap(0,0,2,2,3,2,1,1))
        self.assertFalse(check_overlap(0,0,10,10,40,40,2,2))
        
        # Challenging/Border Cases
        self.assertFalse(check_overlap(0,0,2,2,0,0,2,2))
        self.assertFalse(check_overlap(0,0,2,2,1,2,1,1))
        self.assertFalse(check_overlap(0,0,2,2,2,2,1,1))
        self.assertFalse(check_overlap(2,2,3,6,1,1,1,3))
    
    def test_initialize_matrices(self):
        total_markers = 10
        fname = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/Biclustering/Test_Input/TenByTen.txt'
        N, n, inter_matrix, upper_tri = initialize_matrices(total_markers, fname)
        
        inter_matrix_exp = dict()
        i, j = 0, 0
        while i < 10:
            j = 0
            while j < 10:
                inter_matrix_exp[(i,j)] = 0
                j += 1
            i += 1
            
        inter_matrix_exp[(2,8)] = 1
        inter_matrix_exp[(5,8)] = 1
        inter_matrix_exp[(5,9)] = 1
        
        self.assertEqual(N, 45)
        self.assertEqual(n, 3)
        self.assertTrue(upper_tri)
        
        i, j = 0, 0
        while i < 10:
            j = 0
            while j < 10:
                self.assertEqual(inter_matrix[(i,j)], inter_matrix_exp[(i,j)])
                j += 1
            i += 1
        
        N, n, inter_matrix, upper_tri = initialize_matrices(total_markers, fname, False)
            
        self.assertEqual(N, 100)
        self.assertEqual(n, 3)
        self.assertFalse(upper_tri)
        
        i, j = 0, 0
        while i < 10:
            j = 0
            while j < 10:
                self.assertEqual(inter_matrix[(i,j)], inter_matrix_exp[(i,j)])
                j += 1
            i += 1  
            
    def test_initialize_matrices2(self):
        indir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/Biclustering/Test_Input/'
        bim_file = indir + 'Test_Data.bim'
        inter_file = indir + 'Test_Interaction.anno'
        chrom1 = '1'
        chrom2 = '1'
        
        N, n, inter_array, upper_tri = initialize_matrices2(bim_file, [inter_file], chrom1, chrom2, 1e-5)
        self.assertEqual(N, 435)
        self.assertEqual(n, 1)
        self.assertTrue(upper_tri)
        
        expected_inter_array = np.zeros((inter_array.shape[0],inter_array.shape[1]))
        expected_inter_array[12][19] = 1
        
        i = 0 
        while i < inter_array.shape[0]:
            j = 0
            while j < inter_array.shape[1]:
                self.assertEqual(inter_array[i][j], expected_inter_array[i][j])
                j += 1
            i += 1
                
        indir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/Biclustering/Test_Input/'
        bim_file = indir + 'Test_Data.bim'
        inter_file = indir + 'Test_Interaction.anno'
        chrom1 = '2'
        chrom2 = '3'
        
        N, n, inter_array, upper_tri = initialize_matrices2(bim_file, [inter_file], chrom1, chrom2, 1e-5)
        
        self.assertEqual(N, 900)
        self.assertEqual(n, 5)
        self.assertFalse(upper_tri)
        
        expected_inter_array = np.zeros((inter_array.shape[0],inter_array.shape[1]))
        for k in [(1,0), (8,11), (10,11), (12,14), (14,14)]:
            expected_inter_array[k[0]][k[1]] = 1
        
        i = 0 
        while i < inter_array.shape[0]:
            j = 0
            while j < inter_array.shape[1]:
                self.assertEqual(inter_array[i][j], expected_inter_array[i][j])
                j += 1
            i += 1
                
        N, n, inter_array, upper_tri = initialize_matrices2(bim_file, [inter_file], chrom1, chrom2, 1e-6)
        
        self.assertEqual(N, 900)
        self.assertEqual(n, 1)
        self.assertFalse(upper_tri)
        
        expected_inter_array = np.zeros((inter_array.shape[0],inter_array.shape[1]))
        for k in [(1,0)]:
            expected_inter_array[k[0]][k[1]] = 1
        
        i = 0 
        while i < inter_array.shape[0]:
            j = 0
            while j < inter_array.shape[1]:
                self.assertEqual(inter_array[i][j], expected_inter_array[i][j])
                j += 1
            i += 1
                
        chrom1 = '3'
        chrom2 = '2'
        N, n, inter_array, upper_tri = initialize_matrices2(bim_file, [inter_file], chrom1, chrom2, 1e-5)
        
        self.assertEqual(N, 900)
        self.assertEqual(n, 5)
        self.assertFalse(upper_tri)
        
        expected_inter_array = np.zeros((inter_array.shape[0],inter_array.shape[1]))
        for k in [(0,1), (11,8), (11,10), (14,12), (14,14)]:
            expected_inter_array[k[0]][k[1]] = 1
        
        i = 0 
        while i < inter_array.shape[0]:
            j = 0
            while j < inter_array.shape[1]:
                self.assertEqual(inter_array[i][j], expected_inter_array[i][j])
                j += 1
            i += 1
                
        chrom1 = '2'
        chrom2 = '3'
        inter_file2 = indir + 'Test_Interaction2.anno'
        N, n, inter_array, upper_tri = initialize_matrices2(bim_file, [inter_file, inter_file2], chrom1, chrom2, 1e-5)
        self.assertEqual(N, 900)
        self.assertEqual(n, 6)
        self.assertFalse(upper_tri)
                
        expected_inter_array = np.zeros((inter_array.shape[0],inter_array.shape[1]))
        for k in [(1,0), (8,11), (10,11), (12,14), (14,14), (27, 25)]:
            expected_inter_array[k[0]][k[1]] = 1
        
        i = 0 
        while i < inter_array.shape[0]:
            j = 0
            while j < inter_array.shape[1]:
                self.assertEqual(inter_array[i][j], expected_inter_array[i][j])
                j += 1
            i += 1
                
                
        N, n, inter_array, upper_tri = initialize_matrices2(bim_file, [inter_file, inter_file2], chrom1, chrom2, 1e-8)
        self.assertEqual(N, 900)
        self.assertEqual(n, 2)
        self.assertFalse(upper_tri)

        expected_inter_array = np.zeros((inter_array.shape[0],inter_array.shape[1]))
        for k in [(8,11), (27, 25)]:
            expected_inter_array[k[0]][k[1]] = 1
        
        i = 0 
        while i < inter_array.shape[0]:
            j = 0
            while j < inter_array.shape[1]:
                self.assertEqual(inter_array[i][j], expected_inter_array[i][j])
                j += 1
            i += 1
        
        
    def test_comp_alt(self):
        inter_matrix = dict()
        i, j = 0, 0
        while i < 10:
            j = 0
            while j < 10:
                inter_matrix[(i,j)] = 0
                j += 1
            i += 1
            
        inter_matrix[(2,8)] = 1
        inter_matrix[(5,8)] = 1
        inter_matrix[(5,9)] = 1
        
        inter_array = convert_dict_to_2D_numpy_array(inter_matrix)
        
        #print(inter_array)
        
        self.assertEqual(comp_alt(0,0,2,2, inter_array, upper_triangular=True), (0,1))
        self.assertEqual(comp_alt(0,0,2,2, inter_array, upper_triangular=False), (0,4))
        
        self.assertEqual(comp_alt(0,0,2,4, inter_array, upper_triangular=True), (0,5))
        self.assertEqual(comp_alt(0,0,2,4, inter_array, upper_triangular=False), (0,8))
        
        self.assertEqual(comp_alt(0,0,3,10, inter_array, upper_triangular=True), (1,24))
        self.assertEqual(comp_alt(0,0,3,10, inter_array, upper_triangular=False), (1,30))
        
        self.assertEqual(comp_alt(5,0,4,4, inter_array, upper_triangular=True), (0,0))
        self.assertEqual(comp_alt(5,3,4,4, inter_array, upper_triangular=True), (0,1))
        self.assertEqual(comp_alt(5,6,2,4, inter_array, upper_triangular=True), (2,7))
        self.assertEqual(comp_alt(5,8,1,1, inter_array, upper_triangular=True), (1,1))
        self.assertEqual(comp_alt(5,8,1,2, inter_array, upper_triangular=True), (2,2))
        self.assertEqual(comp_alt(5,8,2,1, inter_array, upper_triangular=True), (1,2))
        
        # Performance testing
        '''
        inter_matrix = dict()

        i, j = 0, 0
        while i < 2000:
            j = 0
            while j < 2000:
                inter_matrix[(i,j)] = 0
                j += 1
            i += 1
            
        inter_matrix[(2,8)] = 1
        inter_matrix[(5,8)] = 1
        inter_matrix[(5,9)] = 1
        
        inter_array = convert_dict_to_2D_numpy_array(inter_matrix)
        
        start = time.time()
        self.assertEqual(comp_alt(0,0,1850,1900, inter_array, upper_triangular = False), (3, 3515000))
        print(time.time()-start)
        '''
        
    def test_compute_k_m_parallel(self):
        
        inter_array = np.zeros((10,10))
        for k in [(2,8), (5,8), (5,9)]:
            inter_array[k[0]][k[1]] = 1
        
        # Note, technically, the last parameter should be 3, but we enter 0 
        # to force all k,m pairs to be written.
        km_results = compute_k_m_parallel(4, inter_array, True, 0, 9, 45, 0)
        
        self.assertEqual(km_results[(0,5,4,4)], (1,16))
        self.assertEqual(km_results[(0,8,3,1)], (1,3))
        self.assertEqual(km_results[(1,7,2,2)], (1,4))
        self.assertEqual(km_results[(2,5,4,4)], (2,15))
        self.assertEqual(km_results[(2,8,4,2)], (3,8))
        self.assertEqual(km_results[(4,7,3,3)], (2,9))
        
        self.assertEqual(km_results[(5,5,4,4)], (1,6))
        self.assertEqual(km_results[(5,5,3,4)], (1,6))
        self.assertEqual(km_results[(5,5,2,4)], (1,5))
        self.assertEqual(km_results[(5,5,1,4)], (1,3))
        self.assertEqual(km_results[(5,6,4,4)], (2,10))
        self.assertEqual(km_results[(5,6,4,3)], (1,6))
        self.assertEqual(km_results[(5,6,3,4)], (2,9))
        self.assertEqual(km_results[(5,6,3,4)], (2,9))
        self.assertEqual(km_results[(5,6,3,4)], (2,9))
        self.assertEqual(km_results[(5,6,3,3)], (1,6))
        self.assertEqual(km_results[(5,6,2,4)], (2,7))
        
        self.assertEqual(km_results[(5,6,3,4)], (2,9))
        self.assertEqual(km_results[(5,6,3,4)], (2,9))
        self.assertEqual(km_results[(5,6,3,4)], (2,9))
        self.assertEqual(km_results[(5,6,3,4)], (2,9))
        self.assertEqual(km_results[(5,6,3,4)], (2,9))
        self.assertEqual(km_results[(5,6,3,4)], (2,9))
        self.assertEqual(km_results[(5,6,3,4)], (2,9))
        
        self.assertEqual(km_results[(5,6,2,4)], (2,7))
        self.assertEqual(km_results[(5,8,1,1)], (1,1))
        self.assertEqual(km_results[(5,8,1,2)], (2,2))
        self.assertEqual(km_results[(5,8,2,1)], (1,2))
        
        inter_array = np.zeros((10,10))
        for k in [(2,8)]:
            inter_array[k[0]][k[1]] = 1
        
        # Note, technically, the last parameter should be 1, but we enter 0 
        # to force all k,m pairs to be written.
        km_results = compute_k_m_parallel(4, inter_array, True, 0, 9, 45, 0)
        self.assertEqual(len(km_results), 63)
      
    def test_compute_interval_pval_parallel(self): 
        inter_array = np.zeros((10,10))
        for k in [(2,8), (5,8), (5,9)]:
            inter_array[k[0]][k[1]] = 1
        
        # Note, technically, the last parameter should be 3, but we enter 0 
        # to force all k,m pairs to be written.
        km_results = compute_k_m_parallel(4, inter_array, True, 0, 9, 45, 0)
        out_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/Biclustering/Test_Output/'
        pval_results = compute_interval_pval_parallel(km_results, 45, 3, 0, 9, 4, out_dir)
        #print(pval_results)
        
        self.assertAlmostEqual(pval_results[(0,5,4,4)], 0.7424947145877379)
        self.assertAlmostEqual(pval_results[(0,8,3,1)], 0.19097956307258637)
        self.assertAlmostEqual(pval_results[(1,7,2,2)], 0.24876673713883019)
        self.assertAlmostEqual(pval_results[(2,5,4,4)], 0.25405214940098664)
        self.assertAlmostEqual(pval_results[(2,8,4,2)], 0.0039464411557434895)
        self.assertAlmostEqual(pval_results[(4,7,3,3)], 0.09725158562367894)
        
        self.assertAlmostEqual(pval_results[(5,5,4,4)], 0.35595489781536405)
        self.assertAlmostEqual(pval_results[(5,5,3,4)], 0.35595489781536405)
        self.assertAlmostEqual(pval_results[(5,5,2,4)], 0.30373502466525815)
        self.assertAlmostEqual(pval_results[(5,5,1,4)], 0.19097956307258632)
        self.assertAlmostEqual(pval_results[(5,6,4,4)], 0.11945031712473582)
        self.assertAlmostEqual(pval_results[(5,6,4,3)], 0.35595489781536405)
        self.assertAlmostEqual(pval_results[(5,6,3,4)], 0.09725158562367894)
        self.assertAlmostEqual(pval_results[(5,6,3,4)], 0.09725158562367894)
        self.assertAlmostEqual(pval_results[(5,6,3,4)], 0.09725158562367894)
        self.assertAlmostEqual(pval_results[(5,6,3,3)], 0.35595489781536405)
        self.assertAlmostEqual(pval_results[(5,6,2,4)], 0.058703312191684544)
        
        self.assertAlmostEqual(pval_results[(5,6,3,4)], 0.09725158562367894)
        self.assertAlmostEqual(pval_results[(5,6,3,4)], 0.09725158562367894)
        self.assertAlmostEqual(pval_results[(5,6,3,4)], 0.09725158562367894)
        self.assertAlmostEqual(pval_results[(5,6,3,4)], 0.09725158562367894)
        self.assertAlmostEqual(pval_results[(5,6,3,4)], 0.09725158562367894)
        self.assertAlmostEqual(pval_results[(5,6,3,4)], 0.09725158562367894)
        self.assertAlmostEqual(pval_results[(5,6,3,4)], 0.09725158562367894)
        
        self.assertAlmostEqual(pval_results[(5,6,2,4)], 0.058703312191684544)
        self.assertAlmostEqual(pval_results[(5,8,1,1)], 0.06666666666666664)
        self.assertAlmostEqual(pval_results[(5,8,1,2)], 0.003030303030303028)
        self.assertAlmostEqual(pval_results[(5,8,2,1)], 0.13030303030303028)
        
        
    def test_trim_intervals(self):
        in_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/Biclustering/Test_Input/'
        filename = in_dir + 'pval_results_0_1.csv'
        
        pval_results = dict()
        with open(filename) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            for row in csv_reader:
                pval_results[(int(row[0]),int(row[1]),int(row[2]),int(row[3]))] = float(row[6])
        s_inter_pairs = sorted(pval_results.items(), key = lambda x: abs(x[1]), reverse = False)
        #print(len(s_inter_pairs))
        #print(s_inter_pairs)
        trimmed_inter_pairs = trim_intervals(s_inter_pairs)
        #print(trimmed_inter_pairs)
        
        #print(len(trimmed_inter_pairs))
        #print(trimmed_inter_pairs)
        expected_trimmed_inter_pairs = [((0, 816, 48, 8), 2.5317372178516744e-24),
                                        ((0, 999, 35, 3), 8.785313396735072e-16),
                                        ((0, 494, 35, 8), 3.520342278584121e-15)]
        
        i = 0
        for inter_pair in trimmed_inter_pairs:
            self.assertEqual(inter_pair, expected_trimmed_inter_pairs[i])
            i += 1
            
        # COMPLETE THIS TEST CASE
        filename = in_dir + 'chr1_chr2_pval_results.csv'
        pval_results = dict()
        with open(filename) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            for row in csv_reader:
                pval_results[(int(row[0]),int(row[1]),int(row[2]),int(row[3]))] = float(row[6])
        s_inter_pairs = sorted(pval_results.items(), key = lambda x: abs(x[1]), reverse = False)
        trimmed_inter_pairs = trim_intervals(s_inter_pairs)
        #print(trimmed_inter_pairs)
        
        expected_trimmed_inter_pairs = [((1, 41, 57, 59), 1.5052159828764328e-29),
                                        ((65, 53, 35, 47), 5.109419105661788e-22),
                                        ((45, 1, 55, 33), 1.5051851637750398e-21)]
        
        i = 0
        for inter_pair in trimmed_inter_pairs:
            self.assertEqual(inter_pair, expected_trimmed_inter_pairs[i])
            i += 1
     
    def test_convert_dict_to_2D_numpy_array(self):
        tested_inter_matrix = dict()
        i, j = 0, 0
        while i < 10:
            j = 0
            while j < 10:
                if (j > i):
                    tested_inter_matrix[(i,j)] = 1
                else:
                    tested_inter_matrix[(i,j)] = 0
                j += 1
            i += 1
        
        array = convert_dict_to_2D_numpy_array(tested_inter_matrix)
        #print(array)
        #for x in np.nditer(array):
        #    print(x)
     
def test_suite():
    # Unit Tests
    unit_test_suite = unittest.TestSuite()
    unit_test_suite.addTest(TestBiClusterCodeBase('test_initialize_matrices'))
    unit_test_suite.addTest(TestBiClusterCodeBase('test_check_overlap'))
    unit_test_suite.addTest(TestBiClusterCodeBase('test_comp_alt'))
    unit_test_suite.addTest(TestBiClusterCodeBase('test_convert_dict_to_2D_numpy_array'))
    unit_test_suite.addTest(TestBiClusterCodeBase('test_compute_k_m_parallel'))
    unit_test_suite.addTest(TestBiClusterCodeBase('test_compute_interval_pval_parallel'))
    unit_test_suite.addTest(TestBiClusterCodeBase('test_generate_gene_blocks'))
    unit_test_suite.addTest(TestBiClusterCodeBase('test_trim_intervals'))
    unit_test_suite.addTest(TestBiClusterCodeBase('test_initialize_matrices2'))
    unit_test_suite.addTest(TestBiClusterCodeBase('test_parse_biclustering_results'))
    
    unit_test_suite.addTest(TestBiClusterCodeBase('test_get_number_of_interacting_snps_in_interval'))
    unit_test_suite.addTest(TestBiClusterCodeBase('test_parse_remma_interaction_data'))
    unit_test_suite.addTest(TestBiClusterCodeBase('test_get_msigdb_enrichment'))
    
    runner = unittest.TextTestRunner()
    runner.run(unit_test_suite)
    

