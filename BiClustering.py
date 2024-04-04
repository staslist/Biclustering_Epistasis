# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 02:54:48 2024

@author: staslist
"""

import csv
import unittest
import numpy as np
import time
from scipy.stats import hypergeom

# Bi-Clustering Algorithm Code

# Assume I have a file that lists all indices in the matrix that are 1.

# GLOBAL VARIABLES

# For tested interaction matrix, assume that every interaction has been tested (exhaustive algorithm)

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

def generate_multicpu_files(total_markers:int, num_cpus:int, out_dir:str, chrom1:str, chrom2:str,
                            pval_cutoff:float):
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
        fname_out_py += '_' + i_s + '_' + i_e + '.py'
        with open(fname_out_py, 'w') as writer:
            writer.write('from BiClustering import *\n')
            writer.write('import sys\n')
            writer.write('if __name__ == "__main__":\n')
            writer.write('\tassert sys.version_info.major == 3\n')
            writer.write('\tassert sys.version_info.minor >= 7\n')
            writer.write("\tindir = '/gpfs/group/home/slistopad/BiClustering/'\n")
            writer.write("\tmarker_file1 = indir + 'epiAA_aa1_approx_parallel_merged_NA3_Combined_Strict_Set_score900_db_exp_1_and_genehancer_score10_regions_maf001_qc3.anno'\n")
            writer.write("\tmarker_file2 = indir + 'epiAA_aa2_approx_parallel_merged_NA3_Combined_Strict_Set_score900_db_exp_1_and_genehancer_score10_regions_maf001_qc3.anno'\n")
            writer.write("\tmarker_file3 = indir + 'epiAD_ad_approx_parallel_merged_NA3_Combined_Strict_Set_score900_db_exp_1_and_genehancer_score10_regions_maf001_qc3.anno'\n")
            writer.write("\tmarker_file4 = indir + 'epiDD_dd_approx_parallel_merged_NA3_Combined_Strict_Set_score900_db_exp_1_and_genehancer_score10_regions_maf001_qc3.anno'\n")
            writer.write('\tmarker_files = [marker_file1, marker_file2, marker_file3, marker_file4]\n')
            writer.write("\tbim_file = indir + 'NA3_Combined_Strict_Set_score900_db_exp_1_and_genehancer_score10_regions_maf001_qc3.bim'\n")
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
        
        fname_out_sh = out_dir + 'BiClustering_Launcher_chr' + chrom1 + '_chr' + chrom2 + '_' + i_s + '_' + i_e + '.sh'
        with open(fname_out_sh, 'w') as writer:
            writer.write("#!/bin/sh\n")
            writer.write("#SBATCH --job-name=SL_BiClustering_Launcher_chr" + chrom1 + '_chr' + chrom2 + '_' + i_s + "_" + i_e + "\n")
            writer.write("#SBATCH --nodes=1\n")
            writer.write("#SBATCH --ntasks=1\n")
            writer.write("#SBATCH --cpus-per-task=1\n")
            writer.write("#SBATCH --mem=32gb\n")
            writer.write("#SBATCH --time=240:00:00\n")
            writer.write("#SBATCH --partition=shared\n\n")
            writer.write("cd $SLURM_SUBMIT_DIR\n")
            writer.write("module load use.own\n")
            writer.write("module load python/3.8.3\n")
            writer.write("export PYTHONPATH=/gpfs/home/slistopad/.local/lib/python3.8/site-packages:$PYTHONPATH\n")
            writer.write("python3 " + 'BiClustering_Launcher_chr' + chrom1 + '_chr' + chrom2 + '_' + i_s + '_' + i_e + '.py')
        counter += 1
    
def parse_biclustering_results(biclustering_file:str, bim_file:str, chrom1:str, chrom2:str, 
                               out_dir:str):
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
            
    filename = out_dir + 'biclustering_results_chr' + chrom1 + '_chr' + chrom2 + '_annotated.txt'
    with open(filename, 'w') as writer:
        i = 0
        num_interval_pairs = len(chrom1_ranges)
        while i < num_interval_pairs:
            #print(chrom1_ranges[i][0])
            #print(chrom1_ranges[i][1])
            #print(chrom1_snps[chrom1_ranges[i][0]:chrom1_ranges[i][1]])
            writer.write(str(chrom1_snps[chrom1_ranges[i][0]:chrom1_ranges[i][1]]) + '\t')
            writer.write(str(chrom2_snps[chrom2_ranges[i][0]:chrom2_ranges[i][1]]) + '\t')
            writer.write(p_values[i] + '\n')
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
    inter_matrix = dict()

    total_markers1 = len(chrom1_snps)
    total_markers2 = len(chrom2_snps)
    
    if(upper_triangular):
        i = 0
        while i < total_markers1:
            N += i
            i += 1
    else:
        N = total_markers1 * total_markers2
    
    i = 0
    while i < total_markers1:
        j = 0
        while j < total_markers2:
            if(upper_triangular):
                if(i < j):
                    j+=1
                    continue
            if((chrom1_snps[i], chrom2_snps[j]) in interactions or
               (chrom2_snps[j], chrom1_snps[i]) in interactions):
                #print(i,j)
                inter_matrix[(i,j)] = 1
                n += 1
            else:
                inter_matrix[(i,j)] = 0 
                
            j += 1
            
        i += 1
                
    return N, n, inter_matrix, upper_triangular
    
  
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
    
    if(compute_k):
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

def compute_k_m_parallel(max_inter_length:int, inter_matrix:dict,
                         upper_triangular:bool, i_start:int, i_end:int, N:int, n:int):
    
    # Our only knowledge about inter_matrix is that it consists of 1s and 0s 
    # and that 1s are typically very sparce.
    
    # We have stronger assumptions about tested_inter_matrix which we deem to be either 
    # upper trinagular (with 1s and 0s) or completely filled with 1s. 
    # Given this we do not need to actually parse tested_inter_matrix to compute m.
    
    # Convert inter_matrix and tested_inter_matrix into numpy arrays.
    inter_array = convert_dict_to_2D_numpy_array(inter_matrix)
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
                        # Do not store the k,m pair if k < expected value (m*n/N)
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
    
    
    return km_results


def compute_interval_pval_parallel(km_results:dict, N:int, n:int, i_start:int, i_end:int,
                                   max_inter_length:int, out_dir:str, chrom1:str = '-1',
                                   chrom2:str = '-1'):
    # Assume that all k values are bigger than expected k-value of (m*n/N)
    
    computed_pvals = dict()
    #correction = 0.05 / ((N**2) * (max_inter_length**2) * 0.5)
    correction = 0.05/((N**2)*0.5)
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
                if(pval_inner <= pval_outer):
                    indeces_to_delete.add(outer_count)
                else:
                    indeces_to_delete.add(inner_count)
                    
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
    
    def test_generate_gene_blocks(self):
        fname = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/Biclustering/Test_Input/'
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
        fname = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/Biclustering/Test_Input/TenByTen.txt'
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
        indir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/Biclustering/Test_Input/'
        bim_file = indir + 'Test_Data.bim'
        inter_file = indir + 'Test_Interaction.anno'
        chrom1 = '1'
        chrom2 = '1'
        
        N, n, inter_matrix, upper_tri = initialize_matrices2(bim_file, [inter_file], chrom1, chrom2, 1e-5)
        self.assertEqual(N, 435)
        self.assertEqual(n, 1)
        self.assertTrue(upper_tri)
        
        for k,v in inter_matrix.items():
            if(k == (19,12)):
                self.assertEqual(v, 1)
            else:
                self.assertEqual(v, 0)
                
        indir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/Biclustering/Test_Input/'
        bim_file = indir + 'Test_Data.bim'
        inter_file = indir + 'Test_Interaction.anno'
        chrom1 = '2'
        chrom2 = '3'
        
        N, n, inter_matrix, upper_tri = initialize_matrices2(bim_file, [inter_file], chrom1, chrom2, 1e-5)
        
        self.assertEqual(N, 900)
        self.assertEqual(n, 5)
        self.assertFalse(upper_tri)
        for k,v in inter_matrix.items():
            if(k in [(1,0), (8,11), (10,11), (12,14), (14,14)]):
                self.assertEqual(v, 1)
            else:
                self.assertEqual(v, 0)
                
        N, n, inter_matrix, upper_tri = initialize_matrices2(bim_file, [inter_file], chrom1, chrom2, 1e-6)
        
        self.assertEqual(N, 900)
        self.assertEqual(n, 1)
        self.assertFalse(upper_tri)
        for k,v in inter_matrix.items():
            if(k in [(1,0)]):
                self.assertEqual(v, 1)
            else:
                self.assertEqual(v, 0)
                
        chrom1 = '3'
        chrom2 = '2'
        N, n, inter_matrix, upper_tri = initialize_matrices2(bim_file, [inter_file], chrom1, chrom2, 1e-5)
        
        self.assertEqual(N, 900)
        self.assertEqual(n, 5)
        self.assertFalse(upper_tri)
        for k,v in inter_matrix.items():
            if(k in [(0,1), (11,8), (11,10), (14,12), (14,14)]):
                self.assertEqual(v, 1)
            else:
                self.assertEqual(v, 0)
                
        chrom1 = '2'
        chrom2 = '3'
        inter_file2 = indir + 'Test_Interaction2.anno'
        N, n, inter_matrix, upper_tri = initialize_matrices2(bim_file, [inter_file, inter_file2], chrom1, chrom2, 1e-5)
        self.assertEqual(N, 900)
        self.assertEqual(n, 6)
        self.assertFalse(upper_tri)
        for k,v in inter_matrix.items():
            if(k in [(1,0), (8,11), (10,11), (12,14), (14,14), (27, 25)]):
                self.assertEqual(v, 1)
            else:
                self.assertEqual(v, 0)
                
                
        N, n, inter_matrix, upper_tri = initialize_matrices2(bim_file, [inter_file, inter_file2], chrom1, chrom2, 1e-8)
        self.assertEqual(N, 900)
        self.assertEqual(n, 2)
        self.assertFalse(upper_tri)
        for k,v in inter_matrix.items():
            if(k in [(8,11), (27, 25)]):
                self.assertEqual(v, 1)
            else:
                self.assertEqual(v, 0)
        
        
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
        
        # Note, technically, the last parameter should be 3, but we enter 0 
        # to force all k,m pairs to be written.
        km_results = compute_k_m_parallel(4, inter_matrix, True, 0, 9, 45, 0)
        
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
        
        inter_matrix = dict()
        i, j = 0, 0
        while i < 10:
            j = 0
            while j < 10:
                inter_matrix[(i,j)] = 0
                j += 1
            i += 1
            
        inter_matrix[(2,8)] = 1
        
        # Note, technically, the last parameter should be 1, but we enter 0 
        # to force all k,m pairs to be written.
        km_results = compute_k_m_parallel(4, inter_matrix, True, 0, 9, 45, 0)
        self.assertEqual(len(km_results), 63)
      
    def test_compute_interval_pval_parallel(self):
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
        
        # Note, technically, the last parameter should be 3, but we enter 0 
        # to force all k,m pairs to be written.
        km_results = compute_k_m_parallel(4, inter_matrix, True, 0, 9, 45, 0)
        out_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/Biclustering/Test_Output/'
        pval_results = compute_interval_pval_parallel(km_results, 45, 3, 0, 9, 4, out_dir)
        #print(pval_results)
        
        self.assertEqual(pval_results[(0,5,4,4)], 0.7424947145877375)
        self.assertEqual(pval_results[(0,8,3,1)], 0.19097956307258632)
        self.assertEqual(pval_results[(1,7,2,2)], 0.24876673713883074)
        self.assertEqual(pval_results[(2,5,4,4)], 0.2540521494009871)
        self.assertEqual(pval_results[(2,8,4,2)], 0.0039464411557434895)
        self.assertEqual(pval_results[(4,7,3,3)], 0.09725158562367894)
        
        self.assertEqual(pval_results[(5,5,4,4)], 0.35595489781536405)
        self.assertEqual(pval_results[(5,5,3,4)], 0.35595489781536405)
        self.assertEqual(pval_results[(5,5,2,4)], 0.30373502466525815)
        self.assertEqual(pval_results[(5,5,1,4)], 0.19097956307258632)
        self.assertEqual(pval_results[(5,6,4,4)], 0.11945031712473582)
        self.assertEqual(pval_results[(5,6,4,3)], 0.35595489781536405)
        self.assertEqual(pval_results[(5,6,3,4)], 0.09725158562367894)
        self.assertEqual(pval_results[(5,6,3,4)], 0.09725158562367894)
        self.assertEqual(pval_results[(5,6,3,4)], 0.09725158562367894)
        self.assertEqual(pval_results[(5,6,3,3)], 0.35595489781536405)
        self.assertEqual(pval_results[(5,6,2,4)], 0.058703312191684544)
        
        self.assertEqual(pval_results[(5,6,3,4)], 0.09725158562367894)
        self.assertEqual(pval_results[(5,6,3,4)], 0.09725158562367894)
        self.assertEqual(pval_results[(5,6,3,4)], 0.09725158562367894)
        self.assertEqual(pval_results[(5,6,3,4)], 0.09725158562367894)
        self.assertEqual(pval_results[(5,6,3,4)], 0.09725158562367894)
        self.assertEqual(pval_results[(5,6,3,4)], 0.09725158562367894)
        self.assertEqual(pval_results[(5,6,3,4)], 0.09725158562367894)
        
        self.assertEqual(pval_results[(5,6,2,4)], 0.058703312191684544)
        self.assertEqual(pval_results[(5,8,1,1)], 0.06666666666666664)
        self.assertEqual(pval_results[(5,8,1,2)], 0.003030303030303028)
        self.assertEqual(pval_results[(5,8,2,1)], 0.13030303030303028)
        
        
    def test_trim_intervals(self):
        in_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/Biclustering/Test_Input/'
        filename = in_dir + 'pval_results_0_1.csv'
        pval_results = dict()
        with open(filename) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            for row in csv_reader:
                pval_results[(int(row[0]),int(row[1]),int(row[2]),int(row[3]))] = float(row[6])
        s_inter_pairs = sorted(pval_results.items(), key = lambda x: abs(x[1]), reverse = False)
        #print(len(s_inter_pairs))
        trimmed_inter_pairs = trim_intervals(s_inter_pairs)
        
        #print(len(trimmed_inter_pairs))
        #print(trimmed_inter_pairs)
        expected_trimmed_inter_pairs = [((0, 816, 48, 8), 2.5317372178516744e-24),
                                        ((0, 999, 35, 3), 8.785313396735072e-16),
                                        ((0, 494, 35, 8), 3.520342278584121e-15)]
        
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
    
    runner = unittest.TextTestRunner()
    runner.run(unit_test_suite)
    

