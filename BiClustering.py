# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 02:54:48 2024

@author: staslist
"""

import csv
from scipy.stats import hypergeom

# Bi-Clustering Algorithm Code

# Assume I have a file that lists all indices in the matrix that are 1.

# GLOBAL VARIABLES

# For tested interaction matrix, assume that every interaction has been tested (exhaustive algorithm)



def initialize_matrices(total_markers:int):
    N = 0
    n = 0
    inter_matrix = dict()
    tested_inter_matrix = dict()
    
    i = 0
    while i < total_markers:
        N += i
        i += 1
    print("N = ", N)

    i = 0
    j = 0
    while i < total_markers:
        while j < total_markers:
            inter_matrix[(i,j)] = 0
            if(i >= j):
                tested_inter_matrix[(i,j)] = 0
            else:
                tested_inter_matrix[(i,j)] = 1
            j += 1
              
        j = 0
        i += 1
    
    indir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/Biclustering/'
    with open(indir + 'biclustering_pub_marker_pairs.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            i,j = int(row[0]), int(row[1])
            inter_matrix[(i,j)] = 1
            n += 1
    
    return N, n, inter_matrix, tested_inter_matrix

def check_overlap(i:int, j:int, a:int, b:int, i2:int, j2:int, a2:int, b2:int):
    # Note if two matrices are identical we treat them as non-overlapping. 
    return ( ((i2 > i) and (i2 < (i+a))) or ((j2 > j) and (j2 < (j+b))) or 
            ((i > i2) and (i < (i2+a2))) or ((j > j2) and (j < (j2+b2))))

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
        
    if( (a == 0 or b == 0) or ((i + a) > 5) or ((j + b) > 5)):
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

def compute_k_m(total_markers:int, max_inter_length:int, inter_matrix:dict, tested_inter_matrix:dict):

    i = 0
    j = 0
    a = 1
    b = 1
    
    k_results = dict()
    m_results = dict()
    km_results = dict()
    pval_results = dict()
    
    
    while i < total_markers:
        while j < total_markers:
            while a <= max_inter_length:
                while b <= max_inter_length:
                    if((i + a) <= total_markers and (j+b) <= total_markers):
                        km_results[(i,j,a,b)] = (comp(i,j,a,b,'k', inter_matrix, tested_inter_matrix),
                                                 comp(i,j,a,b,'m', inter_matrix, tested_inter_matrix))
                    b += 1
                b = 1
                a += 1
            a = 1
            b = 1
            j += 1
        a = 1
        b = 1
        j = 0
        i += 1
        print(i)
        
    return km_results

def compute_interval_pval(km_results:dict, N:int, n:int):
    pval_results = []
    for key,val in km_results.items():
        k,m = val[0], val[1]
        max_k = min(m, n)
        min_k = max(0, m-(N-n))
        pval_results[key] = 0
        if(k > (m*n/N)):
            while k <= max_k:
                pval_results[key] += hypergeom.pmf(k,N,n,m)
                k += 1
        else:
            while k >= min_k:
                pval_results[key] += hypergeom.pmf(k, N, n, m)
                k -= 1
                
    # Now go through interval pairs startin from most significant, and remove all overlapping 
    # intervals for the current interval
    
    s_inter_pairs = sorted(pval_results.items(), key = lambda x: abs(x[1]), reverse = False)
    #print(s_inter_pairs)
    
    return s_inter_pairs
   
def trim_intervals(sorted_intervals:list):
    outer_count = 0
    indeces_to_delete = []
    for tup in sorted_intervals:
        # Skip this interval pair if it is already marked for delition
        if(outer_count in indeces_to_delete):
            outer_count += 1
            continue
        inner_count = 0
        curr_tup = tup
        i,j,a,b = curr_tup[0][0],curr_tup[0][1],curr_tup[0][2],curr_tup[0][3]
        pval_outer = curr_tup[1]
        for tup in sorted_intervals:
            if(inner_count in indeces_to_delete):
                inner_count += 1
                continue
            i2,j2,a2,b2 = tup[0][0],tup[0][1],tup[0][2],tup[0][3]
            pval_inner = tup[1]
            if(check_overlap(i,j,a,b,i2,j2,a2,b2)):
                if(pval_inner < pval_outer):
                    indeces_to_delete.append(outer_count)
                else:
                    indeces_to_delete.append(inner_count)
                    
            inner_count += 1
            
        outer_count += 1
            
    culled_inter_pairs = []
    i = 0
    while i < len(s_inter_pairs):
        if(i not in indeces_to_delete):
            culled_inter_pairs.append(s_inter_pairs[i])
        i += 1
        
    return culled_inter_pairs
    

if __name__ == "__main__":
    N,n,inter_matrix,tested_inter_matrix = initialize_matrices(1211)
    print("Initialized matrices.")
    km_results = compute_k_m(1211, 60, inter_matrix, tested_inter_matrix)
    print("Computed k and m values.")
    s_inter_pairs = compute_interval_pval(km_results, N, n)
    print("Computed p-values for interval pairs.")
    trimmed_inter_pairs = trim_intervals(s_inter_pairs)
    print("Trimmed interval pairs.")




        

    
