# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 01:14:51 2024

@author: staslist
"""

from BiClustering import *
import sys

if __name__ == "__main__":
    assert sys.version_info.major == 3
    assert sys.version_info.minor >= 7  
    
    test_suite()
    
    # out_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/Biclustering/'
    # generate_multicpu_files(1211, 64, out_dir)
    
    # SINGLE CPU VERSION
    
    marker_file = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/Biclustering/biclustering_pub_marker_pairs_reduced.csv'
    N,n,inter_matrix,up_tri = initialize_matrices(1211,marker_file)
    out_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/Biclustering/Test_Output/'
    print("Initialized matrices.")
    km_results = compute_k_m_parallel(1211, 60, inter_matrix, up_tri, 0, 1211)
    print("Computed k and m values.")
    s_inter_pairs = compute_interval_pval_parallel(km_results, N, n, 0, 1211, out_dir)
    print("Computed p-values for interval pairs.")
    filename = out_dir + 'pval_results_0_1211.csv'
    pval_results = []
    with open(filename) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            pval_results[(row[0],row[1],row[2],row[3])] = row[6]
    s_inter_pairs = sorted(pval_results.items(), key = lambda x: abs(x[1]), reverse = False)
    trimmed_inter_pairs = trim_intervals(s_inter_pairs)
    print("Trimmed interval pairs.")
    
    
    
    # MULTIPLE CPU VERSION
    '''
    # Initialize matrices on each CPU (computationally inexpensive)
    marker_file = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/Biclustering/biclustering_pub_marker_pairs.csv'
    N,n,inter_matrix,up_tri = initialize_matrices(1211,marker_file)
    out_dir = 'C:/'
    # Divide total number of markers by number of CPUs to be used to attain target i-interval
    # Each CPU involved will have a different i_start and i_end
    i_start, i_end = 0, 10
    km_results = compute_k_m_parallel(1211, 60, inter_matrix, up_tri, i_start,
                                      i_end)
    # Now compute p-value only for k,m counts computed on this CPU
    compute_interval_pval_parallel(km_results, N, n, i_start, i_end, out_dir)
    # Now on a single CPU read in all the p-value results, sort them, and trim the intervals
    # Assume that all the p-value results can be stored in memory and stored on a single CPU
    pval_results = []
    i_intervals = [(0,10), (11, 20)]
    for i_pair in i_intervals:
        filename = out_dir + 'pval_results_' + str(i_pair[0]) + '_' + str(i_pair[1]) + '.csv'
        with open(filename) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            for row in csv_reader:
                pval_results[(row[0],row[1],row[2],row[3])] = row[6]
                
    s_inter_pairs = sorted(pval_results.items(), key = lambda x: abs(x[1]), reverse = False)
    trimmed_inter_pairs = trim_intervals(s_inter_pairs)
    filename = out_dir + 'biclustering_results.txt'
    with open(filename, 'w') as writer:
        for inter_pair in trimmed_inter_pairs:
            writer.write(str(inter_pair[0][0])+','+str(inter_pair[0][1])+','+str(inter_pair[0][2])
                         +','+str(inter_pair[0][3])+','+str(inter_pair[1])+'\n')
    '''