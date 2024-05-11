# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 01:14:51 2024

@author: staslist
"""

from BiClustering import *
import sys

if __name__ == "__main__":
    # To ensure that dictionaries are ordered. 
    assert sys.version_info.major == 3
    assert sys.version_info.minor >= 7  
    
    test_suite()
    '''
    indir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/REMMA_Results/Bio_Filter_Results/'
    marker_file1 = indir + 'epiAA_aa1_approx_parallel_merged_NA3_Combined_Strict_Set_score900_db_exp_1_and_genehancer_score10_regions_maf001_qc3.anno'
    marker_file2 = indir + 'epiAA_aa2_approx_parallel_merged_NA3_Combined_Strict_Set_score900_db_exp_1_and_genehancer_score10_regions_maf001_qc3.anno'
    marker_file3 = indir + 'epiAD_ad_approx_parallel_merged_NA3_Combined_Strict_Set_score900_db_exp_1_and_genehancer_score10_regions_maf001_qc3.anno'
    marker_file4 = indir + 'epiDD_dd_approx_parallel_merged_NA3_Combined_Strict_Set_score900_db_exp_1_and_genehancer_score10_regions_maf001_qc3.anno'
    marker_files = [marker_file1, marker_file2, marker_file3, marker_file4]
    bim_file = indir + 'NA3_Combined_Strict_Set_score900_db_exp_1_and_genehancer_score10_regions_maf001_qc3.bim'
    chrom1,chrom2 = '6','6'
    N,n,inter_matrix,up_tri = initialize_matrices2(bim_file,marker_files, chrom1, chrom2, 1e-06)
    out_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/REMMA_Results/Bio_Filter_Results/PostBiclustering/'
    i_start, i_end = 4332, 4333
    km_results = compute_k_m_parallel(60, inter_matrix, up_tri, i_start, i_end, N, n)
    compute_interval_pval_parallel(km_results, N, n, i_start, i_end, 60, out_dir, chrom1=chrom1, chrom2=chrom2)
    '''
    
    # fname = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/AUD_Resources/Start_Sets/'
    # fname += 'Expanded_Ranges/Combined_Strict_Set_score900_db_exp_1_and_genehancer_score_10_ranges.txt'
    # blocks = generate_gene_blocks(fname, 100000000)
    # print(len(blocks))
    # print(blocks)
    
    # out_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/Biclustering/'
    # chrom_list = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,19,20,22]
    # for chrom in chrom_list:
    # generate_multicpu_files(1487, 10, out_dir, '19', '19', 1e-6)
    
    # SINGLE CPU VERSION REMMA DATA
    '''
    indir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/REMMA_Results/Bio_Filter_Results/'
    marker_file1 = indir + 'epiAA_aa1_approx_parallel_merged_NA3_Combined_Strict_Set_score900_db_exp_1_and_genehancer_score10_regions_maf001_qc3.anno'
    marker_file2 = indir + 'epiAA_aa2_approx_parallel_merged_NA3_Combined_Strict_Set_score900_db_exp_1_and_genehancer_score10_regions_maf001_qc3.anno'
    marker_file3 = indir + 'epiAD_ad_approx_parallel_merged_NA3_Combined_Strict_Set_score900_db_exp_1_and_genehancer_score10_regions_maf001_qc3.anno'
    marker_file4 = indir + 'epiDD_dd_approx_parallel_merged_NA3_Combined_Strict_Set_score900_db_exp_1_and_genehancer_score10_regions_maf001_qc3.anno'
    marker_files = [marker_file1, marker_file2, marker_file3, marker_file4]
    bim_file = indir + 'NA3_Combined_Strict_Set_score900_db_exp_1_and_genehancer_score10_regions_maf001_qc3.bim'
    
    
    out_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/REMMA_Results/Bio_Filter_Results/PostBiclustering/'
    i = 1
    elements = set()
    while i <= 22:
        j = i
        while j <= 22:
            chrom1,chrom2 = str(i), str(j)
            try:
                parse_biclustering_results(out_dir + 'biclustering_results_chr'+chrom1+'_chr'+chrom2+'.txt',
                                            bim_file, chrom1, chrom2, out_dir)
                fname = out_dir + 'biclustering_results_chr'+chrom1+'_chr'+chrom2+'_annotated.txt'
                with open(fname) as csv_file:
                    csv_reader = csv.reader(csv_file, delimiter='\t')
                    for row in csv_reader:
                        interval0 = row[0][1:-1]
                        interval1 = row[1][1:-1]
                        if(',' in interval0):
                            row0_elements = interval0.split(', ')
                            for ele in row0_elements:
                                elements.add(ele)
                        else:
                            elements.add(interval0)
                        
                        if(',' in interval1):
                            row1_elements = interval1.split(', ')
                            for ele in row1_elements:
                                elements.add(ele)
                        else:
                            elements.add(interval1)
                            
            except FileNotFoundError:
                pass
            j += 1
        i += 1
    '''
    '''
    chrom_list = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,19,20,22]
    for chrom in chrom_list:
        N,n,inter_matrix,up_tri = initialize_matrices2(bim_file,marker_files, str(chrom), str(chrom), 1e-6)
        out_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/Biclustering/Test_Output/'
        print("Initialized matrices.")
        print("N = ", N)
        print("n = ", n)
        print(convert_dict_to_2D_numpy_array(inter_matrix).shape)
    '''
    '''
    #print('Inter Matrix Shape = ', convert_dict_to_2D_numpy_array(inter_matrix).shape)
    km_results = compute_k_m_parallel(60, inter_matrix, up_tri, 152, 160, N, n)
    print("Computed k and m values.")
    print(len(km_results))
    s_inter_pairs = compute_interval_pval_parallel(km_results, N, n, 152, 160, 60, out_dir,
                                                   chrom1=chrom1, chrom2 = chrom2)
    print("Computed p-values for interval pairs.")
    filename = out_dir + 'chr' + chrom1 + '_chr' + chrom2 + '_pval_results_0_1.csv'
    pval_results = dict()
    with open(filename) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            pval_results[(int(row[0]),int(row[1]),int(row[2]),int(row[3]))] = float(row[6])
    s_inter_pairs = sorted(pval_results.items(), key = lambda x: abs(x[1]), reverse = False)
    trimmed_inter_pairs = trim_intervals(s_inter_pairs)
    print("Trimmed interval pairs.")
    '''
    
    # SINGLE CPU VERSION PUBLICATION DATA
    '''
    marker_file = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/Biclustering/biclustering_pub_marker_pairs.csv'
    N,n,inter_matrix,up_tri = initialize_matrices(1211,marker_file)
    out_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/Biclustering/Test_Output/'
    print("Initialized matrices.")
    print("N = ", N)
    print("n = ", n)
    #km_results = compute_k_m_parallel(60, inter_matrix, up_tri, 0, 1211, N, n)
    km_results = compute_k_m_parallel(60, inter_matrix, up_tri, 0, 1, N, n)
    print("Computed k and m values.")
    print(len(km_results))
    #s_inter_pairs = compute_interval_pval_parallel(km_results, N, n, 0, 1211, 60, out_dir)
    s_inter_pairs = compute_interval_pval_parallel(km_results, N, n, 0, 1, 60, out_dir)
    print("Computed p-values for interval pairs.")
    #filename = out_dir + 'pval_results_0_1211.csv'
    filename = out_dir + 'pval_results_0_1.csv'
    pval_results = dict()
    with open(filename) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            pval_results[(int(row[0]),int(row[1]),int(row[2]),int(row[3]))] = float(row[6])
    s_inter_pairs = sorted(pval_results.items(), key = lambda x: abs(x[1]), reverse = False)
    trimmed_inter_pairs = trim_intervals(s_inter_pairs)
    print("Trimmed interval pairs.")
    '''
    
    
    # MULTIPLE CPU VERSION PUBLICATION DATA
    '''
    # COMPUTE k,m,p-values SAME AS FOR SINGLE CPU VERSION PUBLICATION DATA
    
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