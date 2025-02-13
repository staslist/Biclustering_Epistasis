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
    
    indir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/REMMA_Results/'
    marker_file1 = 'epiAA_aa1_approx_parallel_merged_NA3_Combined_Fixed_Set_score900_db001_and_exp001_1_and_genehancer_score25_maf001_region_qc3.anno'
    marker_file2 = 'epiAA_aa2_approx_parallel_merged_NA3_Combined_Fixed_Set_score900_db001_and_exp001_1_and_genehancer_score25_maf001_region_qc3.anno'
    marker_file3 = 'epiAA_aa3_approx_parallel_merged_NA3_Combined_Fixed_Set_score900_db001_and_exp001_1_and_genehancer_score25_maf001_region_qc3.anno'
    marker_file4 = 'epiAD_ad_approx_parallel_merged_NA3_Combined_Fixed_Set_score900_db001_and_exp001_1_and_genehancer_score25_maf001_region_qc3.anno'
    marker_file5 = 'epiDD_dd_approx_parallel_merged_NA3_Combined_Fixed_Set_score900_db001_and_exp001_1_and_genehancer_score25_maf001_region_qc3.anno'
    
    marker_files = [indir + marker_file1, indir + marker_file2, indir + marker_file3, indir + marker_file4, indir + marker_file5]
    marker_files2 = [marker_file1, marker_file2, marker_file3, marker_file4, marker_file5]
    
    bim_file = 'NA3_Combined_Set_score900_db001_and_exp001_1_and_genehancer_score25_maf001_region_qc3.bim'
    
    
    out_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/Biclustering/'
    '''
    chrom_pairs = [(1,2), (1,3), (1,4), (1,6), (1,10), (1,13), (1,16), (2,2), (2,3), (2,5), (2,6), (2,7),
                   (2,8), (2,11), (2,13), (2,15), (2,16), (2,17), (2,18), (2,19), (3,4), (3,5), (3,6),
                   (3,7), (3,8), (3,11), (3,12), (3,15), (3,17), (3,19), (4,4), (4,5), (4,7), 
                   (4,8), (4,9), (4,10), (4,14), (4,15), (4,16), (4,17), (4,18), (5,6), (5,7),
                   (5,8), (5,11), (5,16), (5,21), (6,7), (6,8), (6,10), (6,11), (6,12), (6,15),
                   (6,16), (6,20), (7,7), (7,8), (7,11), (7,14), (7,16), (7,21), (8,9), (8,11), 
                   (8,16), (8,18), (8,20), (9,13), (9,15), (10,16), (10,17), (10,20), (11,17), 
                   (11,21), (12,14), (12,17), (12,19), (15,18), (15,19), (15,21), (16,16), (16,17),
                   (16, 19), (17,19), (17,21)]
    cum_N = 0
    for chrom_pair in chrom_pairs:
        chrom1 = chrom_pair[0]
        chrom2 = chrom_pair[1]    
        
        if(chrom2 < chrom1):
            continue
    
        N,n,inter_array,up_tri = initialize_matrices2(indir + bim_file, marker_files, str(chrom1), str(chrom2), 1e-8)
        out_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/Biclustering/'
        #print("Initialized matrices.")
        print("N = ", N)
        print("n = ", n)
        print(inter_array.shape)
        
        cum_N += N
        
        # if(n > 1):
        #     generate_multicpu_files(marker_files2, bim_file, inter_array.shape[0], inter_array.shape[0]//1000, out_dir, str(chrom1), str(chrom2), 1e-8)
    '''
    
    '''
    bim_file = 'NA3_Combined_Set_score900_db001_and_exp001_1_and_genehancer_score25_maf001_region_qc3.bim'
    bim_file = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/REMMA_Results/' + bim_file
    out_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/REMMA_Results/Combined_Fixed2_Set_Biclustering/'
    gene_location_file = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/AUD_Resources/Start_Sets/Expanded_Ranges/'
    gene_location_file += 'Combined_Set_score900_db001_and_exp001_1_and_genehancer_score25_ranges.txt' 
    i = 1
    pairs = set()
    while i <= 22:
        j = i
        while j <= 22:
            chrom1,chrom2 = str(i), str(j)
            
        # chrom1,chrom2 = str(chrom_pair[0]), str(chrom_pair[1])
            try:
                bicluster_file = out_dir + 'biclustering_results_chr'+chrom1+'_chr'+chrom2+'.txt'
                parse_biclustering_results(bicluster_file, bim_file, gene_location_file,
                                            chrom1, chrom2, marker_files, 1e-8, out_dir, True)
                
                print('chrom1: ', chrom1, '; chrom2: ', chrom2)
                # print(get_number_of_interacting_snps_in_interval(bicluster_file, bim_file, gene_location_file,
                #                                                  chrom1, chrom2, marker_files, 1e-8))
                
                # fname = out_dir + 'biclustering_results_chr'+chrom1+'_chr'+chrom2+'_annotated.txt'
                # with open(fname) as csv_file:
                #     csv_reader = csv.reader(csv_file, delimiter='\t')
                #     for row in csv_reader:
                        
                #         elements = set()
                #         elements2 = set()
                        
                #         interval0 = row[0][1:-1]
                #         interval1 = row[1][1:-1]
                #         if(',' in interval0):
                #             row0_elements = interval0.split(', ')
                #             for ele in row0_elements:
                #                 index = ele.find('_ALT')
                #                 if(index != -1):
                #                     ele = ele[0:index]
                #                 elements.add(ele.replace("'", ""))
                #         else:
                #             index = interval0.find('_ALT')
                #             if(index != -1):
                #                 interval0 = interval0[0:index]
                #             elements.add(interval0.replace("'", ""))
                        
                #         if(',' in interval1):
                #             row1_elements = interval1.split(', ')
                #             for ele in row1_elements:
                #                 index = ele.find('_ALT')
                #                 if(index != -1):
                #                     ele = ele[0:index]
                #                 elements2.add(ele.replace("'", ""))
                #         else:
                #             index = interval1.find('_ALT')
                #             if(index != -1):
                #                 interval1 = interval1[0:index]
                #             elements2.add(interval1.replace("'", ""))
                        
                #         for ele in elements:
                #             for ele2 in elements2:
                #                 pairs.add((ele, ele2))
                #                 print((ele, ele2))
                              
            except FileNotFoundError:
                pass
            j += 1
        i += 1
    '''
    
    # gene_pairs = set()
    # out_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/REMMA_Results/Combined_Fixed2_Set_Biclustering/'
    # fname = out_dir + 'biclustering_regulatory_gene_pairs_temporary.txt'
    # with open(fname) as csv_file:
    #     csv_reader = csv.reader(csv_file, delimiter=',')
    #     for row in csv_reader:
    #         gene_pairs.add((row[0], row[1]))
    # get_msigdb_enrichment(gene_pairs, out_dir)
    
    # fname_out = out_dir + 'biclustering_results_regulatory_elements.txt'
    # with open(fname_out, 'w') as writer:
    #     for element in elements:
    #         writer.write(element.replace("'", "") + '\n')
    
    
    '''
    out_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/REMMA_Results/Combined_Fixed2_Set_Biclustering/'
    fname = out_dir + 'biclustering_regulatory_gene_pairs_temporary.txt'
    to_write = []
    with open(fname) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            to_write.append(row[0])
            to_write.append(row[1])
            
    fname_out = out_dir + 'biclustering_genes.txt'
    with open(fname_out, 'w') as writer:
        for gene in to_write:
            writer.write(gene.replace("'", "") + '\n')
    '''        
            
            
    
    '''
    fname_out = 'combined_set_background_genes.txt'
    background_genes = []
    indir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/AUD_Resources/Start_Sets/Expanded_Ranges/'
    fname = indir + 'Combined_Set_score900_db001_and_exp001_1_and_genehancer_score25_ranges.txt'
    with open(fname) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        for row in csv_reader:
            background_genes.append(row[3])
            
    with open(indir + fname_out, 'w') as writer:
        for gene in background_genes:
            writer.write(gene.replace("'", "") + '\n')
    '''
    
    # TEMPORARY CODE
    # indir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/REMMA_Results/Combined_Fixed_Set_Biclustering/'
    # fname = 'biclustering_regulatory_gene_pairs_temporary.txt'
    # gene_networks = generate_gene_networks(indir + fname)
    # print(gene_networks)
    
    