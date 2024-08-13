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
    marker_file1 = 'epiAA_aa1_approx_parallel_merged_NA3_Combined_Set_score900_db001_and_exp001_1_and_genehancer_score25_maf001_region_qc3.anno'
    marker_file2 = 'epiAA_aa2_approx_parallel_merged_NA3_Combined_Set_score900_db001_and_exp001_1_and_genehancer_score25_maf001_region_qc3.anno'
    marker_file3 = 'epiAA_aa3_approx_parallel_merged_NA3_Combined_Set_score900_db001_and_exp001_1_and_genehancer_score25_maf001_region_qc3.anno'
    marker_file4 = 'epiAD_ad_approx_parallel_merged_NA3_Combined_Set_score900_db001_and_exp001_1_and_genehancer_score25_maf001_region_qc3.anno'
    marker_file5 = 'epiDD_dd_approx_parallel_merged_NA3_Combined_Set_score900_db001_and_exp001_1_and_genehancer_score25_maf001_region_qc3.anno'
    
    marker_files = [indir + marker_file1, indir + marker_file2, indir + marker_file3, indir + marker_file4, indir + marker_file5]
    marker_files2 = [marker_file1, marker_file2, marker_file3, marker_file4, marker_file5]
    '''
    bim_file = 'NA3_Combined_Set_score900_db001_and_exp001_1_and_genehancer_score25_maf001_region_qc3.bim'
    
    
    out_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/Biclustering/'
    chrom_list = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
    #chrom_list = [1,2]
    for chrom1 in chrom_list:
        for chrom2 in chrom_list:
            
            if(chrom2 < chrom1):
                continue
        
            N,n,inter_array,up_tri = initialize_matrices2(indir + bim_file, marker_files, str(chrom1), str(chrom2), 1e-8)
            out_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/Biclustering/'
            #print("Initialized matrices.")
            print("N = ", N)
            print("n = ", n)
            print(inter_array.shape)
            
            if(n > 1):
                generate_multicpu_files(marker_files2, bim_file, inter_array.shape[0], inter_array.shape[0]//1000, out_dir, str(chrom1), str(chrom2), 1e-8)
    '''

    '''
    bim_file = 'NA3_Combined_Set_score900_db001_and_exp001_1_and_genehancer_score25_maf001_region_qc3.bim'
    bim_file = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/REMMA_Results/' + bim_file
    out_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/REMMA_Results/Combined_Set_Biclustering/'
    gene_location_file = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/AUD_Resources/Start_Sets/Expanded_Ranges/'
    gene_location_file += 'Combined_Set_score900_db001_and_exp001_1_and_genehancer_score25_ranges.txt' 
    i = 1
    pairs = set()
    while i <= 22:
        j = i
        #j = 6
        while j <= 22:
            chrom1,chrom2 = str(i), str(j)
            
        # chrom1,chrom2 = str(chrom_pair[0]), str(chrom_pair[1])
            try:
                bicluster_file = out_dir + 'biclustering_results_chr'+chrom1+'_chr'+chrom2+'.txt'
                # parse_biclustering_results(bicluster_file, bim_file, gene_location_file,
                #                            chrom1, chrom2, marker_files, 1e-8, out_dir, True)
                
                print('chrom1: ', chrom1, '; chrom2: ', chrom2)
                print(get_number_of_interacting_snps_in_interval(bicluster_file, bim_file, gene_location_file,
                                                                 chrom1, chrom2, marker_files, 1e-8))
                
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
    # out_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/REMMA_Results/Combined_Set_Biclustering/'
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