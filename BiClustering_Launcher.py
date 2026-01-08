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
    
    #test_suite()
    
    indir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/REMMA_Results/'
    marker_file1 = 'epiAA_aa1_approx_parallel_merged_NA3_Combined_Fixed_Set_score900_db001_and_exp001_1_and_genehancer_score25_maf001_region_qc3.anno'
    marker_file2 = 'epiAA_aa2_approx_parallel_merged_NA3_Combined_Fixed_Set_score900_db001_and_exp001_1_and_genehancer_score25_maf001_region_qc3.anno'
    marker_file3 = 'epiAA_aa3_approx_parallel_merged_NA3_Combined_Fixed_Set_score900_db001_and_exp001_1_and_genehancer_score25_maf001_region_qc3.anno'
    marker_file4 = 'epiAD_ad_approx_parallel_merged_NA3_Combined_Fixed_Set_score900_db001_and_exp001_1_and_genehancer_score25_maf001_region_qc3.anno'
    marker_file5 = 'epiDD_dd_approx_parallel_merged_NA3_Combined_Fixed_Set_score900_db001_and_exp001_1_and_genehancer_score25_maf001_region_qc3.anno'
    
    marker_files = [indir + marker_file1, indir + marker_file2, indir + marker_file3, indir + marker_file4, indir + marker_file5]
    marker_files2 = [marker_file1, marker_file2, marker_file3, marker_file4, marker_file5]
    
    bim_file = 'NA3_Combined_Set_score900_db001_and_exp001_1_and_genehancer_score25_maf001_region_qc3.bim'
    bim_file = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/REMMA_Results/' + bim_file
    out_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/REMMA_Results/Combined_Fixed2_Set_Biclustering/'
    gene_location_file = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/Gene_Annotations/'
    gene_location_file += 'ucsc-hg19.txt'
    #gene_location_file = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/AUD_Resources/Start_Sets/Expanded_Ranges/'
    #gene_location_file += 'Combined_Set_score900_db001_and_exp001_1_and_genehancer_score25_ranges.txt' 
    gene_hancer_file = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/AUD_Resources/Genehancer/'
    gene_hancer_file += 'GeneHancer_AnnotSV_hg19_v5.18.txt'
    
    
    element_pairs_dict_all = dict()
    element_pairs_annotated_all = []
    interacting_snps_all = set()
    interacting_snps_freqs_all = dict()
    interacting_snp_pairs_ld = dict()
    interacting_snp_pairs = set()
    reg_elements_all = set()
    i = 1
    while i <= 22:
        j = i
        while j <= 22:
            
            chrom1,chrom2 = str(i), str(j)
            #chrom1,chrom2 = '1','2'
            try:
                # bicluster_file = out_dir + 'biclustering_results_chr'+chrom1+'_chr'+chrom2+'.txt'
                # parse_biclustering_results(bicluster_file, bim_file, gene_location_file, gene_hancer_file,
                #                          chrom1, chrom2, marker_files, 1e-8, out_dir, False, False, True)
                
                
                bicluster_file_annot = out_dir + 'biclustering_results_chr'+chrom1+'_chr'+chrom2+'_annotated.txt'
                ele_pairs_dict,ele_pairs_annotated,unique_inter_snps,inter_snp_pairs_all,ele_pairs_per_interval,reg_elements = parse_biclustering_annotation(bicluster_file_annot)
                reg_elements_all = reg_elements_all | reg_elements
                
                for k,v in ele_pairs_dict.items():
                    element_pairs_dict_all[k] = v
                    
                element_pairs_annotated_all += ele_pairs_annotated
                
                # Code to compute LD Rsq
                # compute_average_LD_Rsq_of_inter_snp_pairs(out_dir, chrom1, chrom2, interacting_snps_all,
                #                                           interacting_snp_pairs_ld, inter_snp_pairs_all, 
                #                                           ele_pairs_per_interval)
                
                
                        
                # CODE TO COMPUTE AVERAGE MAF OF INTERACTING SNPs PER INTERVAL
                # compute_average_maf_of_inter_snps(out_dir, chrom1, chrom2, interacting_snps_all, 
                #                                   interacting_snps_freqs_all, unique_inter_snps)
                
            except FileNotFoundError:
                pass
                
            j += 1
        i += 1
        print(i)
    
    def compute_average_LD_Rsq_of_inter_snp_pairs(in_dir:str, chrom1:str, chrom2:str, interacting_snps_all:set,
                                                  interacting_snp_pairs_ld:dict, inter_snp_pairs_all:list, 
                                                  ele_pairs_per_interval:list):
        bicluster_file_snplist = in_dir + 'biclustering_results_chr'+chrom1+'_chr'+chrom2+'_snplist_interacting_pairs.txt'
        with open(bicluster_file_snplist) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ')
            for row in csv_reader:
                snp1 = row[0]
                snp2 = row[2]
                interacting_snps_all.add(snp1)
                interacting_snps_all.add(snp2)
                
                snp1_split = snp1.split(':')
                snp2_split = snp2.split(':')
                snp1_chrom,snp1_loci = snp1_split[0],snp1_split[1]
                snp2_chrom,snp2_loci = snp2_split[0],snp2_split[1]
                
                ld_fname = 'biclustering_results_ld_' + snp1_chrom + '_' + snp1_loci
                ld_fname += '_' + snp2_chrom + '_' + snp2_loci + '.log'
                ld_filepath = in_dir + 'Interacting_SNP_LD/' + ld_fname
                #print(ld_filepath)
                try: 
                    with open(ld_filepath) as csv_file:
                        csv_reader = csv.reader(csv_file, delimiter=' ')
                        for row in csv_reader:
                            if('R-sq' in row):
                                equal_index = row.index('=')
                                R_sq = row[equal_index+1]
                                if((snp1, snp2) in interacting_snp_pairs_ld):
                                    # If there are two R^2 values average them
                                    # This is not an ideal solution, but it will do for now
                                    interacting_snp_pairs_ld[(snp1, snp2)] = (interacting_snp_pairs_ld[(snp1, snp2)] + float(R_sq))/2
                                else:
                                    interacting_snp_pairs_ld[(snp1, snp2)] = float(R_sq)
                        # if we did not find R-sq row then one of the snps is monomorphic
                        if((snp1, snp2) not in interacting_snp_pairs_ld):
                            interacting_snp_pairs_ld[(snp1, snp2)] = 0
                                
                                
                except FileNotFoundError:
                    ld_fname = 'biclustering_results_ld_' + snp2_chrom + '_' + snp2_loci
                    ld_fname += '_' + snp1_chrom + '_' + snp1_loci + '.log'
                    ld_filepath = in_dir + 'Interacting_SNP_LD/' + ld_fname
                    try: 
                        with open(ld_filepath) as csv_file:
                            csv_reader = csv.reader(csv_file, delimiter=' ')
                            for row in csv_reader:
                                if('R-sq' in row):
                                    equal_index = row.index('=')
                                    R_sq = row[equal_index+1]
                                    if((snp1, snp2) in interacting_snp_pairs_ld):
                                        # If there are two R^2 values average them
                                        # This is not an ideal solution, but it will do for now
                                        interacting_snp_pairs_ld[(snp1, snp2)] = (interacting_snp_pairs_ld[(snp1, snp2)] + float(R_sq))/2
                                    else:
                                        interacting_snp_pairs_ld[(snp1, snp2)] = float(R_sq)
                            # if we did not find R-sq row then one of the snps is monomorphic
                            if((snp1, snp2) not in interacting_snp_pairs_ld):
                                interacting_snp_pairs_ld[(snp1, snp2)] = 0
                    except FileNotFoundError:
                        print('Could not find file for: ', snp1, ', ', snp2)
                        
        z = 1
        for inter_snp_pairs_in_interval in inter_snp_pairs_all:
            avg_LD_Rsq = 0
            for inter_snp_pair in inter_snp_pairs_in_interval:
                try:
                    avg_LD_Rsq += interacting_snp_pairs_ld[inter_snp_pair]
                except KeyError:
                    try:
                        reverse_pair = (inter_snp_pair[1], inter_snp_pair[0])
                        avg_LD_Rsq += interacting_snp_pairs_ld[reverse_pair]
                    except KeyError:
                        print("Missing LD R^2 data for: ", inter_snp_pair)
                    
            avg_LD_Rsq = avg_LD_Rsq / len(inter_snp_pairs_in_interval)
            print('In biclustering_results_chr'+chrom1+'_chr'+chrom2 + ', interacting interval #' + str(z))
            print('Gene Pairs: ',  ele_pairs_per_interval[z-1])
            print(" the average LD R^2 of all interacting snp pairs is: " + str(avg_LD_Rsq))
            print()
            z += 1
    
    
    def compute_average_maf_of_inter_snps(in_dir:str, chrom1:str, chrom2:str, interacting_snps_all:set, 
                                          interacting_snps_freqs_all:dict, unique_inter_snps:list):
        bicluster_file_freq = in_dir + 'Interacting_SNP_Freqs/' + 'biclustering_results_chr'+chrom1+'_chr'+chrom2+'_maf.frq'
        with open(bicluster_file_freq) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter='\t')
            header = True
            for row in csv_reader:
                if(header):
                    header = False
                    continue
                split_row = row[0].split()
                snp = split_row[1]
                maf = split_row[4]
                if(snp in interacting_snps_all):
                    interacting_snps_freqs_all[snp] = float(maf)
                    
        z = 1
        for unique_inter_snps_in_interval in unique_inter_snps:
            avg_maf = 0
            for unique_inter_snp in unique_inter_snps_in_interval:
                #print(unique_inter_snps_in_interval)
                avg_maf += interacting_snps_freqs_all[unique_inter_snp]
            avg_maf = avg_maf / len(unique_inter_snps_in_interval)
            print('In biclustering_results_chr'+chrom1+'_chr'+chrom2 + ', interacting interval #' + str(z), end='')
            print(" the average maf of all interacting snps is: " + str(avg_maf))
            z += 1
    
    
    
    output_file = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/REMMA_Results/Combined_Fixed2_Set_Biclustering/'
    output_file += 'biclustering_results_all_chr_annotated.txt'
    with open(output_file, "w") as f:
        for element_pair,pval in element_pairs_dict_all.items():
            f.write(element_pair + '\t' + str(pval) + '\n')
            
            
    
    element_pairs_annotated2 = []
    
    bicluster_file_annot = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/'
    bicluster_file_annot += 'AoU_Plink_Replication/hg38_Max_Interacting_Intervals_AI_AUD_pval_1e-08/'
    bicluster_file_annot += 'biclustering_of_hg38_AI_AUD_AoU_Max_Interacting_Intervals_Resulting_Annotated_1e_08.txt'
    with open(bicluster_file_annot) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            interval_pair = row[0]
            pval = row[1]
            interval_pair_split = interval_pair.split("},")
            gene_set1 = interval_pair_split[0]
            gene_set2 = interval_pair_split[1]
            
            gene_set1 = gene_set1.split(',')
            gene_set1_formatted = set()
            for gene in gene_set1:
                gene = gene.strip("({'}) ")
                gene_set1_formatted.add(gene)
                
            gene_set2 = gene_set2.split(',')
            gene_set2_formatted = set()
            for gene in gene_set2:
                gene = gene.strip("({'}) ")
                gene_set2_formatted.add(gene)
                
            element_pair_tuple = (gene_set1_formatted, gene_set2_formatted)
            if(element_pair_tuple not in element_pairs_annotated2):
                element_pairs_annotated2.append( element_pair_tuple )
    
    # lets compare the hg19 biclustering annotated genes to hg38 biclustering annotated genes, they should 
    # be similar

    element_pairs_annotated_flattened_AI_AUD = []
    elements_annotated_flattened_AI_AUD = set()
    for tup in element_pairs_annotated_all:
        element_set1 = tup[0]
        element_set2 = tup[1]
        for element1 in element_set1:
            elements_annotated_flattened_AI_AUD.add(element1)
            for element2 in element_set2:
                element_pairs_annotated_flattened_AI_AUD.append( (element1, element2) )
                elements_annotated_flattened_AI_AUD.add(element2)
                
    element_pairs_annotated_flattened_AoU_AI_AUD = []
    elements_annotated_flattened_AoU_AI_AUD = set()
    for tup in element_pairs_annotated2:
        element_set1 = tup[0]
        element_set2 = tup[1]
        for element1 in element_set1:
            elements_annotated_flattened_AoU_AI_AUD.add(element1)
            for element2 in element_set2:
                element_pairs_annotated_flattened_AoU_AI_AUD.append( (element1, element2) )
                elements_annotated_flattened_AoU_AI_AUD.add(element2)
              
    num_replicated = 0
    for gene_pair in element_pairs_annotated_flattened_AI_AUD:
        gene_pair_reverse = (gene_pair[1], gene_pair[0])
        if(gene_pair in element_pairs_annotated_flattened_AoU_AI_AUD or 
           gene_pair_reverse in element_pairs_annotated_flattened_AoU_AI_AUD):
            print(gene_pair)
            num_replicated += 1
            
    print(num_replicated)
     
    #out_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/REMMA_Results/Combined_Fixed2_Set_Biclustering/'
    #get_msigdb_enrichment(element_pairs_annotated_flattened_AI_AUD, out_dir)
    
    regulated_genes_all = set()
    
    indir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/AUD_Resources/Genehancer/'
    fname = indir + 'GeneHancer_v5.18.gff'
    with open(fname) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        header = True
        for row in csv_reader:
            if(header):
                header = False
                continue
            info = row[8]
            info_split = row[8].split(';')
            connected_genes = info_split[1:]
            genehancer_id = info_split[0].split('=')[1]
            if(genehancer_id in reg_elements_all):
                i = 0
                while i < len(connected_genes):
                    connected_gene = connected_genes[i].split('=')[1]
                    connected_score = float(connected_genes[i+1].split('=')[1])
                    if(connected_score >= 25):
                        regulated_genes_all.add(connected_gene)
                    i += 2
                
    
    
    '''
    indir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/REMMA_Results/Combined_Fixed2_Set_Biclustering/'
    hg38_annot_gene_file = indir + 'biclustering_results_hg38_interacting_intervals_AI_AUD.txt'
    
    hg38_biclustering_annotated_genes = set()
    
    with open(hg38_annot_gene_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        for row in csv_reader:
            if('GH' in row[3] and len(row[3]) == 11):
                pass
            elif('_ALT' in row[3]):
                pass
            else:
                hg38_biclustering_annotated_genes.add(row[3])
                
    print(len(hg38_biclustering_annotated_genes))
    print(len(hg38_biclustering_annotated_genes & biclustering_annotated_genes))
    print(hg38_biclustering_annotated_genes - biclustering_annotated_genes)
    print(biclustering_annotated_genes - hg38_biclustering_annotated_genes)
    '''
    
    '''
    chrom_pairs = [('15','21')]
    for chrom_pair in chrom_pairs:
        chrom1 = chrom_pair[0]
        chrom2 = chrom_pair[1]  
        
        bicluster_file = out_dir + 'biclustering_results_chr'+chrom1+'_chr'+chrom2+'.txt'
        parse_biclustering_results(bicluster_file, bim_file, gene_location_file,gene_hancer_file,
                                   chrom1, chrom2, marker_files, 1e-8, out_dir, False, False, True)
    '''
    
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
    
    