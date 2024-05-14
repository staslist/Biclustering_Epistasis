# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 07:34:11 2023

@author: staslist
"""

import csv
import unittest
import requests ## python -m pip install requests

def process_remma_results(full_dir:str, fname:str, generate_annotation:bool = True):
    # Process REMMA Results
    file = full_dir + fname
    inter_pairs = {}
    with open(file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        line_num = 0
        for row in csv_reader:
            # print(row)
            if(line_num > 1):
                inter_pairs[row[2] + '_' + row[9]] = float(row[18])
            line_num += 1
            
    full_dir2 = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/AUD_Resources/Start_Sets/Expanded_Ranges/'
    file2 = full_dir2 + 'Combined_Strict_Set_score900_db_exp_1_and_genehancer_score_10_ranges.txt'
    full_dir3 = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/AUD_Resources/Genehancer/'
    file3 = full_dir3 + 'GeneHancer_AnnotSV_hg19_v5.18.txt'
    gene_locations = []
    reg_locations = []

    with open(file2) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        for row in csv_reader:
            gene_locations.append((row[0],row[1],row[2],row[3]))
            
    with open(file3) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        line_num = 0
        for row in csv_reader:
            if(line_num == 0):
                line_num +=1 
                continue
            reg_locations.append((row[0][3:],row[1],row[2],row[3]))
            line_num += 1
          
    #print(gene_locations)
            
    #print(inter_pairs)
    s_inter_pairs = sorted(inter_pairs.items(), key = lambda x: abs(x[1]), reverse = False)
    # print(s_inter_pairs)

    s_inter_pairs2 = []
    for pair in s_inter_pairs:
        if(pair[1] < 1e-7):
            s_inter_pairs2.append(pair)
    
    if(not generate_annotation):
        return s_inter_pairs2

    s_gene_pairs = []
    s_reg_pairs = []

    num_pairs = len(s_inter_pairs2)
    i = 0
    #print(num_pairs)
    for pair in s_inter_pairs2:
        loci = pair[0].split('_')
        temp1 = loci[0].split(':')
        chrom1 = temp1[0]
        location1 = int(temp1[1])
        temp2 = loci[1].split(':')
        chrom2 = temp2[0]
        location2 = int(temp2[1])
        p_val = pair[1]
        
        gene1 = 'Intergenic'
        gene2 = 'Intergenic'
        
        flag1, flag2 = False, False
        for gene_location in gene_locations:
            if(gene_location[0] == chrom1 and ((location1 > int(gene_location[1])
                                                     and location1 < int(gene_location[2])))):
                gene1 = gene_location[3]
                flag1 = True
            if(gene_location[0] == chrom2 and ((location2 > int(gene_location[1])
                                                     and location2 < int(gene_location[2])))):
                gene2 = gene_location[3]
                flag2 = True
                
            if(flag1 and flag2):
                break
        
        s_gene_pairs.append((gene1, chrom1 + ':' + str(location1), gene2,
                             chrom2 + ':' + str(location2), p_val))
        
        
        reg1, reg2 = 'NOT-FOUND', 'NOT-FOUND'
        
        flag1, flag2 = False, False
        for reg_location in reg_locations:
            if(reg_location[0] == chrom1 and ((location1 > int(reg_location[1])
                                                      and location1 < int(reg_location[2])))):
                reg1 = reg_location[3]
                flag1 = True
            if(reg_location[0] == chrom2 and ((location2 > int(reg_location[1])
                                                      and location2 < int(reg_location[2])))):
                reg2 = reg_location[3]
                flag2 = True
                
            if(flag1 and flag2):
                break
        
        s_reg_pairs.append((reg1, chrom1 + ':' + str(location1), reg2,
                              chrom2 + ':' + str(location2), p_val))
        
        # print(i)
        i += 1


    fname_out = full_dir + 'epiAA_aa1_approx_parallel_merged_NA3_Combined_Strict_Set_score900_db_exp_1_and_genehancer_score_10_regions_maf001_qc3_multibed_gene_results.txt'

    with open(fname_out, 'w') as writer:
        for tup in s_gene_pairs:
            writer.write(tup[0] + ' ' + str(tup[1]) + ' ' + str(tup[2]) +
                         ' ' + str(tup[3]) + ' ' + str(tup[4]))
            writer.write('\n')
              
    fname_out = full_dir + 'epiAA_aa1_approx_parallel_merged_NA3_Combined_Strict_Set_score900_db_exp_1_and_genehancer_score_10_regions_maf001_qc3_multibed_reg_results.txt'

    with open(fname_out, 'w') as writer:
        for tup in s_reg_pairs:
            writer.write(tup[0] + ' ' + str(tup[1]) + ' ' + str(tup[2]) +
                          ' ' + str(tup[3]) + ' ' + str(tup[4]))
            writer.write('\n')

def generate_remma_scripts(out_dir:str, filter_type:str, comp_type:str = 'approx',
                           grm_type:str = 'aa1', num_jobs:int = 2, p_cut_off:float = 1.0e-05):
 
    assert(comp_type in ['approx', 'exact'])
    assert(grm_type in ['aa1', 'aa2', 'aa3', 'ad', 'dd'])
    assert(num_jobs > 1)
    
    i = 1
    while i <= num_jobs:
        fname_base = 'REMMA_' + grm_type + '_' + comp_type + '_parallel' + str(i) + '_'
        fname_base += 'NA3_' + filter_type + '_region_qc3' 
        
        fname_base2 = grm_type + '_' + comp_type + '_parallel' + str(i) + '_NA3_' + filter_type + '_region_qc3' 
        fname_base2_par1 = grm_type + '_' + comp_type + '_parallel1_NA3_' + filter_type + '_region_qc3' 
        fname_base3 = 'NA3_' + filter_type + '_region_qc3' 
        fname_py = fname_base + '.py'
        fname_sh = fname_base + '.sh'
        
        with open(out_dir + fname_py, 'w') as writer:
            writer.write('import logging\n')
            writer.write('logging.basicConfig(level=logging.INFO)\n')
            writer.write('import numpy as np\n')
            writer.write('import os\n')
            writer.write('from gmat.gmatrix import agmat\n')
            writer.write('from gmat.gmatrix import dgmat_as\n')
            writer.write('from gmat.uvlmm.uvlmm_varcom import wemai_multi_gmat\n')
            writer.write('from gmat.remma import annotation_snp_pos\n')
            
            if(comp_type == 'exact'):
                if(grm_type in ['aa1', 'aa2', 'aa3']):
                    writer.write('from gmat.remma.remma_epiAA import remma_epiAA_parallel\n')
                elif(grm_type == 'ad'):
                    writer.write('from gmat.remma.remma_epiAD import remma_epiAD_parallel\n')
                elif(grm_type == 'dd'):
                    writer.write('from gmat.remma.remma_epiDD import remma_epiDD_parallel\n')
            elif(comp_type == 'approx'):
                if(grm_type in ['aa1', 'aa2', 'aa3']):
                    writer.write('from gmat.remma.remma_epiAA import remma_epiAA_approx_parallel\n')
                elif(grm_type == 'ad'):
                    writer.write('from gmat.remma.remma_epiAD import remma_epiAD_approx_parallel\n')
                elif(grm_type == 'dd'):
                    writer.write('from gmat.remma.remma_epiDD import remma_epiDD_approx_parallel\n')
            
            writer.write('# Step 1: Calculate the genomic relationship matrix\n')
            writer.write("home_dir = '/gpfs/group/home/slistopad/REMMA/data/native_american/'\n")
            writer.write("bed_file = home_dir + '" + fname_base3 + "'\n")
            if(i == 1):
                writer.write("agmat(bed_file, home_dir + 'add_genom_rel_matrix_" + fname_base2_par1 + "')\n")
                writer.write("dgmat_as(bed_file, home_dir + 'dom_genom_rel_matrix_" + fname_base2_par1 + "')\n\n")
            
            writer.write("# Step 2: Estimate the variances\n")
            writer.write("add_genom_rel_matrix = home_dir + 'add_genom_rel_matrix_" + fname_base2_par1 + ".agrm.mat_fmt'\n")
            writer.write("dom_genom_rel_matrix = home_dir + 'dom_genom_rel_matrix_" + fname_base2_par1 + ".dgrm_as.mat_fmt'\n")
            writer.write("pheno_file = home_dir + 'pheno_impute_median.ped'  # phenotypic file\n")
            writer.write("ag = np.loadtxt(add_genom_rel_matrix)  # load the additive genomic relationship matrix\n")
            writer.write("dg = np.loadtxt(dom_genom_rel_matrix)\n")
            if(grm_type in ['aa1']):
                writer.write("gmat_lst = [ag, ag*ag]\n")
            elif(grm_type == 'aa2'):
                writer.write("gmat_lst = [ag, dg, ag*ag]\n")
            elif(grm_type == 'aa3'):
                writer.write("gmat_lst = [ag, dg, ag*ag, ag*dg, dg*dg]\n")
            elif(grm_type == 'ad'):
                writer.write("gmat_lst = [ag, dg, ag*ag, ag*dg, dg*dg]\n")
            elif(grm_type == 'dd'):
                writer.write("gmat_lst = [ag, dg, ag*ag, ag*dg, dg*dg]\n")
            
            if(i == 1):
                writer.write("wemai_multi_gmat(pheno_file, bed_file, gmat_lst, out_file= home_dir + 'var_" + fname_base2_par1 + ".txt')\n\n")
            
            writer.write("# Step 3: Test\n")
            writer.write("var_com = np.loadtxt(home_dir + 'var_" + fname_base2_par1 + ".txt')\n")
            
            if(comp_type == 'exact'):
                if(grm_type in ['aa1', 'aa2', 'aa3']):
                    writer.write("remma_epiAA_parallel(pheno_file, bed_file, gmat_lst, var_com, parallel=[" + str(num_jobs) + ", " + str(i) + "],p_cut=" +str(p_cut_off) + ", out_file=home_dir + 'epiAA_" + fname_base2 + "')\n")
                elif(grm_type == 'ad'):
                    writer.write("remma_epiAD_parallel(pheno_file, bed_file, gmat_lst, var_com, parallel=[" + str(num_jobs) + ", " + str(i) + "],p_cut=" +str(p_cut_off) + ", out_file=home_dir + 'epiAD_" + fname_base2 + "')\n")
                elif(grm_type == 'dd'):
                    writer.write("remma_epiDD_parallel(pheno_file, bed_file, gmat_lst, var_com, parallel=[" + str(num_jobs) + ", " + str(i) + "],p_cut=" +str(p_cut_off) + ", out_file=home_dir + 'epiDD_" + fname_base2 + "')\n")
            elif(comp_type == 'approx'):
                if(grm_type in ['aa1', 'aa2', 'aa3']):
                    writer.write("remma_epiAA_approx_parallel(pheno_file, bed_file, gmat_lst, var_com, parallel=[" + str(num_jobs) + ", " + str(i) + "],p_cut=" +str(p_cut_off) + ", out_file=home_dir + 'epiAA_" + fname_base2 + "')\n")
                elif(grm_type == 'ad'):
                    writer.write("remma_epiAD_approx_parallel(pheno_file, bed_file, gmat_lst, var_com, parallel=[" + str(num_jobs) + ", " + str(i) + "],p_cut=" +str(p_cut_off) + ", out_file=home_dir + 'epiAD_" + fname_base2 + "')\n")
                elif(grm_type == 'dd'):
                    writer.write("remma_epiDD_approx_parallel(pheno_file, bed_file, gmat_lst, var_com, parallel=[" + str(num_jobs) + ", " + str(i) + "],p_cut=" +str(p_cut_off) + ", out_file=home_dir + 'epiDD_" + fname_base2 + "')\n")
            
            
            if(i == 1):
                if(comp_type == 'exact'):
                    if(grm_type in ['aa1', 'aa2', 'aa3']):
                        writer.write("prefix = home_dir + 'epiAA_" + grm_type + "_parallel'\n")
                    elif(grm_type == 'ad'):
                        writer.write("prefix = home_dir + 'epiAD_" + grm_type + "_parallel'\n")
                    elif(grm_type == 'dd'):
                        writer.write("prefix = home_dir + 'epiDD_" + grm_type + "_parallel'\n")
                elif(comp_type == 'approx'):
                    if(grm_type in ['aa1', 'aa2', 'aa3']):
                        writer.write("prefix = home_dir + 'epiAA_" + grm_type + "_approx_parallel'\n")
                    elif(grm_type == 'ad'):
                        writer.write("prefix = home_dir + 'epiAD_" + grm_type + "_approx_parallel'\n")
                    elif(grm_type == 'dd'):
                        writer.write("prefix = home_dir + 'epiDD_" + grm_type + "_approx_parallel'\n")
                    
                writer.write("parallel_num = " + str(num_jobs) + "  # the number of parallels\n")
                writer.write("with open(prefix + '_merged_" + fname_base3 + "', 'w') as fout:\n")
                writer.write("\twith open(prefix + '1_" + fname_base3 + ".1') as fin:\n")
                writer.write("\t\thead_line = fin.readline()\n")
                writer.write("\t\tfout.write(head_line)\n")
                writer.write("\tfor i in range(1, parallel_num+1):\n")
                writer.write("\t\twith open(prefix + str(i) + '_" + fname_base3 + ".' + str(i)) as fin:\n")
                writer.write("\t\t\thead_line = fin.readline()\n")
                writer.write("\t\t\tfor line in fin:\n")
                writer.write("\t\t\t\tfout.write(line)\n")
                writer.write("\t\tos.remove(prefix + str(i) + '_" + fname_base3 + ".' + str(i))\n")
                
                writer.write("# Step 4: Select top SNPs and add the SNP position\n")
                
                if(comp_type == 'exact'):
                    if(grm_type in ['aa1', 'aa2', 'aa3']):
                        writer.write("res_file = home_dir + 'epiAA_" + grm_type + "_parallel_merged_" + fname_base3 + "'\n")
                    elif(grm_type == 'ad'):
                        writer.write("res_file = home_dir + 'epiAD_" + grm_type + "_parallel_merged_" + fname_base3 + "'\n")
                    elif(grm_type == 'dd'):
                        writer.write("res_file = home_dir + 'epiDD_" + grm_type + "_parallel_merged_" + fname_base3 + "'\n")
                elif(comp_type == 'approx'):
                    if(grm_type in ['aa1', 'aa2', 'aa3']):
                        writer.write("res_file = home_dir + 'epiAA_" + grm_type + "_approx_parallel_merged_" + fname_base3 + "'\n")
                    elif(grm_type == 'ad'):
                        writer.write("res_file = home_dir + 'epiAD_" + grm_type + "_approx_parallel_merged_" + fname_base3 + "'\n")
                    elif(grm_type == 'dd'):
                        writer.write("res_file = home_dir + 'epiDD_" + grm_type + "_approx_parallel_merged_" + fname_base3 + "'\n")
                
                writer.write("annotation_snp_pos(res_file, bed_file, p_cut=" + str(p_cut_off) + ", dis=0)")
        
        with open(out_dir + fname_sh, 'w') as writer2:
            writer2.write("#!/bin/sh\n")
            writer2.write("#SBATCH --job-name=SL_" + fname_base + "\n")
            writer2.write("#SBATCH --nodes=1\n")
            writer2.write("#SBATCH --ntasks=1\n")
            writer2.write("#SBATCH --cpus-per-task=4\n")
            writer2.write("#SBATCH --mem=48gb\n")
            writer2.write("#SBATCH --time=240:00:00\n")
            writer2.write("#SBATCH --partition=shared\n\n")
            
            writer2.write("cd $SLURM_SUBMIT_DIR\n")
            writer2.write("module load use.own\n")
            writer2.write("module load python/3.8.3\n")
            writer2.write("export PYTHONPATH=/gpfs/home/slistopad/.local/lib/python3.8/site-packages:$PYTHONPATH\n")
            writer2.write("python3 " + fname_py + "\n")
            
        i += 1

def convert_plink_raw_to_simpleM(file:str):
    '''Ignore first six columns. Takes each column starting with 7th, read it in fully 
    starting from line 2, and then write it out as a row into another file with same name,
    but change .raw to _simpleM.txt'''
    rows = dict()
    with open(file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        j = 0
        for row in csv_reader:
            if(j == 0):
                j += 1
                continue
            
            i = 6
            while i < len(row):
                try:
                    rows[i-6].append(row[i])
                except KeyError:
                    rows[i-6] = [row[i]]
                i += 1
            
            j += 1
            print(j)
    
    fname_out = file[0:-4] + '_simpleM.txt'
    with open(fname_out, 'w') as writer:
        for k,v in rows.items():
            i = 1
            limit = len(v)
            for ele in v:
                writer.write(ele)
                if(i < limit):
                    writer.write(' ')
                else:
                    writer.write('\n')
                i += 1
    
def reformat_gene_range_file(file:str, fname_out:str):
    '''Remove 'chr' from each entry in 1st column. Read in the file.
    Then write out a copy without 'chr'. '''
    
    to_write = []
    with open(file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        for row in csv_reader:
            chr_num = row[0][3:]
            line = chr_num + ' ' + row[1] + ' ' + row[2] + ' ' + row[3]
            to_write.append(line)
            
    with open(fname_out, 'w') as writer:
        for line in to_write:
            writer.write(line + '\n')
            

def compare_genome_references(file1:str, file2:str):
    # assume chromosomes match by default
    map1 = dict()
    map2 = dict()
    
    with open(file1) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            map1[row[3]] = [row[1], row[2]]
            
    with open(file2) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            map2[row[3]] = [row[1], row[2]]
            
    mismatched_genes = []            
            
    for k,v in map1.items():
        if(k not in map2):
            continue
        if((v[0] != map2[k][0]) or (v[1] != map2[k][1])):
            mismatched_genes.append(k)
            
    print(len(mismatched_genes))

def verify_gene_ranges(gene_range_file:str):
    incorrect_genes = []
    incorrect_genes_ranges = 0
    suspicious_gene_ranges = 0
    
    with open(gene_range_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        line_num = 0
        for row in csv_reader:
            loc1 = int(row[1])
            loc2 = int(row[2])
            length = loc2 - loc1
            # print(length)
            if(length > 2500000):
                #print("Incorrect length")
                #print(row[3])
                incorrect_genes_ranges += 1
                incorrect_genes.append(row[3])
            # elif(length > 500000):
            #     print("Concerning length")
            #     print(row[3])
            #     suspicious_gene_ranges += 1
            
    return incorrect_genes

def generate_all_assosciated_transcript_ranges(gene_list:list, gene_transc_map_file:str, 
                                               transc_range_file:str, out_dir:str):
    transcript_gene_map = dict()
    to_write = []
    
    with open(gene_transc_map_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            #print(row)
            transcript = row[0]
            gene = row[1]
            if(gene in gene_list):
                transcript_gene_map[transcript] = gene
                
    with open(transc_range_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            if(row[0] in transcript_gene_map):
                to_write.append(row[1] + ' ' + row[2] + ' ' + row[3] + ' ' + transcript_gene_map[row[0]])
                
    # fname_out = 'ucsc-hg19_converted_supplement_helper.txt'
    # with open(out_dir + fname_out, 'w') as writer:
    #     for line in to_write:
    #         writer.write(line + '\n')
    
    alt_gene_ranges = []
    for gene in gene_list:
        current_gene_ranges = set()
        for line in to_write:
            split_line = line.split(' ')
            chromosome = split_line[0]
            #print(split_line)
            if(split_line[3] == gene):
                current_gene_ranges.add((int(split_line[1]), int(split_line[2])))
        
        #print(current_gene_ranges)
        current_gene_ranges = list(current_gene_ranges)
        
        post_clustering = set()
        while (len(current_gene_ranges) > 0):
            indeces_to_remove = [0]
            current_gene_range = current_gene_ranges[0]
            merge_is_not_complete = True
            while(merge_is_not_complete):
                merges = 0
                i = 0
                while i < len(current_gene_ranges):
                    if((current_gene_range == current_gene_ranges[i]) or (i in indeces_to_remove)):
                        i += 1
                        continue
                    elif(_check_overlap_and_merge(current_gene_range, current_gene_ranges[i])):
                        current_gene_range = _check_overlap_and_merge(current_gene_range, current_gene_ranges[i])
                        indeces_to_remove.append(i)
                        
                        #print(current_gene_range)
                        #print(indeces_to_remove)
                        
                        merges += 1
                        
                    i += 1
                # If no merges have been done we are finished.
                if(merges == 0):
                    merge_is_not_complete = False
                
            post_clustering.add(current_gene_range)
            
            i = 0
            new_current_gene_ranges = []
            while i < len(current_gene_ranges):
                if(i not in indeces_to_remove):
                    new_current_gene_ranges.append(current_gene_ranges[i])
                i += 1
                
            current_gene_ranges = new_current_gene_ranges
            #print("Updated current gene ranges: ", current_gene_ranges)
            
        num_clusters = 0
        for cluster in post_clustering:
            if(num_clusters == 0):
                alt_gene_ranges.append( (chromosome, (cluster), gene) )
            elif(num_clusters == 1):
                alt_gene_ranges.append( (chromosome, (cluster), gene + '_ALT') )
            else:
                alt_gene_ranges.append( (chromosome, (cluster), gene + '_ALT' + str(num_clusters)) )
            num_clusters += 1
                     
    return alt_gene_ranges
                    
                    

def _check_overlap_and_merge(gene_range1:tuple, gene_range2:tuple):
    # First check for complete overlap (one within another), then for partial overlap
    if(gene_range2[0] > gene_range1[0] and gene_range2[1] < gene_range1[1]):
        return gene_range1
    
    if(gene_range1[0] > gene_range2[0] and gene_range1[1] < gene_range2[1]):
        return gene_range2
    
    if(gene_range2[0] < gene_range1[1] and gene_range2[0] > gene_range1[0]):
        return (gene_range1[0], gene_range2[1])
    
    if(gene_range2[1] < gene_range1[1] and gene_range2[1] > gene_range1[0]):
        return (gene_range2[0], gene_range1[1])
    
    if(gene_range1[0] < gene_range2[1] and gene_range1[0] > gene_range2[0]):
        return (gene_range2[0], gene_range1[1])
        
    if(gene_range1[1] < gene_range2[1] and gene_range1[1] > gene_range2[0]):
        return (gene_range1[0], gene_range2[1])

    return False

def build_ucsc_gene_ranges_file(ucsc_file:str, ucsc_ref_file:str):
    transcript_to_gene = dict()
    with open(ucsc_ref_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            transcript_to_gene[row[0]] = row[1]
            
    gene_to_loci_map = dict()
    with open(ucsc_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            #print(row)
            gene_name = transcript_to_gene[row[0]]
            #print(gene_name)
            if(gene_name not in gene_to_loci_map):
                gene_to_loci_map[gene_name] = [row[1], row[2], row[3]]
            else:
                if(int(row[2]) < int(gene_to_loci_map[gene_name][1])):
                    gene_to_loci_map[gene_name][1] = row[2]
                elif(int(row[3]) > int(gene_to_loci_map[gene_name][2])):
                    gene_to_loci_map[gene_name][2] = row[3]
    
    fname_out = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/'
    fname_out += 'ucsc-hg19_converted.txt'
    with open(fname_out, 'w') as writer:
        for k,v in gene_to_loci_map.items():
            writer.write(v[0][3:] + '\t' + v[1] + '\t' + v[2] + '\t' + k + '\n')

def build_ncbi_refseq_gene_ranges_file(ncbi_refseq_file:str):
    gene_to_loci_map = dict()
    with open(ncbi_refseq_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            if(row[3] not in gene_to_loci_map):
                gene_to_loci_map[row[3]] = [row[0], row[1], row[2]]
            else:
                if(int(row[1]) < int(gene_to_loci_map[row[3]][1])):
                    gene_to_loci_map[row[3]][1] = row[1]
                elif(int(row[2]) > int(gene_to_loci_map[row[3]][2])):
                    gene_to_loci_map[row[3]][2] = row[2]
    
    fname_out = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/'
    fname_out += 'ncbi_refseq-hg19_converted.txt'
    with open(fname_out, 'w') as writer:
        for k,v in gene_to_loci_map.items():
            writer.write(v[0] + '\t' + v[1] + '\t' + v[2] + '\t' + k + '\n')

def combine_gene_sets(set_names:list):
    # This function takes the listed gene sets, combines them into one (union),
    # and writes out the result into a file.
    combined_gene_set = set([])
    for set_name in set_names:
    
        start_set_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/AUD_Resources/Start_Sets/'
        full_dir = start_set_dir + set_name + '.txt'
        
        with open(full_dir) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter='\t')
            for row in csv_reader:
                combined_gene_set.add(row[0])
                
    name_out = 'combined_set.txt' 
    fname_out = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/'
    fname_out += 'AUD_Resources/Start_Sets/'
    fname_out += name_out
    with open(fname_out, 'w') as writer:
        for ele in combined_gene_set:
            writer.write(str(ele) + '\n')
            
def combine_gene_sets_intersection(set_names:list):
    # This function takes the listed gene sets and returns their intersection.
    # Assume that two sets are provided.
    assert(len(set_names) == 2)
    genesets = []
    for set_name in set_names:
    
        start_set_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/AUD_Resources/Start_Sets/'
        full_dir = start_set_dir + set_name + '.txt'
        current_set = set()
        
        with open(full_dir) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter='\t')
            for row in csv_reader:
                current_set.add(row[0])
                
        genesets.append(current_set)
        
    return genesets[0] & genesets[1]

def gene_hancer_expansion2(reg_elements_fname:str, score_filter:float = 0):
    # This function expands a list of regulatory elements into a list of genes 
    # that these elements regulate. 
    expanded_set_ids = set([])
    
    
    reg_list = []
    with open(reg_elements_fname) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            reg_list.append(row[0])
            
    full_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/AUD_Resources/Genehancer/'
    full_dir += 'GeneHancer_v5.18.gff'
    with open(full_dir) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        first_row = False
        for row in csv_reader:
            if(not first_row):
                first_row += True
                continue
            temp = row[8]
            temp = temp.split(';')
            genehancer_id = temp[0][14:]
            
            if(genehancer_id in reg_list):
            
                num_items = len(temp)
                i = 0
                while i < (num_items-1)/2:
                    gene_name = temp[i*2+1][15:]
                    score = temp[i*2+2][6:]
                    
                    #print(gene_name)
                    #print(score)
                    
                    if(float(score) > score_filter):
                        expanded_set_ids.add(gene_name)
                    i += 1
                    
    
    print(expanded_set_ids)

def gene_hancer_expansion(set_name:str, score_filter:float = 0):
    # This function returns all the genes that regulate the expression of 
    # genes present within the gene_set. 
    expanded_set_ids = set([])
    
    # Read in the gene names from the file.
    start_set_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/AUD_Resources/Start_Sets/'
    full_dir = start_set_dir + set_name + '.txt'
    
    gene_list = []
    
    with open(full_dir) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            gene_list.append(row[0])
    
    full_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/AUD_Resources/Genehancer/'
    full_dir += 'GeneHancer_v5.18.gff'
    with open(full_dir) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        row_num = 0
        for row in csv_reader:
            if(row_num == 0):
                row_num += 1
                continue
            temp = row[8]
            temp = temp.split(';')
            genehancer_id = temp[0][14:]
            
            gene_name = temp[1][15:]
            score = temp[2][6:]
            if(gene_name in gene_list and float(score) > score_filter):
                expanded_set_ids.add(genehancer_id)
            row_num += 1
    
    # Now we need to convert the genehancer ids into loci
    to_write = []
    full_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/AUD_Resources/Genehancer/'
    full_dir += 'GeneHancer_AnnotSV_hg19_v5.18.txt'
    with open(full_dir) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            if(row[3] in expanded_set_ids):
                to_write.append(row)
    
    # Write out the loci of regulatory elements
    name_out = set_name + '_genehancer_ranges_score' + str(score_filter) + '.txt' 
    fname_out = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/'
    fname_out += 'AUD_Resources/Start_Sets/Expanded_Ranges/'
    fname_out += name_out
    with open(fname_out, 'w') as writer:
        for line in to_write:
            for ele in line:
                writer.write(str(ele) + ' ')
            writer.write('\n')
                

def bio_filter_stringdb(set_name:str, score_filter:int, exp_score_filter:float = 0,
                        db_score_filter:float = 0, bool_oper:str = 'or'):
    
    # Note we are only expanding the original set of genes/proteins using their direct 
    # neighbors. Currently there is a database, experimental, and overall score requirement 
    # to be considered a valid neighbor. 
    
    string_api_url = "https://version-12-0.string-db.org/api"
    output_format = "tsv-no-header"
    method = "get_string_ids"
    
    # Read in the gene names from the file.
    start_set_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/AUD_Resources/Start_Sets/'
    full_dir = start_set_dir + set_name + '.txt'
    
    gene_list = []
    
    with open(full_dir) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            gene_list.append(row[0])
    
    ##
    ## Set parameters
    ##
    
    params = {
    
        "identifiers" : "\r".join(gene_list), # your protein list
        "species" : 9606, # species NCBI identifier 
        "limit" : 1, # only one (best) identifier per input protein
        "echo_query" : 1, # see your input identifiers in the output
        "caller_identity" : "stanislav_listopad" # your app name
    
    }
    
    ##
    ## Construct URL
    ##
    
    
    request_url = "/".join([string_api_url, output_format, method])
    
    ##
    ## Call STRING
    ##
    
    results = requests.post(request_url, data=params)
    
    ##
    ## Read and parse the results
    ##
    
    string_start_set = []
    start_set_echo = []
    
    for line in results.text.strip().split("\n"):
        l = line.split("\t")
        input_identifier, string_identifier = l[0], l[2]
        #print("Input:", input_identifier, "STRING:", string_identifier, sep="\t")
        string_start_set.append(string_identifier)
        start_set_echo.append(input_identifier)
        
    # Check that the name mapping worked correctly.
    # Note that non-coding genes will not appear in STRING. These can be simply ignored for 
    # purposes of protein-protein interaction expansion. 
    
    if(set_name == 'Kegg_AUD_Genes'):
        if(len(string_start_set) != 6):
            raise ValueError('Error.')
    elif(set_name == 'Diseases_AUD_Genes'):
        if(len(string_start_set) != 15):
            raise ValueError('Error.')
    elif(set_name == 'MVP_AUD_AUDIT_Genes'):
        if(len(string_start_set) != 24):
            raise ValueError('Error.')
    elif(set_name == 'Malaa_AUD_Genes'):
        if(len(string_start_set) != 53):
            raise ValueError('Error.')
    elif(set_name == 'Disgenet_AUD_Genes'):
        if(len(string_start_set) != 208):
            raise ValueError('Error.')
    elif(set_name == 'Combined_Strict_Set'):
        if(len(string_start_set) != 46):
            raise ValueError('Error.')
            
    print("Number of genes in the starting gene set: ", len(gene_list))
    print("Number of genes that mapped within STRING database: ", print(len(string_start_set)))
    
    
    output_format = "tsv-no-header"
    method = "interaction_partners"
    
    ##
    ## Construct the request
    ##
    
    request_url = "/".join([string_api_url, output_format, method])
    
    ##
    ## Set parameters
    ##
    
    params = {
    
        "identifiers" : "%0d".join(string_start_set), # your protein
        "species" : 9606, # species NCBI identifier 
        "limit" : 1000,
        "required_score" : score_filter,
        "caller_identity" : "stanislav_listopad" # your app name
    
    }
    
    
    ##
    ## Call STRING
    ##
    
    response = requests.post(request_url, data=params)
    
    ##
    ## Read and parse the results
    ##
    
    expanded_set = set(gene_list)
    
    for line in response.text.strip().split("\n"):
    
        l = line.strip().split("\t")
        query_ensp = l[0]
        query_name = l[2]
        partner_name = l[3]
        combined_score = l[5]
        coexpr_score = l[9]
        exp_score = l[10]
        db_score = l[11]
    
        ## print
        if(bool_oper == 'or'):
            cond = (float(db_score) > db_score_filter or float(exp_score) > exp_score_filter)
        elif(bool_oper == 'and'):
            cond = (float(db_score) > db_score_filter and float(exp_score) > exp_score_filter)
        if(cond):
            #print("\t".join([query_name, partner_name, exp_score,db_score, combined_score]))
            expanded_set.add(partner_name)
    
    print(len(expanded_set))
    for gene in gene_list:
        if(gene not in expanded_set):
            print("ERROR!")
            
    # Now filter the glist-hg19.txt file to only include the expanded_set genes.
    # Now write the above out into a new file. Name it after original starting gene set 
    # and key expansion settings. 
    
    to_write = []
    to_write_genes = []
    
    window = 2 #unit of measure is kilobytes
    # the purpose of the window is to account for some of regulatory elements located 
    # near the protein coding genes
    
    full_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/ucsc-hg19_converted.txt'
    with open(full_dir) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            if(row[3] in expanded_set):
                if(window == 0):
                    to_write.append(row)
                else:
                    chr_id = row[0]
                    chr_num = row[3:]
                    temp = [row[0]]
                    temp.append(int(row[1]) - window * 1000)
                    temp.append(int(row[2]) + window * 1000)
                    temp.append(row[3])
                    to_write.append(temp)
                to_write_genes.append(row[3])
                
    # Check how many genes from expanded set are missing in the hg19 ranges file.
    missing_genes = []            
    for gene in expanded_set:
        if (gene not in to_write_genes):
            missing_genes.append(gene)
            
    print(missing_genes)
            
    num_missing_genes = len(missing_genes)
    print("Ratio of missing genes/proteins to all of genes/proteins in expanded set:")
    ratio = num_missing_genes / len(expanded_set)
    print(ratio)
    if(ratio > 0.1):
        print(missing_genes)
        raise ValueError("The number of missing genes/proteins is larger than 10% of expanded set.")
    
    
    # Output start set name, required overall score, all sub_scores used, and number 
    # of expansions performed for each starting gene. 
    if(db_score_filter == 0 and exp_score_filter == 0):
        name_out = set_name + '_' + 'score' + str(score_filter) + '_db_exp_1_ranges.txt' 
    else:
        name_out = set_name + '_' + 'score' + str(score_filter) + '_db'
        name_out += str(db_score_filter) + '_' + bool_oper + '_exp' + str(exp_score_filter) +'_1_ranges.txt' 
    fname_out = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/'
    fname_out += 'AUD_Resources/Start_Sets/Expanded_Ranges/'
    fname_out += name_out
    with open(fname_out, 'w') as writer:
        for line in to_write:
            for ele in line:
                writer.write(str(ele) + ' ')
            writer.write('\n')


class TestBioFilterCodeBase(unittest.TestCase):  
    
    def test_generate_all_assosciated_transcript_ranges(self):
        full_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/AUD_Resources/Start_Sets/Expanded_Ranges/'
        gene_range_file = full_dir + 'Combined_Set_score900_db001_and_exp001_1_and_genehancer_score25_ranges.txt'
        incorrect_genes = verify_gene_ranges(gene_range_file)
        main_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/'
        gene_transc_map = main_dir + 'ucsc-hg19-kgXref.txt'
        transc_range_file = main_dir + 'ucsc-hg19.txt'
        alt_gene_ranges = generate_all_assosciated_transcript_ranges(['MICB'], gene_transc_map, transc_range_file, main_dir)
        
def test_suite():
    # Unit Tests
    unit_test_suite = unittest.TestSuite()
    unit_test_suite.addTest(TestBioFilterCodeBase('test_generate_all_assosciated_transcript_ranges'))
    
    runner = unittest.TextTestRunner()
    runner.run(unit_test_suite)


# Generate remma scripts
# out_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/REMMA_Results/Scripts/'
# generate_remma_scripts(out_dir, 'Disgenet_score900_db_exp_1', num_jobs = 96, 
#                        p_cut_off = 1e-11)
    
# fname = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/REMMA_Results/Bio_Filter_Results/'
# fname += 'NA3_Combined_Strict_Set_score900_db09_exp04_1_regions_maf001_qc3_reg_merged.txt'
# gene_hancer_expansion2(fname, 10.0)          

# print(len(combine_gene_sets_intersection(['Combined_Strict_Set', 'Combined_Set_ALL'])))
# combine_gene_sets(['AI_AUD_Genes', 'Diseases_AUD_Genes','Disgenet_AUD_Genes_Trimmed_0_3_GDA_Score',
#                    'Disgenet_AD_Genes_Trimmed_0_3_GDA_Score', 'Kegg_AD_Genes', 'Mala_AD_Genes',
#                    'Mala_AUD_Genes', 'MVP_AUD_AUDIT_2019_Genes', 'MVP_AUD_AUDIT_2023_Genes',
#                    'Molecular_Targets_Ethanol_Broad'])

# bio_filter_stringdb('Combined_Set', 900, 0.01, 0.01, 'and')
# gene_hancer_expansion('Combined_Set', 25)
# full_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/AUD_Resources/Start_Sets/Expanded_Ranges/'
# reformat_gene_range_file(full_dir + 'Combined_Set_genehancer_ranges_score25.txt',
#                          full_dir + 'Combined_Set_genehancer_ranges_score25_formatted.txt')


# Compare remma results
'''
full_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/REMMA_Results/'
full_dir1 = full_dir + 'Bio_Filter_Results/Additional_Results/'
full_dir2 = full_dir + 'Bio_Filter_Results/'

results1 = process_remma_results(full_dir1, 'epiAA_aa1_approx_parallel_merged_NA3_Combined_Strict_Set_score900_db_exp_1_and_genehancer_score10_regions_maf001_qc3_multibed.anno', False)
results2 = process_remma_results(full_dir2, 'epiAA_aa1_approx_parallel_merged_NA3_Combined_Strict_Set_score900_db_exp_1_and_genehancer_score10_regions_maf001_qc3.anno', False)

inter_counter = 0
pos_differences = []

counter = 0
pos1 = 1
for tup in results1:
    # if(counter == 400):
    #     break
    snp_pair = tup[0]
    pos2 = 1
    for tup2 in results2:
        if(snp_pair == tup2[0]):
            inter_counter += 1
            pos_differences.append(pos1-pos2)
            pos2 = 1
            break
        pos2 += 1
    pos1 += 1
    counter += 1
        
print(inter_counter)

abs_avg_pos_diff = 0
accum = 0
counter = 1
for pos_diff in pos_differences:
    accum += abs(pos_diff)
    counter += 1
    # if(counter > 300):
    #     break
abs_avg_pos_diff = accum/len(pos_differences)
print(abs_avg_pos_diff)
'''

full_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/AUD_Resources/Start_Sets/Expanded_Ranges/'
gene_range_file = full_dir + 'Combined_Set_score900_db001_and_exp001_1_and_genehancer_score25_ranges.txt'
incorrect_genes = verify_gene_ranges(gene_range_file)
main_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/'
gene_transc_map = main_dir + 'ucsc-hg19-kgXref.txt'
transc_range_file = main_dir + 'ucsc-hg19.txt'
generate_all_assosciated_transcript_ranges(['MICB'], gene_transc_map, transc_range_file, main_dir)

# full_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/REMMA_Results/Bio_Filter_Results/'
# convert_plink_raw_to_simpleM(full_dir + 'NA3_Combined_Strict_Set_score900_db_exp_1_and_genehancer_score10_regions_qc3.raw')
