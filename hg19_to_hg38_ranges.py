# -*- coding: utf-8 -*-
"""
Created on Sun Oct 12 08:52:47 2025

@author: slistopad
"""

import csv

# Use the hg38 ranges for these genes and regulatory elements to create a hg38 
# combined expanded set ranges file. Add 2kb window around genes only.

combined_set_expanded_hg38_ranges = dict()
genes_hg38_ranges = dict()
regulatory_hg38_ranges = dict()

genehancer_hg38_file = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/AUD_Resources/Genehancer/GeneHancer_AnnotSV_elements_v5.25.txt'
ncbi_hg38_file = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/Gene_Annotations/ncbi_hg38.tsv'
gencode_hg38_file = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/Gene_Annotations/gencode.v38.chr_patch_hapl_scaff.annotation.gtf'

with open(genehancer_hg38_file) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter='\t')
    first_line = True
    for row in csv_reader:
        if(first_line):
            first_line = False
            continue
        regulatory_hg38_ranges[row[3]] = [row[0], (int(row[1]),int(row[2]))]

'''        
# This annotation does not contain alternative gene locations
with open(ncbi_hg38_file) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter='\t')
    first_line = True
    for row in csv_reader:
        if(first_line):
            first_line = False
            continue
        genes_hg38_ranges[row[6]] = (row[3],int(row[1])-2000,int(row[2])+2000)
'''

with open(gencode_hg38_file) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter='\t')
    i = 0
    for row in csv_reader:
        if(i <= 4):
            i += 1
            continue
        if(row[2] == 'gene'):
            add_info = row[8].split(' ')
            gene_name = add_info[5][1:-2]
            if(gene_name in genes_hg38_ranges.keys()):
                genes_hg38_ranges[gene_name].append( (int(row[3]),int(row[4])) )
            else:
                genes_hg38_ranges[gene_name] = [row[0], (int(row[3]),int(row[4]))]
        i += 1


annotated_intervals_genes_hg38 = dict()
annotated_genes = set()
annotated_intervals_reg_hg38 = dict()
input_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/REMMA_Results/Combined_Fixed2_Set_Biclustering/'
i = 1 
while i <= 22:
    print(i)
    j = i
    while j <= 22:
        #print(i)
        #print(j)
        # note that for any given sniplist interval file the same interval may show up multiple times,
        # because only INTERVAL PAIRS have to be unique
        # so keep track of intervals already processed and do not annotate the same one multiple times
        intervals_already_annotated = []
        try:
            fname = input_dir + 'biclustering_results_chr' + str(i) + '_chr' + str(j) + '_snplist_interval_hg38.txt'
            with open(fname) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter='-')
                for row in csv_reader:
                    # If the file is empty move onto the next file
                    if(len(row) == 0):
                        break
                    
                    interval_annotated = False
                    
                    chrom = row[0]
                    l_boundary = int(row[1])
                    u_boundary = int(row[2])
                    
                    if((chrom, l_boundary, u_boundary) in intervals_already_annotated):
                        continue
                    
                    for k,v in genes_hg38_ranges.items():
                        for gene_range in v[1:]:
                            if(chrom == v[0]):
                                if((l_boundary > gene_range[0] and l_boundary < gene_range[1]) or (u_boundary > gene_range[0] and u_boundary < gene_range[1])
                                   or (l_boundary < gene_range[0] and u_boundary > gene_range[1])):
                                    annotated_genes.add(k)
                                    interval_annotated = True
                                    if(k in annotated_intervals_genes_hg38.keys()):
                                        annotated_intervals_genes_hg38[k].append((chrom, l_boundary, u_boundary))
                                    else:
                                        annotated_intervals_genes_hg38[k] = [(chrom, l_boundary, u_boundary)]
                                
                    for k,v in regulatory_hg38_ranges.items():
                        for gene_range in v[1:]:
                            if(chrom[3:] == v[0]):
                                if((l_boundary > gene_range[0] and l_boundary < gene_range[1]) or (u_boundary > gene_range[0] and u_boundary < gene_range[1])
                                   or (l_boundary < gene_range[0] and u_boundary > gene_range[1])):
                                    interval_annotated = True
                                    if(k in annotated_intervals_reg_hg38.keys()):
                                        annotated_intervals_reg_hg38[k].append((chrom, l_boundary, u_boundary))
                                    else:
                                        annotated_intervals_reg_hg38[k] = [(chrom, l_boundary, u_boundary)]
                                        
                    #if(not interval_annotated):
                        #print(chrom, l_boundary, u_boundary)
                        #raise ValueError("Could not map the interval to a gene or a regulatory element.")
                    intervals_already_annotated.append((chrom, l_boundary, u_boundary))
            
        except FileNotFoundError:
            j += 1
            continue
                            
        j += 1
    i += 1

# Now that we have mapped the each gene/reg element to all the overlapping interacting intervals
# we need to use those interacting intervals to select 1/8th of the genes length, such 
# as the selected 1/8th of the gene is centered around each interacting interval evenly.
# For regulatory element, simply select the entire length of the element, since they are short. 

# So the algorithm for the genes using an example:
# Gene ABC: chr1 1 100,000 ABC
# Two interacting intervals, one at 25,000-25,100 and one at 60,000-60,500
# So need to select 1/8 of 100,000 = 12,500 bases centered on these two regions.
# That is 6,250 bases for each intervals. 
# The mean of interval one is 25,050 and for interval two it is 60,250.
# Therefore, the selected portions are (25,050-3,125; 25,050+3,125) and (60,250-3,125; 60,250+3,125)
# The only other possible case of interest is where the interval extends past the end or beginning of the gene.
# Example: for same gene imagine there is interval (99,800 - 101,000) and it is the only interval.
# Then we simply select (100,000-12,500; 100,000) as our included gene portion. 

# Another note, if gene has multiple alternative loci, use the loci within which the intervals 
# are located. 

# Change from 1/8th to 1/10th of gene length due to larger number of genes that will be 
# included in the analysis. 

input_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/REMMA_Results/Combined_Fixed2_Set_Biclustering/'
fname = input_dir + 'biclustering_results_hg38_all_annotated_intervals.txt'
gene_to_interval_centers = dict()
gene_to_interval_boundaries = dict()
annotated_reg_elements = set()
with open(fname) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter='\t')
    for row in csv_reader:
        gene_reg_name = row[0]
        if('GH' in gene_reg_name and len(gene_reg_name) == 11):
            annotated_reg_elements.add(gene_reg_name)
            continue
        intervals = row[1]
        num_intervals = intervals.count('(')
        # now we need to compute the middle of each interval
        intervals_split = intervals.split(',')
        i = 0 
        interval_centers = []
        interval_boundaries = []
        while i < num_intervals:
            l_boundary = int(intervals_split[i*3 + 1].strip(' )]'))
            u_boundary = int(intervals_split[i*3 + 2].strip(' )]'))
            center = int((l_boundary + u_boundary)/2)
            interval_centers.append(center)
            interval_boundaries.append((l_boundary,u_boundary))
            i += 1
        gene_to_interval_centers[gene_reg_name] = interval_centers
        gene_to_interval_boundaries[gene_reg_name] = interval_boundaries
        
        
#test_dict = dict()
#test_dict['C2orf42'] = gene_to_interval_centers['C2orf42']
#print('C2orf42 loci are: ' + str(genes_hg38_ranges['C2orf42']))
#print('C2orf42 annotated interval boundaries are: ' + str(gene_to_interval_boundaries['C2orf42']))
#print('C2orf42 annotated interval centers are: ' + str(gene_to_interval_centers['C2orf42']))

selected_gene_ranges = dict()
for k,v in gene_to_interval_centers.items():
    gene_loci = genes_hg38_ranges[k]
    chrom = gene_loci[0]
    gene_ranges = []
    num_interval_centers = len(v)
    i = 0
    interval_boundaries = gene_to_interval_boundaries[k]
    for interval_center in v:
        boundaries = interval_boundaries[i]
        for alt_loci in gene_loci[1:]:
            if(interval_center >= alt_loci[0] and interval_center <= alt_loci[1]):
                gene_length = alt_loci[1] - alt_loci[0]
                div_factor = 1.5 * num_interval_centers
                if( ((interval_center - (gene_length/div_factor)) >= alt_loci[0]) and
                   ((interval_center + (gene_length/div_factor)) <= alt_loci[1]) ):
                    gene_ranges.append((interval_center - int(gene_length/div_factor), interval_center + int(gene_length/div_factor)))
                elif( (interval_center - (gene_length/div_factor)) < alt_loci[0] ):
                    gene_ranges.append((alt_loci[0], alt_loci[0] + int(gene_length/div_factor)))
                elif( (interval_center + (gene_length/div_factor)) > alt_loci[1] ):
                    gene_ranges.append((alt_loci[1] - int(gene_length/div_factor), alt_loci[1]))
                else:
                    raise ValueError()
            # if lower boundary is within the gene, but the center isnt, that means the upper boundary is 
            # to the right of the gene, so select the right portion of the gene
            elif(boundaries[0] >= alt_loci[0] and boundaries[0] <= alt_loci[1]):
                gene_ranges.append((alt_loci[1] - int(gene_length/div_factor), alt_loci[1]))
            # if upper boundary is within the gene, but the center isnt, that means the lower boundary is 
            # to the left of the gene, so select the left portion of the gene
            elif(boundaries[1] >= alt_loci[0] and boundaries[1] <= alt_loci[1]):
                gene_ranges.append((alt_loci[0], alt_loci[0] + int(gene_length/div_factor)))
            # case where the interval is bigger than the gene, in that case include the entire gene
            elif(boundaries[0] < alt_loci[0] and boundaries[1] > alt_loci[1]):
                gene_ranges.append((alt_loci[0], alt_loci[1]))
            else:
                pass
                #print(k)
                #print(alt_loci)
                #print(interval_center)
                #print(interval_boundaries)
                
        i += 1
    selected_gene_ranges[k] = gene_ranges

# Now we are ready to write out the SNP-ranges, for genes we mapped each gene to the portions 
# of that gene that we want to include in the analysis. If a gene maps to multiple gene ranges,
# use the _alt# suffix with the gene name for each entry.
# Ex: gene ABCD maps to 4 intervals 
# chr1 100 1000 ABCD
# chr1 12000 13000 ABCD_ALT
# chr1 14000 15000 ABCD_ALT2
# chr1 22000 22500 ABCD_ALT3


gene_ranges_to_write = []
for k,v in selected_gene_ranges.items():
    chrom = genes_hg38_ranges[k][0]
    i = 0
    for gene_range in v:
        gene_name = k
        if(i == 1):
            gene_name = k + '_ALT'
        elif(i > 1):
            gene_name = k + '_ALT' + str(i)
        gene_ranges_to_write.append( (chrom, gene_range[0], gene_range[1], gene_name) )
        i += 1
        
for reg_element in annotated_reg_elements:
    chrom = regulatory_hg38_ranges[reg_element][0]
    loci = regulatory_hg38_ranges[reg_element][1]
    gene_ranges_to_write.append( (chrom, loci[0], loci[1], reg_element) )
    
output_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/REMMA_Results/Combined_Fixed2_Set_Biclustering/'
output_fname = output_dir + 'biclustering_results_hg38_verylarge_interacting_intervals_AI_AUD.txt'

# Some of the gene ranges in this file overlap, if that is an issue for plink will need to 
# merge them, there is also a small number of duplicate regions
with open(output_fname, 'w') as writer:
    for entry in gene_ranges_to_write:
        writer.write(str(entry[0])+' '+str(entry[1])+' '+str(entry[2])+' '+str(entry[3])+'\n')


'''
hg19_genes = set()
fname = input_dir + 'biclustering_genes.txt'
with open(fname) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter='\t')
    for row in csv_reader:   
        hg19_genes.add(row[0])
        
print(len(annotated_genes))
print(len(hg19_genes))
print(len(annotated_genes & hg19_genes))
print(hg19_genes - annotated_genes)
print(annotated_genes - hg19_genes)
'''

'''
output_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/REMMA_Results/Combined_Fixed2_Set_Biclustering/'
output_fname = output_dir + 'biclustering_results_hg38_all_annotated_intervals.txt'

with open(output_fname, 'w') as writer:
    for k,v in annotated_intervals_genes_hg38.items():
        writer.write(k + '\t' + str(v) + '\n')
    for k,v in annotated_intervals_reg_hg38.items():
        writer.write(k + '\t' + str(v) + '\n')
'''
       
'''
# Gather all the interacting genes/regulatory elements from REMMA/Biclustering analysis
interacting_elements = set()
input_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/REMMA_Results/Combined_Fixed2_Set_Biclustering/'
i = 1 
while i <= 22:
    j = i
    while j <= 22:
        #print(i)
        #print(j)
        try:
            fname = input_dir + 'biclustering_results_chr' + str(i) + '_chr' + str(j) + '_annotated.txt'
            with open(fname) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter='\t')
                for row in csv_reader:
                    # If the file is empty move onto the next file
                    if(len(row) == 0):
                        break
                    
                    interacting_genes = row[0]
                    # print(interacting_genes)
                    recording_gene_name = False
                    current_gene = ''
                    for character in interacting_genes:
                        if(character == "'"):
                            if(recording_gene_name):
                                current_gene_stripped = current_gene.strip()
                                interacting_elements.add(current_gene_stripped)
                                current_gene = ''
                            recording_gene_name = not recording_gene_name
                        else:
                            if(recording_gene_name):
                                current_gene += character
            #print(interacting_elements)
        except FileNotFoundError:
            j += 1
            continue
                            
        j += 1
    i += 1
              
def parse_gene_regulatory_element_name(element:str, interacting_elements_cleaned:set):
    if('_ALT' in element):
        starting_index = element.find('_ALT')
        #print(element)
        cleaned_element = element[0:starting_index]
        interacting_elements_cleaned.add(cleaned_element)
    else:
        interacting_elements_cleaned.add(element)
        
    return interacting_elements_cleaned
    
interacting_elements_cleaned = set()
for element in interacting_elements:
    
    multiple_elements = element.split(' ')
    if(len(multiple_elements) > 1):
        for ele in multiple_elements:
            interacting_elements_cleaned = parse_gene_regulatory_element_name(ele, interacting_elements_cleaned)
    else:
        interacting_elements_cleaned = parse_gene_regulatory_element_name(element, interacting_elements_cleaned)
        
output_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/REMMA_Results/Combined_Fixed2_Set_Biclustering/'
output_fname = output_dir + 'biclustering_results_hg38_AI_AUD.txt'
'''
'''
with open(output_fname, 'w') as writer:
    for element in interacting_elements_cleaned:
        try:
            element_loci = genes_hg38_ranges[element]
        except KeyError:
            try:
                element_loci = regulatory_hg38_ranges[element]
            except KeyError:
                print("Element: " + element + " not found in ncbi or genehancer data.")
        writer.write(str(element_loci[0]) + ' ' + str(element_loci[1]) + ' ' + str(element_loci[2]) + ' ' + element + '\n')
'''

'''
# Some genes/regulatory elements need to be lifted over manually, as they are absent 
# in the hg38 NCBI/Genehancer data files
combined_set_expanded_file = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/AUD_Resources/Start_Sets/Expanded_Ranges/'
combined_set_expanded_file += 'Combined_Set_score900_db001_and_exp001_1_and_genehancer_score25_ranges.txt'
combined_set_expanded_elements = []

with open(combined_set_expanded_file) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=' ')
    for row in csv_reader:
        combined_set_expanded_elements.append(row[3])
      
combined_set_expanded_file_hg38 = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/AUD_Resources/Start_Sets/Expanded_Ranges/'
combined_set_expanded_file_hg38 += 'Combined_Set_score900_db001_and_exp001_1_and_genehancer_score25_ranges_hg38.txt'

with open(combined_set_expanded_file_hg38, 'w') as writer:
    for element in combined_set_expanded_elements:
        try:
            element_loci = genes_hg38_ranges[element]
        except KeyError:
            try:
                element_loci = regulatory_hg38_ranges[element]
            except KeyError:
                print("Element: " + element + " not found in ncbi or genehancer data.")
    writer.write(str(element_loci[0]) + ' ' + str(element_loci[1]) + ' ' + str(element_loci[2]) + ' ' + element + '\n')
'''

def map_snp_to_gene(chrom:str, loci:int, gene_database:dict):
    gene = ''
    for k,v in gene_database.items():
        gene_chrom = v[0]
        if(chrom == gene_chrom):
            for gene_range in v[1:]:
                if(loci > gene_range[0] and loci < gene_range[1]):
                    gene = k
                    return gene
    return gene
'''
input_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_AI_AUD/AoU_Plink_Replication/'
fname = input_dir + 'AI_AUD_Interacting_Intervals_Epi_Upper_Boundary_ACAF.tsv'
gene_pairs = set()

with open(fname) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter='\t')
    for row in csv_reader:
        snp1_info = row[0].split(':')
        chr1,snp1_loci = snp1_info[0],snp1_info[1]
        snp2_info = row[1].split(':')
        chr2,snp2_loci = snp2_info[0],snp2_info[1]
        p_val = float(row[2])
        
        gene1 = map_snp_to_gene(chr1, int(snp1_loci), genes_hg38_ranges)
        gene2 = map_snp_to_gene(chr2, int(snp2_loci), genes_hg38_ranges)
        if( (gene2, gene1) not in gene_pairs):
            gene_pairs.add( (gene1, gene2, p_val) )
            
sorted_results = sorted(gene_pairs, key=lambda pair: pair[2])

output_fname = input_dir + 'AI_AUD_Interacting_Intervals_Epi_Upper_Boundary_ACAF_Genes.tsv'
with open(output_fname, 'w') as writer:
    for gene_pair in sorted_results:
        writer.write(gene_pair[0] + '\t' + gene_pair[1] + '\t' + str(gene_pair[2]) + '\n')
'''            