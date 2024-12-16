# -*- coding: utf-8 -*-
"""
Created on Thu May 16 00:15:40 2024

@author: staslist
"""

# Circos Plotting Code
import pycircos
import collections
import matplotlib.pyplot as plt
import numpy as np

Garc    = pycircos.Garc
Gcircle = pycircos.Gcircle
'''
colors = ['darkred', 'red', 'coral', 'peru', 'darkorange', 'orange', 'goldenrod', 'gold', 'olive', 
          'lawngreen', 'forestgreen', 'seagreen', 'lightseagreen', 'teal', 'cyan', 'deepskyblue', 
          'dodgerblue', 'darkblue', 'blueviolet', 'indigo', 'purple', 'magenta', 'deeppink', 'crimson']
'''
colors = ['lightcoral', 'red', 'coral', 'peru', 'darkorange', 'orange', 'goldenrod', 'gold', 'yellow', 
          'lawngreen', 'forestgreen', 'seagreen', 'lightseagreen', 'paleturquoise', 'cyan', 'deepskyblue', 
          'dodgerblue', 'cornflowerblue', 'blueviolet', 'mediumorchid', 'hotpink', 'magenta', 'deeppink', 'crimson']
indir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/Circos/example_notebooks/sample_data/'
#Set chromosomes
circle = Gcircle() 
with open(indir + "example_data_chromosome_general.csv") as f:
    f.readline()
    i = 0
    for line in f:
        line   = line.rstrip().split(",") 
        name   = line[0]
        length = int(line[-1]) 
        arc    = Garc(arc_id=name, size=length, interspace=3, raxis_range=(950,1000), labelposition=60, label_visible=True,
                      facecolor = colors[i])
        circle.add_garc(arc) 
        
        i += 1
circle.set_garcs() 

'''
color_dict = {"gneg":"#FFFFFF00", "gpos25":"#EEEEEE", "gpos50":"#BBBBBB", "gpos75":"#777777",
              "gpos100":"#000000", "gvar":"#FFFFFF00", "stalk":"#C01E27", 
              "acen":"#D82322"}

arcdata_dict = collections.defaultdict(dict)
with open(indir + "example_data_chromosome_cytoband.csv") as f:
    f.readline()
    for line in f:
        line  = line.rstrip().split(",")
        name  = line[0]     
        start = int(line[1])-1 
        width = int(line[2])-(int(line[1])-1) 
        if name not in arcdata_dict:
            arcdata_dict[name]["positions"] = []
            arcdata_dict[name]["widths"]    = [] 
            arcdata_dict[name]["colors"]    = [] 
        arcdata_dict[name]["positions"].append(start) 
        arcdata_dict[name]["widths"].append(width)
        arcdata_dict[name]["colors"].append(color_dict[line[-1]])

for key in arcdata_dict:
    circle.barplot(key, data=[1]*len(arcdata_dict[key]["positions"]), positions=arcdata_dict[key]["positions"], 
                   width=arcdata_dict[key]["widths"], raxis_range=[950,1000], facecolor=arcdata_dict[key]["colors"]) 

circle.figure
'''

#linkplot
indir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/REMMA_Results/Combined_Fixed2_Set_Biclustering/'
values_all   = [] 
arcdata_dict = collections.defaultdict(dict)
with open(indir + "biclustering_results_all.txt") as f:
    #f.readline()
    i = 0
    for line in f:

        line  = line.rstrip().split(",")
        name1  = line[0]     
        start1 = int(line[1])-1
        end1   = int(line[2])
        name2  = line[3]     
        start2 = int(line[4])-1
        end2   = int(line[5])
        source = (name1, start1, end1, 950)
        destination = (name2, start2, end2, 950)
        
        circle.chord_plot(source, destination, facecolor=circle.garc_dict[name1].facecolor,
                          edgecolor=circle.garc_dict[name1].facecolor, linewidth=0.15)
        i += 1
 
    
gene_loci_all = []
indir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch5_NA_Cohort/AUD_Resources/Start_Sets/Expanded_Ranges/'
with open(indir + "Combined_Set_score900_db001_and_exp001_1_and_genehancer_score25_ranges.txt") as f2:
    i = 0
    for line in f2:
        
        line  = line.rstrip().split(" ")
        name = line[0]
        start1 = int(line[1])
        end1 = int(line[2])
        gene_loci_all.append([name,start1,end1])
        
        i += 1
        
#print(gene_loci_all)
    
i = 1
while i <= 22:
    gene_loci_current = []
    for gene_loci in gene_loci_all:
        if(gene_loci[0] == str(i)):
            gene_loci_current.append(gene_loci)
            
           
    if(len(gene_loci_current) == 0):
        i += 1
        continue
    arc_name = 'chr' + str(i)
    data = list()
    positions = list()
    widths = list()        
    for gene_loci in gene_loci_current:
        data.append(1)
        positions.append(gene_loci[1])
        widths.append(gene_loci[2] - gene_loci[1])
        
    
    circle.barplot(arc_name, data, positions, widths, (950,1000),
                   facecolor = 'black', edgecolor = 'black', linewidth = 0.08)
        
    i += 1

circle.figure
plt.savefig('filename.png', dpi=300, bbox_inches="tight")