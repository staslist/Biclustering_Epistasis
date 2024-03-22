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

def generate_multicpu_files(total_markers:int, num_cpus:int, out_dir:str):
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
        fname_out_py = out_dir + 'BiClustering_Launcher' + i_s + '_' + i_e + '.py'
        with open(fname_out_py, 'w') as writer:
            writer.write('from BiClustering import *\n')
            writer.write('if __name__ == "__main__":\n')
            writer.write("\tmarker_file = '/gpfs/group/home/slistopad/BiClustering/biclustering_pub_marker_pairs.csv'\n")
            writer.write('\tN,n,inter_matrix,up_tri = initialize_matrices(' + tot_marks + ', marker_file)\n')
            writer.write("\tout_dir = '/gpfs/group/home/slistopad/BiClustering/'\n")
            writer.write('\ti_start, i_end = '+i_s+', '+i_e+'\n')
            writer.write('\tkm_results = compute_k_m_parallel('+ tot_marks + ',60, inter_matrix, up_tri, i_start, i_end)\n')
            writer.write('\tcompute_interval_pval_parallel(km_results, N, n, i_start, i_end, out_dir)\n')
            
            if(counter == 0):
                writer.write('\tpval_results = []\n')
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
                writer.write("\t\tfilename = out_dir + 'pval_results_' + str(i_pair[0]) + '_' + str(i_pair[1]) + '.csv'\n")
                writer.write('\t\twith open(filename) as csv_file:\n')
                writer.write("\t\t\tcsv_reader = csv.reader(csv_file, delimiter=',')\n")
                writer.write("\t\t\tfor row in csv_reader:\n")
                writer.write("\t\t\t\tpval_results[(row[0],row[1],row[2],row[3])] = row[6]\n")
                        
                writer.write("\ts_inter_pairs = sorted(pval_results.items(), key = lambda x: abs(x[1]), reverse = False)\n")            
                writer.write("\ttrimmed_inter_pairs = trim_intervals(s_inter_pairs)\n")
                writer.write("\tfilename = out_dir + 'biclustering_results.txt'\n")
                writer.write("\twith open(filename, 'w') as writer:\n")
                writer.write("\t\tfor inter_pair in trimmed_inter_pairs:\n")
                writer.write("\t\t\twriter.write(str(inter_pair[0][0])+','+str(inter_pair[0][1])+','+str(inter_pair[0][2])+','+str(inter_pair[0][3])+','+str(inter_pair[1]))\n")
        
        fname_out_sh = out_dir + 'BiClustering_Launcher' + i_s + '_' + i_e + '.sh'
        with open(fname_out_sh, 'w') as writer:
            writer.write("#!/bin/sh\n")
            writer.write("#SBATCH --job-name=SL_BiClustering_Launcher" + i_s + "_" + i_e + "\n")
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
            writer.write("python3 " + 'BiClustering_Launcher' + i_s + '_' + i_e + '.py')
        counter += 1
         
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
    # Note if two matrices are identical we treat them as non-overlapping. 
    # Note that to overlap one matrix's origin point must be INSIDE of another matrix (not on border).
    if(i == i2 and j == j2 and a==a2 and b==b2):
        return False
    else:
        return ( ((i+a) > i2 > i) and ((j+b) > j2 > j) ) or ( ((i2+a2) > i > i2) and ((j2+b2) > j > j2) )

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

def compute_k_m_parallel(total_markers:int, max_inter_length:int, inter_matrix:dict,
                         upper_triangular:bool, i_start:int, i_end:int):
    
    # Our only knowledge about inter_matrix is that it consists of 1s and 0s 
    # and that 1s are typically very sparce.
    
    # We have stronger assumptions about tested_inter_matrix which we deem to be either 
    # upper trinagular (with 1s and 0s) or completely filled with 1s. 
    # Given this we do not need to actually parse tested_inter_matrix to compute m.
    
    # Convert inter_matrix and tested_inter_matrix into numpy arrays.
    inter_array = convert_dict_to_2D_numpy_array(inter_matrix)
    
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
        while j < total_markers:
            # If we the biggest matrix at (i,j) has k = 0, then all of the matrices within it 
            # will also have k = 0. 
            k_is_zero = False
            # a_thresh and b_thresh represent biggest submatrix that is all zeros
            a_thresh = 0
            b_thresh = 0
            while a >= 1:
                while b >= 1:
                    if((i + a) <= total_markers and (j+b) <= total_markers):
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
                        km_results[(i,j,a,b)] = (k,m)
                        
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
                                   out_dir:str):
    filename = out_dir + 'pval_results_' + str(i_start) + '_' + str(i_end) + '.csv'
    with open(filename, 'w') as writer:
    
        pval_results = dict()
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
            
            if(pval_results[key] != 1.0):
                # Do not write out the p-value for this interval pair if its 1.0
                # Assumming that number of interacting marker pairs is sparce this drastically 
                # cuts down on number of writes. 
                writer.write(str(key[0])+','+str(key[1])+','+str(key[2])+','+str(key[3])+',')
                writer.write(str(val[0]) + ',' + str(m) + ',' + str(pval_results[key]) + '\n')
    
    return

def trim_intervals(sorted_intervals:list):
    # Now go through interval pairs starting from most significant, and remove all overlapping 
    # intervals for the current interval
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
    while i < len(sorted_intervals):
        if(i not in indeces_to_delete):
            culled_inter_pairs.append(sorted_intervals[i])
        i += 1
        
    return culled_inter_pairs

class TestBiClusterCodeBase(unittest.TestCase):  
    
    def test_check_overlap(self):
        self.assertFalse(check_overlap(0,0,2,2,0,0,2,2))
        self.assertTrue(check_overlap(0,0,2,2,1,1,2,2))
        self.assertFalse(check_overlap(0,0,2,2,1,2,1,1))
        self.assertFalse(check_overlap(0,0,2,2,2,2,1,1))
        self.assertFalse(check_overlap(0,0,2,2,3,2,1,1))
        
        self.assertTrue(check_overlap(2,2,3,6,1,1,3,3))
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
        
        km_results = compute_k_m_parallel(10, 4, inter_matrix, True, 5, 6)

        # print(km_results)
        
        self.assertEqual(km_results[(5,0,4,4)], (0,0))
        self.assertEqual(km_results[(5,3,4,4)], (0,1))
        self.assertEqual(km_results[(5,6,2,4)], (2,7))
        self.assertEqual(km_results[(5,8,1,1)], (1,1))
        self.assertEqual(km_results[(5,8,1,2)], (2,2))
        self.assertEqual(km_results[(5,8,2,1)], (1,2))
        
     
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
    
    runner = unittest.TextTestRunner()
    runner.run(unit_test_suite)
    

