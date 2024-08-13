# -*- coding: utf-8 -*-
"""
Created on Sat Jun 22 03:22:53 2024

@author: staslist
"""

import os

current_dir = '/gpfs/group/home/slistopad/BiClustering/'

all_names = os.listdir()
for name in all_names:
    #print(name)
    #print(name[-3:])
    if(name[-3:] == 'out' or name[-3:] == 'csv'):
        path = os.path.join(current_dir, name)
        os.remove(path)