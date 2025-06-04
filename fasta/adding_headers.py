#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  4 14:10:36 2025
This script was created to add a header to a bunch of sequences. The header
is just the name query followed by a counter
@author: danielcm
"""

import os
path = "/Users/danielcm/Desktop/"
os.chdir(path)
i = 1
f_out = open("mimotopes_headers.fasta","w")

with open("mimotope_info.fasta","r") as file:
    lines = file.readlines()
    for line in lines:
        f_out.write(">query"+str(i)+"\n")
        f_out.write(line)
        i = i+1
        
        
