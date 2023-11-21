#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 12:21:48 2023

@author: guidodavoli
"""

import os
import oroglobo_parameters as oropar


# if output folders do not exist, create them

outfolders=oropar.paths_out

for x in outfolders:
    
    outfolder=outfolders[x]
    if not os.path.exists(outfolder):
        
        os.makedirs(outfolder)
        
# if working folders do not exist, create them

outfolders=oropar.paths_work

for x in outfolders:
    
    outfolder=outfolders[x]
    if not os.path.exists(outfolder):
        
        os.makedirs(outfolder)