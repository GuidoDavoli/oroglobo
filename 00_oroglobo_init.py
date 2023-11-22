#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 12:21:48 2023

@author: guidodavoli
"""

import os
import oroglobo_parameters as oropar


gridname=oropar.model_grid["GRIDNAME"]

# if output folders do not exist, create them

folders=oropar.paths_out

for x in folders:
    
    folder=folders[x].replace("*GRIDNAME*", gridname) 
    if not os.path.exists(folder):
        
        os.makedirs(folder)
        
# if working folders do not exist, create them

folders=oropar.paths_work

for x in folders:
    
    folder=folders[x].replace("*GRIDNAME*", gridname)
    if not os.path.exists(folder):
        
        os.makedirs(folder)