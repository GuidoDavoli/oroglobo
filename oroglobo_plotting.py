#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Guido Davoli - CNR ISAC

"""

import matplotlib.pyplot as plt
import numpy as np

def orography_plot(data,img,dpi):
    
    plt.figure()
    cmap = plt.get_cmap('terrain')
    cmap.set_bad('white') # points=nan will be displayed in this color
    dataplot = np.ma.masked_equal(data, 0) # set to nan where data=0
    plt.imshow(dataplot,origin='lower',cmap=cmap,vmin=0,vmax=8500)
    plt.colorbar(location='bottom')
    plt.savefig(img,dpi=dpi)
    plt.show()
    
    return
    
    