#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 15:42:02 2018

@author: mohit
"""

import matplotlib.pyplot as plt
import numpy as np
from os.path import expanduser
home = expanduser("~")
#import matplotlib.colors as colors
#import matplotlib.cm as cm
sam="100k"
tag="v1_"+sam
directory=home+'/CLUT/code/data/'+tag

with open(directory+ "/results.csv") as f:
    d = {}
    next(f)    # skip first header line
    for line in f:
#        print(line.split(','))
        nurses,numSol,numSam,TP2,FN,trConsReject,time = line.split(',')
        d[numSol] = d.get(numSol, []) + [(int(TP2)*100/int(numSam))]

# sum and average of column Value by ID
ax=np.zeros(len(d))
y_recall=np.zeros(len(d))
yerr_recall=np.zeros(len(d))
i=0
for id_ in d:
    ax[i]=int(id_)
    average = sum(d[id_]) / float(len(d[id_]))
    stdev = np.std([d[id_]])/np.sqrt(len(d))
    y_recall[i]=average
    yerr_recall[i]=stdev
    
    
    i+=1
    print('{}: avg = {:.2f}, std = {:.2f}'.format(id_, average, stdev))
#print(ax)
#print(y_recall)
ax=ax[0:5]
y_recall=y_recall[0:5]
yerr_recall=yerr_recall[0:5]
plt.figure(1)
plt.plot(ax,y_recall)
plt.fill_between(ax, y_recall - yerr_recall, y_recall + yerr_recall, linewidth=0, alpha=0.35)
plt.xlabel('Number of Solutions')
plt.ylabel('Recall')
plt.title("Recall vs # of Sol for "+sam+" samples")
plt.savefig(directory+'/'+tag+'_recall_first5.png', bbox_inches='tight', pad_inches=0, dpi=120)
plt.show()

#plt.show()


























