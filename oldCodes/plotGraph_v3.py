#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 15:32:35 2018

@author: mohit
"""

import matplotlib.pyplot as plt
import numpy as np
#import matplotlib.colors as colors
#import matplotlib.cm as cm
sam="100k"
tag="v3_"+sam
directory='/home/mohit/CLUT/code/data/'+tag

with open(directory+ "/results.csv") as f:
    d = {}
    next(f)    # skip first header line
    for line in f:
#        print(line.split(','))
        nurses,numSol,numSam,TP1,TP2,FP,FN,leConsReject,trConsReject,realObj,learnObj,time = line.split(',')
        print(nurses,numSol,numSam,TP1,TP2,FP,FN,leConsReject,trConsReject,realObj,learnObj,time)
        d[numSol] = d.get(numSol, []) + [(int(TP1)*100/int(numSam), int(TP2)*100/int(numSam), (int(float(learnObj))-3)*10)]

# sum and average of column Value by ID
ax=np.zeros(len(d))
y_recall=np.zeros(len(d))
yerr_recall=np.zeros(len(d))
y_prec=np.zeros(len(d))
yerr_prec=np.zeros(len(d))
y_dist=np.zeros(len(d))
yerr_dist=np.zeros(len(d))
i=0
for id_ in d:
    ax[i]=int(id_)
    average = sum(t[1] for t in d[id_]) / float(len(d[id_]))
    stdev = np.std([t[1] for t in d[id_]])/np.sqrt(len(d))
    y_recall[i]=average
    yerr_recall[i]=stdev
    
    average = sum(t[0] for t in d[id_]) / float(len(d[id_]))
    stdev = np.std([t[0] for t in d[id_]])/np.sqrt(len(d))
    y_prec[i]=average
    yerr_prec[i]=stdev
    
    average = sum(t[2] for t in d[id_]) / float(len(d[id_]))
    stdev = np.std([t[2] for t in d[id_]])/np.sqrt(len(d))
    y_dist[i]=average
    yerr_dist[i]=stdev
    
    i+=1
    print('{}: avg = {:.2f}, std = {:.2f}'.format(id_, average, stdev))

plt.figure(1)
plt.plot(ax,y_recall)
plt.fill_between(ax, y_recall - yerr_recall, y_recall + yerr_recall, linewidth=0, alpha=0.35)
plt.xlabel('Number of Solutions')
plt.ylabel('Recall')
plt.title("Recall vs # of Sol for "+sam+" samples")
plt.savefig(directory+'/'+tag+'_recall.png', bbox_inches='tight', pad_inches=0, dpi=120)
#plt.show()
plt.figure(2)
plt.plot(ax,y_prec)
plt.fill_between(ax, y_prec - yerr_prec, y_prec + yerr_prec, linewidth=0, alpha=0.35)
#plt.axis([0, 51, 0, 102])
plt.xlabel('Number of Solutions')
plt.ylabel('Precision')
plt.title("Precision vs # of Sol for "+sam+" samples")
plt.savefig(directory+'/'+tag+'_precision.png', bbox_inches='tight', pad_inches=0, dpi=120)
#plt.show()

plt.figure(3)
plt.plot(ax,y_dist)
plt.fill_between(ax, y_dist - yerr_dist, y_dist + yerr_dist, linewidth=0, alpha=0.35)
plt.axis([0, 200, 0, 50])
plt.xlabel('Number of Solutions')
plt.ylabel('% Distance')
plt.title("% Distance vs # of Sol for "+sam+" samples")
plt.savefig(directory+'/'+tag+'_distance.png', bbox_inches='tight', pad_inches=0, dpi=120)