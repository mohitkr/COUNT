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
def plot(sam):
    tag="v22_"+sam
    directory=home+'/COUNT/code/data/'+tag
    ax=np.zeros(5)
    y_recall=np.zeros(5)
    yerr_recall=np.zeros(5)
    y_precision=np.zeros(5)
    yerr_precision=np.zeros(5)
    times=np.zeros(5)
    times_err=np.zeros(5)
    i=0
    with open(directory+ "/results.csv") as f:
        d = {}
        next(f)    # skip first header line
        for line in f:
    #        print(line.split(','))
            nurses,numSol,numSam,TP1,TP1_err,TP2,TP2_err,FP,FP_err,FN,FN_err,time,time_err = line.split(',')
            ax[i]=float(numSol)
            y_recall[i] = (float(TP2)*100)/int(numSam)
            yerr_recall[i] = (float(TP2_err)*100)/int(numSam)
            y_precision[i] = (float(TP1)*100)/int(numSam)
            yerr_precision[i] = (float(TP1_err)*100)/int(numSam)
            times[i] = float(time)
            times_err[i] = float(time_err)
            i+=1
    
    #print(ax)
    print(y_recall)
    print(yerr_recall)
    #ax=ax[0:5]
    plt.figure(1)
    plt.style.use('ggplot')
    plt.plot(ax,y_recall)
    plt.fill_between(ax, y_recall - yerr_recall, y_recall + yerr_recall, linewidth=0, alpha=0.35)
    plt.xlabel('Number of Examples')
    plt.ylabel('Recall')
    plt.axis([0, 200, 0, 100])
    
    plt.figure(2)
    plt.style.use('ggplot')
    plt.plot(ax,y_precision)
    plt.fill_between(ax, y_precision - yerr_precision, y_precision + yerr_precision, linewidth=0, alpha=0.35)
    plt.xlabel('Number of Examples')
    plt.ylabel('Precision')
    plt.axis([0, 200, 0, 100])
    
    plt.figure(3)
    plt.style.use('ggplot')
    plt.plot(ax,times)
    plt.fill_between(ax, times - times_err, times + times_err, linewidth=0, alpha=0.35)
    plt.xlabel('Number of Examples')
    plt.ylabel('Time Taken (in sec)')

plot("10000_0")
plot("10000_1")
plot("10000_2")


directory=home+'/COUNT/code/data/'
plt.figure(1)
plt.legend(['Small','Medium','Large'], loc='lower right', prop={'size':15})
plt.savefig(directory+'/recall.png', bbox_inches='tight', pad_inches=0, dpi=120)
#plt.show()


plt.figure(2)
plt.legend(['Small','Medium','Large'], loc='lower right', prop={'size':15})
plt.savefig(directory+'/precision.png', bbox_inches='tight', pad_inches=0, dpi=120)
#plt.show()

plt.figure(3)
plt.legend(['Small','Medium','Large'], loc='upper left', prop={'size':15})
plt.savefig(directory+'/time.png', bbox_inches='tight', pad_inches=0, dpi=120)
plt.show()
























