#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 11:16:20 2018

@author: mohit
"""
import generateSol as gs
import generate_v1 as gfl
import readData as rd
import numpy as np
import glob
import os
import csv
from scipy import stats
import time
import sys
import shutil
from os.path import expanduser
home = expanduser("~")

numSam=int(sys.argv[1])
numFiles=numSam

hospSize=sys.argv[2]

tag="v1_"+str(numSam)+"_"+hospSize

directory=home+'/COUNT/code/data/'+tag
#directory='/home/mohit/CLUT/code/data/'+tag
num_constrType=12
constrList=[[(0,),(1,)],[(0,),(2,)],[(0,),(1,2)],[(1,),(0,)],[(1,),(2,)],[(1,),(0,2)],[(2,),(0,)],[(2,),(1,)],[(2,),(0,1)],[(0,1),(2,)],[(0,2),(1,)],[(1,2),(0,)]]
if int(hospSize)==0:
    num_nurses=10
    num_days=28
    num_shifts=4   
    orderingNotImp=[2]
    num_constr=6
    bounds=np.zeros([num_constrType,num_constr])
    bounds[2,0]=4
    bounds[2,1]=6
    bounds[6,0]=9
    bounds[6,1]=16
    bounds[6,2]=1
    bounds[6,3]=7
    bounds[6,4]=2
    bounds[6,5]=8   
    bounds[9,0]=1
    bounds[9,1]=2
    bounds[10,1]=1
#    bounds[11,4]=1
#    bounds[11,5]=4

if int(hospSize)==1:
    num_nurses=31
    num_days=28
    num_shifts=4   
    orderingNotImp=[2]
    num_constr=6
    bounds=np.zeros([num_constrType,num_constr])
    bounds[2,0]=16
    bounds[2,1]=24
    bounds[6,0]=6
    bounds[6,1]=16
    bounds[6,2]=1
    bounds[6,3]=7
    bounds[6,4]=2
    bounds[6,5]=8   
    bounds[9,0]=2
    bounds[9,1]=9
    bounds[10,1]=1
#    bounds[11,4]=1
#    bounds[11,5]=4
    
if int(hospSize)==2:
    num_nurses=49
    num_days=28
    num_shifts=4   
    orderingNotImp=[2]
    num_constr=6
    bounds=np.zeros([num_constrType,num_constr])
    bounds[2,0]=20
    bounds[2,1]=29
    bounds[6,0]=6
    bounds[6,1]=16
    bounds[6,2]=1
    bounds[6,3]=5
    bounds[6,4]=1
    bounds[6,5]=7   
    bounds[9,0]=1
    bounds[9,1]=8
    bounds[10,1]=1
#    bounds[11,4]=1
#    bounds[11,5]=4
bounds=bounds.astype(np.int64)

if not os.path.exists(directory):
    os.makedirs(directory)
if not os.path.exists(directory+"/solutions"):
    os.makedirs(directory+"/solutions")
    
my_csv = open(directory+ "/results.csv" ,"w+")
detail_csv = open(directory+ "/detail_results.csv" ,"w+")
csvWriter = csv.writer(my_csv,delimiter=',')
#row=['Nurses','# of Sol', '# of Sample', 'TP1', 'TP2','FP','FN','LeConsReject', 'TrConsReject', 'Time Taken']
row=['Nurses','# of Sol', '# of Sample', 'TP2', 'TP2_err','FN','FN_err', 'Time Taken', 'Time Taken_err']
csvWriter.writerow(row)  

#for num_nurses in range(6,16,3):
dataDir=directory+"/solutions"+"/*N"+str(num_nurses)+"*.csv"
for fl in glob.glob(dataDir): 
      os.remove(fl)       
gfl.generateSolution(num_nurses,num_days,num_shifts,orderingNotImp,numFiles,directory+"/solutions",constrList,bounds)
bounds_prev=np.zeros([num_constrType,num_constr])
for numSol in [10]:
    print("########################## Number of Nurses:",num_nurses," NumSol: ",numSol," ##########################")
    try:
        os.remove(directory+'/result_N'+str(num_nurses)+'.csv')
    except OSError:
        pass
    numSeed=5
    tot_tp=np.zeros(numSeed)
    tot_fn=np.zeros(numSeed)
    tot_time=np.zeros(numSeed)
    for seed in range(numSeed):
        detail_csv = open(directory+ "/detail_results.csv" ,"a")
        csvWriter_detail = csv.writer(detail_csv,delimiter=',')
        try:
            os.remove(directory+'/result_N'+str(num_nurses)+'_seed'+str(seed)+'.csv')
        except OSError:
            pass
        start=time.clock()
        rd.learnConstraints(dataDir,numSol,num_nurses,1,directory,numFiles,seed)
        end=time.clock()
        data=rd.readCSV(directory+'/result_N'+str(num_nurses)+'_seed'+str(seed)+'.csv')

        data_transpose=list(zip(*data))
        data_int=np.zeros([len(data_transpose),len(data_transpose[0])-1])
        for i in range(len(data_transpose)):
            for j in range(1,len(data_transpose[i])):
                if data_transpose[i][j]!='':
                    data_int[i,j-1]=int(data_transpose[i][j])
        
        bounds_learned=np.zeros([num_constrType,num_constr])
        k=0
        
        for i in range(len(data_transpose)):
            if (i+1)%7 != 0:
                if (k%6)%2==0:
                    if stats.mode(data_int[i], axis=None)[0][0] == 0:
                        bounds_learned[int(k/6),k%6]=0
                    else:
                        bounds_learned[int(k/6),k%6]=np.min(data_int[i])
                if (k%6)%2!=0:
                    if stats.mode(data_int[i], axis=None)[0][0] == 0:
                        bounds_learned[int(k/6),k%6]=0
                    else:
                        bounds_learned[int(k/6),k%6]=np.max(data_int[i])
                k+=1
        bounds_learned=bounds_learned.astype(np.int64)
#        print(bounds_learned)
        if not np.array_equal(bounds_learned,bounds_prev):
            totSam2,fn,TrConsReject,constrRej3,constrRej4=gfl.generateSample(num_nurses,num_days,num_shifts,orderingNotImp,numSam,num_constrType,constrList,bounds,bounds_learned)
        else:
            totSam2,fn,TrConsReject,constrRej3,constrRej4=totSam2_prev,fn_prev,TrConsReject_prev,constrRej3_prev,constrRej4_prev
        totSam2_prev,fn_prev,TrConsReject_prev,constrRej3_prev,constrRej4_prev=totSam2,fn,TrConsReject,constrRej3,constrRej4
        bounds_prev=bounds_learned
        
        row=[]
        row.extend([num_nurses])
        row.extend([numSol])
        row.extend([numSam])
        row.extend([totSam2])
        row.extend([fn])
        row.extend([TrConsReject])
        row.extend([end-start])
        csvWriter_detail.writerow(row)
        detail_csv.close()
        print(fn)
        tot_tp[seed]=totSam2
        tot_fn[seed]=fn
        tot_time[seed]=(end-start)
    
    row=[]
    row.extend([num_nurses])
    row.extend([numSol])
    row.extend([numSam])
    row.extend([sum(tot_tp)/numSeed])
    row.extend([np.std(tot_tp)/np.sqrt(numSeed)])
    row.extend([sum(tot_fn)/numSeed])
    row.extend([np.std(tot_fn)/np.sqrt(numSeed)])
    row.extend([sum(tot_time)/numSeed])
    row.extend([np.std(tot_time)/np.sqrt(numSeed)])
    csvWriter.writerow(row)
    print(row)
    print(constrRej3)
    print(constrRej4)
    
shutil.rmtree(directory+"/solutions")
os.makedirs(directory+"/solutions")
    
my_csv.close()   






