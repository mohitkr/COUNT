#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 10:52:49 2018

@author: mohit
"""

import generateSol as gs
import generate_v4 as gfl
import readData_bk as rd
import numpy as np
import glob
import os
import csv
from scipy import stats
import random
import time
import sys
import shutil
from os.path import expanduser
home = expanduser("~")

numSam=int(sys.argv[1])
hospSize=sys.argv[2]

numFiles=numSam

tag="v4_"+str(numSam)+"_"+str(hospSize)

directory=home+'/COUNT/code/data/'+tag

num_nurses=15
num_days=7
num_shifts=3   
orderingNotImp=[2]

num_constrType=12
num_constr=6
constrList=[[(0,),(1,)],[(0,),(2,)],[(0,),(1,2)],[(1,),(0,)],[(1,),(2,)],[(1,),(0,2)],[(2,),(0,)],[(2,),(1,)],[(2,),(0,1)],[(0,1),(2,)],[(0,2),(1,)],[(1,2),(0,)]]

bounds=np.zeros([num_constrType,num_constr])
subset1_bounds=np.zeros([num_constrType,num_constr])
subset2_bounds=np.zeros([num_constrType,num_constr])
if int(hospSize)==0:
    num_nurses=10
    num_days=28
    num_shifts=4   
    bounds[2,0]=4
    bounds[2,1]=6
#    bounds[6,0]=9
#    bounds[6,1]=16
    bounds[6,2]=1
    bounds[6,3]=7
    bounds[6,4]=2
#    bounds[6,5]=8   
    bounds[9,0]=1
    bounds[9,1]=2
    bounds[10,1]=1
    
    subset1_bounds[6,0]=9
    subset1_bounds[6,1]=16
    subset1_bounds[6,5]=8   
    
    subset2_bounds[6,0]=4
    subset2_bounds[6,1]=9
    subset2_bounds[6,5]=4 
    
if int(hospSize)==1:
    num_nurses=31
    num_days=28
    num_shifts=4   
    bounds[2,0]=16
    bounds[2,1]=24
#    bounds[6,0]=9
#    bounds[6,1]=16
    bounds[6,2]=1
#    bounds[6,3]=7
    bounds[6,4]=2
#    bounds[6,5]=8   
    bounds[9,0]=2
    bounds[9,1]=9
    bounds[10,1]=1
    
    subset1_bounds[6,0]=4
    subset1_bounds[6,1]=13
    subset1_bounds[6,3]=8
    subset1_bounds[6,5]=4   
    
    subset2_bounds[6,0]=6
    subset2_bounds[6,1]=16
    subset2_bounds[6,3]=7
    subset2_bounds[6,5]=8
    
if int(hospSize)==2:
    num_nurses=49
    num_days=28
    num_shifts=4   
    bounds[2,0]=20
    bounds[2,1]=29
#    bounds[6,0]=9
#    bounds[6,1]=16
    bounds[6,2]=1
#    bounds[6,3]=5
#    bounds[6,4]=1
#    bounds[6,5]=8   
    bounds[9,0]=1
    bounds[9,1]=8
    bounds[10,1]=1
    
    subset1_bounds[6,0]=6
    subset1_bounds[6,1]=16
    subset1_bounds[6,3]=5
    subset1_bounds[6,4]=1
    subset1_bounds[6,5]=7   
    
    subset2_bounds[6,0]=3
    subset2_bounds[6,1]=10
    subset2_bounds[6,3]=8
    subset2_bounds[6,4]=2
    subset2_bounds[6,5]=5 
    
bounds=bounds.astype(np.int64)
subset1_bounds=subset1_bounds.astype(np.int64)
subset2_bounds=subset2_bounds.astype(np.int64)

if not os.path.exists(directory):
    os.makedirs(directory)
if not os.path.exists(directory+"/solutions"):
    os.makedirs(directory+"/solutions")
    
my_csv = open(directory+ "/results.csv" ,"w+")
csvWriter = csv.writer(my_csv,delimiter=',')
detail_csv = open(directory+ "/detail_results.csv" ,"w+")
row=['Nurses','# of Sol', '# of Sample','TP','TP_err','FN','FN_err', 'Time Taken', 'Time Taken_err']
csvWriter.writerow(row)  


dataDir=directory+"/solutions"+"/*N"+str(num_nurses)+"*.csv"
for fl in glob.glob(dataDir): 
      os.remove(fl) 
            
nurse_skill=np.zeros(num_nurses)
for i in range(num_nurses):
    random.seed(i)
    nurse_skill[i]=random.randint(0,1)
print(nurse_skill)
high_nurse=[i for i, x in enumerate(nurse_skill) if x]
low_nurse=[i for i, x in enumerate(nurse_skill) if not x]
gfl.generateSolution(num_nurses,num_days,num_shifts,orderingNotImp,numFiles,directory+"/solutions",constrList,bounds,subset1_bounds,subset2_bounds,nurse_skill,high_nurse,low_nurse)
#    
try:
    os.remove(directory+'/result_N'+str(num_nurses)+"_"+"0"+'_forAll.csv')
    os.remove(directory+'/result_N'+str(num_nurses)+"_"+"00"+'_forAll.csv')
    os.remove(directory+'/result_N'+str(num_nurses)+"_"+"01"+'_forAll.csv')
except OSError:
    pass
rd.learnConstraintsForAll(dataDir,num_nurses,1,directory,nurse_skill)

bounds_prev=np.zeros([num_constrType,num_constr])
bounds_prev0=np.zeros([num_constrType,num_constr])
bounds_prev1=np.zeros([num_constrType,num_constr])

data=rd.readCSV(directory+'/result_N'+str(num_nurses)+"_"+"0"+'_forAll.csv')
data_transpose=list(zip(*data))
data_int=np.zeros([len(data_transpose),len(data_transpose[0])-1])
for i in range(len(data_transpose)):
    for j in range(1,len(data_transpose[i])):
        if data_transpose[i][j]!='':
            data_int[i,j-1]=int(data_transpose[i][j])

bounds_tr=np.zeros([len(data_int[0]),num_constrType,num_constr])
for j in range(len(data_int[0])):
    k=0
    for i in range(83):
        if (i+1)%7 != 0:
            bounds_tr[j,int(k/6),k%6]=data_int[i,j]
            k+=1
bounds_tr=bounds_tr.astype(np.int64)

data=rd.readCSV(directory+'/result_N'+str(num_nurses)+"_"+"00"+'_forAll.csv')
data_transpose=list(zip(*data))
data_int=np.zeros([len(data_transpose),len(data_transpose[0])-1])
for i in range(len(data_transpose)):
    for j in range(1,len(data_transpose[i])):
        if data_transpose[i][j]!='':
            data_int[i,j-1]=int(data_transpose[i][j])

bounds_tr0=np.zeros([len(data_int[0]),num_constrType,num_constr])
for j in range(len(data_int[0])):
    k=0
    for i in range(83):
        if (i+1)%7 != 0:
            bounds_tr0[j,int(k/6),k%6]=data_int[i,j]
            k+=1
bounds_tr0=bounds_tr0.astype(np.int64)

data=rd.readCSV(directory+'/result_N'+str(num_nurses)+"_"+"01"+'_forAll.csv')
data_transpose=list(zip(*data))
data_int=np.zeros([len(data_transpose),len(data_transpose[0])-1])
for i in range(len(data_transpose)):
    for j in range(1,len(data_transpose[i])):
        if data_transpose[i][j]!='':
            data_int[i,j-1]=int(data_transpose[i][j])

bounds_tr1=np.zeros([len(data_int[0]),num_constrType,num_constr])
for j in range(len(data_int[0])):
    k=0
    for i in range(83):
        if (i+1)%7 != 0:
            bounds_tr1[j,int(k/6),k%6]=data_int[i,j]
            k+=1
bounds_tr1=bounds_tr1.astype(np.int64)

for numSol in [1,10,50,100,200]:
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
            os.remove(directory+'/result_N'+str(num_nurses)+'_seed'+str(seed)+"_0"+'.csv')
            os.remove(directory+'/result_N'+str(num_nurses)+'_seed'+str(seed)+"_00"+'.csv')
            os.remove(directory+'/result_N'+str(num_nurses)+'_seed'+str(seed)+"_01"+'.csv')
        except OSError:
            pass
        start=time.clock()
        rd.learnConstraints(dataDir,numSol,num_nurses,1,directory,numFiles,seed,nurse_skill)
        end=time.clock()
        
        data=rd.readCSV(directory+'/result_N'+str(num_nurses)+'_seed'+str(seed)+"_0"+'.csv')
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
        
        
        data=rd.readCSV(directory+'/result_N'+str(num_nurses)+'_seed'+str(seed)+"_00"+'.csv')
        data_transpose=list(zip(*data))
        data_int=np.zeros([len(data_transpose),len(data_transpose[0])-1])
        for i in range(len(data_transpose)):
            for j in range(1,len(data_transpose[i])):
                if data_transpose[i][j]!='':
                    data_int[i,j-1]=int(data_transpose[i][j])
        bounds_learned0=np.zeros([num_constrType,num_constr])
        k=0
        for i in range(len(data_transpose)):
            if (i+1)%7 != 0:
                if (k%6)%2==0:
                    if stats.mode(data_int[i], axis=None)[0][0] == 0:
                        bounds_learned0[int(k/6),k%6]=0
                    else:
                        bounds_learned0[int(k/6),k%6]=np.min(data_int[i])
                if (k%6)%2!=0:
                    if stats.mode(data_int[i], axis=None)[0][0] == 0:
                        bounds_learned0[int(k/6),k%6]=0
                    else:
                        bounds_learned0[int(k/6),k%6]=np.max(data_int[i])
                k+=1
        bounds_learned0=bounds_learned0.astype(np.int64)
        
        
        data=rd.readCSV(directory+'/result_N'+str(num_nurses)+'_seed'+str(seed)+"_01"+'.csv')
        data_transpose=list(zip(*data))
        data_int=np.zeros([len(data_transpose),len(data_transpose[0])-1])
        for i in range(len(data_transpose)):
            for j in range(1,len(data_transpose[i])):
                if data_transpose[i][j]!='':
                    data_int[i,j-1]=int(data_transpose[i][j])
        bounds_learned1=np.zeros([num_constrType,num_constr])
        k=0
        for i in range(len(data_transpose)):
            if (i+1)%7 != 0:
                if (k%6)%2==0:
                    if stats.mode(data_int[i], axis=None)[0][0] == 0:
                        bounds_learned1[int(k/6),k%6]=0
                    else:
                        bounds_learned1[int(k/6),k%6]=np.min(data_int[i])
                if (k%6)%2!=0:
                    if stats.mode(data_int[i], axis=None)[0][0] == 0:
                        bounds_learned1[int(k/6),k%6]=0
                    else:
                        bounds_learned1[int(k/6),k%6]=np.max(data_int[i])
                k+=1
        bounds_learned1=bounds_learned1.astype(np.int64)
#        print(bounds_learned)
#        print(bounds_tr)
        if not (np.array_equal(bounds_learned,bounds_prev) and np.array_equal(bounds_learned0,bounds_prev0) and np.array_equal(bounds_learned1,bounds_prev1)):
            fn=0
            for i in range(len(bounds_tr)):
                accept=1
                for j in range(num_constrType):
                    for k in range(num_constr):
                        if bounds_learned[j,k] != 0:
#                            print(j,k)
                            if k%2==0 and bounds_learned[j,k]>bounds_tr[i,j,k]:
                                accept=0
                                break
                            if k%2==1 and bounds_learned[j,k]<bounds_tr[i,j,k]:
                                accept=0
                                break
                        if bounds_learned0[j,k] != 0:
#                            print(j,k)
                            if k%2==0 and bounds_learned0[j,k]>bounds_tr0[i,j,k]:
                                accept=0
                                break
                            if k%2==1 and bounds_learned0[j,k]<bounds_tr0[i,j,k]:
                                accept=0
                                break
                        if bounds_learned1[j,k] != 0:
#                            print(j,k)
                            if k%2==0 and bounds_learned1[j,k]>bounds_tr1[i,j,k]:
                                accept=0
                                break
                            if k%2==1 and bounds_learned1[j,k]<bounds_tr1[i,j,k]:
                                accept=0
                                break
                    if accept==0:
                        break
                if accept==0:
                    fn+=1
        else:
            fn=fn_prev
        fn_prev=fn
        bounds_prev=bounds_learned
        bounds_prev0=bounds_learned0
        bounds_prev1=bounds_learned1
        
        row=[]
        row.extend([num_nurses])
        row.extend([numSol])
        row.extend([numSam])
        row.extend([fn])
        row.extend([end-start])
        csvWriter_detail.writerow(row)
        detail_csv.close()
        print(fn)
        tot_tp[seed]=numSam-fn
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
#    print(constrRej3)
#    print(constrRej4)
    
shutil.rmtree(directory+"/solutions")
os.makedirs(directory+"/solutions")
    
my_csv.close()   
#
#for numSol in [1,10,50,100]:
#    print("########################## Number of Nurses:",num_nurses," ##########################")
#    try:
#        os.remove(directory+'/result_N'+str(num_nurses)+'.csv')
#    except OSError:
#        pass
#    
#    start=time.clock()
#    rd.learnConstraints(dataDir,numSol,num_nurses,1,directory,numFiles)
#    end=time.clock()
#    
#    data=rd.readCSV(directory+'/result_N'+str(num_nurses)+'.csv')
#    data_transpose=list(zip(*data))
#    data_int=np.zeros([len(data_transpose),len(data_transpose[0])-1])
#    for i in range(len(data_transpose)):
#        for j in range(1,len(data_transpose[i])):
#            if data_transpose[i][j]!='':
#                data_int[i,j-1]=int(data_transpose[i][j])
#    
#    bounds_learned=np.zeros([num_constrType,num_constr])
#    
#    k=0
#    for i in range(len(data_transpose)):
#        if (i+1)%7 != 0:
#            if (k%6)%2==0:
#                if stats.mode(data_int[i], axis=None)[0][0] == 0:
#                    bounds_learned[int(k/6),k%6]=0
#                else:
#                    bounds_learned[int(k/6),k%6]=np.min(data_int[i])
#            if (k%6)%2!=0:
#                if stats.mode(data_int[i], axis=None)[0][0] == 0:
#                    bounds_learned[int(k/6),k%6]=0
#                else:
#                    bounds_learned[int(k/6),k%6]=np.max(data_int[i])
#            k+=1
#    bounds_learned=bounds_learned.astype(np.int64)
#    
#    
#    
#    data=rd.readCSV(directory+'/subset_result_N'+str(num_nurses)+'.csv')
#    data_transpose=list(zip(*data))
#    data_int=np.zeros([len(data_transpose),len(data_transpose[0])-1])
#    for i in range(len(data_transpose)):
#        for j in range(1,len(data_transpose[i])):
#            if data_transpose[i][j]!='':
#                data_int[i,j-1]=int(data_transpose[i][j])
#    
#    subset_bounds_learned=np.zeros([num_constrType,num_constr])
#    
#    k=0
#    for i in range(len(data_transpose)):
#        if (i+1)%7 != 0:
#            if (k%6)%2==0:
#                if stats.mode(data_int[i], axis=None)[0][0] == 0:
#                    subset_bounds_learned[int(k/6),k%6]=0
#                else:
#                    subset_bounds_learned[int(k/6),k%6]=np.min(data_int[i])
#            if (k%6)%2!=0:
#                if stats.mode(data_int[i], axis=None)[0][0] == 0:
#                    subset_bounds_learned[int(k/6),k%6]=0
#                else:
#                    subset_bounds_learned[int(k/6),k%6]=np.max(data_int[i])
#            k+=1
#    subset_bounds_learned=subset_bounds_learned.astype(np.int64)
##        totSam1,fp,LeConsReject,constrRej1,constrRej2=gfl.generateSample(num_nurses,num_days,num_shifts,orderingNotImp,numSam,num_constrType,constrList,bounds_learned,bounds)
#
#    if not np.array_equal(bounds_learned,bounds_prev):
#        totSam2,fn,TrConsReject,constrRej3,constrRej4=gfl.generateSample(num_nurses,num_days,num_shifts,orderingNotImp,numSam,num_constrType,constrList,bounds,bounds_learned,subset_bounds_learned)
#    else:
#        totSam2,fn,TrConsReject,constrRej3,constrRej4=totSam2_prev,fn_prev,TrConsReject_prev,constrRej3_prev,constrRej4_prev
#    totSam2_prev,fn_prev,TrConsReject_prev,constrRej3_prev,constrRej4_prev=totSam2,fn,TrConsReject,constrRej3,constrRej4
#    bounds_prev=bounds_learned
#    
#    row=[]
#    row.extend([num_nurses])
#    row.extend([numSol])
#    row.extend([numSam])
##        row.extend([totSam1])
#    row.extend([totSam2])
##        row.extend([fp])
#    row.extend([fn])
##        row.extend([LeConsReject])
#    row.extend([TrConsReject])
#    row.extend([end-start])
#    csvWriter.writerow(row)
#    print(row)
##        print(constrRej1)
##        print(constrRej2)
#    print(constrRej3)
#    print(constrRej4)
#        
#shutil.rmtree(directory+"/solutions")
#os.makedirs(directory+"/solutions")
#    
#my_csv.close()   






