#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 11:53:31 2018

@author: mohit
"""

import generateSol as gs
import generate_v3 as gfl
import readData as rd
import numpy as np
import glob
import os
import csv
import random
from scipy import stats
import time
import sys
import shutil
from os.path import expanduser
home = expanduser("~")

numSam=int(sys.argv[1])
numFiles=numSam

hospSize=sys.argv[2]

tag="v33_"+str(numSam)+"_"+hospSize

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
    bounds[6,2]=2
    bounds[6,3]=2
    bounds[6,4]=4
    bounds[6,5]=4   
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
bounds=bounds.astype(np.int64)
if not os.path.exists(directory):
    os.makedirs(directory)
if not os.path.exists(directory+"/solutions"):
    os.makedirs(directory+"/solutions")
    
my_csv = open(directory+ "/results.csv" ,"w+")
csvWriter = csv.writer(my_csv,delimiter=',')
detail_csv = open(directory+ "/detail_results.csv" ,"w+")
row=['Nurses','# of Sol', '# of Sample', 'FP','FP_err','FN','FN_err', 'Time Taken', 'Time Taken_err', 'Obj1','minObj1', 'Obj2', 'Obj2_err']
csvWriter.writerow(row)   



dataDir=directory+"/solutions"+"/*N"+str(num_nurses)+"*.csv"
for fl in glob.glob(dataDir): 
      os.remove(fl)
      
#nurse_preference={}
#n=round(num_nurses*0.2)
#for i in range(n):
#    random.seed(i)
#    nurse_preference[i]=(random.randint(0,num_days-1),random.randint(0,num_shifts-1)) 
    
gfl.generateSolution(num_nurses,num_days,num_shifts,orderingNotImp,numFiles,directory+"/solutions",constrList,bounds)

try:
    os.remove(directory+'/result_N'+str(num_nurses)+'_forAll.csv')
except OSError:
    pass
rd.learnConstraintsForAll(dataDir,num_nurses,1,directory)

bounds_prev=np.zeros([num_constrType,num_constr])

objVal=0
minObjVal=0
i=0
for fl in glob.glob(directory+"/solutions"+"/*objVal"+str(num_nurses)+"*.csv"):
    tmp=float(rd.readCSV(fl)[0][0])
    if i==0:
        minObjVal=tmp
    else:
        minObjVal=min(minObjVal,tmp)
    objVal += tmp
    i+=1
objVal=objVal/numFiles

data=rd.readCSV(directory+'/result_N'+str(num_nurses)+'_forAll.csv')
data_transpose=list(zip(*data))
data_int=np.zeros([len(data_transpose),len(data_transpose[0])-1])
for i in range(len(data_transpose)):
    for j in range(1,len(data_transpose[i])):
        if data_transpose[i][j]!='':
            data_int[i,j-1]=int(data_transpose[i][j])

bounds_tr=np.zeros([len(data_int[0]),num_constrType,num_constr])
#print(data_transpose)
#print(data_int)
for j in range(len(data_int[0])):
    k=0
    for i in range(83):
        if (i+1)%7 != 0:
            bounds_tr[j,int(k/6),k%6]=data_int[i,j]
            k+=1
bounds_tr=bounds_tr.astype(np.int64)

#for numSol in [100]:
for numSol in [1]:
    print("########################## Number of Nurses:",num_nurses," NumSol: ",numSol," ##########################")
    try:
        os.remove(directory+'/result_N'+str(num_nurses)+'.csv')
    except OSError:
        pass
#        gfl.generateSolution(num_nurses,num_days,num_shifts,orderingNotImp,numSol,directory+"/solutions",constrList,bounds)
    numSeed=1
    tot_tp1=np.zeros(numSeed)
    tot_tp2=np.zeros(numSeed)
    tot_fn=np.zeros(numSeed)
    tot_fp=np.zeros(numSeed)
    tot_time=np.zeros(numSeed)
    tot_objVal1=np.zeros(numSeed)
    tot_objVal2=np.zeros(numSeed)
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
        
        if not np.array_equal(bounds_learned,bounds_prev):
            totSam1,fp,LeConsReject,constrRej1,constrRej2,objVal2=gfl.generateSample2(num_nurses,num_days,num_shifts,orderingNotImp,numSam,num_constrType,constrList,bounds_learned,bounds)
#            totSam2,fn,TrConsReject,constrRej3,constrRej4,objVal1=gfl.generateSample1(num_nurses,num_days,num_shifts,orderingNotImp,numSam,num_constrType,constrList,bounds,bounds_learned)
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
                    if accept==0:
                        break
                if accept==0:
                    fn+=1
        else:
            totSam1,fp,LeConsReject,constrRej1,constrRej2,objVal2=totSam1_prev,fp_prev,LeConsReject_prev,constrRej1_prev,constrRej2_prev,objVal2_prev
            fn=fn_prev
        fn_prev=fn
        totSam1_prev,fp_prev,LeConsReject_prev,constrRej1_prev,constrRej2_prev,objVal2_prev=totSam1,fp,LeConsReject,constrRej1,constrRej2,objVal2
        bounds_prev=bounds_learned
        
        row=[]
        row.extend([num_nurses])
#        row.extend([len(nurse_preference)])
        row.extend([numSol])
        row.extend([numSam])
#        row.extend([totSam1])
#        row.extend([totSam2])
        row.extend([fp])
        row.extend([fn])
#        row.extend([LeConsReject])
#        row.extend([TrConsReject])
#        row.extend([objVal1])
        row.extend([objVal2])
        row.extend([end-start])
        csvWriter_detail.writerow(row)
        detail_csv.close()
        
#        tot_tp1[seed]=totSam1
#        tot_tp2[seed]=totSam2
        tot_fp[seed]=fp
        tot_fn[seed]=fn
#        tot_objVal1[seed]=objVal1
        tot_objVal2[seed]=objVal2
        tot_time[seed]=(end-start)
        print(row)
    
    row=[]
    row.extend([num_nurses])
    row.extend([numSol])
    row.extend([numSam])
#    row.extend([sum(tot_tp1)/numSeed])
#    row.extend([np.std(tot_tp1)/np.sqrt(numSeed)])
#    row.extend([sum(tot_tp2)/numSeed])
#    row.extend([np.std(tot_tp2)/np.sqrt(numSeed)])
    row.extend([sum(tot_fp)/numSeed])
    row.extend([np.std(tot_fp)/np.sqrt(numSeed)])
    row.extend([sum(tot_fn)/numSeed])
    row.extend([np.std(tot_fn)/np.sqrt(numSeed)])
    row.extend([sum(tot_time)/numSeed])
    row.extend([np.std(tot_time)/np.sqrt(numSeed)])
    row.extend([objVal])
    row.extend([minObjVal])
    row.extend([sum(tot_objVal2)/numSeed])
    row.extend([np.std(tot_objVal2)/np.sqrt(numSeed)])
    csvWriter.writerow(row)
    print(row)
    
shutil.rmtree(directory+"/solutions")
os.makedirs(directory+"/solutions")
        
my_csv.close()    