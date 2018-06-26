#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 11:45:50 2018

@author: mohit
"""

import sampler
import countor
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

def readBounds(file,num_constrType,num_constr):
    data=rd.readCSV(file)
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
    return bounds_tr.astype(np.int64)

def aggrBounds(selbounds,num_constrType,num_constr):
    bounds_learned=np.zeros([num_constrType,num_constr])
    for i in range(num_constrType):
        for j in range(num_constr):
            if j%2==0:
                if stats.mode(selbounds[:,i,j])[0][0] == 0:
                    bounds_learned[int((i*num_constr+j)/6),(i*num_constr+j)%6]=0
                else:
                    bounds_learned[int((i*num_constr+j)/6),(i*num_constr+j)%6]=np.min(selbounds[:,i,j])
            if j%2!=0:
                if stats.mode(selbounds[:,i,j])[0][0] == 0:
                    bounds_learned[int((i*num_constr+j)/6),(i*num_constr+j)%6]=0
                else:
                    bounds_learned[int((i*num_constr+j)/6),(i*num_constr+j)%6]=np.max(selbounds[:,i,j])
    return bounds_learned.astype(np.int64)

###########checks if bound2 is more constrained than bound1##################
def moreConstrained(bound1,bound2,num_constrType,num_constr):
    output=1
    for i in range(num_constrType):
        for j in range(num_constr):
            if j%2==0 and bound2[i,j]<bound1[i,j]:
                output=0
                break
            if j%2==1 and bound2[i,j]>bound1[i,j]:
                output=0
                break
        if output==0:
            break
    return output

home = expanduser("~")

numSam=int(sys.argv[1])
bk=int(sys.argv[2])
mt=int(sys.argv[3])
hs=int(sys.argv[4])

extraConstPerc = 10

numFiles=numSam

tag=str(bk)+str(mt)+str(hs)+"_"+str(numSam)

directory=os.getcwd()+'/data/'+tag

num_nurses=15
num_days=7
num_shifts=3   
orderingNotImp=[2]

num_constrType=12
num_constr=6
constrList=[[(0,),(1,)],[(0,),(2,)],[(0,),(1,2)],[(1,),(0,)],[(1,),(2,)],[(1,),(0,2)],[(2,),(0,)],[(2,),(1,)],[(2,),(0,1)],[(0,1),(2,)],[(0,2),(1,)],[(1,2),(0,)]]

tbounds=np.zeros([num_constrType,num_constr])
tbounds0=np.zeros([num_constrType,num_constr])
tbounds1=np.zeros([num_constrType,num_constr])
if hs==0:
    num_nurses=10
    num_days=28
    num_shifts=4   
    tbounds[2,0]=4
    tbounds[2,1]=6
    tbounds[6,2]=1
    tbounds[6,3]=7
    tbounds[6,4]=2  
    tbounds[9,0]=1
    tbounds[9,1]=2
    tbounds[10,1]=1 
    
    if bk==0:
        tbounds[6,0]=9
        tbounds[6,1]=16
        tbounds[6,5]=8 
    
    if bk==1:
        tbounds0[6,0]=9
        tbounds0[6,1]=16
        tbounds0[6,5]=8   
        
        tbounds1[6,0]=4
        tbounds1[6,1]=9
        tbounds1[6,5]=4 
    
if hs==1:
    num_nurses=31
    num_days=28
    num_shifts=4   
    tbounds[2,0]=16
    tbounds[2,1]=24
    tbounds[6,2]=1
    tbounds[6,4]=2 
    tbounds[9,0]=2
    tbounds[9,1]=9
    tbounds[10,1]=1
    
    if bk==0:
        tbounds[6,0]=6
        tbounds[6,1]=16
        tbounds[6,3]=7
        tbounds[6,5]=8  
    
    if bk==1:
        tbounds0[6,0]=4
        tbounds0[6,1]=13
        tbounds0[6,3]=8
        tbounds0[6,5]=4   
        
        tbounds1[6,0]=6
        tbounds1[6,1]=16
        tbounds1[6,3]=7
        tbounds1[6,5]=8
    
if hs==2:
    num_nurses=49
    num_days=28
    num_shifts=4   
    tbounds[2,0]=20
    tbounds[2,1]=29
    tbounds[6,2]=1  
    tbounds[9,0]=1
    tbounds[9,1]=8
    tbounds[10,1]=1
    
    if bk==0:
        tbounds[6,0]=6
        tbounds[6,1]=16
        tbounds[6,3]=5
        tbounds[6,4]=1
        tbounds[6,5]=7 
    
    if bk==1:
        tbounds0[6,0]=6
        tbounds0[6,1]=16
        tbounds0[6,3]=5
        tbounds0[6,4]=1
        tbounds0[6,5]=7   
        
        tbounds1[6,0]=3
        tbounds1[6,1]=10
        tbounds1[6,3]=8
        tbounds1[6,4]=2
        tbounds1[6,5]=5 
    
tbounds=tbounds.astype(np.int64)
if bk==1:
    tbounds0=tbounds0.astype(np.int64)
    tbounds1=tbounds1.astype(np.int64)

soln=directory+"/solutions"
result=directory+"/results"

if not os.path.exists(directory):
    os.makedirs(directory)
if not os.path.exists(soln):
    os.makedirs(soln)
if not os.path.exists(result):
    os.makedirs(result)
   
nurse_skill=np.zeros(num_nurses)
if bk==1:            
    for i in range(num_nurses):
        random.seed(i)
        nurse_skill[i]=random.randint(0,1)
#    print(nurse_skill)

for fl in glob.glob(soln+"/*.csv"): 
      os.remove(fl) 
sampler.generateSample(num_nurses,num_days,num_shifts,numSam,extraConstPerc,nurse_skill,tbounds,tbounds0,tbounds1,soln,bk,mt)

for fl in glob.glob(result+"/*.csv"): 
      os.remove(fl) 
countor.learnConstraintsForAll(directory,num_nurses,nurse_skill,bk,mt,hs)

tag=str(bk)+str(mt)+str(hs)
file=result+"/learnedBounds"+"_"+tag+"0.csv"
lbounds=readBounds(file,num_constrType,num_constr)

if bk==1:
    file=result+"/learnedBounds"+"_"+tag+"00.csv"
    lbounds0=readBounds(file,num_constrType,num_constr)
    
    file=result+"/learnedBounds"+"_"+tag+"01.csv"
    lbounds1=readBounds(file,num_constrType,num_constr)


bounds_prev=np.zeros([num_constrType,num_constr])
bounds_prev0=np.zeros([num_constrType,num_constr])
bounds_prev1=np.zeros([num_constrType,num_constr])

for numSol in [1, 10, 50, 100]:
    print("########################## Number of Nurses:",num_nurses," NumSol: ",numSol," ##########################")
    
    recall,fn,precision,fp=0,0,0,0
    numSeed=4
    tot_rec=np.zeros(numSeed)
    tot_pre=np.zeros(numSeed)
    tot_fn=np.zeros(numSeed)
    tot_fp=np.zeros(numSeed)
    tot_time=np.zeros(numSeed)
    
    for seed in range(numSeed):
        random.seed(seed)
        selRows=random.sample(range(0,numSam),numSol)
        selbounds=np.array([lbounds[i] for i in selRows])
        bounds_learned=aggrBounds(selbounds,num_constrType,num_constr)
        
        bounds_learned0=np.zeros([num_constrType,num_constr])
        bounds_learned1=np.zeros([num_constrType,num_constr])
        if bk==1:
            selbounds0=np.array([lbounds0[i] for i in selRows])
            selbounds1=np.array([lbounds1[i] for i in selRows])
            bounds_learned0=aggrBounds(selbounds0,num_constrType,num_constr)
            bounds_learned1=aggrBounds(selbounds1,num_constrType,num_constr)
        
        selbounds=np.array([lbounds[i] for i in range(len(lbounds)) if i not in selRows])
        for i in range(len(selbounds)):
            accept=0
            accept=moreConstrained(bounds_learned,selbounds[i],num_constrType,num_constr)
            if accept==0 and bk==1:
                accept=moreConstrained(bounds_learned0,selbounds0[i],num_constrType,num_constr)
                if accept==0:
                    accept=moreConstrained(bounds_learned1,selbounds1[i],num_constrType,num_constr)
            recall+=accept
        
        tmpDir=directory+"/tmp"
        if not os.path.exists(tmpDir):
            os.makedirs(tmpDir)
        for fl in glob.glob(tmpDir+"/*.csv"): 
            os.remove(fl) 
        sampler.generateSample(num_nurses,num_days,num_shifts,numSam,extraConstPerc,nurse_skill,bounds_learned,bounds_learned0,bounds_learned1,tmpDir,bk,0)
        
        countor.learnConstraintsForAll(tmpDir,num_nurses,nurse_skill,bk,0,hs)
        tag=str(bk)+str(0)+str(hs)
        file=tmpDir+"/results"+"/learnedBounds"+"_"+tag+"0.csv"
        tmpBounds=readBounds(file,num_constrType,num_constr)
        if bk==1:
            file=tmpDir+"/results"+"/learnedBounds"+"_"+tag+"00.csv"
            tmpBounds0=readBounds(file,num_constrType,num_constr)
            
            file=tmpDir+"/results"+"/learnedBounds"+"_"+tag+"01.csv"
            tmpBounds1=readBounds(file,num_constrType,num_constr)
        
        for i in range(len(tmpBounds)):
            accept=0
            accept=moreConstrained(tbounds,tmpBounds[i],num_constrType,num_constr)
            if accept==0 and bk==1:
                accept=moreConstrained(tbounds0,tmpBounds0[i],num_constrType,num_constr)
                if accept==0:
                    accept=moreConstrained(tbounds1,tmpBounds1[i],num_constrType,num_constr)
            precision+=accept







    
    try:
        os.remove(directory+'/result_N'+str(num_nurses)+'.csv')
    except OSError:
        pass
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
            fp=0
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
            fp=gfl.generateSample(num_nurses,num_days,num_shifts,orderingNotImp,numSol,constrList,bounds_learned,bounds_learned0,bounds_learned1,bounds,subset1_bounds,subset2_bounds,nurse_skill,high_nurse,low_nurse,nurse_preference,directory,seed)
        
        else:
            fn=fn_prev
            fp=fp_prev
        fn_prev=fn
        fp_prev=fp
        bounds_prev=bounds_learned
        bounds_prev0=bounds_learned0
        bounds_prev1=bounds_learned1
        
        row=[]
        row.extend([num_nurses])
        row.extend([numSol])
        row.extend([numSam])
        row.extend([fn])
        row.extend([fp])
        row.extend([end-start])
        csvWriter_detail.writerow(row)
        detail_csv.close()
        print(fn)
        tot_tp[seed]=numSam-fn
        tot_fn[seed]=fn
        tot_fp[seed]=fp
        tot_time[seed]=(end-start)
    
    row=[]
    row.extend([num_nurses])
    row.extend([numSol])
    row.extend([numSam])
    row.extend([sum(tot_tp)/numSeed])
    row.extend([np.std(tot_tp)/np.sqrt(numSeed)])
    row.extend([sum(tot_fn)/numSeed])
    row.extend([np.std(tot_fn)/np.sqrt(numSeed)])
    row.extend([sum(tot_fp)/numSeed])
    row.extend([np.std(tot_fp)/np.sqrt(numSeed)])
    row.extend([sum(tot_time)/numSeed])
    row.extend([np.std(tot_time)/np.sqrt(numSeed)])
    csvWriter.writerow(row)
    print(row)
#    print(constrRej3)
#    print(constrRej4)
    
shutil.rmtree(directory+"/solutions")
os.makedirs(directory+"/solutions")
    
my_csv.close()   