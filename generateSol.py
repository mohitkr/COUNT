#!/home/mohit/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  14 12:20:40 2018

@author: mohit
"""
import csv
import numpy as np
import pandas as pd
from gurobipy import *
import os.path
import itertools as it
import constraintFormulation as cf

    
def generateSam(num_nurses,num_days,num_shifts,orderingNotImp,numSol,direc,num_constrType,constrList,bounds,bounds_learned):
    
    N=list(range(num_nurses))
    D=list(range(num_days))
    S=list(range(num_shifts))
    variables=[D,S,N]
    
    #Forbidden Shift Successions
    F=[(num_shifts-1,0)]
        
    #Weekends
    W=[(5,6)]
    
    #minimum number of nurses required for covering a shift s on day d
    C_day_min={}
    for d in D:
        C_day_min[d]=3
                
    #Maximum number of nurses required for covering a shift s on day d
    C_day_max={}
    for d in D:
        C_day_max[d]=8
            
    #minimum number of nurses required for covering a shift s on day d
    C={}
    for d in D:
        for s in S:
            C[d,s]=1
                
    #Maximum number of nurses required for covering a shift s on day d
    C_max={}
    for d in D:
        for s in S:
            C_max[d,s]=3
    
    #Min Consecutive Working Days
    CW_min={}
    for n in N:
        CW_min[n]=1
        
    #Max Consecutive Working Days
    CW_max={}
    for n in N:
        CW_max[n]=5
       
    #Min Consecutive Working Days in shift s
    CWs_min={}
    for n in N:
        for s in S:
            CWs_min[n,s]=1
        
    #Max Consecutive Working Days in shift s
    CWs_max={}
    for n in N:
        for s in S:
            CWs_max[n,s]=4
     
    #Min Consecutive Free Days
    CF_min={}
    for n in N:
        CF_min[n]=2
        
    #Max Consecutive Free Days
    CF_max={}
    for n in N:
        CF_max[n]=4
    
    #Max Working Days
    WD_max={}
    for n in N:
        WD_max[n]=6
        
    #Min Working Days
    WD_min={}
    for n in N:
        WD_min[n]=2

    rSam=0
    tSam=0
    lSam=0
        
    try:
        m=Model("nspSolver")
        m.setParam(GRB.Param.OutputFlag,0)
        ########### Decision Variables #############
        o = m.addVars(N,D,S, vtype=GRB.BINARY, name="o")
        p = m.addVars(N,D, vtype=GRB.BINARY, name="p")
        S1=m.addVars(D,S,vtype=GRB.CONTINUOUS, lb=-100 , name="S1")
        S11=m.addVars(D,S,vtype=GRB.CONTINUOUS, name="S11")
        tw = m.addVars(N,D,D, vtype=GRB.CONTINUOUS, name="tw")
        sw = m.addVars(N,D,D, vtype=GRB.CONTINUOUS, name="sw")
        tw1 = m.addVars(N,D,S,D, vtype=GRB.CONTINUOUS, name="tw1")
        sw1 = m.addVars(N,D,S,D, vtype=GRB.CONTINUOUS, name="sw1")
        
        tf = m.addVars(N,D,D, vtype=GRB.CONTINUOUS, name="tf")
        sf = m.addVars(N,D,D, vtype=GRB.CONTINUOUS, name="sf")
        
        ########### Required Constraints #############
        m.addConstrs(
                (o.sum(n,d,'*') == p[n,d] 
                for n in N for d in D),"po")
        
        
        ########### Hard Constraints #############
        m.addConstrs(
                (p[n,d]<=1 for n in N for d in D),"singleAssignmentPerDay")
        m.addConstrs(
                (p.sum('*',d) >= C_day_min[d] for d in D),"minNurses")
        m.addConstrs(
                (p.sum('*',d) <= C_day_max[d] for d in D),"maxNurses")
        m.addConstrs(
                (o.sum('*',d,s) >= C[d,s]
                for d in D for s in S),"underStaffing")
        m.addConstrs(
                (o.sum('*',d,s) <= C_max[d,s]
                for d in D for s in S),"overStaffing")
        m.addConstrs(
                (o[n,d,f1]+o[n,d+1,f2]<=1 
                 for n in N for d in D[:-1] for f1,f2 in F),"shiftTypeSuccession")
        m.addConstrs(
                (quicksum(p[n,d] for d in D) <= WD_max[n] for n in N),"maxWorkingDay")
        m.addConstrs(
                (quicksum(p[n,d] for d in D) >= WD_min[n] for n in N),"minWorkingDay")
        
        m.addConstrs(( tw[n,0,0] == p[n,0] for n in N), "MaxConsWork")
        m.addConstrs(( tw[n,d1+1,0] <= p[n,d1+1]
                for n in N for d1 in D if d1<len(D)-1), "MaxConsWork")
        m.addConstrs(( tw[n,d1+1,0] <= 1-p[n,d1]
                for n in N for d1 in D if d1<len(D)-1), "MaxConsWork")
        m.addConstrs(( tw[n,d1+1,0] >= p[n,d1+1]-p[n,d1]
                for n in N for d1 in D if d1<len(D)-1), "MaxConsWork")
        
        m.addConstrs(( tw[n,0,d2] == 0
                for n in N for d2 in D if d2>0), "MaxConsWork")
        m.addConstrs(( tw[n,d1,d2] <= tw[n,d1-1,d2-1]
                for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsWork")
        m.addConstrs(( tw[n,d1,d2] <= p[n,d1]
                for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsWork")
        m.addConstrs(( tw[n,d1,d2] >= p[n,d1]+tw[n,d1-1,d2-1]-1
                for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsWork")
        m.addConstr(
                (quicksum(tw[n,d1,d2] for n in N for d1 in D for d2 in range(CW_max[n],len(D)))==0),"maxconswork"
                )
        
        m.addConstrs(( sw[n,0,d2] == 0
                for n in N for d2 in D), "MinConsWork")
        m.addConstrs(( sw[n,d1,d2] <= tw[n,d1-1,d2]
                for n in N for d1 in D for d2 in D if d1>0), "MinConsWork")
        m.addConstrs(( sw[n,d1,d2] <= 1-p[n,d1]
                for n in N for d1 in D for d2 in D if d1>0), "MinConsWork")
        m.addConstrs(( sw[n,d1,d2] >= tw[n,d1-1,d2]-p[n,d1]
                for n in N for d1 in D for d2 in D if d1>0), "MinConsWork")
        m.addConstr(
                (quicksum(sw[n,d1,d2]*(CW_min[n]-d2) for n in N for d1 in D for d2 in range(CW_min[n]))==0),"minconswork"
                )
        
        m.addConstrs(( tw1[n,0,s,0] == o[n,0,s] for n in N for s in S), "MaxConsSameShift")
        m.addConstrs(( tw1[n,d1+1,s,0] <= o[n,d1+1,s]
                for s in S for n in N for d1 in D if d1<len(D)-1), "MaxConsSameShift")
        m.addConstrs(( tw1[n,d1+1,s,0] <= 1-o[n,d1,s]
                for s in S for n in N for d1 in D if d1<len(D)-1), "MaxConsSameShift")
        m.addConstrs(( tw1[n,d1+1,s,0] >= o[n,d1+1,s]-o[n,d1,s]
                for s in S for n in N for d1 in D if d1<len(D)-1), "MaxConsSameShift")
        
        m.addConstrs(( tw1[n,0,s,d2] == 0
                for s in S for n in N for d2 in D if d2>0), "MaxConsSameShift")
        m.addConstrs(( tw1[n,d1,s,d2] <= tw1[n,d1-1,s,d2-1]
                for s in S for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsSameShift")
        m.addConstrs(( tw1[n,d1,s,d2] <= o[n,d1,s]
                for s in S for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsSameShift")
        m.addConstrs(( tw1[n,d1,s,d2] >= o[n,d1,s]+tw1[n,d1-1,s,d2-1]-1
                for s in S for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsSameShift")
        m.addConstr(
                (quicksum(tw1[n,d1,s,d2] for s in S for n in N for d1 in D for d2 in range(CWs_max[n,s],len(D)))==0),"maxconssameshift"
                )
        
        m.addConstrs(( sw1[n,0,s,d2] == 0
                for s in S for n in N for d2 in D), "MinConsSameShift")
        m.addConstrs(( sw1[n,d1,s,d2] <= tw1[n,d1-1,s,d2]
                for s in S for n in N for d1 in D for d2 in D if d1>0), "MinConsSameShift")
        m.addConstrs(( sw1[n,d1,s,d2] <= 1-o[n,d1,s]
                for s in S for n in N for d1 in D for d2 in D if d1>0), "MinConsSameShift")
        m.addConstrs(( sw1[n,d1,s,d2] >= tw1[n,d1-1,s,d2]-o[n,d1,s]
                for s in S for n in N for d1 in D for d2 in D if d1>0), "MinConsSameShift")
        m.addConstr(
                (quicksum(sw1[n,d1,s,d2]*(CWs_min[n,s]-d2) for s in S for n in N for d1 in D for d2 in range(CWs_min[n,s]))==0),"minconssameshift"
                )
        
        m.addConstrs(( tf[n,0,0] == 1-p[n,0] for n in N), "MaxConsFree")
        m.addConstrs(( tf[n,d1+1,0] <= p[n,d1]
                for n in N for d1 in D if d1<len(D)-1), "MaxConsFree")
        m.addConstrs(( tf[n,d1+1,0] <= 1-p[n,d1+1]
                for n in N for d1 in D if d1<len(D)-1), "MaxConsFree")
        m.addConstrs(( tf[n,d1+1,0] >= p[n,d1]-p[n,d1+1]
                for n in N for d1 in D if d1<len(D)-1), "MaxConsFree")
        
        m.addConstrs(( tf[n,0,d2] == 0
                for n in N for d2 in D if d2>0), "MaxConsFree")
        m.addConstrs(( tf[n,d1,d2] <= tf[n,d1-1,d2-1]
                for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsFree")
        m.addConstrs(( tf[n,d1,d2] <= 1-p[n,d1]
                for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsFree")
        m.addConstrs(( tf[n,d1,d2] >= tf[n,d1-1,d2-1]-p[n,d1]
                for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsFree")
        m.addConstr(
                (quicksum(tf[n,d1,d2] for n in N for d1 in D for d2 in range(CF_max[n],len(D)))==0),"maxconsfree"
                )
        
        m.addConstrs(( sf[n,0,d2] == 0
                for n in N for d2 in D), "MinConsFree")
        m.addConstrs(( sf[n,d1,d2] <= tf[n,d1-1,d2]
                for n in N for d1 in D for d2 in D if d1>0), "MinConsFree")
        m.addConstrs(( sf[n,d1,d2] <= p[n,d1]
                for n in N for d1 in D for d2 in D if d1>0), "MinConsFree")
        m.addConstrs(( sf[n,d1,d2] >= tf[n,d1-1,d2]+p[n,d1]-1
                for n in N for d1 in D for d2 in D if d1>0), "MinConsFree")
        m.addConstr(
                (quicksum(sf[n,d1,d2]*(CF_min[n]-d2) for n in N for d1 in D for d2 in range(CF_min[n]))==0),"minconsfree"
                )
        
        m.setParam(GRB.Param.PoolSolutions, numSol)
        m.setParam(GRB.Param.PoolSearchMode, 2)
        m.optimize()
        nSolutions = m.SolCount
        print('Number of solutions found: ' + str(nSolutions))
        if (m.status==GRB.Status.INFEASIBLE):
            m.computeIIS()
            print('\nThe following constraint(s) cannot be satisfied:')
            for c in m.getConstrs():
                if c.IISConstr:
                    print('%s' % c.constrName)
        if m.status==GRB.Status.OPTIMAL:
            m.write('m.sol')
        
        constrRej1=np.zeros(num_constrType)
        constrRej2=np.zeros(num_constrType)
        for i in range(nSolutions):
            m.setParam(GRB.Param.SolutionNumber,i)
#            print(m.ObjVal)
            solution=m.getAttr('xn', o)
            tmp=np.zeros([num_nurses,num_days,num_shifts])
            for key in solution:
                tmp[key]=round(solution[key])
            print(tmp)
            tSample=np.swapaxes(np.swapaxes(tmp,0,1),1,2)
            accept=1
            for c in range(num_constrType):
                if sum(bounds[c,:])>0:
                    subset=constrList[c]
                    newset=subset[0]+subset[1]
                    idTensor=cf.tensorIndicator(tSample,newset, variables)
                    sumSet = range(len(subset[0]),len(newset))
                    sumTensor_max,sumTensor_min=cf.tensorSum(idTensor,sumSet, np.array(variables)[list(newset)],0)
                    
                    if (bounds[c][0]!=0 and bounds[c][0]>sumTensor_min) or (bounds[c][1]!=0 and bounds[c][1]<sumTensor_max):
                        accept=0
                        constrRej1[c]+=1
                        break
                    if len(set(subset[1]))==1 and len(set(orderingNotImp) & set(subset[1]))==0:
                        minConsZero,maxConsZero,minConsNonZero,maxConsNonZero = cf.tensorConsZero(idTensor,sumSet, np.array(variables)[list(newset)])
                        if (bounds[c][2]!=0 and bounds[c][2]>minConsZero) or (bounds[c][3]!=0 and bounds[c][3]<maxConsZero) and (bounds[c][4]!=0 and bounds[c][4]>minConsNonZero) or (bounds[c][5]!=0 and bounds[c][5]<maxConsNonZero):
                            accept=0
                            constrRej1[c]+=1
                            break         
            if accept==0:
                rSam+=1
            else:
                tSam+=1
            
            accept=1
            for c in range(num_constrType):
                if sum(bounds_learned[c,:])>0:
                    subset=constrList[c]
                    newset=subset[0]+subset[1]
                    idTensor=cf.tensorIndicator(tSample,newset, variables)
                    sumSet = range(len(subset[0]),len(newset))
                    sumTensor_max,sumTensor_min=cf.tensorSum(idTensor,sumSet, np.array(variables)[list(newset)],0)
    #                print(bounds_learned)
                    if (bounds_learned[c][0]!=0 and bounds_learned[c][0]>sumTensor_min) or (bounds_learned[c][1]!=0 and bounds_learned[c][1]<sumTensor_max):
#                        print("here:",sumTensor_min,sumTensor_max)
                        accept=0
                        constrRej2[c]+=1
                        break
                    if len(set(subset[1]))==1 and len(set(orderingNotImp) & set(subset[1]))==0:
                        minConsZero,maxConsZero,minConsNonZero,maxConsNonZero = cf.tensorConsZero(idTensor,sumSet, np.array(variables)[list(newset)])
                        if (bounds_learned[c][2]!=0 and bounds_learned[c][2]>minConsZero) or (bounds_learned[c][3]!=0 and bounds_learned[c][3]<maxConsZero) and (bounds_learned[c][4]!=0 and bounds_learned[c][4]>minConsNonZero) or (bounds_learned[c][5]!=0 and bounds_learned[c][5]<maxConsNonZero):
#                            print("there:",minConsZero,maxConsZero,minConsNonZero,maxConsNonZero)
                            accept=0
                            constrRej2[c]+=1
                            break
            if accept==0:
                lSam+=1
        return tSam,rSam,lSam,constrRej1,constrRej2
            
    except GurobiError as e:
        print('Error code ' + str(e.errno) + ": " + str(e))
    
    except AttributeError:
        print('Encountered an attribute error')
  

def generateSolution(num_nurses,num_days,num_shifts,orderingNotImp,numSol,direc):
    
    N=list(range(num_nurses))
    D=list(range(num_days))
    S=list(range(num_shifts))
    variables=[D,S,N]
    
    #Forbidden Shift Successions
    F=[(num_shifts-1,0)]
        
    #Weekends
    W=[(5,6)]
    
    #minimum number of nurses required for covering a shift s on day d
    C_day_min={}
    for d in D:
        C_day_min[d]=3
                
    #Maximum number of nurses required for covering a shift s on day d
    C_day_max={}
    for d in D:
        C_day_max[d]=8
            
    #minimum number of nurses required for covering a shift s on day d
    C={}
    for d in D:
        for s in S:
            C[d,s]=1
                
    #Maximum number of nurses required for covering a shift s on day d
    C_max={}
    for d in D:
        for s in S:
            C_max[d,s]=3
    
    #Min Consecutive Working Days
    CW_min={}
    for n in N:
        CW_min[n]=1
        
    #Max Consecutive Working Days
    CW_max={}
    for n in N:
        CW_max[n]=5
       
    #Min Consecutive Working Days in shift s
    CWs_min={}
    for n in N:
        for s in S:
            CWs_min[n,s]=1
        
    #Max Consecutive Working Days in shift s
    CWs_max={}
    for n in N:
        for s in S:
            CWs_max[n,s]=4
     
    #Min Consecutive Free Days
    CF_min={}
    for n in N:
        CF_min[n]=2
        
    #Max Consecutive Free Days
    CF_max={}
    for n in N:
        CF_max[n]=4
    
    #Max Working Days
    WD_max={}
    for n in N:
        WD_max[n]=6
        
    #Min Working Days
    WD_min={}
    for n in N:
        WD_min[n]=2
        
    try:
        m=Model("nspSolver")
        m.setParam(GRB.Param.OutputFlag,0)
        ########### Decision Variables #############
        o = m.addVars(N,D,S, vtype=GRB.BINARY, name="o")
        p = m.addVars(N,D, vtype=GRB.BINARY, name="p")
        S1=m.addVars(D,S,vtype=GRB.CONTINUOUS, lb=-100 , name="S1")
        S11=m.addVars(D,S,vtype=GRB.CONTINUOUS, name="S11")
        tw = m.addVars(N,D,D, vtype=GRB.CONTINUOUS, name="tw")
        sw = m.addVars(N,D,D, vtype=GRB.CONTINUOUS, name="sw")
        tw1 = m.addVars(N,D,S,D, vtype=GRB.CONTINUOUS, name="tw1")
        sw1 = m.addVars(N,D,S,D, vtype=GRB.CONTINUOUS, name="sw1")
        
        tf = m.addVars(N,D,D, vtype=GRB.CONTINUOUS, name="tf")
        sf = m.addVars(N,D,D, vtype=GRB.CONTINUOUS, name="sf")
        
        ########### Required Constraints #############
        m.addConstrs(
                (o.sum(n,d,'*') == p[n,d] 
                for n in N for d in D),"po")
        
        
        ########### Hard Constraints #############
        m.addConstrs(
                (p[n,d]<=1 for n in N for d in D),"singleAssignmentPerDay")
        m.addConstrs(
                (p.sum('*',d) >= C_day_min[d] for d in D),"minNurses")
        m.addConstrs(
                (p.sum('*',d) <= C_day_max[d] for d in D),"maxNurses")
        m.addConstrs(
                (o.sum('*',d,s) >= C[d,s]
                for d in D for s in S),"underStaffing")
        m.addConstrs(
                (o.sum('*',d,s) <= C_max[d,s]
                for d in D for s in S),"overStaffing")
        m.addConstrs(
                (o[n,d,f1]+o[n,d+1,f2]<=1 
                 for n in N for d in D[:-1] for f1,f2 in F),"shiftTypeSuccession")
        m.addConstrs(
                (quicksum(p[n,d] for d in D) <= WD_max[n] for n in N),"maxWorkingDay")
        m.addConstrs(
                (quicksum(p[n,d] for d in D) >= WD_min[n] for n in N),"minWorkingDay")
        
        m.addConstrs(( tw[n,0,0] == p[n,0] for n in N), "MaxConsWork")
        m.addConstrs(( tw[n,d1+1,0] <= p[n,d1+1]
                for n in N for d1 in D if d1<len(D)-1), "MaxConsWork")
        m.addConstrs(( tw[n,d1+1,0] <= 1-p[n,d1]
                for n in N for d1 in D if d1<len(D)-1), "MaxConsWork")
        m.addConstrs(( tw[n,d1+1,0] >= p[n,d1+1]-p[n,d1]
                for n in N for d1 in D if d1<len(D)-1), "MaxConsWork")
        
        m.addConstrs(( tw[n,0,d2] == 0
                for n in N for d2 in D if d2>0), "MaxConsWork")
        m.addConstrs(( tw[n,d1,d2] <= tw[n,d1-1,d2-1]
                for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsWork")
        m.addConstrs(( tw[n,d1,d2] <= p[n,d1]
                for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsWork")
        m.addConstrs(( tw[n,d1,d2] >= p[n,d1]+tw[n,d1-1,d2-1]-1
                for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsWork")
        m.addConstr(
                (quicksum(tw[n,d1,d2] for n in N for d1 in D for d2 in range(CW_max[n],len(D)))==0),"maxconswork"
                )
        
        m.addConstrs(( sw[n,0,d2] == 0
                for n in N for d2 in D), "MinConsWork")
        m.addConstrs(( sw[n,d1,d2] <= tw[n,d1-1,d2]
                for n in N for d1 in D for d2 in D if d1>0), "MinConsWork")
        m.addConstrs(( sw[n,d1,d2] <= 1-p[n,d1]
                for n in N for d1 in D for d2 in D if d1>0), "MinConsWork")
        m.addConstrs(( sw[n,d1,d2] >= tw[n,d1-1,d2]-p[n,d1]
                for n in N for d1 in D for d2 in D if d1>0), "MinConsWork")
        m.addConstr(
                (quicksum(sw[n,d1,d2]*(CW_min[n]-1-d2) for n in N for d1 in D for d2 in range(CW_min[n]-1))==0),"minconswork"
                )
        
        m.addConstrs(( tw1[n,0,s,0] == o[n,0,s] for n in N for s in S), "MaxConsSameShift")
        m.addConstrs(( tw1[n,d1+1,s,0] <= o[n,d1+1,s]
                for s in S for n in N for d1 in D if d1<len(D)-1), "MaxConsSameShift")
        m.addConstrs(( tw1[n,d1+1,s,0] <= 1-o[n,d1,s]
                for s in S for n in N for d1 in D if d1<len(D)-1), "MaxConsSameShift")
        m.addConstrs(( tw1[n,d1+1,s,0] >= o[n,d1+1,s]-o[n,d1,s]
                for s in S for n in N for d1 in D if d1<len(D)-1), "MaxConsSameShift")
        
        m.addConstrs(( tw1[n,0,s,d2] == 0
                for s in S for n in N for d2 in D if d2>0), "MaxConsSameShift")
        m.addConstrs(( tw1[n,d1,s,d2] <= tw1[n,d1-1,s,d2-1]
                for s in S for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsSameShift")
        m.addConstrs(( tw1[n,d1,s,d2] <= o[n,d1,s]
                for s in S for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsSameShift")
        m.addConstrs(( tw1[n,d1,s,d2] >= o[n,d1,s]+tw1[n,d1-1,s,d2-1]-1
                for s in S for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsSameShift")
        m.addConstr(
                (quicksum(tw1[n,d1,s,d2] for s in S for n in N for d1 in D for d2 in range(CWs_max[n,s],len(D)))==0),"maxconssameshift"
                )
        
        m.addConstrs(( sw1[n,0,s,d2] == 0
                for s in S for n in N for d2 in D), "MinConsSameShift")
        m.addConstrs(( sw1[n,d1,s,d2] <= tw1[n,d1-1,s,d2]
                for s in S for n in N for d1 in D for d2 in D if d1>0), "MinConsSameShift")
        m.addConstrs(( sw1[n,d1,s,d2] <= 1-o[n,d1,s]
                for s in S for n in N for d1 in D for d2 in D if d1>0), "MinConsSameShift")
        m.addConstrs(( sw1[n,d1,s,d2] >= tw1[n,d1-1,s,d2]-o[n,d1,s]
                for s in S for n in N for d1 in D for d2 in D if d1>0), "MinConsSameShift")
        m.addConstr(
                (quicksum(sw1[n,d1,s,d2]*(CWs_min[n,s]-1-d2) for s in S for n in N for d1 in D for d2 in range(CWs_min[n,s]-1))==0),"minconssameshift"
                )
        
        m.addConstrs(( tf[n,0,0] == 1-p[n,0] for n in N), "MaxConsFree")
        m.addConstrs(( tf[n,d1+1,0] <= p[n,d1]
                for n in N for d1 in D if d1<len(D)-1), "MaxConsFree")
        m.addConstrs(( tf[n,d1+1,0] <= 1-p[n,d1+1]
                for n in N for d1 in D if d1<len(D)-1), "MaxConsFree")
        m.addConstrs(( tf[n,d1+1,0] >= p[n,d1]-p[n,d1+1]
                for n in N for d1 in D if d1<len(D)-1), "MaxConsFree")
        
        m.addConstrs(( tf[n,0,d2] == 0
                for n in N for d2 in D if d2>0), "MaxConsFree")
        m.addConstrs(( tf[n,d1,d2] <= tf[n,d1-1,d2-1]
                for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsFree")
        m.addConstrs(( tf[n,d1,d2] <= 1-p[n,d1]
                for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsFree")
        m.addConstrs(( tf[n,d1,d2] >= tf[n,d1-1,d2-1]-p[n,d1]
                for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsFree")
        m.addConstr(
                (quicksum(tf[n,d1,d2] for n in N for d1 in D for d2 in range(CF_max[n],len(D)))==0),"maxconsfree"
                )
        
        m.addConstrs(( sf[n,0,d2] == 0
                for n in N for d2 in D), "MinConsFree")
        m.addConstrs(( sf[n,d1,d2] <= tf[n,d1-1,d2]
                for n in N for d1 in D for d2 in D if d1>0), "MinConsFree")
        m.addConstrs(( sf[n,d1,d2] <= p[n,d1]
                for n in N for d1 in D for d2 in D if d1>0), "MinConsFree")
        m.addConstrs(( sf[n,d1,d2] >= tf[n,d1-1,d2]+p[n,d1]-1
                for n in N for d1 in D for d2 in D if d1>0), "MinConsFree")
        m.addConstr(
                (quicksum(sf[n,d1,d2]*(CF_min[n]-1-d2) for n in N for d1 in D for d2 in range(CF_min[n]-1))==0),"minconsfree"
                )
        
        m.setParam(GRB.Param.PoolSolutions, numSol)
        m.setParam(GRB.Param.PoolSearchMode, 2)
        m.optimize()
        nSolutions = m.SolCount
        print('Number of solutions found: ' + str(nSolutions))
        if (m.status==GRB.Status.INFEASIBLE):
            m.computeIIS()
            print('\nThe following constraint(s) cannot be satisfied:')
            for c in m.getConstrs():
                if c.IISConstr:
                    print('%s' % c.constrName)
        if m.status==GRB.Status.OPTIMAL:
            m.write('m.sol')
        
        for i in range(nSolutions):
            m.setParam(GRB.Param.SolutionNumber,i)
#            print(m.ObjVal)
            p1=m.getAttr('xn', p)
            solution=m.getAttr('xn', o)
            tmp_sol=np.zeros([num_nurses,num_days,num_shifts])
            for key in solution:
                tmp_sol[key]=round(solution[key])
            tmp_sol=tmp_sol.astype(np.int64)
            with open(os.path.join(direc, "sol_wo_skill_N"+str(num_nurses)+"_hard"+str(i)+".csv") ,"w+") as my_csv:
                csvWriter = csv.writer(my_csv,delimiter=',')
                
                row=['']
                for j in range(num_days):
                    row.extend(['D'+str(j)]*num_shifts)
                csvWriter.writerow(row)
                
                row=[]
                for j in range(num_shifts):
                    row.extend(['S'+str(j)])
                csvWriter.writerow([''] + row*num_days)
                
                tmp_sol.astype(int)
                for j in range(num_nurses):
                    row=['N'+str(j)]
                    row.extend(tmp_sol[j].flatten())
                    csvWriter.writerow(row)
            
    except GurobiError as e:
        print('Error code ' + str(e.errno) + ": " + str(e))
    
    except AttributeError:
        print('Encountered an attribute error')      
        


    
def generateSam_v2(num_nurses,num_days,num_shifts,orderingNotImp,numSol,direc,num_constrType,constrList,bounds,bounds_learned,nurse_preference):
    
    N=list(range(num_nurses))
    D=list(range(num_days))
    S=list(range(num_shifts))
    variables=[D,S,N]
    
    #Forbidden Shift Successions
    F=[(num_shifts-1,0)]
        
    #Weekends
    W=[(5,6)]
    
    #minimum number of nurses required for covering a shift s on day d
    C_day_min={}
    for d in D:
        C_day_min[d]=3
                
    #Maximum number of nurses required for covering a shift s on day d
    C_day_max={}
    for d in D:
        C_day_max[d]=8
            
    #minimum number of nurses required for covering a shift s on day d
    C={}
    for d in D:
        for s in S:
            C[d,s]=1
                
    #Maximum number of nurses required for covering a shift s on day d
    C_max={}
    for d in D:
        for s in S:
            C_max[d,s]=3
    
    #Min Consecutive Working Days
    CW_min={}
    for n in N:
        CW_min[n]=1
        
    #Max Consecutive Working Days
    CW_max={}
    for n in N:
        CW_max[n]=5
       
    #Min Consecutive Working Days in shift s
    CWs_min={}
    for n in N:
        for s in S:
            CWs_min[n,s]=1
        
    #Max Consecutive Working Days in shift s
    CWs_max={}
    for n in N:
        for s in S:
            CWs_max[n,s]=4
     
    #Min Consecutive Free Days
    CF_min={}
    for n in N:
        CF_min[n]=2
        
    #Max Consecutive Free Days
    CF_max={}
    for n in N:
        CF_max[n]=4
    
    #Max Working Days
    WD_max={}
    for n in N:
        WD_max[n]=6
        
    #Min Working Days
    WD_min={}
    for n in N:
        WD_min[n]=2

    rSam=0
    tSam=0
    lSam=0
        
    try:
        m=Model("nspSolver")
        m.setParam(GRB.Param.OutputFlag,0)
        ########### Decision Variables #############
        o = m.addVars(N,D,S, vtype=GRB.BINARY, name="o")
        p = m.addVars(N,D, vtype=GRB.BINARY, name="p")
        S1=m.addVars(D,S,vtype=GRB.CONTINUOUS, lb=-100 , name="S1")
        S11=m.addVars(D,S,vtype=GRB.CONTINUOUS, name="S11")
        tw = m.addVars(N,D,D, vtype=GRB.CONTINUOUS, name="tw")
        sw = m.addVars(N,D,D, vtype=GRB.CONTINUOUS, name="sw")
        tw1 = m.addVars(N,D,S,D, vtype=GRB.CONTINUOUS, name="tw1")
        sw1 = m.addVars(N,D,S,D, vtype=GRB.CONTINUOUS, name="sw1")
        
        tf = m.addVars(N,D,D, vtype=GRB.CONTINUOUS, name="tf")
        sf = m.addVars(N,D,D, vtype=GRB.CONTINUOUS, name="sf")
        
        ########### Required Constraints #############
        m.addConstrs(
                (o.sum(n,d,'*') == p[n,d] 
                for n in N for d in D),"po")
        
        
        ########### Hard Constraints #############
        m.addConstrs(
                (p[n,d]<=1 for n in N for d in D),"singleAssignmentPerDay")
        m.addConstrs(
                (p.sum('*',d) >= C_day_min[d] for d in D),"minNurses")
        m.addConstrs(
                (p.sum('*',d) <= C_day_max[d] for d in D),"maxNurses")
        m.addConstrs(
                (o.sum('*',d,s) >= C[d,s]
                for d in D for s in S),"underStaffing")
        m.addConstrs(
                (o.sum('*',d,s) <= C_max[d,s]
                for d in D for s in S),"overStaffing")
        m.addConstrs(
                (o[n,d,f1]+o[n,d+1,f2]<=1 
                 for n in N for d in D[:-1] for f1,f2 in F),"shiftTypeSuccession")
        m.addConstrs(
                (quicksum(p[n,d] for d in D) <= WD_max[n] for n in N),"maxWorkingDay")
        m.addConstrs(
                (quicksum(p[n,d] for d in D) >= WD_min[n] for n in N),"minWorkingDay")
        
        m.addConstrs(( tw[n,0,0] == p[n,0] for n in N), "MaxConsWork")
        m.addConstrs(( tw[n,d1+1,0] <= p[n,d1+1]
                for n in N for d1 in D if d1<len(D)-1), "MaxConsWork")
        m.addConstrs(( tw[n,d1+1,0] <= 1-p[n,d1]
                for n in N for d1 in D if d1<len(D)-1), "MaxConsWork")
        m.addConstrs(( tw[n,d1+1,0] >= p[n,d1+1]-p[n,d1]
                for n in N for d1 in D if d1<len(D)-1), "MaxConsWork")
        
        m.addConstrs(( tw[n,0,d2] == 0
                for n in N for d2 in D if d2>0), "MaxConsWork")
        m.addConstrs(( tw[n,d1,d2] <= tw[n,d1-1,d2-1]
                for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsWork")
        m.addConstrs(( tw[n,d1,d2] <= p[n,d1]
                for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsWork")
        m.addConstrs(( tw[n,d1,d2] >= p[n,d1]+tw[n,d1-1,d2-1]-1
                for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsWork")
        m.addConstr(
                (quicksum(tw[n,d1,d2] for n in N for d1 in D for d2 in range(CW_max[n],len(D)))==0),"maxconswork"
                )
        
        m.addConstrs(( sw[n,0,d2] == 0
                for n in N for d2 in D), "MinConsWork")
        m.addConstrs(( sw[n,d1,d2] <= tw[n,d1-1,d2]
                for n in N for d1 in D for d2 in D if d1>0), "MinConsWork")
        m.addConstrs(( sw[n,d1,d2] <= 1-p[n,d1]
                for n in N for d1 in D for d2 in D if d1>0), "MinConsWork")
        m.addConstrs(( sw[n,d1,d2] >= tw[n,d1-1,d2]-p[n,d1]
                for n in N for d1 in D for d2 in D if d1>0), "MinConsWork")
        m.addConstr(
                (quicksum(sw[n,d1,d2]*(CW_min[n]-1-d2) for n in N for d1 in D for d2 in range(CW_min[n]-1))==0),"minconswork"
                )
        
        m.addConstrs(( tw1[n,0,s,0] == o[n,0,s] for n in N for s in S), "MaxConsSameShift")
        m.addConstrs(( tw1[n,d1+1,s,0] <= o[n,d1+1,s]
                for s in S for n in N for d1 in D if d1<len(D)-1), "MaxConsSameShift")
        m.addConstrs(( tw1[n,d1+1,s,0] <= 1-o[n,d1,s]
                for s in S for n in N for d1 in D if d1<len(D)-1), "MaxConsSameShift")
        m.addConstrs(( tw1[n,d1+1,s,0] >= o[n,d1+1,s]-o[n,d1,s]
                for s in S for n in N for d1 in D if d1<len(D)-1), "MaxConsSameShift")
        
        m.addConstrs(( tw1[n,0,s,d2] == 0
                for s in S for n in N for d2 in D if d2>0), "MaxConsSameShift")
        m.addConstrs(( tw1[n,d1,s,d2] <= tw1[n,d1-1,s,d2-1]
                for s in S for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsSameShift")
        m.addConstrs(( tw1[n,d1,s,d2] <= o[n,d1,s]
                for s in S for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsSameShift")
        m.addConstrs(( tw1[n,d1,s,d2] >= o[n,d1,s]+tw1[n,d1-1,s,d2-1]-1
                for s in S for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsSameShift")
        m.addConstr(
                (quicksum(tw1[n,d1,s,d2] for s in S for n in N for d1 in D for d2 in range(CWs_max[n,s],len(D)))==0),"maxconssameshift"
                )
        
        m.addConstrs(( sw1[n,0,s,d2] == 0
                for s in S for n in N for d2 in D), "MinConsSameShift")
        m.addConstrs(( sw1[n,d1,s,d2] <= tw1[n,d1-1,s,d2]
                for s in S for n in N for d1 in D for d2 in D if d1>0), "MinConsSameShift")
        m.addConstrs(( sw1[n,d1,s,d2] <= 1-o[n,d1,s]
                for s in S for n in N for d1 in D for d2 in D if d1>0), "MinConsSameShift")
        m.addConstrs(( sw1[n,d1,s,d2] >= tw1[n,d1-1,s,d2]-o[n,d1,s]
                for s in S for n in N for d1 in D for d2 in D if d1>0), "MinConsSameShift")
        m.addConstr(
                (quicksum(sw1[n,d1,s,d2]*(CWs_min[n,s]-1-d2) for s in S for n in N for d1 in D for d2 in range(CWs_min[n,s]-1))==0),"minconssameshift"
                )
        
        m.addConstrs(( tf[n,0,0] == 1-p[n,0] for n in N), "MaxConsFree")
        m.addConstrs(( tf[n,d1+1,0] <= p[n,d1]
                for n in N for d1 in D if d1<len(D)-1), "MaxConsFree")
        m.addConstrs(( tf[n,d1+1,0] <= 1-p[n,d1+1]
                for n in N for d1 in D if d1<len(D)-1), "MaxConsFree")
        m.addConstrs(( tf[n,d1+1,0] >= p[n,d1]-p[n,d1+1]
                for n in N for d1 in D if d1<len(D)-1), "MaxConsFree")
        
        m.addConstrs(( tf[n,0,d2] == 0
                for n in N for d2 in D if d2>0), "MaxConsFree")
        m.addConstrs(( tf[n,d1,d2] <= tf[n,d1-1,d2-1]
                for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsFree")
        m.addConstrs(( tf[n,d1,d2] <= 1-p[n,d1]
                for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsFree")
        m.addConstrs(( tf[n,d1,d2] >= tf[n,d1-1,d2-1]-p[n,d1]
                for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsFree")
        m.addConstr(
                (quicksum(tf[n,d1,d2] for n in N for d1 in D for d2 in range(CF_max[n],len(D)))==0),"maxconsfree"
                )
        
        m.addConstrs(( sf[n,0,d2] == 0
                for n in N for d2 in D), "MinConsFree")
        m.addConstrs(( sf[n,d1,d2] <= tf[n,d1-1,d2]
                for n in N for d1 in D for d2 in D if d1>0), "MinConsFree")
        m.addConstrs(( sf[n,d1,d2] <= p[n,d1]
                for n in N for d1 in D for d2 in D if d1>0), "MinConsFree")
        m.addConstrs(( sf[n,d1,d2] >= tf[n,d1-1,d2]+p[n,d1]-1
                for n in N for d1 in D for d2 in D if d1>0), "MinConsFree")
        m.addConstr(
                (quicksum(sf[n,d1,d2]*(CF_min[n]-1-d2) for n in N for d1 in D for d2 in range(CF_min[n]-1))==0),"minconsfree"
                )
        
#        for i in range(len(nurse_preference)):
#            m.addConstr((o[i,nurse_preference[i][0],nurse_preference[i][1]] == 0),"nursePref")
        
        m.setParam(GRB.Param.PoolSolutions, numSol)
        m.setParam(GRB.Param.PoolSearchMode, 2)
        m.optimize()
        nSolutions = m.SolCount
        print('Number of solutions found: ' + str(nSolutions))
        if (m.status==GRB.Status.INFEASIBLE):
            m.computeIIS()
            print('\nThe following constraint(s) cannot be satisfied:')
            for c in m.getConstrs():
                if c.IISConstr:
                    print('%s' % c.constrName)
        if m.status==GRB.Status.OPTIMAL:
            m.write('m.sol')
        
        constrRej1=np.zeros(num_constrType)
        constrRej2=np.zeros(num_constrType)
        for i in range(nSolutions):
            m.setParam(GRB.Param.SolutionNumber,i)
#            print(m.ObjVal)
            solution=m.getAttr('xn', o)
            tmp=np.zeros([num_nurses,num_days,num_shifts])
            for key in solution:
                tmp[key]=round(solution[key])
            tSample=np.swapaxes(np.swapaxes(tmp,0,1),1,2)
            accept=1
            for c in range(num_constrType):
                if sum(bounds[c,:])>0:
                    subset=constrList[c]
                    newset=subset[0]+subset[1]
                    idTensor=cf.tensorIndicator(tSample,newset, variables)
                    sumSet = range(len(subset[0]),len(newset))
                    sumTensor_max,sumTensor_min=cf.tensorSum(idTensor,sumSet, np.array(variables)[list(newset)],0)
                    
                    if (bounds[c][0]!=0 and bounds[c][0]>sumTensor_min) or (bounds[c][1]!=0 and bounds[c][1]<sumTensor_max):
                        accept=0
                        constrRej1[c]+=1
                        break
                    if len(set(subset[1]))==1 and len(set(orderingNotImp) & set(subset[1]))==0:
                        minConsZero,maxConsZero,minConsNonZero,maxConsNonZero = cf.tensorConsZero(idTensor,sumSet, np.array(variables)[list(newset)])
                        if (bounds[c][2]!=0 and bounds[c][2]>minConsZero) or (bounds[c][3]!=0 and bounds[c][3]<maxConsZero) and (bounds[c][4]!=0 and bounds[c][4]>minConsNonZero) or (bounds[c][5]!=0 and bounds[c][5]<maxConsNonZero):
                            accept=0
                            constrRej1[c]+=1
                            break         
            if accept==0:
                rSam+=1
            else:
                tSam+=1
            
            accept=1
            for c in range(num_constrType):
                if sum(bounds_learned[c,:])>0:
                    subset=constrList[c]
                    newset=subset[0]+subset[1]
                    idTensor=cf.tensorIndicator(tSample,newset, variables)
                    sumSet = range(len(subset[0]),len(newset))
                    sumTensor_max,sumTensor_min=cf.tensorSum(idTensor,sumSet, np.array(variables)[list(newset)],0)
    #                print(bounds_learned)
                    if (bounds_learned[c][0]!=0 and bounds_learned[c][0]>sumTensor_min) or (bounds_learned[c][1]!=0 and bounds_learned[c][1]<sumTensor_max):
#                        print("here:",sumTensor_min,sumTensor_max)
                        accept=0
                        constrRej2[c]+=1
                        break
                    if len(set(subset[1]))==1 and len(set(orderingNotImp) & set(subset[1]))==0:
                        minConsZero,maxConsZero,minConsNonZero,maxConsNonZero = cf.tensorConsZero(idTensor,sumSet, np.array(variables)[list(newset)])
                        if (bounds_learned[c][2]!=0 and bounds_learned[c][2]>minConsZero) or (bounds_learned[c][3]!=0 and bounds_learned[c][3]<maxConsZero) and (bounds_learned[c][4]!=0 and bounds_learned[c][4]>minConsNonZero) or (bounds_learned[c][5]!=0 and bounds_learned[c][5]<maxConsNonZero):
#                            print("there:",minConsZero,maxConsZero,minConsNonZero,maxConsNonZero)
                            accept=0
                            constrRej2[c]+=1
                            break
            if accept==0:
                lSam+=1
        return tSam,rSam,lSam,constrRej1,constrRej2
            
    except GurobiError as e:
        print('Error code ' + str(e.errno) + ": " + str(e))
    
    except AttributeError:
        print('Encountered an attribute error')
  

def generateSolution_v2(num_nurses,num_days,num_shifts,orderingNotImp,numSol,direc,nurse_preference):
    
    N=list(range(num_nurses))
    D=list(range(num_days))
    S=list(range(num_shifts))
    variables=[D,S,N]
    
    #Forbidden Shift Successions
    F=[(num_shifts-1,0)]
        
    #Weekends
    W=[(5,6)]
    
    #minimum number of nurses required for covering a shift s on day d
    C_day_min={}
    for d in D:
        C_day_min[d]=3
                
    #Maximum number of nurses required for covering a shift s on day d
    C_day_max={}
    for d in D:
        C_day_max[d]=8
            
    #minimum number of nurses required for covering a shift s on day d
    C={}
    for d in D:
        for s in S:
            C[d,s]=1
                
    #Maximum number of nurses required for covering a shift s on day d
    C_max={}
    for d in D:
        for s in S:
            C_max[d,s]=3
    
    #Min Consecutive Working Days
    CW_min={}
    for n in N:
        CW_min[n]=1
        
    #Max Consecutive Working Days
    CW_max={}
    for n in N:
        CW_max[n]=5
       
    #Min Consecutive Working Days in shift s
    CWs_min={}
    for n in N:
        for s in S:
            CWs_min[n,s]=1
        
    #Max Consecutive Working Days in shift s
    CWs_max={}
    for n in N:
        for s in S:
            CWs_max[n,s]=4
     
    #Min Consecutive Free Days
    CF_min={}
    for n in N:
        CF_min[n]=2
        
    #Max Consecutive Free Days
    CF_max={}
    for n in N:
        CF_max[n]=4
    
    #Max Working Days
    WD_max={}
    for n in N:
        WD_max[n]=6
        
    #Min Working Days
    WD_min={}
    for n in N:
        WD_min[n]=2
        
    try:
        m=Model("nspSolver")
        m.setParam(GRB.Param.OutputFlag,0)
        ########### Decision Variables #############
        o = m.addVars(N,D,S, vtype=GRB.BINARY, name="o")
        p = m.addVars(N,D, vtype=GRB.BINARY, name="p")
        S1=m.addVars(D,S,vtype=GRB.CONTINUOUS, lb=-100 , name="S1")
        S11=m.addVars(D,S,vtype=GRB.CONTINUOUS, name="S11")
        tw = m.addVars(N,D,D, vtype=GRB.CONTINUOUS, name="tw")
        sw = m.addVars(N,D,D, vtype=GRB.CONTINUOUS, name="sw")
        tw1 = m.addVars(N,D,S,D, vtype=GRB.CONTINUOUS, name="tw1")
        sw1 = m.addVars(N,D,S,D, vtype=GRB.CONTINUOUS, name="sw1")
        
        tf = m.addVars(N,D,D, vtype=GRB.CONTINUOUS, name="tf")
        sf = m.addVars(N,D,D, vtype=GRB.CONTINUOUS, name="sf")
        
        ########### Required Constraints #############
        m.addConstrs(
                (o.sum(n,d,'*') == p[n,d] 
                for n in N for d in D),"po")
        
        
        ########### Hard Constraints #############
        m.addConstrs(
                (p[n,d]<=1 for n in N for d in D),"singleAssignmentPerDay")
        m.addConstrs(
                (p.sum('*',d) >= C_day_min[d] for d in D),"minNurses")
        m.addConstrs(
                (p.sum('*',d) <= C_day_max[d] for d in D),"maxNurses")
        m.addConstrs(
                (o.sum('*',d,s) >= C[d,s]
                for d in D for s in S),"underStaffing")
        m.addConstrs(
                (o.sum('*',d,s) <= C_max[d,s]
                for d in D for s in S),"overStaffing")
        m.addConstrs(
                (o[n,d,f1]+o[n,d+1,f2]<=1 
                 for n in N for d in D[:-1] for f1,f2 in F),"shiftTypeSuccession")
        m.addConstrs(
                (quicksum(p[n,d] for d in D) <= WD_max[n] for n in N),"maxWorkingDay")
        m.addConstrs(
                (quicksum(p[n,d] for d in D) >= WD_min[n] for n in N),"minWorkingDay")
        
        m.addConstrs(( tw[n,0,0] == p[n,0] for n in N), "MaxConsWork")
        m.addConstrs(( tw[n,d1+1,0] <= p[n,d1+1]
                for n in N for d1 in D if d1<len(D)-1), "MaxConsWork")
        m.addConstrs(( tw[n,d1+1,0] <= 1-p[n,d1]
                for n in N for d1 in D if d1<len(D)-1), "MaxConsWork")
        m.addConstrs(( tw[n,d1+1,0] >= p[n,d1+1]-p[n,d1]
                for n in N for d1 in D if d1<len(D)-1), "MaxConsWork")
        
        m.addConstrs(( tw[n,0,d2] == 0
                for n in N for d2 in D if d2>0), "MaxConsWork")
        m.addConstrs(( tw[n,d1,d2] <= tw[n,d1-1,d2-1]
                for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsWork")
        m.addConstrs(( tw[n,d1,d2] <= p[n,d1]
                for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsWork")
        m.addConstrs(( tw[n,d1,d2] >= p[n,d1]+tw[n,d1-1,d2-1]-1
                for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsWork")
        m.addConstr(
                (quicksum(tw[n,d1,d2] for n in N for d1 in D for d2 in range(CW_max[n],len(D)))==0),"maxconswork"
                )
        
        m.addConstrs(( sw[n,0,d2] == 0
                for n in N for d2 in D), "MinConsWork")
        m.addConstrs(( sw[n,d1,d2] <= tw[n,d1-1,d2]
                for n in N for d1 in D for d2 in D if d1>0), "MinConsWork")
        m.addConstrs(( sw[n,d1,d2] <= 1-p[n,d1]
                for n in N for d1 in D for d2 in D if d1>0), "MinConsWork")
        m.addConstrs(( sw[n,d1,d2] >= tw[n,d1-1,d2]-p[n,d1]
                for n in N for d1 in D for d2 in D if d1>0), "MinConsWork")
        m.addConstr(
                (quicksum(sw[n,d1,d2]*(CW_min[n]-d2) for n in N for d1 in D for d2 in range(CW_min[n]))==0),"minconswork"
                )
        
        m.addConstrs(( tw1[n,0,s,0] == o[n,0,s] for n in N for s in S), "MaxConsSameShift")
        m.addConstrs(( tw1[n,d1+1,s,0] <= o[n,d1+1,s]
                for s in S for n in N for d1 in D if d1<len(D)-1), "MaxConsSameShift")
        m.addConstrs(( tw1[n,d1+1,s,0] <= 1-o[n,d1,s]
                for s in S for n in N for d1 in D if d1<len(D)-1), "MaxConsSameShift")
        m.addConstrs(( tw1[n,d1+1,s,0] >= o[n,d1+1,s]-o[n,d1,s]
                for s in S for n in N for d1 in D if d1<len(D)-1), "MaxConsSameShift")
        
        m.addConstrs(( tw1[n,0,s,d2] == 0
                for s in S for n in N for d2 in D if d2>0), "MaxConsSameShift")
        m.addConstrs(( tw1[n,d1,s,d2] <= tw1[n,d1-1,s,d2-1]
                for s in S for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsSameShift")
        m.addConstrs(( tw1[n,d1,s,d2] <= o[n,d1,s]
                for s in S for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsSameShift")
        m.addConstrs(( tw1[n,d1,s,d2] >= o[n,d1,s]+tw1[n,d1-1,s,d2-1]-1
                for s in S for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsSameShift")
        m.addConstr(
                (quicksum(tw1[n,d1,s,d2] for s in S for n in N for d1 in D for d2 in range(CWs_max[n,s],len(D)))==0),"maxconssameshift"
                )
        
        m.addConstrs(( sw1[n,0,s,d2] == 0
                for s in S for n in N for d2 in D), "MinConsSameShift")
        m.addConstrs(( sw1[n,d1,s,d2] <= tw1[n,d1-1,s,d2]
                for s in S for n in N for d1 in D for d2 in D if d1>0), "MinConsSameShift")
        m.addConstrs(( sw1[n,d1,s,d2] <= 1-o[n,d1,s]
                for s in S for n in N for d1 in D for d2 in D if d1>0), "MinConsSameShift")
        m.addConstrs(( sw1[n,d1,s,d2] >= tw1[n,d1-1,s,d2]-o[n,d1,s]
                for s in S for n in N for d1 in D for d2 in D if d1>0), "MinConsSameShift")
        m.addConstr(
                (quicksum(sw1[n,d1,s,d2]*(CWs_min[n,s]-d2) for s in S for n in N for d1 in D for d2 in range(CWs_min[n,s]))==0),"minconssameshift"
                )
        
        m.addConstrs(( tf[n,0,0] == 1-p[n,0] for n in N), "MaxConsFree")
        m.addConstrs(( tf[n,d1+1,0] <= p[n,d1]
                for n in N for d1 in D if d1<len(D)-1), "MaxConsFree")
        m.addConstrs(( tf[n,d1+1,0] <= 1-p[n,d1+1]
                for n in N for d1 in D if d1<len(D)-1), "MaxConsFree")
        m.addConstrs(( tf[n,d1+1,0] >= p[n,d1]-p[n,d1+1]
                for n in N for d1 in D if d1<len(D)-1), "MaxConsFree")
        
        m.addConstrs(( tf[n,0,d2] == 0
                for n in N for d2 in D if d2>0), "MaxConsFree")
        m.addConstrs(( tf[n,d1,d2] <= tf[n,d1-1,d2-1]
                for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsFree")
        m.addConstrs(( tf[n,d1,d2] <= 1-p[n,d1]
                for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsFree")
        m.addConstrs(( tf[n,d1,d2] >= tf[n,d1-1,d2-1]-p[n,d1]
                for n in N for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsFree")
        m.addConstr(
                (quicksum(tf[n,d1,d2] for n in N for d1 in D for d2 in range(CF_max[n],len(D)))==0),"maxconsfree"
                )
        
        m.addConstrs(( sf[n,0,d2] == 0
                for n in N for d2 in D), "MinConsFree")
        m.addConstrs(( sf[n,d1,d2] <= tf[n,d1-1,d2]
                for n in N for d1 in D for d2 in D if d1>0), "MinConsFree")
        m.addConstrs(( sf[n,d1,d2] <= p[n,d1]
                for n in N for d1 in D for d2 in D if d1>0), "MinConsFree")
        m.addConstrs(( sf[n,d1,d2] >= tf[n,d1-1,d2]+p[n,d1]-1
                for n in N for d1 in D for d2 in D if d1>0), "MinConsFree")
        m.addConstr(
                (quicksum(sf[n,d1,d2]*(CF_min[n]-d2) for n in N for d1 in D for d2 in range(CF_min[n]))==0),"minconsfree"
                )
        
        for i in range(len(nurse_preference)):
            m.addConstr((o[i,nurse_preference[i][0],nurse_preference[i][1]] == 0),"nursePref")
        
        m.setParam(GRB.Param.PoolSolutions, numSol)
        m.setParam(GRB.Param.PoolSearchMode, 2)
        m.optimize()
        nSolutions = m.SolCount
        print('Number of solutions found: ' + str(nSolutions))
        if (m.status==GRB.Status.INFEASIBLE):
            m.computeIIS()
            print('\nThe following constraint(s) cannot be satisfied:')
            for c in m.getConstrs():
                if c.IISConstr:
                    print('%s' % c.constrName)
        if m.status==GRB.Status.OPTIMAL:
            m.write('m.sol')
        
        for i in range(nSolutions):
            m.setParam(GRB.Param.SolutionNumber,i)
#            print(m.ObjVal)
            p1=m.getAttr('xn', p)
            solution=m.getAttr('xn', o)
            tmp_sol=np.zeros([num_nurses,num_days,num_shifts])
            for key in solution:
                tmp_sol[key]=round(solution[key])
            tmp_sol=tmp_sol.astype(np.int64)
            with open(os.path.join(direc, "sol_nursePref_N"+str(num_nurses)+"_hard"+str(i)+".csv") ,"w+") as my_csv:
                csvWriter = csv.writer(my_csv,delimiter=',')
                
                row=['']
                for j in range(num_days):
                    row.extend(['D'+str(j)]*num_shifts)
                csvWriter.writerow(row)
                
                row=[]
                for j in range(num_shifts):
                    row.extend(['S'+str(j)])
                csvWriter.writerow([''] + row*num_days)
                
                tmp_sol.astype(int)
                for j in range(num_nurses):
                    row=['N'+str(j)]
                    row.extend(tmp_sol[j].flatten())
                    csvWriter.writerow(row)
            
    except GurobiError as e:
        print('Error code ' + str(e.errno) + ": " + str(e))
    
    except AttributeError:
        print('Encountered an attribute error')      