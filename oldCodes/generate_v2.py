#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  6 10:22:17 2018

@author: mohit
"""

import csv
import numpy as np
import pandas as pd
from gurobipy import *
import os.path
import itertools as it
import constraintFormulation as cf
import random
    
def generateSample(num_nurses,num_days,num_shifts,orderingNotImp,numSol,num_constrType,constrList,bounds_learned,bounds_tr,nurse_preference,learner):
#    print(bounds_learned)
    N=list(range(num_nurses))
    D=list(range(num_days))
    Ds=list(range(num_days+1))
    S=list(range(num_shifts))
    Ss=list(range(num_shifts+1))
    variables=[D,S,N]
     
        
    #Forbidden Shift Successions
    F=[(num_shifts-1,0)]
        
    #Weekends
    W=[(5,6)]

    rSam=0
    tSam=0
    lSam=0
        
    try:
        m=Model("nspSolver")
        m.setParam(GRB.Param.OutputFlag,0)
        ########### Decision Variables #############
        o = m.addVars(N,D,S, vtype=GRB.BINARY, name="o")
        p = m.addVars(N,D, vtype=GRB.BINARY, name="p")
        q = m.addVars(N,S, vtype=GRB.BINARY, name="p")
        r = m.addVars(S,D, vtype=GRB.BINARY, name="p")
        
        tw = m.addVars(N,D,D, vtype=GRB.BINARY, name="tw")
        sw = m.addVars(N,Ds,D, vtype=GRB.BINARY, name="sw")
        tw1 = m.addVars(N,D,S,D, vtype=GRB.BINARY, name="tw1")
        sw1 = m.addVars(N,Ds,S,D, vtype=GRB.BINARY, name="sw1")
        
        
        tws = m.addVars(N,S,S, vtype=GRB.BINARY, name="tws")
        sws = m.addVars(N,Ss,S, vtype=GRB.BINARY, name="sws")
        tfs = m.addVars(N,S,S, vtype=GRB.BINARY, name="tfs")
        sfs = m.addVars(N,Ss,S, vtype=GRB.BINARY, name="sfs")
        
        tf = m.addVars(N,D,D, vtype=GRB.BINARY, name="tf")
        sf = m.addVars(N,Ds,D, vtype=GRB.BINARY, name="sf")
        tw1f = m.addVars(N,D,S,D, vtype=GRB.BINARY, name="tw1f")
        sw1f = m.addVars(N,Ds,S,D, vtype=GRB.BINARY, name="sw1f")
        
        ########### Required Constraints #############
        m.addConstrs(
                (o.sum(n,d,'*') == p[n,d] 
                for n in N for d in D),"po")
        m.addConstrs((q[n,s] <= o.sum(n,'*',s) for n in N for s in S ),"qo")
        m.addConstrs((q[n,s]*o.sum(n,'*',s) == o.sum(n,'*',s) for n in N for s in S ),"qo")
        m.addConstrs((r[s,d] <= o.sum('*',d,s) for d in D for s in S ),"qo")
        m.addConstrs((r[s,d]*o.sum('*',d,s) == o.sum('*',d,s) for d in D for s in S ),"qo")
        
        ########### Hard Constraints #############
        
        for i in range(len(bounds_learned)):
            if bounds_learned[i,0]>0:
                if constrList[i]==[(0,),(1,)]:
                    m.addConstrs((r.sum('*',d) >= bounds_learned[i,0] for d in D),"constr")
                    
                elif constrList[i]==[(0,),(2,)]:
                    m.addConstrs((p.sum('*',d) >= bounds_learned[i,0] for d in D),"constr")
                    
                elif constrList[i]==[(0,),(1,2)]:
                    m.addConstrs((o.sum('*',d,'*') >= bounds_learned[i,0] for d in D),"constr")
                    
                elif constrList[i]==[(1,),(0,)]:
                    m.addConstrs((r.sum(s,'*') >= bounds_learned[i,0] for s in S),"constr")
                    
                elif constrList[i]==[(1,),(2,)]:
                    m.addConstrs((q.sum('*',s) >= bounds_learned[i,0] for s in S),"constr")
                    
                elif constrList[i]==[(1,),(0,2)]:
                    m.addConstrs((o.sum('*','*',s) >= bounds_learned[i,0] for s in S),"constr")
                    
                elif constrList[i]==[(2,),(0,)]:
                    m.addConstrs((p.sum(n,'*') >= bounds_learned[i,0] for n in N),"constr")
                    
                elif constrList[i]==[(2,),(1,)]:
                    m.addConstrs((q.sum(n,'*') >= bounds_learned[i,0] for n in N),"constr")
                    
                elif constrList[i]==[(2,),(0,1)]:
                    m.addConstrs((o.sum(n,'*','*') >= bounds_learned[i,0] for n in N),"constr")
                    
                elif constrList[i]==[(0,1),(2,)]:
                    m.addConstrs((o.sum('*',d,s) >= bounds_learned[i,0] for d in D for s in S),"constr")
                    
                elif constrList[i]==[(0,2),(1,)]:
                    m.addConstrs((o.sum(n,d,'*') >= bounds_learned[i,0] for d in D for n in N),"constr")
                    
                elif constrList[i]==[(1,2),(0,)]:
                    m.addConstrs((o.sum(n,'*',s) >= bounds_learned[i,0] for n in N for s in S),"constr")
                    
            if bounds_learned[i,1]>0:
                if constrList[i]==[(0,),(1,)]:
                    m.addConstrs((r.sum('*',d) <= bounds_learned[i,1] for d in D),"constr")
                    
                elif constrList[i]==[(0,),(2,)]:
                    m.addConstrs((p.sum('*',d) <= bounds_learned[i,1] for d in D),"constr")
                    
                elif constrList[i]==[(0,),(1,2)]:
                    m.addConstrs((o.sum('*',d,'*') <= bounds_learned[i,1] for d in D),"constr")
                    
                elif constrList[i]==[(1,),(0,)]:
                    m.addConstrs((r.sum(s,'*') <= bounds_learned[i,1] for s in S),"constr")
                    
                elif constrList[i]==[(1,),(2,)]:
                    m.addConstrs((q.sum('*',s) <= bounds_learned[i,1] for s in S),"constr")
                    
                elif constrList[i]==[(1,),(0,2)]:
                    m.addConstrs((o.sum('*','*',s) <= bounds_learned[i,1] for s in S),"constr")
                    
                elif constrList[i]==[(2,),(0,)]:
                    m.addConstrs((p.sum(n,'*') <= bounds_learned[i,1] for n in N),"constr")
                    
                elif constrList[i]==[(2,),(1,)]:
                    m.addConstrs((q.sum(n,'*') <= bounds_learned[i,1] for n in N),"constr")
                    
                elif constrList[i]==[(2,),(0,1)]:
                    m.addConstrs((o.sum(n,'*','*') <= bounds_learned[i,1] for n in N),"constr")
                    
                elif constrList[i]==[(0,1),(2,)]:
                    m.addConstrs((o.sum('*',d,s) <= bounds_learned[i,1] for d in D for s in S),"constr")
                    
                elif constrList[i]==[(0,2),(1,)]:
                    m.addConstrs((o.sum(n,d,'*') <= bounds_learned[i,1] for d in D for n in N),"constr")
                    
                elif constrList[i]==[(1,2),(0,)]:
                    m.addConstrs((o.sum(n,'*',s) <= bounds_learned[i,1] for n in N for s in S),"constr")
            
        if learner==0:
#            m.addConstrs(
#                    (o[n,d,f1]+o[n,d+1,f2]<=1 
#                     for n in N for d in D[:-1] for f1,f2 in F),"shiftTypeSuccession")
            
            for i in range(len(nurse_preference)):
                m.addConstr((o[i,nurse_preference[i][0],nurse_preference[i][1]] == 0),"nursePref")
        
        if bounds_learned[6,5] + bounds_learned[6,4] >0:
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
            if bounds_learned[6,5]>0:
                m.addConstr(
                        (quicksum(tw[n,d1,d2] for n in N for d1 in D for d2 in range(bounds_learned[6,5],len(D)))==0),"maxconswork"
                        )
            if bounds_learned[6,4]>0:
                m.addConstrs(( sw[n,0,d2] == 0
                        for n in N for d2 in D), "MinConsWork")
                m.addConstrs(( sw[n,d1,d2] <= tw[n,d1-1,d2]
                        for n in N for d1 in D for d2 in D if d1>0), "MinConsWork")
                m.addConstrs(( sw[n,d1,d2] <= 1-p[n,d1]
                        for n in N for d1 in D for d2 in D if d1>0), "MinConsWork")
                m.addConstrs(( sw[n,d1,d2] >= tw[n,d1-1,d2]-p[n,d1]
                        for n in N for d1 in D for d2 in D if d1>0), "MinConsWork")
                
                m.addConstrs(( sw[n,num_days,d2] == tw[n,num_days-1,d2]
                        for n in N for d2 in D), "MinConsWork")
                
                m.addConstr(
                        (quicksum(sw[n,d1,d2]*(bounds_learned[6,4]-1-d2) for n in N for d1 in Ds for d2 in range(bounds_learned[6,4]-1))==0),"minconswork"
                        )
        
        if bounds_learned[6,3] + bounds_learned[6,2]>0:
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
            if bounds_learned[6,3]>0:
                m.addConstr(
                        (quicksum(tf[n,d1,d2] for n in N for d1 in D for d2 in range(bounds_learned[6,3],len(D)))==0),"maxconsfree"
                        )
            
            if bounds_learned[6,2]>0:
                m.addConstrs(( sf[n,0,d2] == 0
                        for n in N for d2 in D), "MinConsFree")
                m.addConstrs(( sf[n,d1,d2] <= tf[n,d1-1,d2]
                        for n in N for d1 in D for d2 in D if d1>0), "MinConsFree")
                m.addConstrs(( sf[n,d1,d2] <= p[n,d1]
                        for n in N for d1 in D for d2 in D if d1>0), "MinConsFree")
                m.addConstrs(( sf[n,d1,d2] >= tf[n,d1-1,d2]+p[n,d1]-1
                        for n in N for d1 in D for d2 in D if d1>0), "MinConsFree")
                
                m.addConstrs(( sf[n,num_days,d2] == tf[n,num_days-1,d2]
                        for n in N for d2 in D), "MinConsFree")
                
                m.addConstr(
                        (quicksum(sf[n,d1,d2]*(bounds_learned[6,2]-1-d2) for n in N for d1 in Ds for d2 in range(bounds_learned[6,2]-1))==0),"minconsfree"
                        )
        if bounds_learned[7,5] + bounds_learned[7,4] > 0:
            m.addConstrs(( tws[n,0,0] == q[n,0] for n in N), "MaxConsWork")
            m.addConstrs(( tws[n,s1+1,0] <= q[n,s1+1]
                    for n in N for s1 in S if s1<len(S)-1), "MaxConsWork")
            m.addConstrs(( tws[n,s1+1,0] <= 1-q[n,s1]
                    for n in N for s1 in S if s1<len(S)-1), "MaxConsWork")
            m.addConstrs(( tws[n,s1+1,0] >= q[n,s1+1]-q[n,s1]
                    for n in N for s1 in S if s1<len(S)-1), "MaxConsWork")
            
            m.addConstrs(( tws[n,0,s2] == 0
                    for n in N for s2 in S if s2>0), "MaxConsWork")
            m.addConstrs(( tws[n,s1,s2] <= tws[n,s1-1,s2-1]
                    for n in N for s1 in S for s2 in S if s1>0 if s2>0), "MaxConsWork")
            m.addConstrs(( tws[n,s1,s2] <= q[n,s1]
                    for n in N for s1 in S for s2 in S if s1>0 if s2>0), "MaxConsWork")
            m.addConstrs(( tws[n,s1,s2] >= q[n,s1]+tws[n,s1-1,s2-1]-1
                    for n in N for s1 in S for s2 in S if s1>0 if s2>0), "MaxConsWork")
            
            if bounds_learned[7,5]>0:
                m.addConstr(
                        (quicksum(tws[n,s1,s2] for n in N for s1 in S for s2 in range(bounds_learned[7,5],len(S)))==0),"maxconswork"
                        )
        
            if bounds_learned[7,4]>0:
                m.addConstrs(( sws[n,0,s2] == 0
                        for n in N for s2 in S), "MinConsWork")
                m.addConstrs(( sws[n,s1,s2] <= tws[n,s1-1,s2]
                        for n in N for s1 in S for s2 in S if s1>0), "MinConsWork")
                m.addConstrs(( sws[n,s1,s2] <= 1-q[n,s1]
                        for n in N for s1 in S for s2 in S if s1>0), "MinConsWork")
                m.addConstrs(( sws[n,s1,s2] >= tws[n,s1-1,s2]-q[n,s1]
                        for n in N for s1 in S for s2 in S if s1>0), "MinConsWork")
                
                m.addConstrs(( sws[n,num_shifts,s2] == tws[n,num_shifts-1,s2]
                        for n in N for s2 in S), "MinConsWork")
                
                m.addConstr(
                        (quicksum(sws[n,s1,s2]*(bounds_learned[7,4]-1-s2) for n in N for s1 in Ss for s2 in range(bounds_learned[7,4]-1))==0),"minconswork"
                        )
        
        if bounds_learned[7,3]+bounds_learned[7,2]>0:
            m.addConstrs(( tfs[n,0,0] == 1-q[n,0] for n in N), "MaxConsFree")
            m.addConstrs(( tfs[n,s1+1,0] <= q[n,s1]
                    for n in N for s1 in S if s1<len(S)-1), "MaxConsFree")
            m.addConstrs(( tfs[n,s1+1,0] <= 1-q[n,s1+1]
                    for n in N for s1 in S if s1<len(S)-1), "MaxConsFree")
            m.addConstrs(( tfs[n,s1+1,0] >= q[n,s1]-q[n,s1+1]
                    for n in N for s1 in S if s1<len(S)-1), "MaxConsFree")
            
            m.addConstrs(( tfs[n,0,s2] == 0
                    for n in N for s2 in S if s2>0), "MaxConsFree")
            m.addConstrs(( tfs[n,s1,s2] <= tfs[n,s1-1,s2-1]
                    for n in N for s1 in S for s2 in S if s1>0 if s2>0), "MaxConsFree")
            m.addConstrs(( tfs[n,s1,s2] <= 1-q[n,s1]
                    for n in N for s1 in S for s2 in S if s1>0 if s2>0), "MaxConsFree")
            m.addConstrs(( tfs[n,s1,s2] >= tfs[n,s1-1,s2-1]-q[n,s1]
                    for n in N for s1 in S for s2 in S if s1>0 if s2>0), "MaxConsFree")
            if bounds_learned[7,3]>0:
                m.addConstr(
                        (quicksum(tfs[n,s1,s2] for n in N for s1 in S for s2 in range(bounds_learned[7,3],len(S)))==0),"maxconsfree"
                        )
            
            if bounds_learned[7,2]>0:
                m.addConstrs(( sfs[n,0,s2] == 0
                        for n in N for s2 in S), "MinConsFree")
                m.addConstrs(( sfs[n,s1,s2] <= tfs[n,s1-1,s2]
                        for n in N for s1 in S for s2 in S if s1>0), "MinConsFree")
                m.addConstrs(( sfs[n,s1,s2] <= q[n,s1]
                        for n in N for s1 in S for s2 in S if s1>0), "MinConsFree")
                m.addConstrs(( sfs[n,s1,s2] >= tfs[n,s1-1,s2]+q[n,s1]-1
                        for n in N for s1 in S for s2 in S if s1>0), "MinConsFree")
                
                m.addConstrs(( sfs[n,num_shifts,s2] == tfs[n,num_shifts-1,s2]
                        for n in N for s2 in S), "MinConsWork")
                
                m.addConstr(
                        (quicksum(sfs[n,s1,s2]*(bounds_learned[7,2]-1-s2) for n in N for s1 in Ss for s2 in range(bounds_learned[7,2]-1))==0),"minconsfree"
                        )
                
        if bounds_learned[11,5]+bounds_learned[11,4]>0:
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
            if bounds_learned[11,5]>0:
                m.addConstr(
                        (quicksum(tw1[n,d1,s,d2] for s in S for n in N for d1 in D for d2 in range(bounds_learned[11,5],len(D)))==0),"maxconssameshift"
                        )
        
            if bounds_learned[11,4]>0:
                m.addConstrs(( sw1[n,0,s,d2] == 0
                        for s in S for n in N for d2 in D), "MinConsSameShift")
                m.addConstrs(( sw1[n,d1,s,d2] <= tw1[n,d1-1,s,d2]
                        for s in S for n in N for d1 in D for d2 in D if d1>0), "MinConsSameShift")
                m.addConstrs(( sw1[n,d1,s,d2] <= 1-o[n,d1,s]
                        for s in S for n in N for d1 in D for d2 in D if d1>0), "MinConsSameShift")
                m.addConstrs(( sw1[n,d1,s,d2] >= tw1[n,d1-1,s,d2]-o[n,d1,s]
                        for s in S for n in N for d1 in D for d2 in D if d1>0), "MinConsSameShift")
                
                m.addConstrs(( sw1[n,num_days,s,d2] == tw1[n,num_days-1,s,d2]
                        for n in N for s in S for d2 in D), "MinConsWork")
                
                m.addConstr(
                        (quicksum(sw1[n,d1,s,d2]*(bounds_learned[11,4]-1-d2) for s in S for n in N for d1 in Ds for d2 in range(bounds_learned[11,4]-1))==0),"minconssameshift"
                        )
                
        if bounds_learned[11,3] + bounds_learned[11,2]>0:
            m.addConstrs(( tw1f[n,0,s,0] == 1-o[n,0,s] for n in N for s in S), "MaxConsFree")
            m.addConstrs(( tw1f[n,d1+1,s,0] <= o[n,d1,s]
                    for n in N for s in S for d1 in D if d1<len(D)-1), "MaxConsFree")
            m.addConstrs(( tw1f[n,d1+1,s,0] <= 1-o[n,d1+1,s]
                    for n in N for s in S for d1 in D if d1<len(D)-1), "MaxConsFree")
            m.addConstrs(( tw1f[n,d1+1,s,0] >= o[n,d1,s]-o[n,d1+1,s]
                    for n in N for s in S for d1 in D if d1<len(D)-1), "MaxConsFree")
            
            m.addConstrs(( tw1f[n,0,s,d2] == 0
                    for n in N for s in S for d2 in D if d2>0), "MaxConsFree")
            m.addConstrs(( tw1f[n,d1,s,d2] <= tw1f[n,d1-1,s,d2-1]
                    for n in N for s in S for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsFree")
            m.addConstrs(( tw1f[n,d1,s,d2] <= 1-o[n,d1,s]
                    for n in N for s in S for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsFree")
            m.addConstrs(( tw1f[n,d1,s,d2] >= tw1f[n,d1-1,s,d2-1]-o[n,d1,s]
                    for n in N for s in S for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsFree")
            if bounds_learned[11,3]>0:
                m.addConstr(
                        (quicksum(tw1f[n,d1,s,d2] for n in N for s in S for d1 in D for d2 in range(bounds_learned[11,3],len(D)))==0),"maxconsfree"
                        )
            
            if bounds_learned[11,2]>0:
                m.addConstrs(( sw1f[n,0,s,d2] == 0
                        for n in N for s in S for d2 in D), "MinConsw1free")
                m.addConstrs(( sw1f[n,d1,s,d2] <= tw1f[n,d1-1,s,d2]
                        for n in N for s in S for d1 in D for d2 in D if d1>0), "MinConsw1free")
                m.addConstrs(( sw1f[n,d1,s,d2] <= o[n,d1,s]
                        for n in N for s in S for d1 in D for d2 in D if d1>0), "MinConsw1free")
                m.addConstrs(( sw1f[n,d1,s,d2] >= tw1f[n,d1-1,s,d2]+o[n,d1,s]-1
                        for n in N for s in S for d1 in D for d2 in D if d1>0), "MinConsw1free")
                
                m.addConstrs(( sw1f[n,num_days,s,d2] == tw1f[n,num_days-1,s,d2]
                        for n in N for s in S for d2 in D), "MinConsWork")
                
                m.addConstr(
                        (quicksum(sw1f[n,d1,s,d2]*(bounds_learned[11,2]-1-d2) for n in N for s in S for d1 in Ds for d2 in range(bounds_learned[11,2]-1))==0),"minconsw1free"
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
            solution=m.getAttr('xn', o)
            tmp=np.zeros([num_nurses,num_days,num_shifts])
            for key in solution:
                tmp[key]=round(solution[key])
            tSample=np.swapaxes(np.swapaxes(tmp,0,1),1,2)
#            print(tSample)
            accept=1
            if learner==1:
                for i in range(len(nurse_preference)):
                    if tSample[nurse_preference[i][0],nurse_preference[i][1],i] != 0:
                        accept=0
                        break
            if accept==1:
                for c in range(num_constrType):
                    if sum(bounds_tr[c,:])>0:
                        subset=constrList[c]
                        newset=subset[0]+subset[1]
                        idTensor=cf.tensorIndicator(tSample,newset, variables)
                        sumSet = range(len(subset[0]),len(newset))
                        sumTensor_max,sumTensor_min=cf.tensorSum(idTensor,sumSet, np.array(variables)[list(newset)],0)
                        
                        if (bounds_tr[c][0]!=0 and bounds_tr[c][0]>sumTensor_min) or (bounds_tr[c][1]!=0 and bounds_tr[c][1]<sumTensor_max):
                            accept=0
                            constrRej1[c]+=1
    #                        print(sumTensor_min,sumTensor_max)
    #                        print(tSample)
    #                        print(bounds_tr[c])
                            break
                        if len(set(subset[1]))==1 and len(set(orderingNotImp) & set(subset[1]))==0:
                            minConsZero,maxConsZero,minConsNonZero,maxConsNonZero = cf.tensorConsZero(idTensor,sumSet, np.array(variables)[list(newset)])
                            if (bounds_tr[c][2]!=0 and bounds_tr[c][2]>minConsZero) or (bounds_tr[c][3]!=0 and bounds_tr[c][3]<maxConsZero) and (bounds_tr[c][4]!=0 and bounds_tr[c][4]>minConsNonZero) or (bounds_tr[c][5]!=0 and bounds_tr[c][5]<maxConsNonZero):
                                accept=0
    #                            print(sumTensor_min,sumTensor_max,minConsZero,maxConsZero,minConsNonZero,maxConsNonZero)
    #                            print(bounds_tr[c])
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
                        if (bounds_learned[c][2]!=0 and minConsZero!=0 and bounds_learned[c][2]>minConsZero) or (bounds_learned[c][3]!=0 and maxConsZero!=0 and bounds_learned[c][3]<maxConsZero) and (bounds_learned[c][4]!=0 and minConsNonZero!=0 and bounds_learned[c][4]>minConsNonZero) or (bounds_learned[c][5]!=0 and maxConsNonZero!=0 and bounds_learned[c][5]<maxConsNonZero):
#                            print(c)
#                            print("there:",minConsZero,maxConsZero,minConsNonZero,maxConsNonZero)
#                            print(bounds_learned[c])
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



def generateSolution(num_nurses,num_days,num_shifts,orderingNotImp,numSol,direc,constrList,bounds_learned,nurse_preference):
    
    N=list(range(num_nurses))
    D=list(range(num_days))
    Ds=list(range(num_days+1))
    S=list(range(num_shifts))
    Ss=list(range(num_shifts+1))
    variables=[D,S,N]
    
    #Forbidden Shift Successions
    F=[(num_shifts-1,0)]
        
    #Weekends
    W=[(5,6)]
    
    try:
        m=Model("nspSolver")
        m.setParam(GRB.Param.OutputFlag,0)
        ########### Decision Variables #############
        o = m.addVars(N,D,S, vtype=GRB.BINARY, name="o")
        p = m.addVars(N,D, vtype=GRB.BINARY, name="p")
        q = m.addVars(N,S, vtype=GRB.BINARY, name="p")
        r = m.addVars(S,D, vtype=GRB.BINARY, name="p")
#        S1=m.addVars(D,S,vtype=GRB.CONTINUOUS, lb=-100 , name="S1")
#        S11=m.addVars(D,S,vtype=GRB.CONTINUOUS, name="S11")
        tw = m.addVars(N,D,D, vtype=GRB.BINARY, name="tw")
        sw = m.addVars(N,Ds,D, vtype=GRB.BINARY, name="sw")
        tw1 = m.addVars(N,D,S,D, vtype=GRB.BINARY, name="tw1")
        sw1 = m.addVars(N,Ds,S,D, vtype=GRB.BINARY, name="sw1")
        
        
        tws = m.addVars(N,S,S, vtype=GRB.BINARY, name="tws")
        sws = m.addVars(N,Ss,S, vtype=GRB.BINARY, name="sws")
        tfs = m.addVars(N,S,S, vtype=GRB.BINARY, name="tfs")
        sfs = m.addVars(N,Ss,S, vtype=GRB.BINARY, name="sfs")
        
        tf = m.addVars(N,D,D, vtype=GRB.BINARY, name="tf")
        sf = m.addVars(N,Ds,D, vtype=GRB.BINARY, name="sf")
        tw1f = m.addVars(N,D,S,D, vtype=GRB.BINARY, name="tw1f")
        sw1f = m.addVars(N,Ds,S,D, vtype=GRB.BINARY, name="sw1f")
        
        ########### Required Constraints #############
        m.addConstrs(
                (o.sum(n,d,'*') == p[n,d] 
                for n in N for d in D),"po")
        m.addConstrs((q[n,s] <= o.sum(n,'*',s) for n in N for s in S ),"qo")
        m.addConstrs((q[n,s]*o.sum(n,'*',s) == o.sum(n,'*',s) for n in N for s in S ),"qo")
        m.addConstrs((r[s,d] <= o.sum('*',d,s) for d in D for s in S ),"ro")
        m.addConstrs((r[s,d]*o.sum('*',d,s) == o.sum('*',d,s) for d in D for s in S ),"ro")
        
        ########### Hard Constraints #############
#        print(bounds_learned)
        for i in range(len(bounds_learned)):
            if bounds_learned[i,0]>0:
#                print(bounds_learned[i,0])
                if constrList[i]==[(0,),(1,)]:
                    m.addConstrs((r.sum('*',d) >= bounds_learned[i,0] for d in D),"constr")
                    
                elif constrList[i]==[(0,),(2,)]:
                    m.addConstrs((p.sum('*',d) >= bounds_learned[i,0] for d in D),"constr")
                    
                elif constrList[i]==[(0,),(1,2)]:
                    m.addConstrs((o.sum('*',d,'*') >= bounds_learned[i,0] for d in D),"constr")
                    
                elif constrList[i]==[(1,),(0,)]:
                    m.addConstrs((r.sum(s,'*') >= bounds_learned[i,0] for s in S),"constr")
                    
                elif constrList[i]==[(1,),(2,)]:
                    m.addConstrs((q.sum('*',s) >= bounds_learned[i,0] for s in S),"constr")
                    
                elif constrList[i]==[(1,),(0,2)]:
                    m.addConstrs((o.sum('*','*',s) >= bounds_learned[i,0] for s in S),"constr")
                    
                elif constrList[i]==[(2,),(0,)]:
                    m.addConstrs((p.sum(n,'*') >= bounds_learned[i,0] for n in N),"constr")
                    
                elif constrList[i]==[(2,),(1,)]:
                    m.addConstrs((q.sum(n,'*') >= bounds_learned[i,0] for n in N),"constr")
                    
                elif constrList[i]==[(2,),(0,1)]:
                    m.addConstrs((o.sum(n,'*','*') >= bounds_learned[i,0] for n in N),"constr")
                    
                elif constrList[i]==[(0,1),(2,)]:
                    m.addConstrs((o.sum('*',d,s) >= bounds_learned[i,0] for d in D for s in S),"constr")
                    
                elif constrList[i]==[(0,2),(1,)]:
                    m.addConstrs((o.sum(n,d,'*') >= bounds_learned[i,0] for d in D for n in N),"constr")
                    
                elif constrList[i]==[(1,2),(0,)]:
                    m.addConstrs((o.sum(n,'*',s) >= bounds_learned[i,0] for n in N for s in S),"constr")
                    
            if bounds_learned[i,1]>0:
#                print(bounds_learned[i,1])
                if constrList[i]==[(0,),(1,)]:
                    m.addConstrs((r.sum('*',d) <= bounds_learned[i,1] for d in D),"constr")
                    
                elif constrList[i]==[(0,),(2,)]:
                    m.addConstrs((p.sum('*',d) <= bounds_learned[i,1] for d in D),"constr")
                    
                elif constrList[i]==[(0,),(1,2)]:
                    m.addConstrs((o.sum('*',d,'*') <= bounds_learned[i,1] for d in D),"constr")
                    
                elif constrList[i]==[(1,),(0,)]:
                    m.addConstrs((r.sum(s,'*') <= bounds_learned[i,1] for s in S),"constr")
                    
                elif constrList[i]==[(1,),(2,)]:
                    m.addConstrs((q.sum('*',s) <= bounds_learned[i,1] for s in S),"constr")
                    
                elif constrList[i]==[(1,),(0,2)]:
                    m.addConstrs((o.sum('*','*',s) <= bounds_learned[i,1] for s in S),"constr")
                    
                elif constrList[i]==[(2,),(0,)]:
                    m.addConstrs((p.sum(n,'*') <= bounds_learned[i,1] for n in N),"constr")
                    
                elif constrList[i]==[(2,),(1,)]:
                    m.addConstrs((q.sum(n,'*') <= bounds_learned[i,1] for n in N),"constr")
                    
                elif constrList[i]==[(2,),(0,1)]:
                    m.addConstrs((o.sum(n,'*','*') <= bounds_learned[i,1] for n in N),"constr")
                    
                elif constrList[i]==[(0,1),(2,)]:
                    m.addConstrs((o.sum('*',d,s) <= bounds_learned[i,1] for d in D for s in S),"constr")
                    
                elif constrList[i]==[(0,2),(1,)]:
                    m.addConstrs((o.sum(n,d,'*') <= bounds_learned[i,1] for d in D for n in N),"constr")
                    
                elif constrList[i]==[(1,2),(0,)]:
                    m.addConstrs((o.sum(n,'*',s) <= bounds_learned[i,1] for n in N for s in S),"constr")
            
        
        m.addConstrs(
                (o[n,d,f1]+o[n,d+1,f2]<=1 
                 for n in N for d in D[:-1] for f1,f2 in F),"shiftTypeSuccession")
        
        for i in range(len(nurse_preference)):
            m.addConstr((o[i,nurse_preference[i][0],nurse_preference[i][1]] == 0),"nursePref")
            
        if bounds_learned[6,5] + bounds_learned[6,4] >0:
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
            if bounds_learned[6,5]>0:
                m.addConstr(
                        (quicksum(tw[n,d1,d2] for n in N for d1 in D for d2 in range(bounds_learned[6,5],len(D)))==0),"maxconswork"
                        )
            if bounds_learned[6,4]>0:
                m.addConstrs(( sw[n,0,d2] == 0
                        for n in N for d2 in D), "MinConsWork")
                m.addConstrs(( sw[n,d1,d2] <= tw[n,d1-1,d2]
                        for n in N for d1 in D for d2 in D if d1>0), "MinConsWork")
                m.addConstrs(( sw[n,d1,d2] <= 1-p[n,d1]
                        for n in N for d1 in D for d2 in D if d1>0), "MinConsWork")
                m.addConstrs(( sw[n,d1,d2] >= tw[n,d1-1,d2]-p[n,d1]
                        for n in N for d1 in D for d2 in D if d1>0), "MinConsWork")
                
                m.addConstrs(( sw[n,num_days,d2] == tw[n,num_days-1,d2]
                        for n in N for d2 in D), "MinConsWork")
                
                m.addConstr(
                        (quicksum(sw[n,d1,d2]*(bounds_learned[6,4]-1-d2) for n in N for d1 in Ds for d2 in range(bounds_learned[6,4]-1))==0),"minconswork"
                        )
        
            
        if bounds_learned[6,3] + bounds_learned[6,2]>0:
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
            if bounds_learned[6,3]>0:
                m.addConstr(
                        (quicksum(tf[n,d1,d2] for n in N for d1 in D for d2 in range(bounds_learned[6,3],len(D)))==0),"maxconsfree"
                        )
            
            if bounds_learned[6,2]>0:
                m.addConstrs(( sf[n,0,d2] == 0
                        for n in N for d2 in D), "MinConsFree")
                m.addConstrs(( sf[n,d1,d2] <= tf[n,d1-1,d2]
                        for n in N for d1 in D for d2 in D if d1>0), "MinConsFree")
                m.addConstrs(( sf[n,d1,d2] <= p[n,d1]
                        for n in N for d1 in D for d2 in D if d1>0), "MinConsFree")
                m.addConstrs(( sf[n,d1,d2] >= tf[n,d1-1,d2]+p[n,d1]-1
                        for n in N for d1 in D for d2 in D if d1>0), "MinConsFree")
                
                m.addConstrs(( sf[n,num_days,d2] == tf[n,num_days-1,d2]
                        for n in N for d2 in D), "MinConsFree")
                
                m.addConstr(
                        (quicksum(sf[n,d1,d2]*(bounds_learned[6,2]-1-d2) for n in N for d1 in Ds for d2 in range(bounds_learned[6,2]-1))==0),"minconsfree"
                        )   
        
        if bounds_learned[7,5] + bounds_learned[7,4] > 0:
            m.addConstrs(( tws[n,0,0] == q[n,0] for n in N), "MaxConsWork")
            m.addConstrs(( tws[n,s1+1,0] <= q[n,s1+1]
                    for n in N for s1 in S if s1<len(S)-1), "MaxConsWork")
            m.addConstrs(( tws[n,s1+1,0] <= 1-q[n,s1]
                    for n in N for s1 in S if s1<len(S)-1), "MaxConsWork")
            m.addConstrs(( tws[n,s1+1,0] >= q[n,s1+1]-q[n,s1]
                    for n in N for s1 in S if s1<len(S)-1), "MaxConsWork")
            
            m.addConstrs(( tws[n,0,s2] == 0
                    for n in N for s2 in S if s2>0), "MaxConsWork")
            m.addConstrs(( tws[n,s1,s2] <= tws[n,s1-1,s2-1]
                    for n in N for s1 in S for s2 in S if s1>0 if s2>0), "MaxConsWork")
            m.addConstrs(( tws[n,s1,s2] <= q[n,s1]
                    for n in N for s1 in S for s2 in S if s1>0 if s2>0), "MaxConsWork")
            m.addConstrs(( tws[n,s1,s2] >= q[n,s1]+tws[n,s1-1,s2-1]-1
                    for n in N for s1 in S for s2 in S if s1>0 if s2>0), "MaxConsWork")
            
            if bounds_learned[7,5]>0:
                m.addConstr(
                        (quicksum(tws[n,s1,s2] for n in N for s1 in S for s2 in range(bounds_learned[7,5],len(S)))==0),"maxconswork"
                        )
        
            if bounds_learned[7,4]>0:
                m.addConstrs(( sws[n,0,s2] == 0
                        for n in N for s2 in S), "MinConsWork")
                m.addConstrs(( sws[n,s1,s2] <= tws[n,s1-1,s2]
                        for n in N for s1 in S for s2 in S if s1>0), "MinConsWork")
                m.addConstrs(( sws[n,s1,s2] <= 1-q[n,s1]
                        for n in N for s1 in S for s2 in S if s1>0), "MinConsWork")
                m.addConstrs(( sws[n,s1,s2] >= tws[n,s1-1,s2]-q[n,s1]
                        for n in N for s1 in S for s2 in S if s1>0), "MinConsWork")
                
                m.addConstrs(( sws[n,num_shifts,s2] == tws[n,num_shifts-1,s2]
                        for n in N for s2 in S), "MinConsWork")
                
                m.addConstr(
                        (quicksum(sws[n,s1,s2]*(bounds_learned[7,4]-1-s2) for n in N for s1 in Ss for s2 in range(bounds_learned[7,4]-1))==0),"minconswork"
                        )
        
        if bounds_learned[7,3]+bounds_learned[7,2]>0:
            m.addConstrs(( tfs[n,0,0] == 1-q[n,0] for n in N), "MaxConsFree")
            m.addConstrs(( tfs[n,s1+1,0] <= q[n,s1]
                    for n in N for s1 in S if s1<len(S)-1), "MaxConsFree")
            m.addConstrs(( tfs[n,s1+1,0] <= 1-q[n,s1+1]
                    for n in N for s1 in S if s1<len(S)-1), "MaxConsFree")
            m.addConstrs(( tfs[n,s1+1,0] >= q[n,s1]-q[n,s1+1]
                    for n in N for s1 in S if s1<len(S)-1), "MaxConsFree")
            
            m.addConstrs(( tfs[n,0,s2] == 0
                    for n in N for s2 in S if s2>0), "MaxConsFree")
            m.addConstrs(( tfs[n,s1,s2] <= tfs[n,s1-1,s2-1]
                    for n in N for s1 in S for s2 in S if s1>0 if s2>0), "MaxConsFree")
            m.addConstrs(( tfs[n,s1,s2] <= 1-q[n,s1]
                    for n in N for s1 in S for s2 in S if s1>0 if s2>0), "MaxConsFree")
            m.addConstrs(( tfs[n,s1,s2] >= tfs[n,s1-1,s2-1]-q[n,s1]
                    for n in N for s1 in S for s2 in S if s1>0 if s2>0), "MaxConsFree")
            if bounds_learned[7,3]>0:
                m.addConstr(
                        (quicksum(tfs[n,s1,s2] for n in N for s1 in S for s2 in range(bounds_learned[7,3],len(S)))==0),"maxconsfree"
                        )
            
            if bounds_learned[7,2]>0:
                m.addConstrs(( sfs[n,0,s2] == 0
                        for n in N for s2 in S), "MinConsFree")
                m.addConstrs(( sfs[n,s1,s2] <= tfs[n,s1-1,s2]
                        for n in N for s1 in S for s2 in S if s1>0), "MinConsFree")
                m.addConstrs(( sfs[n,s1,s2] <= q[n,s1]
                        for n in N for s1 in S for s2 in S if s1>0), "MinConsFree")
                m.addConstrs(( sfs[n,s1,s2] >= tfs[n,s1-1,s2]+q[n,s1]-1
                        for n in N for s1 in S for s2 in S if s1>0), "MinConsFree")
                
                m.addConstrs(( sfs[n,num_shifts,s2] == tfs[n,num_shifts-1,s2]
                        for n in N for s2 in S), "MinConsWork")
                
                m.addConstr(
                        (quicksum(sfs[n,s1,s2]*(bounds_learned[7,2]-1-s2) for n in N for s1 in Ss for s2 in range(bounds_learned[7,2]-1))==0),"minconsfree"
                        )
                
        if bounds_learned[11,5]+bounds_learned[11,4]>0:
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
            if bounds_learned[11,5]>0:
                m.addConstr(
                        (quicksum(tw1[n,d1,s,d2] for s in S for n in N for d1 in D for d2 in range(bounds_learned[11,5],len(D)))==0),"maxconssameshift"
                        )
        
            if bounds_learned[11,4]>0:
                m.addConstrs(( sw1[n,0,s,d2] == 0
                        for s in S for n in N for d2 in D), "MinConsSameShift")
                m.addConstrs(( sw1[n,d1,s,d2] <= tw1[n,d1-1,s,d2]
                        for s in S for n in N for d1 in D for d2 in D if d1>0), "MinConsSameShift")
                m.addConstrs(( sw1[n,d1,s,d2] <= 1-o[n,d1,s]
                        for s in S for n in N for d1 in D for d2 in D if d1>0), "MinConsSameShift")
                m.addConstrs(( sw1[n,d1,s,d2] >= tw1[n,d1-1,s,d2]-o[n,d1,s]
                        for s in S for n in N for d1 in D for d2 in D if d1>0), "MinConsSameShift")
                
                m.addConstrs(( sw1[n,num_days,s,d2] == tw1[n,num_days-1,s,d2]
                        for n in N for s in S for d2 in D), "MinConsWork")
                
                m.addConstr(
                        (quicksum(sw1[n,d1,s,d2]*(bounds_learned[11,4]-1-d2) for s in S for n in N for d1 in Ds for d2 in range(bounds_learned[11,4]-1))==0),"minconssameshift"
                        )
                
        if bounds_learned[11,3] + bounds_learned[11,2]>0:
            m.addConstrs(( tw1f[n,0,s,0] == 1-o[n,0,s] for n in N for s in S), "MaxConsFree")
            m.addConstrs(( tw1f[n,d1+1,s,0] <= o[n,d1,s]
                    for n in N for s in S for d1 in D if d1<len(D)-1), "MaxConsFree")
            m.addConstrs(( tw1f[n,d1+1,s,0] <= 1-o[n,d1+1,s]
                    for n in N for s in S for d1 in D if d1<len(D)-1), "MaxConsFree")
            m.addConstrs(( tw1f[n,d1+1,s,0] >= o[n,d1,s]-o[n,d1+1,s]
                    for n in N for s in S for d1 in D if d1<len(D)-1), "MaxConsFree")
            
            m.addConstrs(( tw1f[n,0,s,d2] == 0
                    for n in N for s in S for d2 in D if d2>0), "MaxConsFree")
            m.addConstrs(( tw1f[n,d1,s,d2] <= tw1f[n,d1-1,s,d2-1]
                    for n in N for s in S for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsFree")
            m.addConstrs(( tw1f[n,d1,s,d2] <= 1-o[n,d1,s]
                    for n in N for s in S for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsFree")
            m.addConstrs(( tw1f[n,d1,s,d2] >= tw1f[n,d1-1,s,d2-1]-o[n,d1,s]
                    for n in N for s in S for d1 in D for d2 in D if d1>0 if d2>0), "MaxConsFree")
            if bounds_learned[11,3]>0:
                m.addConstr(
                        (quicksum(tw1f[n,d1,s,d2] for n in N for s in S for d1 in D for d2 in range(bounds_learned[11,3],len(D)))==0),"maxconsfree"
                        )
            
            if bounds_learned[11,2]>0:
                m.addConstrs(( sw1f[n,0,s,d2] == 0
                        for n in N for s in S for d2 in D), "MinConsw1free")
                m.addConstrs(( sw1f[n,d1,s,d2] <= tw1f[n,d1-1,s,d2]
                        for n in N for s in S for d1 in D for d2 in D if d1>0), "MinConsw1free")
                m.addConstrs(( sw1f[n,d1,s,d2] <= o[n,d1,s]
                        for n in N for s in S for d1 in D for d2 in D if d1>0), "MinConsw1free")
                m.addConstrs(( sw1f[n,d1,s,d2] >= tw1f[n,d1-1,s,d2]+o[n,d1,s]-1
                        for n in N for s in S for d1 in D for d2 in D if d1>0), "MinConsw1free")
                
                m.addConstrs(( sw1f[n,num_days,s,d2] == tw1f[n,num_days-1,s,d2]
                        for n in N for s in S for d2 in D), "MinConsWork")
                
                m.addConstr(
                        (quicksum(sw1f[n,d1,s,d2]*(bounds_learned[11,2]-1-d2) for n in N for s in S for d1 in Ds for d2 in range(bounds_learned[11,2]-1))==0),"minconsw1free"
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
        