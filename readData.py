import sys
import csv
import numpy as np
from collections import OrderedDict
import constraintFormulation as cf
import operator
from numpy.core.defchararray import index
import time
import glob
import os.path
import random
# import xlsxwriter

def readCSV(fileName):
    with open(fileName, 'rU') as csvFile:
        reader = csv.reader(csvFile, delimiter=',')
        data = list(reader)
        data = np.asarray(data)
        return data
        
def cleanData(data):
    variables=[]
    rows, cols = data.shape
#     finding number of variables
    blankRows=0
    while(not data[blankRows,0]):
        blankRows+=1
    blankCols=0
    while(not data[0,blankCols]):
        blankCols+=1
    for i in range(blankRows):
        variables.append(list(OrderedDict.fromkeys(data[i,blankCols:])))
    for i in range(blankCols):
        variables.append(list(OrderedDict.fromkeys(data[blankRows:,i])))
    variables_mat=np.matrix(variables)

    lengths = []
    for i in range(variables_mat.shape[1]):
        lengths.append(len(variables_mat[0,i]))
    dataTensor=np.zeros(lengths)
    
    for i in range(blankRows, rows):
        for j in range(blankCols, cols):
            if data[i,j].astype(int)==1:
                index=()
                for k in range(blankRows):
                    index=index+(variables[k].index(data[k,j]),)
                for k in range(blankCols):
                    index=index+(variables[blankRows+k].index(data[i,k]),)
                dataTensor[index]=1
    return dataTensor,variables

def findConstraints(dataTensor,variables,orderingNotImp,extraConstraints,repeatDim):
    r=set([v for v in range(len(variables)) if v not in repeatDim])
    subsets=cf.split(r,(),repeatDim)
    for l in range(len(subsets)):
        subset=subsets[l]
        newset=subset[0]+subset[1]
        print(subset)
        print(newset)
        # this value will be used to filter max constraints
        maxPossible=1
        for i in range(len(subset[1])):
            maxPossible*=len(variables[subset[1][i]])   
        idTensor=cf.tensorIndicator(dataTensor,newset, variables)
        sumSet = range(len(subset[0]),len(newset))
        print(sumSet)
        if subset[0] in extraConstraints:
            sumTensor=cf.tensorSum(idTensor,sumSet, np.array(variables)[list(newset)],1)
            print('M:',subset[0],' S:',subset[1])
            print(sumTensor)
        else:
            sumTensor_max,sumTensor_min=cf.tensorSum(idTensor,sumSet, np.array(variables)[list(newset)],0)
            if (sumTensor_min != 0 and sumTensor_min!=maxPossible) or (sumTensor_max != maxPossible and sumTensor_max!=0): 
                print('M:',subset[0],' S:',subset[1])
            if sumTensor_min != 0 and sumTensor_min!=maxPossible: 
                print('min: ',sumTensor_min)
            if sumTensor_max != maxPossible and sumTensor_max!=0: 
                print('max: ',sumTensor_max)      
        
        if len(set(subset[1]))==1 and len(set(orderingNotImp) & set(subset[1]))==0:
            minConsZero,maxConsZero,minConsNonZero,maxConsNonZero = cf.tensorConsZero(idTensor,sumSet, np.array(variables)[list(newset)])
            if (minConsZero!=0 and minConsZero!=maxPossible) or (maxConsZero!=maxPossible and maxConsZero>0) or (minConsNonZero!=0 and minConsNonZero!=maxPossible) or (maxConsNonZero!=maxPossible and maxConsNonZero>0):
                print('M:',subset[0],' S:',subset[1])
            if minConsZero!=0 and minConsZero!=maxPossible:
                print('min cons Zero: ',minConsZero)
            if maxConsZero!=maxPossible and maxConsZero>0:
                print('max cons Zero: ',maxConsZero)
            if minConsNonZero!=0 and minConsNonZero!=maxPossible:
                print('min cons NonZero: ',minConsNonZero)
            if maxConsNonZero!=maxPossible and maxConsNonZero>0:
                print('max cons NonZero: ',maxConsNonZero)

def saveConstraints(dataTensor,variables,orderingNotImp,repeatDim,ind,num_nurses,direc,seed):
    r=set([v for v in range(len(variables)) if v not in repeatDim])
    subsets=cf.split(r,(),repeatDim)
#    direc='/home/mohit/CLUT/code/data'
    with open(os.path.join(direc, "result_N"+str(num_nurses)+'_seed'+str(seed)+".csv") ,"a") as my_csv:
        csvWriter = csv.writer(my_csv,delimiter=',')
        if ind==0:
            row=(['']*2)
            for subset in subsets:
                row.extend(subset)
                row.extend(['']*5)
            csvWriter.writerow(row)
        else:
            row=[]
            for l in range(len(subsets)):
                subset=subsets[l]
                newset=subset[0]+subset[1]
                # this value will be used to filter max constraints
                maxPossible=1
                for i in range(len(subset[1])):
                    maxPossible*=len(variables[subset[1][i]])   
                idTensor=cf.tensorIndicator(dataTensor,newset, variables)
                sumSet = range(len(subset[0]),len(newset))
                
                sumTensor_max,sumTensor_min=cf.tensorSum(idTensor,sumSet, np.array(variables)[list(newset)],0)
                if sumTensor_min != 0 and sumTensor_min!=maxPossible: 
                    row.extend([sumTensor_min])
                else:
                    row.extend([''])
                if sumTensor_max != maxPossible and sumTensor_max!=0: 
                    row.extend([sumTensor_max]) 
                else:
                    row.extend([''])
                if len(set(subset[1]))==1 and len(set(orderingNotImp) & set(subset[1]))==0:
#                    print(subset)
                    minConsZero,maxConsZero,minConsNonZero,maxConsNonZero = cf.tensorConsZero(idTensor,sumSet, np.array(variables)[list(newset)])
                    if minConsZero!=0 and minConsZero!=maxPossible:
                        row.extend([minConsZero])
                    else:
                        row.extend([''])
                    if maxConsZero!=maxPossible and maxConsZero>0:
                        row.extend([maxConsZero])
                    else:
                        row.extend([''])
                    if minConsNonZero!=0 and minConsNonZero!=maxPossible:
                        row.extend([minConsNonZero])
                    else:
                        row.extend([''])
                    if maxConsNonZero!=maxPossible and maxConsNonZero>0:
                        row.extend([maxConsNonZero])
                    else:
                        row.extend([''])
                else:
                    row.extend(['']*4)
                row.extend([''])
            csvWriter.writerow(row)

def saveConstraintsForAll(dataTensor,variables,orderingNotImp,repeatDim,ind,num_nurses,direc):
    r=set([v for v in range(len(variables)) if v not in repeatDim])
    subsets=cf.split(r,(),repeatDim)
#    direc='/home/mohit/CLUT/code/data'
    with open(os.path.join(direc, "result_N"+str(num_nurses)+"_forAll.csv") ,"a") as my_csv:
        csvWriter = csv.writer(my_csv,delimiter=',')
        if ind==0:
            row=(['']*2)
            for subset in subsets:
                row.extend(subset)
                row.extend(['']*5)
            csvWriter.writerow(row)
        else:
            row=[]
            for l in range(len(subsets)):
                subset=subsets[l]
                newset=subset[0]+subset[1]
                # this value will be used to filter max constraints
                maxPossible=1
                for i in range(len(subset[1])):
                    maxPossible*=len(variables[subset[1][i]])   
                idTensor=cf.tensorIndicator(dataTensor,newset, variables)
                sumSet = range(len(subset[0]),len(newset))
                
                sumTensor_max,sumTensor_min=cf.tensorSum(idTensor,sumSet, np.array(variables)[list(newset)],0)
                row.extend([sumTensor_min])
                row.extend([sumTensor_max]) 
                if len(set(subset[1]))==1 and len(set(orderingNotImp) & set(subset[1]))==0:
#                    print(subset)
                    minConsZero,maxConsZero,minConsNonZero,maxConsNonZero = cf.tensorConsZero(idTensor,sumSet, np.array(variables)[list(newset)])
                    row.extend([minConsZero])
                    row.extend([maxConsZero])
                    row.extend([minConsNonZero])
                    row.extend([maxConsNonZero])
                    
                else:
                    row.extend(['']*4)
                row.extend([''])
            csvWriter.writerow(row)

def createList(data):
    output=[]
    values=[]
    vals=[]
    n=data.shape[0]
    for i in range(n):
        if i==0:
            values.append(data[i,0])
            vals.append(data[i,1])
        elif i!=0 and data[i,1]==data[i-1,1] and i!=n-1:
            values.append(data[i,0])
        elif i!=0 and data[i,1]==data[i-1,1] and i==n-1:
            values.append(data[i,0])
            output.append(values)
        else:
            output.append(values)
            values=[]
            values.append(data[i,0])
            vals.append(data[i,1])
    return output, vals
 
def vectorizeExtraInfo(extraInfo,variables):
    dim=int(extraInfo[0,0])
    extraInfo=extraInfo[1:,:]
    sortedData=np.array(sorted(extraInfo, key=operator.itemgetter(1)))
    lst,vals=createList(sortedData)
    variable=variables[dim]
    filterTensors=[]
    for k in range(len(lst)):
        filterTensor = np.zeros((len(lst[k]),len(variable)))
        i=0
        j=0
        while i<len(lst[k]) or j<len(variable):
            if variable[j] in lst[k]:
                filterTensor[i,j]=1
                i+=1
            j+=1
        filterTensors.append(filterTensor)
    return filterTensors,lst,vals

def learnConstraints(dataDir,numSol,num_nurses,numTables,direc,numFiles,seed):
    random.seed(seed)
#    print(numFiles,numSol)
    selectFiles=random.sample(range(0,numFiles),numSol)
#    print(selectFiles)
    start=time.clock()
    ind=0
    for fl in glob.glob(dataDir):  
        if ind==0 or ind in selectFiles:
            data = readCSV(fl)
            dataTensor,variables=cleanData(data)
            lenVar=[]
            for i in range(len(variables)):
                lenVar.append(len(variables[i]))
            extraConstraints = []
            orderingNotImp=[2]
            repeatDim=()
            if ind==0:
                saveConstraints(dataTensor,variables,orderingNotImp,repeatDim,0,num_nurses,direc,seed)
            if ind in selectFiles:
                saveConstraints(dataTensor,variables,orderingNotImp,repeatDim,1,num_nurses,direc,seed)
        ind+=1

def learnConstraintsForAll(dataDir,num_nurses,numTables,direc):
    start=time.clock()
    ind=0
    for fl in glob.glob(dataDir):
        data = readCSV(fl)
        dataTensor,variables=cleanData(data)
        lenVar=[]
        for i in range(len(variables)):
            lenVar.append(len(variables[i]))
        extraConstraints = []
        orderingNotImp=[2]
        repeatDim=()
        if ind==0:
            saveConstraintsForAll(dataTensor,variables,orderingNotImp,repeatDim,0,num_nurses,direc)
        saveConstraintsForAll(dataTensor,variables,orderingNotImp,repeatDim,1,num_nurses,direc)
        ind+=1
    
#    for i in range(numTables-1):
#        fileName = sys.argv[i+3]
#        extraInfo = readCSV(fileName)
#        filterTensors,lst,vals=vectorizeExtraInfo(extraInfo,variables)
#        for i in range(len(vals)):
#            dim=int(extraInfo[0,0])
#            updatedVariables=variables[:]
#            updatedVariables[dim]=[x for x in variables[dim] if x in lst[i]]
#            mat=np.tensordot(filterTensors[i],dataTensor, [1,dim])
#            for j in range(dim):
#                mat=np.swapaxes(mat,j,j+1)
#            print(mat.shape)
#            print("\nConstraints for filter: ", vals[i])
#            findConstraints(mat,updatedVariables,orderingNotImp,extraConstraints,repeatDim)
#            saveConstraints(dataTensor,variables,orderingNotImp,repeatDim,1)
    print("\nTime Taken: ",time.clock()-start,' secs')
















