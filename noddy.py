#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 20:34:39 2017

@author: mintuser
"""


import matplotlib.pyplot as plt
import numpy as np



def noddy(eig_AO_store,AO,d_I,d_M):
#The noddy algorithm for the transfrom of two electron integrals to the MO basis
    #The first index function
    def index1(i, j):
            ij=(i*(i+1)/2)+j
            return int(ij)
        
    #The second index function
    def index2(i, j, k, l):
        if (i>j):
            ij=index1(i,j)
        else:
            ij=index1(j,i)      
        if (k>l):
            kl=index1(k,l)
        else:
            kl=index1(l,k)     
        if (ij> kl):
            ijkl=int(index1(ij,kl))
        else:
            ijkl=int(index1(kl,ij))
        return int(ijkl)
    #The second index function used for the mp2 program
    def index3(i,j):
        if (i>j):
            return index1(i,j)
        else:
            return index1(j,i)


    ijkl=0
    for i in range(AO):
        for j in range(i+1):
            for k in range(i+1):
                lim=j if (i==k) else k 
                for l in range(lim+1):               
                    for p in range(AO):
                        for q in range(AO):
                            for r in range(AO):
                                for s in range(AO):
                                    #print(p,q,r,s)
                                    pqrs=index2(p+1,q+1,r+1,s+1)
                                    #print(ijkl," ",pqrs," ", d_I[pqrs])
                                    d_M[ijkl]+=eig_AO_store[p,i]*eig_AO_store[q,j]*eig_AO_store[r,k]*eig_AO_store[s,l]*d_I[pqrs]
                    ijkl=ijkl+1
    return d_M