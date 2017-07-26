
"""
The Second-Order Moller-Plesset Perturbation Theory (MP2) Energy code
"""
#%matplotlib notebook
import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as LA
import timing
import print_perfect
import hartree
import noddy

#Please check the number of occupied MO in the molecule
OCC=5;
print("The number of occupied MO orbitals in the molecule is: ", OCC)

#Please check the number of atomic orbitals
AO=7;
print("The number of Atomic orbitals used is: ", AO)

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


eig_AO_store=np.zeros((AO,AO))
eps=np.zeros((AO))
d_I =np.zeros(1000)

eig_AO_store,eps,d_I=hartree.hartree(OCC,AO,eig_AO_store,eps,d_I)


d_M=np.zeros((int((AO*(AO+1)/2)*((AO*(AO+1)/2)+1)/2)))

d_M=noddy.noddy(eig_AO_store,AO,d_I,d_M)

#np.set_printoptions(threshold=10000)                       
#for a in range(d_M.size):
#    print(a,"  ", '{:11.8f}'.format(d_M[a]))

#print(eig_AO_store,eps,d_I)

emp2=0.0
for i in range(OCC):
    for a in range(OCC,AO):
        for j in range(OCC):
            for b in range(OCC,AO):
                iajb=index2(i,a, j, b)
                ibja=index2(i,b, j, a)
                emp2+=d_M[iajb]*(2*d_M[iajb]-d_M[ibja])/(eps[i]+eps[j]-eps[a]-eps[b])

print("This is the mp2 energy",emp2)




"""

X=np.zeros((AO,AO))
Y=np.zeros((AO,AO))
TMP=np.zeros((AO*(AO+1)/2,AO*(AO+1)/2))
d_M_f=np.zeros(((AO*(AO+1)/2)*((AO*(AO+1)/2)+1)/2))

ij=0
for i in range(AO):
    for j in range(i+1):
        kl=0
        for k in range(AO):
            for l in range(k+1):
                ijkl=index3(ij,kl)
                X[k,l]=X[l,k]=d_I[ijkl]
                kl=kl+1
        Y=np.zeros((AO,AO))
        Y=eig_AO_store.transpose().dot(X)
        X=np.zeros((AO,AO))
        X=Y.dot(eig_AO_store)
        kl=0
        for k in range(AO):
            for l in range(k+1):
                TMP[kl,ij]=X[k,l]
                kl=kl+1
    ij=ij+1



kl=0
for k in range(AO):
    for l in range(k+1):
        Y=np.zeros((AO,AO))
        X=np.zeros((AO,AO))
        ij=0
        for i in range(AO):
            for j in range (i+1):
                X[i,j]=X[j,i]=TMP[kl][ij]
                ij=ij+1
        Y=np.zeros((AO,AO))
        Y=eig_AO_store.transpose().dot(X)
        X=np.zeros((AO,AO))
        X=Y.dot(eig_AO_store)
        ij=0
        for i in range (AO):
            for j in range(i+1):
                klij=index3(kl,ij)
                print("this is x[i,j]",X[i,j]," ",klij)
                d_M_f[klij]=X[i,j]
                ij=ij+1
    kl=kl+1
        
np.set_printoptions(threshold=10000)                       
for a in range(d_M_f.size):
    print(a,"  ", '{:11.8f}'.format(d_M_f[a]))       

"""

















