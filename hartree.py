#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 20:27:14 2017

@author: mintuser
"""
import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as LA
import timing
import print_perfect




def hartree (OCC,AO,eig_AO_store,eps,d_I):
    
    nu_repul=np.loadtxt('enuc.dat')
    print("\n\nNuclear repulsion energy in Hartree: ",nu_repul,"\n")

    #The overlap matrix
    overlap=np.loadtxt('s.dat')
    over_mat=np.zeros((int(overlap[overlap.shape[0]-1,0]), int(overlap[overlap.shape[0]-1,0])))
    for row in range(overlap.shape[0]):
        i1=int(overlap[row,0])
        i2=int(overlap[row,1])
        val=overlap[row,2]
        over_mat[i1-1,i2-1]=over_mat[i2-1,i1-1]=val
    print ("Overlap Integral Matrix \n \n")
    print_perfect.print_perfect(over_mat)
    
    
    #The kinetic energy matrix
    kinetic=np.loadtxt('t.dat')
    kinetic_mat=np.zeros((int(kinetic[kinetic.shape[0]-1,0]),int(kinetic[kinetic.shape[0]-1,0])))
    for row in range(kinetic.shape[0]):
        i1=int(kinetic[row,0])
        i2=int(kinetic[row,1])
        val=kinetic[row,2]
        kinetic_mat[i1-1,i2-1]=kinetic_mat[i2-1,i1-1]=val
    print ("\n Kinetic Energy Matrix \n \n")
    print_perfect.print_perfect(kinetic_mat)
    
    
    #The nuclear attraction matrix
    nuk_att=np.loadtxt('v.dat')
    nuk_mat=np.zeros((int(nuk_att[nuk_att.shape[0]-1,0]),int(nuk_att[nuk_att.shape[0]-1,0])))
    for row in range(nuk_att.shape[0]):
        i1=int(nuk_att[row,0])
        i2=int(nuk_att[row,1])
        val=nuk_att[row,2]
        nuk_mat[i1-1,i2-1]=nuk_mat[i2-1,i1-1]=val
    print ("\n  Nuclear Att. Integrals are  \n \n")
    print_perfect.print_perfect(kinetic_mat)
        
    
    
    
    #np.set_printoptions(precision=4)
    #The the core Hamiltonian matrix
    H_core=kinetic_mat+nuk_mat
    print ("\n This is the Core Hamiltonian Matrix \n \n")
    print_perfect.print_perfect(H_core)
    
    
    
    """
    Reading in the the two electron integrals
    """
    
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
    
    
    
    
    ele_int=np.loadtxt('eri.dat')
    #np.set_printoptions(threshold=10000)
    
    timing.log('Tracking the element allocation ')
    for row in range(ele_int.shape[0]):
        i=int(ele_int[row,0])
        j=int(ele_int[row,1])
        k=int(ele_int[row,2])
        l=int(ele_int[row,3])
        val=ele_int[row,4]
        ijkl=index2(i,j,k,l)
        d_I[ijkl]=val  
    timing.endlog()   
    #Diagonilizing the matrix
    eigval, eigvec=LA.eigh(over_mat)
    print("The eigenvalues are:", "\n \n",eigval)
    print("\n\n\nThe eigenvectors are: \n \n ")
    print_perfect.print_perfect(eigvec)
    eig_t=eigvec.transpose()
    eigval=eigval**(-1/2)
    eigval_m=np.diag(eigval)
    np.set_printoptions(precision=5)
    print("\n\n\nThe eigenvalues in matrix form:\n\n ",eigval_m)
    s_1by2=eigvec.dot(eigval_m).dot(eig_t)
    print("\n\n\nThe symmetric orthogonalization matrix is:\n\n ")
    print_perfect.print_perfect(s_1by2)
    
    
    
    # Building the initial (guess denstiy matrix)
    
    
    
    i_F=s_1by2.transpose().dot(H_core).dot(s_1by2)
    print("The transformed Fock matrix is :\n \n ")
    print_perfect.print_perfect(i_F)
    eigval1, eigvec1=LA.eigh(i_F)
    AO_eigvec=s_1by2.dot(eigvec1)
    print("\nThe initial MO coeffcient is  :\n \n ")
    print_perfect.print_perfect(AO_eigvec)
    
    
    #The initial density matrix
    AO_eigvec=np.delete(AO_eigvec,range(OCC,AO_eigvec.shape[1]),1)
    print("\n\nThe new MO matrix with only the occupied orbitals is\n\n ")
    print_perfect.print_perfect(AO_eigvec)
    
    
    print("\n\nThe new MO matrix (transpose) is\n\n ")
    print_perfect.print_perfect(AO_eigvec.transpose())
    
    
    i_Den=AO_eigvec.dot(AO_eigvec.transpose())
    print("\nThe initial density matrix is  :\n \n ")
    print_perfect.print_perfect(i_Den)
    
        
    #Computing the Initial SCF Energy
    ene1=0.0
    for row in range(H_core.shape[0]):
        for col in range (H_core.shape[1]):
            ene1+=i_Den[row,col]*(H_core[row,col]+i_F[row,col])
    print("\n The initial energy is (in Hartree) ", ene1)
    print("\n The total energy is (in Hartree) ", ene1+nu_repul)
    
    
    # Computing the new density matrix
    # The loop begins here 
    
    print()
    fock=np.zeros((H_core.shape[0],H_core.shape[1]))
    ene=ene1
    xx=np.zeros((1000)) 
    ene222=np.zeros((1000))
    
    # Learning about the time function as well
    
    timing.log('Tracking the SCF loop')
    print("Iteration"," ", "Electronic energy", "\t Total Energy", "\t\t Delta E", "\t Delta RMS")
    
    
    for x in range(1000):
                for i in range(H_core.shape[0]):
                    for j in range(H_core.shape[0]):
                        fock[i,j]=H_core[i,j]
                        for k in range(H_core.shape[0]):
                            for l in range(H_core.shape[0]):
                                fock[i,j]+=i_Den[k,l]*(2.0*d_I[int(index2(i+1,j+1,k+1,l+1))]-d_I[int(index2(i+1,k+1,j+1,l+1))])
    
                new_fock=s_1by2.transpose().dot(fock).dot(s_1by2)
                eigval2, eigvec2=LA.eigh(new_fock)
                eps=eigval2
                AO_eigvec2=s_1by2.dot(eigvec2)
                eig_AO_store=np.copy(AO_eigvec2)
                AO_eigvec2=np.delete(AO_eigvec2,range(OCC,AO_eigvec2.shape[1]),1)
                den=AO_eigvec2.dot(AO_eigvec2.transpose())
                ene_new=0.0
                for row in range(H_core.shape[0]):
                    for col in range (H_core.shape[1]):
                        #I had made a mistake here I was using the new fock matrix instead of the old one
                        ene_new+=den[row,col]*(H_core[row,col]+fock[row,col])
                rms=0.0
                rms_last=0.0
                for row in range(den.shape[0]):
                    for col in range(den.shape[1]):
                        rms+=(den[row,col]-i_Den[row,col])**2
                rms=rms**0.5
                xx[x]=x
                ene222[x]=ene_new
                print(x,"\t",'{:15.13f}'.format(ene_new),"\t", '{:15.13f}'.format(ene_new+nu_repul),"   ",'{:13.10f}'.format(ene_new-ene)," ",'{:13.10f}'.format(rms_last-rms))
                if ( (abs(ene_new-ene) < 1e-14 and abs(rms_last-rms)<1e-14) ):
                    break
                i_Den=den
                rms_last=rms
                ene=ene_new
    timing.endlog()
    
    #####################################################
    plt.plot(xx,ene222, '--',linewidth=1)
    plt.xlim(0,25)
    plt.ylim(-76,-85)
    #####################################################
    
    
    
    #This implements the extras in Project 3 - The molecular orbital basis to check whether transformation is 
    #actually required or not.
    
    print("\nThe final coefficient matrix obtained form HF is \n")
    print_perfect.print_perfect(eig_AO_store)
    
    timing.log('Tracking the matrix multiplication time')
    F_MO2=eig_AO_store.transpose().dot(fock).dot(eig_AO_store)
    timing.endlog()
    print("\nThe Fock matrix in the MO basis is \n")
    print_perfect.print_perfect(F_MO2)
    

    return eig_AO_store,eps,d_I