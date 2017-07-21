
"""
The Second-Order Moller-Plesset Perturbation Theory (MP2) Energy code
"""

import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as LA
import timing

#Please check the number of occupied MO in the molecule
OCC=5;
print("The number of occupied MO orbitals in the molecule is: ", OCC)

#Please check the number of atomic orbitals
AO=7;
print("The number of Atomic orbitals used is: ", AO)

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
for row in range(over_mat.shape[0]):
    for col in range(over_mat.shape[1]):
        print('{:11.8f}'.format(over_mat[row,col]),end=' ')
    print()


#The kinetic energy matrix
kinetic=np.loadtxt('t.dat')
kinetic_mat=np.zeros((int(kinetic[kinetic.shape[0]-1,0]),int(kinetic[kinetic.shape[0]-1,0])))
for row in range(kinetic.shape[0]):
    i1=int(kinetic[row,0])
    i2=int(kinetic[row,1])
    val=kinetic[row,2]
    kinetic_mat[i1-1,i2-1]=kinetic_mat[i2-1,i1-1]=val
print ("\n Kinetic Energy Matrix \n \n")
for row in range(kinetic_mat.shape[0]):
    for col in range(kinetic_mat.shape[1]):
        print('{:11.8f}'.format(kinetic_mat[row,col]),end=' ')
    print()


#The nuclear attraction matrix
nuk_att=np.loadtxt('v.dat')
nuk_mat=np.zeros((int(nuk_att[nuk_att.shape[0]-1,0]),int(nuk_att[nuk_att.shape[0]-1,0])))
for row in range(nuk_att.shape[0]):
    i1=int(nuk_att[row,0])
    i2=int(nuk_att[row,1])
    val=nuk_att[row,2]
    nuk_mat[i1-1,i2-1]=nuk_mat[i2-1,i1-1]=val
print ("\n  Nuclear Att. Integrals are  \n \n")
for row in range(nuk_mat.shape[0]):
    for col in range(nuk_mat.shape[1]):
        print('{:11.7f}'.format(nuk_mat[row,col]),end=' ')
    print()
    



#np.set_printoptions(precision=4)
#The the core Hamiltonian matrix
H_core=kinetic_mat+nuk_mat
print ("\n This is the Core Hamiltonian Matrix \n \n")
for row in range(H_core.shape[0]):
    for col in range(H_core.shape[1]):
        print('{:11.7f}'.format(H_core[row,col]),end=' ')
    print()



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
    




ele_int=np.loadtxt('eri.dat')
#np.set_printoptions(threshold=10000)

timing.log('Tracking the element allocation ')
d_I =np.zeros(1000)
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
for row in range(eigvec.shape[0]):
    for col in range(eigvec.shape[1]):
        print('{:11.7f}'.format(eigvec[row,col]),end=' ')
    print()
eig_t=eigvec.transpose()
eigval=eigval**(-1/2)
eigval_m=np.diag(eigval)
np.set_printoptions(precision=5)
print("\n\n\nThe eigenvalues in matrix form:\n\n ",eigval_m)
s_1by2=eigvec.dot(eigval_m).dot(eig_t)
print("\n\n\nThe symmetric orthogonalization matrix is:\n\n ")
for row in range(s_1by2.shape[0]):
    for col in range(s_1by2.shape[1]):
        print('{:11.7f}'.format(s_1by2[row,col]),end=' ')
    print()



# Building the initial (guess denstiy matrix)



i_F=s_1by2.transpose().dot(H_core).dot(s_1by2)
print("The transformed Fock matrix is :\n \n ")
for row in range(i_F.shape[0]):
    for col in range(i_F.shape[1]):
        print('{:11.7f}'.format(i_F[row,col]),end=' ')
    print()
eigval1, eigvec1=LA.eigh(i_F)
AO_eigvec=s_1by2.dot(eigvec1)
print("\nThe initial MO coeffcient is  :\n \n ")
for row in range(AO_eigvec.shape[0]):
    for col in range(AO_eigvec.shape[1]):
        print('{:11.7f}'.format(AO_eigvec[row,col]),end=' ')
    print()


#The initial density matrix
AO_eigvec=np.delete(AO_eigvec,range(OCC,AO_eigvec.shape[1]),1)
print("\n\nThe new MO matrix with only the occupied orbitals is\n\n ")
for row in range(AO_eigvec.shape[0]):
    for col in range(AO_eigvec.shape[1]):
        print('{:11.7f}'.format(AO_eigvec[row,col]),end=' ')
    print()

print("\n\nThe new MO matrix (transpose) is\n\n ")
for row in range(AO_eigvec.transpose().shape[0]):
    for col in range(AO_eigvec.transpose().shape[1]):
        print('{:11.7f}'.format(AO_eigvec.transpose()[row,col]),end=' ')
    print()


i_Den=AO_eigvec.dot(AO_eigvec.transpose())
print("\nThe initial density matrix is  :\n \n ")
for row in range(i_Den.shape[0]):
    for col in range(i_Den.shape[1]):
        print('{:11.7f}'.format(i_Den[row,col]),end=' ')
    print()
    
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
counter=0
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
            AO_eigvec2=s_1by2.dot(eigvec2)
            eig_AO_store=AO_eigvec2
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




plt.plot(xx,ene222, '--',linewidth=1)
plt.xlim(0,25)
plt.ylim(-76,-85)




#This implements the extras in Project 3 - The molecular orbital basis to check whether transformation is 
#actually required or not.


print("\nThe final coefficient matrix obtained form HF is \n")

for row in range(eig_AO_store.shape[0]):
    for col in range(eig_AO_store.shape[1]):
        print('{:11.7f}'.format(eig_AO_store[row,col]),end=' ')
    print()
timing.log('Tracking the matrix multiplication time')
F_MO2=eig_AO_store.transpose().dot(fock).dot(eig_AO_store)
timing.endlog()
print("\nThe Fock matrix in the MO basis is \n")
for row in range(F_MO2.shape[0]):
    for col in range(F_MO2.shape[1]):
        print('{:11.7f}'.format(F_MO2[row,col]),end=' ')
    print()


d_M=np.zeros((d_I.size))

ijkl=0
for i in range(AO):
    for j in range(i):
        for k in range(i):
            lim=j if (i==k) else k 
            for l in range(lim):
                for p in range(AO):
                    for q in range(AO):
                        for r in range(AO):
                            for s in range(AO):
                                #print(p,q,r,s)
                                pqrs=index2(p+1,q+1,r+1,s+1)
                                #print(ijkl," ",pqrs," ", d_I[pqrs])
                                d_M[ijkl]+=eig_AO_store[p,i]*eig_AO_store[q,j]*eig_AO_store[r,k]*eig_AO_store[s,l]*d_I[pqrs]
            ijkl=ijkl+1
            
np.set_printoptions(threshold=10000)                       
for a in range(d_M.size):
    print(a,"  ", d_M[a])






























