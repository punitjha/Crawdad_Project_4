//The mass header file

#include "mass.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <vector>
#include <stdio.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <armadillo>
#include <complex>

#include <array>   // to be able to use this library compile using this g++ -std=c++0x FileNames
void geometry(arma::mat &geom, int &n_atoms, arma::vec &z_val);
int index1(int i, int j);
int index_call (int i, int j, int k, int l);
using namespace std;

int main()
{


/******************************************************************/	
//Reading the geometries and nuclear repulsion energy from the file 
/******************************************************************/	


	arma::mat geom;
	arma::vec z_val;
	int n_atoms=0;
	geometry(geom, n_atoms, z_val);
	geom.print("\nThis is the geometrical coordinates of the molecule\n");
	double enuc=0.0;
	ifstream enuc1 ("enuc.dat", ios::in);
	if (enuc1.is_open())
	{
		 cout<<"\nReading from the nuclear repulsion energy file\n";
		 enuc1>>enuc;
		 cout<<"\nThe nuclear repulsion energy is:\n"<<enuc<<endl;
	}

/**************************************************************/	
//Reading the one electron integrals 
/**************************************************************/	

//Reading the overlap matrix elements
	arma::mat S_file;
	S_file.load("s.dat");
	int last_s=S_file(S_file.n_rows-1,0);
	arma::mat S=arma::zeros(last_s,last_s);
	for(int i=0; i<S_file.n_rows; i++)
	{
		int i1=S_file(i,0);
		int i2=S_file(i,1);
		double val=S_file(i,2);
		S(i1-1,i2-1)=val;
		S(i2-1,i1-1)=val;
	}	
	S.print("\nThe overlap intergral elements are\n ");

	arma::mat T_file;
	T_file.load("t.dat");
	int last_t=T_file(T_file.n_rows-1,0);
	arma::mat T=arma::zeros(last_t,last_t);
	for(int i=0; i<T_file.n_rows; i++)
	{
		int i1=T_file(i,0);
		int i2=T_file(i,1);
		double val=T_file(i,2);
		T(i1-1,i2-1)=val;
		T(i2-1,i1-1)=val;
	}	
	T.print("\nThe kinetic energy intergals are\n");

	arma::mat V_file;
	V_file.load("v.dat");
	int last_v=V_file(V_file.n_rows-1,0);
	arma::mat V=arma::zeros(last_v,last_v);
	for(int i=0; i<V_file.n_rows; i++)
	{
		int i1=V_file(i,0);
		int i2=V_file(i,1);
		double val=V_file(i,2);
		V(i1-1,i2-1)=val;
		V(i2-1,i1-1)=val;
	}	
	V.print("\nThe nuclear attaraction intergral elements are\n ");

	arma::mat H_c=T+V;
	real(H_c).print("\nThis is the core Hamiltonian\n");

/**************************************************************/	
//Reading the two electron integrals 
/**************************************************************/	

//here we have to store the 4D matrix involving the 2 electron integrals
//into a one dimensional matrix
	arma::mat d_file;
	d_file.load("eri.dat");
	arma::vec d_I =arma::zeros(1000);	
	for(int i1=0; i1<d_file.n_rows; i1++)
	{
		int i=d_file(i1,0);
		int j=d_file(i1,1);
		int k=d_file(i1,2); 
		int l=d_file(i1,3);
		double val=d_file(i1,4);
		int ij, kl,ijkl;
		if(i>j)
		{
			 ij=index1(i,j);
		}
		else    
		{
			 ij=index1(j,i);
		}
		if (k>l)
		{
			 kl=index1(k,l);
		}
		else
		{
			kl=index1(l,k);
		}
		if(ij > kl)
		{
			 ijkl=ij*(ij+1)/2+(kl);
		}
		else //if (ij < kl)

		{
			 ijkl=kl*(kl+1)/2+(ij);
		}
		d_I(ijkl)=val;
	}	


/**************************************************************/
//Building the orthogonalization matrix
/**************************************************************/	

	arma::vec eigval;
	arma::mat eigvec;
	arma::eig_sym (eigval,eigvec,S);
	real(eigvec).print("\nPrinting the eigenvectors of the overlap matrix\n");
	real(eigval).print("\nPrinting the eigenvalues of the overlap matrix\n");
	arma::mat tran_eigvec=eigvec.t();
	arma::mat prod=tran_eigvec*eigvec;
	real(prod).print("\nThe product of the eigenvec and its transpose\n");
	arma::mat sqrt_eigen=arma::pow(eigval,-0.5);
	real(sqrt_eigen).print("\nThe square root of the eigenvalue matrix\n");
	arma::mat diag_ee=arma::diagmat(sqrt_eigen);
	real(diag_ee).print("\nThe diagonal matrix of the eigenvalues\n");
	arma::mat S_sqrt=eigvec*diag_ee*tran_eigvec;
	arma::mat r_S_sqrt=arma::real(S_sqrt);
	real(S_sqrt).print("\nThe orthogonalized symmetric overlap matrix\n");


/**************************************************************/
//Building the initial guess density matrix 
/**************************************************************/	

	arma::mat I_Gs=S_sqrt.t()*H_c*S_sqrt;
	I_Gs.print("\nThis is the initial Fock matrix(orthogonal basis)\n");
//Note that this contains the inital orbital range
	arma::vec eig11; 
	arma::mat eigenvec;
	arma::eig_sym(eig11,eigenvec,I_Gs);
	real(eig11).print("\nThe eigenvalues (initially)\n");
	real(eigenvec).print("\nThe eigenvectors (initially)\n");
	//eigvec1.print("\nThe initial MO coefficient Matrix is \n");
//Transforming the eigenvectors into the original basis
	arma::mat eig_AO = S_sqrt*eigenvec;
	real(eig_AO).print("\nThe initial MO cofficients\n");
//The density matrix
//Since only the first MO of water are occupied so the lowest five 
//MO are taken into consideration and they are effectively squared.	
	eig_AO.shed_cols(5,6);
	arma::mat den_I = eig_AO*eig_AO.t();	
	real(den_I).print("\nThe initial guess density matrix\n");

/**************************************************************/
//Computing the initial SCF energy matrix and the Fock Matrix
/**************************************************************/	

//Computing the SCF electronic energy using the density matrix

	double energy=0.0;
	for(int i=0; i< H_c.n_rows; i++)
	{
		for(int j=0; j<H_c.n_cols; j++)
		{
			energy+=den_I(i,j)*(H_c(i,j)+I_Gs(i,j));
		}
	}
	
	printf("\n %s\t %s\t %s\t %s\t\t %s\n %d %19.5f \n","Iteration", "Electronic energy","E_total","Delta(E)","RMS(D)",0,energy);

/**************************************************************/
//Computing the new Fock Matrix and starting the SCF loop
/**************************************************************/	
	arma::mat density_initial=den_I;
	double energy_last=energy;
//setting the precision here
	cout.precision(15);
	arma::mat F_mv=H_c;
        arma::mat coeff_car=den_I;	
//the SCF loop starts here
	for(int ii=1; ii<100000; ii++)
	{	
		F_mv=H_c;
		for(int i=0;i<H_c.n_rows;i++)
		{
			for(int j=0; j<H_c.n_cols; j++)
			{
				for(int k=0;k<H_c.n_cols;k++)
				{
					for(int l=0;l<H_c.n_cols;l++)
					{
						F_mv(i,j)+=den_I(k,l)*(2*d_I(index_call(i+1,j+1,k+1,l+1))-d_I(index_call(i+1,k+1,j+1,l+1)));
					}
				}
			}
		}
		F_mv.print("\nThe Fock matrix with 2 electron integrals \n");
		arma::mat F_store=S_sqrt.t()*F_mv*S_sqrt;
		arma::vec eig_store;
		arma::mat eigenvec_store;
		arma::eig_sym(eig_store,eigenvec_store,F_store);
		arma:: mat eig_AO_store=S_sqrt*eigenvec_store;
		coeff_car=eig_AO_store;
		//eig_AO_store.print("\n");
		eig_AO_store.shed_cols(5,6);
		arma:: mat density_store=eig_AO_store*eig_AO_store.t();
		double energy1=0.0;
		for(int i=0; i< H_c.n_rows; i++)
		{
			for(int j=0; j<H_c.n_cols; j++)
			{
				energy1+=density_store(i,j)*(H_c(i,j)+F_mv(i,j));
			}
		}
		double d_ene=energy_last-energy1;
		double rms=0.0;
		double rms_last=0.0;
		for(int i=0; i<density_store.n_rows;i++)
		{
			for (int j=0;j<density_store.n_cols;j++)
			{
				rms+=pow(density_store(i,j)-density_initial(i,j),2);
			}
		}
		rms=pow(rms,0.5);
		den_I=density_store;
		printf("\n %d\t  %10.12f \t %10.12f \t %10.12f \t %10.12f \n ",ii,energy1,energy1+enuc,d_ene,rms);	
		//F_mv.print("\n");
		//den_I.print("\nThe density matrix after quitting the loop\n");
		if (abs(d_ene)<1e-12 && abs(rms_last-rms)< 1e-12 )
			break;
//Changing the stuff here
		density_initial=density_store;
		energy_last=energy1;
		rms_last=rms;
		
	}		
/**************************************************************/
//The MO-basis Fock-Fatrix
/**************************************************************/
	//F_mv.print("\n");
	//coeff_car.print("\n");	
	arma::mat F_MO2=coeff_car.t()*F_mv*coeff_car;
	cout<<"\n\nChecking for the Diagonal form of Fork Matrix in the basis of the Molecular Orbitals\n\n";

	for (int i=0; i<F_MO2.n_rows;i++)
	{
		for(int j=0; j<F_MO2.n_cols;j++)
		{
			printf("%15.7f",F_MO2(i,j));
		}
		cout<<endl;
	}
			
/********************************************************************/
//One Electron Properties- electronic contribution to dipole operator
/*********************************************************************/
//Reading the files
	arma::mat mux_file;
	mux_file.load("mux.dat");
	arma::mat muy_file;
	muy_file.load("muy.dat");
	arma::mat muz_file;
	muz_file.load("muz.dat");
//Constructing the matrix (one dimensional based on index)
	int ele=muz_file(muz_file.n_rows-1,0);
	arma::mat mux=arma::zeros(ele,ele);
	arma::mat muy=arma::zeros(ele,ele);
	arma::mat muz=arma::zeros(ele,ele);
	for(int i=0; i< mux_file.n_rows;i++)
	{
		int x1=mux_file(i,0);
		int x2=mux_file(i,1);
		double val_x=mux_file(i,2);
		mux(x1-1,x2-1)=val_x;

		int y1=muy_file(i,0);
		int y2=muy_file(i,1);
		double val_y=muy_file(i,2);
		muy(y1-1,y2-1)=val_y;

		int z1=muz_file(i,0);
		int z2=muz_file(i,1);
		double val_z=muz_file(i,2);
		muz(z1-1,z2-1)=val_z;

	}
	mux=symmatl(mux);
	muy=symmatl(muy);
	muz=symmatl(muz);
	mux.print("\nThis is the mu_x integrals\n");
	muy.print("\nThis is the mu_y integrals\n");
	muz.print("\nThis is the mu_z integrals\n");

	double mu_x=0.0;
	double mu_y=0.0;
	double mu_z=0.0;
	double N_x=0.0;
	double N_y=0.0;	
	double N_z=0.0;

	for(int i=0; i<mux.n_rows; i++)
	{
		for(int j=0; j<mux.n_cols;j++)
		{
			mu_x+=2*den_I(i,j)*mux(i,j);
			mu_y+=2*den_I(i,j)*muy(i,j);
			mu_z+=2*den_I(i,j)*muz(i,j);
		}
	}
	for(int i=0; i<n_atoms; i++)
	{
		N_x+=z_val(i)*geom(i,0);
		N_y+=z_val(i)*geom(i,1);
		N_z+=z_val(i)*geom(i,2);
	}
	cout.precision(15);
	printf(" \n %s \t %8.5f \n %s \t %8.5f \n %s \t %8.5f \n","The dipole moment operators are: \n\n mu_x",mu_x+N_x,"mu_y",mu_y+N_y,"mu_z",mu_z+N_z);
	printf("\n %s \t %8.5f\n", "The total dipole moment is :", mu_x+N_x+mu_y+N_y+mu_z+N_z);

//*****************************************************************************
//Starting the mp2 calculations here to  verify the python codes
//*****************************************************************************

	arma::vec d_MO=arma::zeros(d_I.n_elem);
	int ijkl=0;
	int pqrs=0;
	for (int i=0; i<7; i++)
	{
		for (int j=0; j<=i; j++)
		{
			for (int k=0; k<=i; k++)
			{
				for ( int l=0; l<= (i==k ? j : k); l++, ijkl++)
				{
					for ( int p=0; p<7; p++)
					{
						for ( int q=0; q<7; q++)
						{
							for(int r=0; r<7; r++)
							{
								for(int s=0; s<7; s++)
								{
									pqrs=index_call(p,q,r,s);
									d_MO(ijkl)+= coeff_car(p,i)*coeff_car(q,j)*coeff_car(r,k)*coeff_car(s,l)*d_I(pqrs);
								}
							}
						}
					}
				}
			}
		}
	}
	for (int i=0; i<d_MO.n_elem; i++)
	{
		cout<<i<<"  "<<d_MO(i)<<endl;
	}


	
}



void geometry(arma::mat &geom, int &n_atoms, arma::vec &z_val)
{
	ifstream myreadfile1("geom.dat");
	if(myreadfile1.is_open())
		{
			cout<<"\nReading the geometry file now\n"<<endl;
			myreadfile1 >> n_atoms;
			geom=arma::zeros(n_atoms,n_atoms);
			z_val=arma::zeros(n_atoms);
			for(int i=0; i<n_atoms; i++)
			{
				myreadfile1 >> z_val(i) >> geom(i,0) >> geom(i,1) >> 
				geom(i,2);
			}
		}
}


int index1(int i, int j)
{
	int ij=i*(i+1)/2+j;
	return ij;	
}


int index_call(int i, int j , int k, int l)
{
	int ij,kl,ijkl;
	if(i>j)
	{
		 ij=index1(i,j);
	}
	else    
	{
		 ij=index1(j,i);
	}
	if (k>l)
	{
		 kl=index1(k,l);
	}
	else
	{
		kl=index1(l,k);
	}
	if(ij > kl)
	{
		 ijkl=ij*(ij+1)/2+(kl);
	}
	else //if (ij < kl)

	{
		 ijkl=kl*(kl+1)/2+(ij);
	}

	return ijkl;
}


