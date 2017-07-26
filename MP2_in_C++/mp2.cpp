//The mass header file

#include "mmult.h"
#include "hartree.h"
#include "4d_eri.h"
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
#include "boost/multi_array.hpp"
#include <cassert>
#include <array>   // to be able to use this library compile using this g++ -std=c++0x FileNames
using namespace std;

int main()
{

//*****************************************************************************
//Starting the mp2 calculations here to  verify the python codes
//*****************************************************************************

	int nao,OCC;
	cout<<"Please enter the number of atomic orbital\n";
	cin>>nao;
	cout<<"Please enter the number of occupied orbital\n";
	cin>>OCC;
	arma::vec d_I=arma::zeros(1000);	
	arma::mat coeff_car=arma::zeros(nao,nao);
	arma::vec eps=arma::zeros(nao);
	hartree(nao, OCC,d_I, coeff_car, eps);

	arma::vec d_MO=arma::zeros((nao*(nao+1)/2)*((nao*(nao+1)/2)+1)/2);
	int ijkl=0;
	int pqrs=0;
	for (int i=0; i<nao; i++)
	{
		for (int j=0; j<=i; j++)
		{
			for (int k=0; k<=i; k++)
			{
				for ( int l=0; l<= (i==k ? j : k); l++, ijkl++)
				{
					for ( int p=0; p<nao; p++)
					{
						for ( int q=0; q<nao; q++)
						{
							for(int r=0; r<nao; r++)
							{
								for(int s=0; s<nao; s++)
								{
									pqrs=index_call(p+1,q+1,r+1,s+1);
									d_MO(ijkl)+= coeff_car(p,i)*coeff_car(q,j)*coeff_car(r,k)*coeff_car(s,l)*d_I(pqrs);
								}
							}
						}
					}
				}
			}
		}
	}
//*****************************************************************************
//This is the implementation of the better transformation algorithm
//*****************************************************************************
//	arma::mat X=arma::zeros(nao,nao);
//	arma::mat Y=arma::zeros(nao,nao);
//	arma::mat TMP=arma::zeros(nao*(nao+1)/2,nao*(nao+1)/2);
//	arma::vec d_MO_f=arma::zeros((nao*(nao+1)/2)*((nao*(nao+1)/2)+1)/2);
//	int i,j,k,l,ij,kl,klij;
//	ijkl=0;
//        for(i=0,ij=0; i<nao; i++)
//		for(j=0; j<=i; j++,ij++){
//			for(k=0,kl=0; k<nao; k++)
//				for(l=0; l<=k;l++,kl++){
//					ijkl=index2(ij,kl);
//					X(k,l)=X(l,k)=d_I(ijkl);
//				}
//			Y=arma::zeros(nao,nao);
//			Y=coeff_car.t()*X;
//			//mmult(coeff_car,1,X,0,Y,nao,nao,nao);
//			Y.print();
//			X=arma::zeros(nao,nao);
//			//X=Y*coeff_car;
//			mmult(Y,0,coeff_car,0,X,nao,nao,nao);
//			for(k=0, kl=0; k<nao; k++)
//				for(l=0; l<=k; l++, kl++)
//					TMP(kl,ij)=X(k,l);
//	}
//	for(k=0, kl=0; k<nao; k++)
//		for(l=0; l<= k; l++, kl++){
//			X=arma::zeros(nao,nao);
//			Y=arma::zeros(nao,nao);
//			for(i=0, ij=0; i<nao; i++)
//				for(j=0; j<=i; j++, ij++)
//					X(i,j)=X(j,i)=TMP(kl,ij);
//			Y=arma::zeros(nao,nao);
//			//Y=coeff_car.t()*X;
//			mmult(coeff_car,1,X,0,Y,nao,nao,nao);
//			X=arma::zeros(nao,nao);
//			//X=Y*coeff_car;
//			mmult(Y,0,coeff_car,0,X,nao,nao,nao);
//			for(i=0, ij=0; i<nao; i++)
//				for(j=0;j<=i;j++,ij++){
//					klij=index2(kl,ij);
//					d_MO_f(klij)=X(i,j);
//				}
//			}
//	d_MO_f.print("\n\nThis is the new MO array \n");



//*****************************************************************************
//This is the implementation of the better transformation algorithm
//*****************************************************************************

//	typedef boost::multi_array<double,4> array_type;
//	typedef array_type::index index;
//	array_type deri(boost::extents[nao][nao][nao][nao]);
//	double *p=deri.data();
//	double (*array)[nao][nao][nao][nao]=(double (*)[nao][nao][nao][nao])p;

	vector<vector<vector<vector<double> > > > deri (nao,vector<vector<vector<double> > >(nao,vector<vector<double> > (nao,vector<double>(nao))));
	d_eri(deri, nao);

//*****************************************************************************
//This is the implementation of the mp2 algorithm 
//*****************************************************************************
	double emp2=0.0;
	for(int i=0; i<OCC; i++)
	{
		for(int j=0;j< OCC; j++)
		{
			for( int b=OCC; b<nao; b++)
			{
				for(int a=OCC; a<nao; a++)
				{
					int iajb=index_call(i,a,j,b);
					int ibja=index_call(i,b,j,a);
					emp2+=d_MO(iajb)*(2*d_MO(iajb)-d_MO(ibja))/(eps(i)+eps(j)-eps(a)-eps(b));
				}
			}
		}
	}
	cout<<"\n\nThis is the mp2 energy:"<<emp2<<endl;
}


