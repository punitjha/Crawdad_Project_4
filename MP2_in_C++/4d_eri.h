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
using namespace std;

void d_eri( vector<vector<vector<vector<double> > > > & deri, int nao) 
{
	arma::mat d_file;
	d_file.load("eri.dat");
	for(int i1=0; i1<d_file.n_rows; i1++)
	{
		int i=d_file(i1,0)-1;
		int j=d_file(i1,1)-1;
		int k=d_file(i1,2)-1; 
		int l=d_file(i1,3)-1;
		double val=d_file(i1,4);
		deri[i][j][k][l]=deri[i][j][l][k]=deri[j][i][k][l]=deri[j][i][l][k]=val;
		deri[k][l][i][j]=deri[l][k][i][j]=deri[k][l][j][i]=deri[l][k][j][i]=val;
	}

	for (int i=0; i<nao; i++)
	{
		for(int j=0; j<nao; j++)
		{
			for(int k=0; k<nao; k++)
			{
				for(int l=0; l<nao; l++)
				{
					cout<<i<<j<<k<<l<<deri[i][j][k][l]<<endl;
				}
			}
		}
	}
}	
