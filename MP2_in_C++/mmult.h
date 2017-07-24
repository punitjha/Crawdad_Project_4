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


void mmult(arma::mat &A, int transa, arma::mat &B, int transb, arma::mat &C, int l, int m, int n)
{
	int i, j, k;
 
	if(!transa && !transb)  
	{
		for(i=0; i < l; i++)  
		{
			for(j=0; j < m; j++) 
			{
				for(k=0; k < n; k++)  
				{
					C(i,j) += A(i,k) * B(k,j);
			       	}
		       	}
		}
	}
	else if(!transa && transb) 
      	{
		for(i=0; i < l; i++)  
		{
			for(j=0; j < m; j++)  
			{
				for(k=0; k < n; k++)  
				{
					C(i,j) += A(i,k) * B(j,k);
			       	}
			}
		}
	}

	else if(transa && !transb)  
	{
		for(i=0; i < l; i++)  
		{
			for(j=0; j < m; j++)  
			{
				for(k=0; k < n; k++)  
				{
				       	C(i,j) += A(k,i) * B(k,j);
				}
			}
		}
	}
	else if(transa && transb)  
	{
		for(i=0; i < l; i++)  
		{
			for(j=0; j < m; j++)  
			{
				for(k=0; k < n; k++)  
				{
					C(i,j) += A(k,i) * B(j,k);
				}
			}
		}
	}
}
