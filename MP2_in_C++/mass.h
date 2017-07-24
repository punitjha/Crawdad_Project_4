//Header for the masses

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
/***********************************************************************************************************/	
using namespace std;

//creating a static array which stores the atomic mass indexed by the array element
//0
//hydrogen
//lithium
double mass[]= {
			0.00000,
			1.007825,
			4.002603,
			6.9415123,
			9.0122,
			10.811,
			12.0107,
			14.0067,
			15.9994,
			18.9984,
			20.1797,
			22.9897,
			24.305,
			26.981,
			28.055,
			30.9738,
			32.065,
			35.4530,
			29.0948,
			39.948,
			40.078
		};


