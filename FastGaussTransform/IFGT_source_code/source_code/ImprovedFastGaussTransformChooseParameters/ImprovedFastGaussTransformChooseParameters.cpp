//-------------------------------------------------------------------
// The code was written by Vikas C. Raykar 
// and is copyrighted under the Lesser GPL: 
//
// Copyright (C) 2006 Vikas C. Raykar
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 or later.
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the GNU Lesser General Public License for more details. 
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, 
// MA 02111-1307, USA.  
//
// The author may be contacted via email at: vikas(at)cs(.)umd(.)edu 
//-------------------------------------------------------------------

//-------------------------------------------------------------------
// File    : ImprovedFastGaussTransformChooseParameters.cpp
// Purpose : Implementation for the 
//           Improved Fast Gauss Transform Chooose Parameters
// Author  : Vikas C. Raykar (vikas@cs.umd.edu)
// Date    : August 22, 2005
//-------------------------------------------------------------------
// Parameter selection for the Improved Fast Gauss Transform (IFGT).
//
// Implementation based on:
//
// Fast computation of sums of Gaussians in high dimensions. 
// Vikas C. Raykar, C. Yang, R. Duraiswami, and N. Gumerov,
// CS-TR-4767, Department of computer science,
// University of Maryland, Collegepark.
//--------------------------------------------------------------------

#include "ImprovedFastGaussTransformChooseParameters.h"
#include <math.h>
#define  min(a,b) (((a)<(b))?(a):(b)) 

//-------------------------------------------------------------------
// Constructor
//-------------------------------------------------------------------

ImprovedFastGaussTransformChooseParameters::ImprovedFastGaussTransformChooseParameters(
			int Dim,
			double Bandwidth,
			double epsilon,
			int MaxNumClusters
		    )		
{	

	d=Dim;
	h=Bandwidth;
	Klimit=MaxNumClusters;
	eps=epsilon;
	R=sqrt((double)d);

	double h_square=h*h;

	r=min(R,h*sqrt((double)log(1.0f/eps)));

	int p_ul=200; // [Upper limit on the truncation number]I have roughly set this to 200

	K=1;
	int p=1;

	double complexity_min=1e16;

	double rx;
	
	for(int i=0;i <Klimit;i++){

		rx=pow((double)i+1,-1.0/(double)d);
		double rx_square=rx*rx;
		double n=min(i+1,pow(r/rx,(double)d));
		double error=1;
	    double temp=1;
		int p=0;
		while((error > eps) & (p <= p_ul)){
			p++;
			double b=min(((rx+sqrt((rx_square)+(2*p*h_square)))/2),rx+r);
			double c=rx-b;
			temp=temp*(((2*rx*b)/h_square)/p);
			error=temp*(exp(-(c*c)/h_square));			
		}	
		double complexity=(i+1)+log((double)i+1)+((1+n)*nchoosek(p-1+d,d));

	   //printf("d=%d r=%f K=%d rx=%f n=%f p=%d terms=%d c=%f\n",d,r,i+1,rx,n,p,nchoosek(p-1+d,d),complexity);
		
		if (complexity < complexity_min ){
			complexity_min=complexity;
			K=i+1;
			p_max=p;		
		}
	}
	
		
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------

ImprovedFastGaussTransformChooseParameters::~ImprovedFastGaussTransformChooseParameters()
{
	
}

//-------------------------------------------------------------------
// Compute the combinatorial number nchoosek.
//-------------------------------------------------------------------

int
ImprovedFastGaussTransformChooseParameters::nchoosek(int n, int k){
	int n_k = n - k;
	
	if (k < n_k)
	{
		k = n_k;
		n_k = n - k;
	}

	int  nchsk = 1; 
	for ( int i = 1; i <= n_k; i++)
	{
		nchsk *= (++k);
		nchsk /= i;
	}

	return nchsk;
}


