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
// File    : ImprovedFastGaussTransformChooseTruncationNumber.cpp
// Purpose : Implementation for the 
//           ImprovedFastGaussTransformChooseTruncationNumber.
// Author  : Vikas C. Raykar (vikas@cs.umd.edu)
// Date    : August 22, 2005
//-------------------------------------------------------------------
// Automatic parameter selection for the
// Improved Fast Gauss Transform (IFGT).
//
// Given the maximum cluster radius this routine computes the
// required truncation number.
//
// Implementation based on:
//
// Fast computation of sums of Gaussians in high dimensions. 
// Vikas C. Raykar, C. Yang, R. Duraiswami, and N. Gumerov,
// CS-TR-4767, Department of computer science,
// University of Maryland, Collegepark.
//--------------------------------------------------------------------

#include "ImprovedFastGaussTransformChooseTruncationNumber.h"
#include <math.h>
#define  min(a,b) (((a)<(b))?(a):(b)) 

ImprovedFastGaussTransformChooseTruncationNumber::ImprovedFastGaussTransformChooseTruncationNumber(
			int Dim,
			double Bandwidth,
			double epsilon,
			double MaxClusterRadius
		    )		
{	

	d=Dim;
	h=Bandwidth;
	rx=MaxClusterRadius;
	eps=epsilon;
	R=sqrt((double)d);

	double h_square=h*h;

	r=min(R,h*sqrt(log(1/eps)));

	int p_ul=300;

	double rx_square=rx*rx;
	
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

	p_max=p;		
		
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------

ImprovedFastGaussTransformChooseTruncationNumber::~ImprovedFastGaussTransformChooseTruncationNumber()
{
	
}

