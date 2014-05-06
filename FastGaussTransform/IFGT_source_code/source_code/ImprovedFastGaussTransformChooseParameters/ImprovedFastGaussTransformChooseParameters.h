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


//---------------------------------------------------------------------
// File    : ImprovedFastGaussTransformChooseParameters.h
// Purpose : Header file for
//          ImprovedFastGaussTransformChooseParameters.cpp
// Author  : Vikas C. Raykar (vikas@cs.umd.edu)
// Date    : August 22, 2005
//--------------------------------------------------------------------
// Parameter selection for the Improved Fast Gauss Transform (IFGT).
//
// Implementation based on:
//
// Fast computation of sums of Gaussians in high dimensions. 
// Vikas C. Raykar, C. Yang, R. Duraiswami, and N. Gumerov,
// CS-TR-4767, Department of computer science,
// University of Maryland, Collegepark.
//--------------------------------------------------------------------

#ifndef IMPROVED_FAST_GAUSS_TRANSFORM_CHOOSE_PARAMETERS_H
#define IMPROVED_FAST_GAUSS_TRANSFORM_CHOOSE_PARAMETERS_H

#include <mex.h>

class ImprovedFastGaussTransformChooseParameters{
	public:
		//constructor 

		ImprovedFastGaussTransformChooseParameters(
			int Dim,
			double Bandwidth,
			double epsilon,
			int MaxNumClusters
		    );		

		//destructor

		~ImprovedFastGaussTransformChooseParameters();				
		
		int K;
		int p_max;
		double r;

	private:
		
		//Parameters

		int d;		
		double h;
		double eps;
		int Klimit;
		double R;

		int nchoosek(int n, int k);

		//MATLAB applications should always call mxMalloc rather than malloc to allocate memory

		void *operator new[] (size_t s){ return mxMalloc(s);}
	    void operator delete[] (void* mem){ mxFree(mem);}
    
};


#endif