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
// File    : ImprovedFastGaussTransformChooseTruncationNumber.h
// Purpose : Header file for ChooseTruncNumber.cpp
// Author  : Vikas C. Raykar (vikas@cs.umd.edu)
// Date    : August 22, 2005
//--------------------------------------------------------------------
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

#ifndef CHOOSE_TRUNC_NUMBER_H
#define CHOOSE_TRUNC_NUMBER_H

#include <mex.h>

class ImprovedFastGaussTransformChooseTruncationNumber{
	public:
		//constructor 
		ImprovedFastGaussTransformChooseTruncationNumber(
			int Dim,
			double Bandwidth,
			double epsilon,
			double MaxClusterRadius
			);		

		//destructor
		~ImprovedFastGaussTransformChooseTruncationNumber();				
		
		int p_max;

	private:

		//Parameters

		int d;		
		double h;
		double eps;
		double R;
		double r;
		double rx;
    
};


#endif