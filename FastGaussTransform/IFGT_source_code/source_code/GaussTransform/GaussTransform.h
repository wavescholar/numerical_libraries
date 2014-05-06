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

//-------------------------------------------------------------
// File    : GaussTransform.h
// Purpose : Header file for GaussTransform.cpp
// Author  : Vikas C. Raykar (vikas@cs.umd.edu)
// Date    : July 08 2005
//-------------------------------------------------------------


#ifndef GAUSS_TRANSFORM_H
#define GAUSS_TRANSFORM_H

class GaussTransform{
	public:
		//constructor 
		GaussTransform(int Dim,
			int NSources,
			int MTargets,
			double *pSources,
			double Bandwidth,
			double *pWeights,
			double *pTargets,
			double *pGaussTransform);

		//destructor
		~GaussTransform();

		//function to evaluate the Gauss Transform.
		void Evaluate();

	private:
		int d;				//dimension of the points.
		int N;				//number of sources.
		int M;				//number of targets.
		double *px;			//pointer to sources, (d*N).
		double  h;			//the source bandwidth.
		double *pq;			//pointer to the weights, (N).
		double *py;		    //pointer to the targets, (d*M).
		double *pG;         //pointer the the evaluated Gauss Transform, (M).

    
};


#endif