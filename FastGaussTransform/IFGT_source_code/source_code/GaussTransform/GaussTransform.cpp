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
// File    : GaussTransform.cpp
// Purpose : Implementation for Gauss Transform.
// Author  : Vikas C. Raykar (vikas@cs.umd.edu)
// Date    : July 08 2005
//-------------------------------------------------------------------

#include "GaussTransform.h"
#include <math.h>

//-------------------------------------------------------------------
// Constructor 
//
// PURPOSE                                                    
// -------   
// Initialize the class. 
// Read the parameters.
//
// PARAMETERS                                                      
// ----------
// Dim			   --> dimension of the points.
// NSources		   --> number of sources.
// MTargets		   --> number of targets.
// pSources		   --> pointer to sources, (d*N).
// Bandwidth	   --> source bandwidth.
// pWeights        --> pointer to the weights, (N).
// pTargets        --> pointer to the targets, (d*M).
// pGaussTransform --> pointer the the evaluated Gauss Transform,(M).
//-------------------------------------------------------------------

GaussTransform::GaussTransform(int Dim,
			int NSources,
			int MTargets,
			double *pSources,
			double Bandwidth,
			double *pWeights,
			double *pTargets,
			double *pGaussTransform)
{	

	d=Dim;
	N=NSources;
	M=MTargets;
	px=pSources;
	h=Bandwidth;
	pq=pWeights;
	py=pTargets;
	pG=pGaussTransform;

}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------

GaussTransform::~GaussTransform()
{
}

//-------------------------------------------------------------------
// Actual function to evaluate the Gauss Transform.
//-------------------------------------------------------------------

void
GaussTransform::Evaluate()
{
	double h_square=h*h;

	for(int j=0; j<M; j++)
	{
		pG[j]=0.0;

		for(int i=0; i<N; i++)
		{
			double norm=0.0;
			for (int k=0; k<d; k++)
			{
				double temp=px[(d*i)+k]-py[(d*j)+k];
				norm = norm + (temp*temp);
			}
		
			pG[j] = pG[j]+(pq[i]*exp(-norm/h_square));			

		}
	}
}

