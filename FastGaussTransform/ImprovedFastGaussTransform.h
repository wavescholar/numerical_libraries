//-------------------------------------------------------------------
// The code was written by Vikas Raykar and Changjiang Yang 
// and is copyrighted under the Lesser GPL: 
//
// Copyright (C) 2006 Vikas Raykar and Changjiang Yang 
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
// The author may be contacted via email at:
// vikas(at)umiacs(.)umd(.)edu, cyang(at)sarnoff(.)com
//-------------------------------------------------------------------


//-------------------------------------------------------------
// File    : ImprovedFastGaussTransform.h
// Purpose : Interface for 
//			 Data Adaptive Improved Fast Gauss Transform.
// Author  : Vikas C. Raykar (vikas@cs.umd.edu)
// Date    : July 15 2005
//-------------------------------------------------------------
// Improved Fast Gauss Transform (IFGT).
// All points have the same truncation number.
// ------------------------------------------------------------
//
// INPUTS [14] 
// ----------------
// Dim			  --> dimension of the points, d.
// NSources		  --> number of sources, N.
// MTargets		  --> number of targets, M.
// pSources		  --> pointer to sources, px(d*N).
// Bandwidth	  --> the source bandwidth, h.
// pWeights       --> pointer to the weights, pq(N).
// pTargets       --> pointer to the targets, py(d*M).
// MaxTruncNumber --> maximum truncation number for the 
//                    Taylor series, p_max.
// NumClusters    --> number of clusters (K).
// pClusterIndex  --> vector of length N where the i th element 
//                    is the cluster number to which the i th point
//					  belongs, pci (K).
//                    pClusterIndex[i] varies between 0 to K-1. 
// pClusterCenter --> pointer to the cluster centers, pcc(d*K).
// pClusterRadii  --> pointer to the cluster radii, pcr(K).
// CutoffRadius   --> source cutoff radius, r.
// epsilon        --> error, eps. 
//
// OUTPUTS [1]
// ----------------
// pGaussTransform --> pointer the the evaluated Gauss Transform, 
//					   pG(M).
//-------------------------------------------------------------------


#ifndef IMPROVED_FAST_GAUSS_TRANSFORM_H
#define IMPROVED_FAST_GAUSS_TRANSFORM_H


class ImprovedFastGaussTransform{
	public:
		//constructor 
		ImprovedFastGaussTransform(int Dim,
			int NSources,
			int MTargets,
			double *pSources,
			double Bandwidth,
			double *pWeights,
			double *pTargets,
			int MaxTruncNumber,
			int NumClusters,
			int *pClusterIndex, 
			double *pClusterCenter,
			double *pClusterRadii,
			double CutoffRadius,
			double epsilon,
		    double *pGaussTransform
			);		

		//destructor
		~ImprovedFastGaussTransform();				

		//function to evaluate the Gauss Transform.
		void Evaluate();

	private:
		//Parameters

		int d;				
		int N;				
		int M;				
		double *px;			
		double  h;			
		double *pq;			
		double *py;	
		int p_max;	
		int K;	
		int *pci; 
		double *pcc;  
		double *pcr;
		double r;
		double eps;


		double *pG;         
	
		//

		int     p_max_total;
		double *constant_series;
		double *source_center_monomials;
		double  source_center_distance_square;
		double *target_center_monomials;
		double  target_center_distance_square;
	    double *dx;
		double *dy;
	    int	   *heads;		
		double *C;
		double h_square;
		double *ry;
		double *ry_square;

		//Functions

		int  nchoosek(int n, int k);
		void compute_constant_series();
		void compute_source_center_monomials(int p);
		void compute_target_center_monomials();
		void compute_C();

		//MATLAB applications should always call mxMalloc rather than malloc to allocate memory

		//void *operator new[] (size_t s){ return mxMalloc(s);}
	 //   void operator delete[] (void* mem){ mxFree(mem);}

    
};


#endif