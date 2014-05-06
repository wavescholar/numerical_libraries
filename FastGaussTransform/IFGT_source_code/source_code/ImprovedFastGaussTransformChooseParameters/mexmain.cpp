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
// File    : mexmain.cpp
// Purpose : Interface between MATLAB and C++
// Author  : Vikas C. Raykar (vikas@cs.umd.edu)
// Date    : August 22, 2005
//-------------------------------------------------------------------
// INPUTS [4] 
//-------------------------------------------------------------------
// Dim			  --> dimension of the points, d.
// Bandwidth	  --> the source bandwidth, h.
// epsilon        --> the desired error, eps.   
// MaxNumClusters --> upper limit on the number of clusters, Klimit.  
//-------------------------------------------------------------------
// OUTPUTS [3]
//-------------------------------------------------------------------
// NumClusters    --> number of clusters, K.
// MaxTruncNumber --> maximum truncation number for the 
//					  Taylor series, p_max.
// CutoffRadius   --> source cutoff radius, r.
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

#include "mex.h"
#include "ImprovedFastGaussTransformChooseParameters.h"

//The gateway function

void mexFunction(int nlhs,				 // Number of left hand side (output) arguments
				 mxArray *plhs[],		 // Array of left hand side arguments
				 int nrhs,               // Number of right hand side (input) arguments
				 const mxArray *prhs[] )  // Array of right hand side arguments
{

  //check for proper number of arguments 
 
  if(nrhs != 4) mexErrMsgTxt("4 inputs required.");
  if(nlhs != 3) mexErrMsgTxt("3 outputs required.");

  //-------------------------------------------------------------------
  // Input arguments
  //-------------------------------------------------------------------
  
  //------ the first input argument: Dim ---------------//

  int argu = 0;

  /* check to make sure the input argument is a scalar */
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu])!=1 ) 
  {
    mexErrMsgTxt("Input 'd' must be a scalar.");
  }

  /*  get the scalar input Dim */
  int Dim = (int) mxGetScalar(prhs[argu]);
  if (Dim <= 0) mexErrMsgTxt("Input 'd' must be a positive number.");

 
  //----- the second input argument: Bandwidth--------------//
  argu = 1;

  /* check to make sure the input argument is a scalar */
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu])!=1 ) 
  {
    mexErrMsgTxt("Input 'h' must be a scalar.");
  }

  /*  get the scalar input CutoffRadius */
  double Bandwidth = (double) mxGetScalar(prhs[argu]);
  if (Bandwidth <= 0.0) mexErrMsgTxt("Input 'h' must be a positive number.");

   //------ the third input argument: epsilon ---------------//

  argu = 2;

  /* check to make sure the input argument is a scalar */
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu])!=1 ) 
  {
    mexErrMsgTxt("Input 'epsilon' must be a scalar.");
  }

  /*  get the scalar input CutoffRadius */
  double epsilon = (double) mxGetScalar(prhs[argu]);
  if (epsilon <= 0.0) mexErrMsgTxt("Input 'epsilon' must be a positive number.");


  //------ the fourth input argument: MaxNumClusters ---------------//

  argu = 3;

   /* check to make sure the input argument is a scalar */
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu])!=1 ) 
  {
    mexErrMsgTxt("Input 'Klimit' must be a scalar.");
  }

  /*  get the scalar input TruncNumber1 */
  int MaxNumClusters = (int) mxGetScalar(prhs[argu]);
  if (MaxNumClusters <= 0) mexErrMsgTxt("Input 'Klimit' must be a positive number.");


  //-------------------------------------------------------------------
  // function call
  //-------------------------------------------------------------------
  
 ImprovedFastGaussTransformChooseParameters* pIFGTCP = 
	  new  ImprovedFastGaussTransformChooseParameters(
	  Dim,
	  Bandwidth,
	  epsilon,
	  MaxNumClusters
	  );

  //-------------------------------------------------------------------
  // Output arguments
  //-------------------------------------------------------------------

  plhs[0]=mxCreateNumericMatrix(1,1,mxUINT32_CLASS,mxREAL);
  int *NumClusters =(int*) mxGetPr(plhs[0]);
  NumClusters[0]=pIFGTCP->K;

  plhs[1]=mxCreateNumericMatrix(1,1,mxUINT32_CLASS,mxREAL);
  int *MaxTruncNumber =(int*) mxGetPr(plhs[1]);
  MaxTruncNumber[0]=pIFGTCP->p_max;

  plhs[2]=mxCreateDoubleMatrix(1,1,mxREAL);
  double * CutoffRadius = mxGetPr(plhs[2]);
  CutoffRadius[0]=pIFGTCP->r;

  delete pIFGTCP ;

  return;
  
}
