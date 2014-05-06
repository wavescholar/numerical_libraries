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


//-------------------------------------------------------------------
// File    : mexmain.cpp
// Purpose : Interface between MATLAB and C++
// Author  : Vikas C. Raykar (vikas@cs.umd.edu)
// Date    : July 15 2005
//-------------------------------------------------------------------
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
// OUTPUTS [2]
// ----------------
// pGaussTransform --> pointer the the evaluated Gauss Transform, 
//					   pG(M).
// pTruncNumber    --> truncation number used for each source, pT(N).
//-------------------------------------------------------------------



#include "mex.h"
#include "DataAdaptiveImprovedFastGaussTransform.h"

//The gateway function

void mexFunction(int nlhs,				// Number of left hand side (output) arguments
				 mxArray *plhs[],		// Array of left hand side arguments
				 int nrhs,              // Number of right hand side (input) arguments
				 const mxArray *prhs[])  // Array of right hand side arguments
{

  //check for proper number of arguments 
 
  if(nrhs != 14) mexErrMsgTxt("14 inputs required.");
  if(nlhs != 2) mexErrMsgTxt("2 outputs required.");

   //////////////////////////////////////////////////////////////
  // Input arguments
  //////////////////////////////////////////////////////////////
  
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

  //------ the second input argument: NSources ---------------//

  argu = 1;

  /* check to make sure the input argument is a scalar */
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu])!=1 ) 
  {
    mexErrMsgTxt("Input 'N' must be a scalar.");
  }

  /*  get the scalar input NSources */
  int NSources = (int) mxGetScalar(prhs[argu]);
  if (NSources <= 0) mexErrMsgTxt("Input 'N' must be a positive number.");

   //------ the third input argument: MTargets ---------------//

  argu = 2;

  /* check to make sure the input argument is a scalar */
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu])!=1 ) 
  {
    mexErrMsgTxt("Input 'M' must be a scalar.");
  }

  /*  get the scalar input NSources */
  int MTargets = (int) mxGetScalar(prhs[argu]);
  if (MTargets <= 0) mexErrMsgTxt("Input 'M' must be a positive number.");

  //----- the fourth input argument: pSources--------------//
  //  The 2D array is column-major: each column represents a point.

  argu = 3;

  /*  create a pointer to the input vector pSources */
  double *pSources = mxGetPr(prhs[argu]);
  
  int mrows = (int) mxGetM(prhs[argu]); //mrows
  int ncols = (int) mxGetN(prhs[argu]); //ncols
  if ( mrows != Dim && ncols != NSources)  mexErrMsgTxt("Input 'X' must be a d x N matrix");

  //----- the fifth input argument: Bandwidth--------------//
  argu = 4;

  /* check to make sure the input argument is a scalar */
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu])!=1 ) 
  {
    mexErrMsgTxt("Input 'h' must be a scalar.");
  }

  /*  get the scalar input CutoffRadius */
  double Bandwidth = (double) mxGetScalar(prhs[argu]);
  if (Bandwidth <= 0.0) mexErrMsgTxt("Input 'h' must be a positive number.");

 
  //----- the sixth input argument: pWeights--------------//
  argu = 5;

  /*  create a pointer to the input vector pWeights */
  double *pWeights = mxGetPr(prhs[argu]);

  mrows = (int)  mxGetM(prhs[argu]); //mrows
  ncols = (int) mxGetN(prhs[argu]); //ncols
  if ( mrows != 1 && ncols != NSources)  mexErrMsgTxt("Input 'q' must be a 1 x N matrix");


   //----- the seventh input argument: pTargets--------------//
  argu = 6;

  /*  create a pointer to the input vector pTargets */
  double *pTargets = mxGetPr(prhs[argu]);
 
  mrows = (int) mxGetM(prhs[argu]); //mrows
  ncols = (int) mxGetN(prhs[argu]); //ncols
  if ( mrows != Dim && ncols != MTargets)  mexErrMsgTxt("Input 'Y' must be a d x M matrix");

   //------ the eighth input argument: p_max ---------------//

  argu = 7;

  /* check to make sure the input argument is a scalar */
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu])!=1 ) 
  {
    mexErrMsgTxt("Input 'p_max' must be a scalar.");
  }

  /*  get the scalar input TruncNumber1 */
  int MaxTruncNumber = (int) mxGetScalar(prhs[argu]);
  if (MaxTruncNumber <= 0) mexErrMsgTxt("Input 'p_max' must be a positive number.");

  
  //------ the ninth input argument: K ---------------//

  argu = 8;

   /* check to make sure the input argument is a scalar */
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu])!=1 ) 
  {
    mexErrMsgTxt("Input 'K' must be a scalar.");
  }

  /*  get the scalar input TruncNumber1 */
  int NumClusters = (int) mxGetScalar(prhs[argu]);
  if (NumClusters <= 0) mexErrMsgTxt("Input 'K' must be a positive number.");

 
  //----- the tenth input argument: pClusterIndex--------------//
  argu = 9;

  /*  create a pointer to the input vector pClusterIndex */
  int *pClusterIndex =(int*) mxGetPr(prhs[argu]);

  mrows = (int) mxGetM(prhs[argu]); //mrows
  ncols = (int) mxGetN(prhs[argu]); //ncols
  if ( mrows != 1 && ncols != NSources)  mexErrMsgTxt("Input 'ClusterIndex' must be a 1 x N matrix");

   //----- the eleventh input argument: pClusterCenter--------------//
  argu = 10;

  /*  create a pointer to the input vector pTargets */
  double * pClusterCenter = mxGetPr(prhs[argu]);
 
  mrows =(int)  mxGetM(prhs[argu]); //mrows
  ncols = (int) mxGetN(prhs[argu]); //ncols
  if ( mrows != Dim && ncols != NumClusters)  mexErrMsgTxt("Input  'ClusterCenter' must be a d x K matrix");

   //----- the twelevth input argument: pClusterRadii--------------//
  argu = 11;

  /*  create a pointer to the input vector pClusterRadii */
  double * pClusterRadii = mxGetPr(prhs[argu]);
 
  mrows = (int) mxGetM(prhs[argu]); //mrows
  ncols = (int) mxGetN(prhs[argu]); //ncols
  if ( mrows != 1 && ncols != NumClusters)  mexErrMsgTxt("Input 'ClusterRadii' must be a 1 x K matrix");

   //------ the thirteenth input argument: r ---------------//

  argu = 12;

  /* check to make sure the input argument is a scalar */
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu])!=1 ) 
  {
    mexErrMsgTxt("Input 'r' must be a scalar.");
  }

  /*  get the scalar input CutoffRadius */
  double CutoffRadius = (double) mxGetScalar(prhs[argu]);
  if (CutoffRadius <= 0.0) mexErrMsgTxt("Input 'r' must be a positive number.");

   //------ the 14th input argument: epsilon ---------------//

  argu = 13;

  /* check to make sure the input argument is a scalar */
  if( !mxIsDouble(prhs[argu]) || mxIsComplex(prhs[argu]) || mxGetN(prhs[argu])*mxGetM(prhs[argu])!=1 ) 
  {
    mexErrMsgTxt("Input 'epsilon' must be a scalar.");
  }

  /*  get the scalar input CutoffRadius */
  double epsilon = (double) mxGetScalar(prhs[argu]);
  if (epsilon <= 0.0) mexErrMsgTxt("Input 'epsilon' must be a positive number.");


  //////////////////////////////////////////////////////////////
  // Output arguments
  //////////////////////////////////////////////////////////////


  plhs[0] = mxCreateDoubleMatrix(1,MTargets,mxREAL);
  double *pGaussTransform = mxGetPr(plhs[0]);

  plhs[1]=mxCreateNumericMatrix(1,NSources,mxUINT32_CLASS,mxREAL);
  int *pTruncNumber =(int*) mxGetPr(plhs[1]);


  //////////////////////////////////////////////////////////////
  // function calls;
  //////////////////////////////////////////////////////////////

  DataAdaptiveImprovedFastGaussTransform* pIFGT = new  DataAdaptiveImprovedFastGaussTransform(
	  Dim,
	  NSources,
	  MTargets,
	  pSources,
	  Bandwidth,
	  pWeights,
	  pTargets,
	  MaxTruncNumber,
	  NumClusters,
	  pClusterIndex,
	  pClusterCenter,
	  pClusterRadii,
	  CutoffRadius,
	  epsilon,
	  pGaussTransform,
	  pTruncNumber
	  );

  pIFGT->Evaluate();

  delete pIFGT;

  return;
  
}
