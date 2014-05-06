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
// Date    : July 08 2005
//-------------------------------------------------------------------

#include "mex.h"
#include "GaussTransform.h"


//The gateway function

void mexFunction(int nlhs,				// Number of left hand side (output) arguments
				 mxArray *plhs[],		// Array of left hand side arguments
				 int nrhs,              // Number of right hand side (input) arguments
				 const mxArray *prhs[])  // Array of right hand side arguments
{

   //check for proper number of arguments 
 
    if(nrhs != 7) mexErrMsgTxt("7 input  arguments required.");
	if(nlhs != 1) mexErrMsgTxt("1 output argument  required.");

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
  
  int mrows = (int)  mxGetM(prhs[argu]); //mrows
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
  double  Bandwidth = (double) mxGetScalar(prhs[argu]);
  if ( Bandwidth <= 0.0) mexErrMsgTxt("Input 'h' must be a positive number.");

 
  //----- the sixth input argument: pWeights--------------//
  argu = 5;

  /*  create a pointer to the input vector pWeights */
  double *pWeights = mxGetPr(prhs[argu]);

  mrows =(int) mxGetM(prhs[argu]); //mrows
  ncols =(int) mxGetN(prhs[argu]); //ncols
  if ( mrows != 1 && ncols != NSources)  mexErrMsgTxt("Input 'q' must be a 1 x N matrix");


   //----- the seventh input argument: pTargets--------------//
  argu = 6;

  /*  create a pointer to the input vector pTargets */
  double *pTargets = mxGetPr(prhs[argu]);
 
  mrows = (int)mxGetM(prhs[argu]); //mrows
  ncols = (int)mxGetN(prhs[argu]); //ncols
  if ( mrows != Dim && ncols != MTargets)  mexErrMsgTxt("Input 'Y' must be a d x M matrix");

  //////////////////////////////////////////////////////////////
  // Output arguments
  //////////////////////////////////////////////////////////////

   /*  set the output pointer to the output result(vector) */
  plhs[0] = mxCreateDoubleMatrix(1,MTargets,mxREAL);
  
  /*  create a C pointer to a copy of the output result(vector)*/
  double *pGaussTransform = mxGetPr(plhs[0]);

  //////////////////////////////////////////////////////////////
  // function calls;
  //////////////////////////////////////////////////////////////

  GaussTransform* pGT = new GaussTransform(
	  Dim,
	  NSources,
	  MTargets,
	  pSources,
	  Bandwidth,
	  pWeights,
	  pTargets,
	  pGaussTransform
	  );

  pGT->Evaluate();

  delete pGT;

  return;
  
}
