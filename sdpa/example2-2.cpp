/* -------------------------------------------------------------

This file is a component of SDPA
Copyright (C) 2004 SDPA Project

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

------------------------------------------------------------- */
/* The beginning of the "example2-2.cpp" */
#include <stdio.h>
#include <stdlib.h>

#include "sdpa-lib.hpp"
#include "sdpa-lib2.hpp"

/*
example2.dat:

*Example 2:
*mDim = 5, nBLOCK = 3, {2,3,-2}
   5  =  mDIM
   3  =  nBLOCK
   2    3   -2   = bLOCKsTRUCT
{1.1, -10, 6.6 , 19 , 4.1}
{
{ { -1.4, -3.2 },
  { -3.2,-28   }   }
{ { 15,  -12,    2.1 },
  {-12,   16,   -3.8 },
  {  2.1, -3.8, 15   }   }
  {  1.8, -4.0 }
}
{
{ {  0.5,  5.2 },
  {  5.2, -5.3 }   }
{ {  7.8, -2.4,  6.0 },
  { -2.4,  4.2,  6.5 },
  {  6.0,  6.5,  2.1 }   }
  { -4.5, -3.5 }
}
{
{ { 1.7,  7.0 },
  { 7.0, -9.3 }   }
{ {-1.9, -0.9, -1.3 },
  {-0.9, -0.8, -2.1 },
  {-1.3, -2.1,  4.0 }   }
  {-0.2, -3.7 }
}
{
{ { 6.3, -7.5 },
  {-7.5, -3.3 }   }
{ { 0.2,  8.8,  5.4 },
  { 8.8,  3.4, -0.4 },
  { 5.4, -0.4,  7.5 }   }
  {-3.3, -4.0 }
}
{
{ { -2.4, -2.5 },
  { -2.5, -2.9 }   }
{ {  3.4, -3.2, -4.5 },
  { -3.2,  3.0, -4.8 },
  { -4.5, -4.8,  3.6 }   }
  {  4.8 , 9.7 }
}
{
{ { -6.5, -5.4 },
  { -5.4, -6.6 }   }
{ {  6.7, -7.2, -3.6 },
  { -7.2,  7.3, -3.0 },
  { -3.6, -3.0, -1.4 }   }
  {  6.1, -1.5 }
}
*/

void PrintMatrix(int nRow, int nCol, double* element, FILE* fpOut);
void PrintVector(int nDim, double* element, FILE* fpOut);

extern "C" int sdpa_EX22 ()
{
	double* 	element;

	SDPA	Problem1;

        Problem1.Method          	= KSH;
	Problem1.InitialPoint	 	= false;
	
	Problem1.pARAM.maxIteration 	= 50;
        Problem1.pARAM.epsilonStar   	= 1.0E-8;
        Problem1.pARAM.lambdaStar    	= 1.0E2;
        Problem1.pARAM.omegaStar     	= 2.0;
        Problem1.pARAM.lowerBound    	= -1.0E5;
        Problem1.pARAM.upperBound    	= 1.0E5;
        Problem1.pARAM.betaStar       	= 0.1;
        Problem1.pARAM.betaBar       	= 0.2;
        Problem1.pARAM.gammaStar     	= 0.9;

	Problem1.DisplayInformation 	= stdout;
	

	SDPA_initialize(Problem1);


//	mDim = 5, nBLOCK = 3, {2,3,-2}
//	5  =  mDIM
//	3  =  nBLOCK
//	2    3   -2   = bLOCKsTRUCT
	Problem1.mDIM	= 5;
	Problem1.nBLOCK	= 3;
	Problem1.bLOCKsTRUCT	= new int [Problem1.nBLOCK];
	Problem1.bLOCKsTRUCT[0] = 2;
	Problem1.bLOCKsTRUCT[1] = 3;
	Problem1.bLOCKsTRUCT[2] = -2;

	
	SDPA_initialize2(Problem1);


//	cVECT = {1.1, -10, 6.6 , 19 , 4.1}
	SDPA_Input_cVECT(Problem1, 1, 1.1);
	SDPA_Input_cVECT(Problem1, 2, -10);
	SDPA_Input_cVECT(Problem1, 3, 6.6);
	SDPA_Input_cVECT(Problem1, 4, 19);
	SDPA_Input_cVECT(Problem1, 5, 4.1);


	for (int i = 0; i <= Problem1.mDIM; i++)
	{
		SDPA_CountUpperTriangle(Problem1, i, 1, 3);
		SDPA_CountUpperTriangle(Problem1, i, 2, 6);
		SDPA_CountUpperTriangle(Problem1, i, 3, 2);
	}

	SDPA_Make_sfMAT(Problem1);

//	Input F_0

//	1st block
//	{ { -1.4, -3.2 },
//	  { -3.2,-28   }   }
	SDPA_InputElement(Problem1, 0, 1, 1, 1, -1.4);
	SDPA_InputElement(Problem1, 0, 1, 1, 2, -3.2);
	SDPA_InputElement(Problem1, 0, 1, 2, 2, -28);

//	2nd block
//	{ { 15,  -12,    2.1 },
//	  {-12,   16,   -3.8 },
//	  {  2.1, -3.8, 15   }   }
	SDPA_InputElement(Problem1, 0, 2, 1, 1, 15);
	SDPA_InputElement(Problem1, 0, 2, 1, 2, -12);
	SDPA_InputElement(Problem1, 0, 2, 1, 3, 2.1);
	SDPA_InputElement(Problem1, 0, 2, 2, 2, 16);
	SDPA_InputElement(Problem1, 0, 2, 2, 3, -3.8);
	SDPA_InputElement(Problem1, 0, 2, 3, 3, 15);

//	3rd block
//	{  1.8, -4.0 }
	SDPA_InputElement(Problem1, 0, 3, 1, 1, 1.8);
	SDPA_InputElement(Problem1, 0, 3, 2, 2, -4.0);


//	Input F_1

//	1st block
//	{ {  0.5,  5.2 },
//	  {  5.2, -5.3 }   }
	SDPA_InputElement(Problem1, 1, 1, 1, 1, 0.5);
	SDPA_InputElement(Problem1, 1, 1, 1, 2, 5.2);
	SDPA_InputElement(Problem1, 1, 1, 2, 2, -5.3);

//	2nd block
//	{ {  7.8, -2.4,  6.0 },
//	  { -2.4,  4.2,  6.5 },
//	  {  6.0,  6.5,  2.1 }   }
	SDPA_InputElement(Problem1, 1, 2, 1, 1, 7.8);
	SDPA_InputElement(Problem1, 1, 2, 1, 2, -2.4);
	SDPA_InputElement(Problem1, 1, 2, 1, 3, 6.0);
	SDPA_InputElement(Problem1, 1, 2, 2, 2, 4.2);
	SDPA_InputElement(Problem1, 1, 2, 2, 3, 6.5);
	SDPA_InputElement(Problem1, 1, 2, 3, 3, 2.1);

//	3rd block
//	{ -4.5, -3.5 }
	SDPA_InputElement(Problem1, 1, 3, 1, 1, -4.5);
	SDPA_InputElement(Problem1, 1, 3, 2, 2, -3.5);


//	Input F_2

//	1st block
//	{ { 1.7,  7.0 },
//	  { 7.0, -9.3 }   }
	SDPA_InputElement(Problem1, 2, 1, 1, 1, 1.7);
	SDPA_InputElement(Problem1, 2, 1, 1, 2, 7.0);
	SDPA_InputElement(Problem1, 2, 1, 2, 2, -9.3);

//	2nd block
//	{ {-1.9, -0.9, -1.3 },
//	  {-0.9, -0.8, -2.1 },
//	  {-1.3, -2.1,  4.0 }   }
	SDPA_InputElement(Problem1, 2, 2, 1, 1, -1.9);
	SDPA_InputElement(Problem1, 2, 2, 1, 2, -0.9);
	SDPA_InputElement(Problem1, 2, 2, 1, 3, -1.3);
	SDPA_InputElement(Problem1, 2, 2, 2, 2, -0.8);
	SDPA_InputElement(Problem1, 2, 2, 2, 3, -2.1);
	SDPA_InputElement(Problem1, 2, 2, 3, 3, 4.0);

//	3rd blcok
//	{-0.2, -3.7 }
	SDPA_InputElement(Problem1, 2, 3, 1, 1, -0.2);
	SDPA_InputElement(Problem1, 2, 3, 2, 2, -3.7);


//	Input F_3

//	1st block
//	{ { 6.3, -7.5 },
//	  {-7.5, -3.3 }   }
	SDPA_InputElement(Problem1, 3, 1, 1, 1, 6.3);
	SDPA_InputElement(Problem1, 3, 1, 1, 2, -7.5);
	SDPA_InputElement(Problem1, 3, 1, 2, 2, -3.3);

//	2nd block
//	{ { 0.2,  8.8,  5.4 },
//	  { 8.8,  3.4, -0.4 },
//	  { 5.4, -0.4,  7.5 }   }
	SDPA_InputElement(Problem1, 3, 2, 1, 1, 0.2);
	SDPA_InputElement(Problem1, 3, 2, 1, 2, 8.8);
	SDPA_InputElement(Problem1, 3, 2, 1, 3, 5.4);
	SDPA_InputElement(Problem1, 3, 2, 2, 2, 3.4);
	SDPA_InputElement(Problem1, 3, 2, 2, 3, -0.4);
	SDPA_InputElement(Problem1, 3, 2, 3, 3, 7.5);

//	3rd blcok
//	{-3.3, -4.0 }
	SDPA_InputElement(Problem1, 3, 3, 1, 1, -3.3);
	SDPA_InputElement(Problem1, 3, 3, 2, 2, -4.0);


//	Input F_4

//	1st block
//	{ { -2.4, -2.5 },
//	  { -2.5, -2.9 }   }
	SDPA_InputElement(Problem1, 4, 1, 1, 1, -2.4);
	SDPA_InputElement(Problem1, 4, 1, 1, 2, -2.5);
	SDPA_InputElement(Problem1, 4, 1, 2, 2, -2.9);

//	2nd block
//	{ {  3.4, -3.2, -4.5 },
//	  { -3.2,  3.0, -4.8 },
//	  { -4.5, -4.8,  3.6 }   }
	SDPA_InputElement(Problem1, 4, 2, 1, 1, 3.4);
	SDPA_InputElement(Problem1, 4, 2, 1, 2, -3.2);
	SDPA_InputElement(Problem1, 4, 2, 1, 3, -4.5);
	SDPA_InputElement(Problem1, 4, 2, 2, 2, 3.0);
	SDPA_InputElement(Problem1, 4, 2, 2, 3, -4.8);
	SDPA_InputElement(Problem1, 4, 2, 3, 3, 3.6);

//	3rd block
//	{  4.8 , 9.7 }
	SDPA_InputElement(Problem1, 4, 3, 1, 1, 4.8);
	SDPA_InputElement(Problem1, 4, 3, 2, 2, 9.7);


//	Input F_5

//	1st block
//	{ { -6.5, -5.4 },
//	  { -5.4, -6.6 }   }
	SDPA_InputElement(Problem1, 5, 1, 1, 1, -6.5);
	SDPA_InputElement(Problem1, 5, 1, 1, 2, -5.4);
	SDPA_InputElement(Problem1, 5, 1, 2, 2, -6.6);

//	2nd block
//	{ {  6.7, -7.2, -3.6 },
//	  { -7.2,  7.3, -3.0 },
//	  { -3.6, -3.0, -1.4 }   }
	SDPA_InputElement(Problem1, 5, 2, 1, 1, 6.7);
	SDPA_InputElement(Problem1, 5, 2, 1, 2, -7.2);
	SDPA_InputElement(Problem1, 5, 2, 1, 3, -3.6);
	SDPA_InputElement(Problem1, 5, 2, 2, 2, 7.3);
	SDPA_InputElement(Problem1, 5, 2, 2, 3, -3.0);
	SDPA_InputElement(Problem1, 5, 2, 3, 3, -1.4);

//	3rd block
//	{  6.1, -1.5 }
	SDPA_InputElement(Problem1, 5, 3, 1, 1, 6.1);
	SDPA_InputElement(Problem1, 5, 3, 2, 2, -1.5);


	if (false == SDPA_Check_sfMAT(Problem1))
		exit(1);	  


	SDPA_Solve(Problem1);


	fprintf(stdout, "\nStop iteration = %d\n", Problem1.Iteration);
       	fprintf(stdout, "objValPrimal   = %10.6e\n", Problem1.PrimalObj);
       	fprintf(stdout, "objValDual     = %10.6e\n", Problem1.DualObj);
       	fprintf(stdout, "p. feas. error = %10.6e\n", Problem1.PrimalError);
       	fprintf(stdout, "d. feas. error = %10.6e\n\n", Problem1.DualError);	
	
	fprintf(stdout, "xVec = \n");
	// Problem1.printResultXVec(stdout);
	element = Problem1.getResultXVec();
	PrintVector(Problem1.mDIM,element,stdout);
	
	fprintf(stdout, "xMat = \n");
	// Problem1.printResultXMat(stdout);
	fprintf(stdout, "{\n");
	for (int l=0; l<Problem1.nBLOCK; ++l) {
	  element = Problem1.getResultXMat(l);
	  int nRow = Problem1.bLOCKsTRUCT[l];
	  int nCol = nRow;
	  if (nRow<0) {
	    nCol = -nCol;
	    nRow = 1;
	  }
	  PrintMatrix(nRow,nCol,element,stdout);
	}
	fprintf(stdout, "}\n");
    
	fprintf(stdout, "yMat = \n");
	// Problem1.printResultYMat(stdout);
	fprintf(stdout, "{\n");
	for (int l=0; l<Problem1.nBLOCK; ++l) {
	  element = Problem1.getResultYMat(l);
	  int nRow = Problem1.bLOCKsTRUCT[l];
	  int nCol = nRow;
	  if (nRow<0) {
	    nCol = -nCol;
	    nRow = 1;
	  }
	  PrintMatrix(nRow,nCol,element,stdout);
	}
	fprintf(stdout, "}\n");
	
	Problem1.Delete();

	return 0;
};	

#define FORMAT "%+8.16e"

void PrintMatrix(int nRow, int nCol, double* element, FILE* fpOut)
{
  fprintf(fpOut,"{");
  for (int i=1; i<=nRow-1; ++i) {
    if (i==1) {
      fprintf(fpOut," ");
    } else {
      fprintf(fpOut,"  ");
    }
    fprintf(fpOut,"{");
    for (int j=1; j<=nCol-1; ++j) {
      fprintf(fpOut, FORMAT",",element[(i-1)+nRow*(j-1)]);
    }
    fprintf(fpOut,FORMAT" },\n",element[(i-1)+nRow*(nCol-1)]);
  }
  if (nRow>1) {
    fprintf(fpOut,"  {");
  }
  for (int j=1; j<=nCol-1; ++j) {
    fprintf(fpOut,FORMAT",",element[(nRow-1)+nRow*(j-1)]);
  }
  fprintf(fpOut,FORMAT" }",element[(nRow-1)+nRow*(nCol-1)]);
  if (nRow>1) {
    fprintf(fpOut,"   }\n");
  } else {
    fprintf(fpOut,"\n");
  }
}

void PrintVector(int nDim, double* element, FILE* fpOut)
{
  fprintf(fpOut,"{");
  for (int j=1; j<=nDim-1; ++j) {
    fprintf(fpOut,FORMAT",",element[j-1]);
  }
  if (nDim>0) {
    fprintf(fpOut,FORMAT"}\n",element[nDim-1]);
  } else {
    fprintf(fpOut,"  }\n");
  }
}
