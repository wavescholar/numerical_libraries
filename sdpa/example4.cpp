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
#include <stdio.h>
#include <stdlib.h>

#include "sdpa-lib.hpp"
#include "sdpa-lib2.hpp"

/*
example1.dat:

"Example 1: mDim = 3, nBLOCK = 1, {2}"
   3  =  mDIM
   1  =  nBOLCK
   2  = bLOCKsTRUCT
{48, -8, 20}
{ {-11,  0}, { 0, 23} }
{ { 10,  4}, { 4,  0} }
{ {  0,  0}, { 0, -8} }
{ {  0, -8}, {-8, -2} }
*/

/*
example1.ini:

{0.0, -4.0, 0.0}
{ {11.0, 0.0}, {0.0, 9.0} }
{ {5.9,  -1.375}, {-1.375, 1.0} }
*/

void PrintMatrix(int nRow, int nCol, double* element, FILE* fpOut);
void PrintVector(int nDim, double* element, FILE* fpOut);

extern "C" int sdpa_EX4 ()
{
	double* 	element;

	SDPA	Problem1;

        Problem1.Method          	= KSH;
	Problem1.InitialPoint	 	= true;
	
	Problem1.pARAM.maxIteration 	= 50;
        Problem1.pARAM.epsilonStar   	= 1.0E-6;
        Problem1.pARAM.lambdaStar    	= 1.0E2;
        Problem1.pARAM.omegaStar     	= 2.0;
        Problem1.pARAM.lowerBound    	= -1.0E5;
        Problem1.pARAM.upperBound    	= 1.0E5;
        Problem1.pARAM.betaStar       	= 0.1;
        Problem1.pARAM.betaBar       	= 0.2;
        Problem1.pARAM.gammaStar     	= 0.9;

	Problem1.DisplayInformation 	= stdout;
	
	SDPA_initialize(Problem1);

	Problem1.mDIM	= 3;
	Problem1.nBLOCK	= 1;
	Problem1.bLOCKsTRUCT	= new int [Problem1.nBLOCK];
	Problem1.bLOCKsTRUCT[0] = 2;
	
	SDPA_initialize2(Problem1);

	SDPA_Input_cVECT(Problem1, 1, 48);
	SDPA_Input_cVECT(Problem1, 2, -8);
	SDPA_Input_cVECT(Problem1, 3, 20);

	SDPA_CountUpperTriangle(Problem1, 0, 1, 2);
	SDPA_CountUpperTriangle(Problem1, 1, 1, 2);
	SDPA_CountUpperTriangle(Problem1, 2, 1, 1);
	SDPA_CountUpperTriangle(Problem1, 3, 1, 2);


	SDPA_Make_sfMAT(Problem1);

	
	SDPA_InputElement(Problem1, 0, 1, 1, 1, -11);
	SDPA_InputElement(Problem1, 0, 1, 2, 2, 23);
	
	SDPA_InputElement(Problem1, 1, 1, 1, 1, 10);
	SDPA_InputElement(Problem1, 1, 1, 1, 2, 4);

	SDPA_InputElement(Problem1, 2, 1, 2, 2, -8);

	SDPA_InputElement(Problem1, 3, 1, 1, 2, -8);
	SDPA_InputElement(Problem1, 3, 1, 2, 2, -2);


	SDPA_Input_IniXMat(Problem1, 1, 1, 1, 11);
	SDPA_Input_IniXMat(Problem1, 1, 2, 2, 9);

	SDPA_Input_InixVec(Problem1, 2, -4);

	SDPA_Input_IniYMat(Problem1, 1, 1, 1, 5.9);
	SDPA_Input_IniYMat(Problem1, 1, 1, 2, -1.375);
	SDPA_Input_IniYMat(Problem1, 1, 2, 2, 1);


	SDPA_Solve(Problem1);


	fprintf(stdout, "\nStop iteration = %d\n", Problem1.Iteration);
       	fprintf(stdout, "objValPrimal   = %10.6e\n", Problem1.PrimalObj);
       	fprintf(stdout, "objValDual     = %10.6e\n", Problem1.DualObj);
       	fprintf(stdout, "p. feas. error = %10.6e\n", Problem1.PrimalError);
       	fprintf(stdout, "d. feas. error = %10.6e\n\n", Problem1.DualError);	


        Problem1.pARAM.epsilonStar   	= 1.0E-8;

	SDPA_Copy_Current_To_Ini(Problem1);

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

//void PrintMatrix(int nRow, int nCol, double* element, FILE* fpOut)
//{
//  fprintf(fpOut,"{");
//  for (int i=1; i<=nRow-1; ++i) {
//    if (i==1) {
//      fprintf(fpOut," ");
//    } else {
//      fprintf(fpOut,"  ");
//    }
//    fprintf(fpOut,"{");
//    for (int j=1; j<=nCol-1; ++j) {
//      fprintf(fpOut, FORMAT",",element[(i-1)+nRow*(j-1)]);
//    }
//    fprintf(fpOut,FORMAT" },\n",element[(i-1)+nRow*(nCol-1)]);
//  }
//  if (nRow>1) {
//    fprintf(fpOut,"  {");
//  }
//  for (int j=1; j<=nCol-1; ++j) {
//    fprintf(fpOut,FORMAT",",element[(nRow-1)+nRow*(j-1)]);
//  }
//  fprintf(fpOut,FORMAT" }",element[(nRow-1)+nRow*(nCol-1)]);
//  if (nRow>1) {
//    fprintf(fpOut,"   }\n");
//  } else {
//    fprintf(fpOut,"\n");
//  }
//}

//void PrintVector(int nDim, double* element, FILE* fpOut)
//{
//  fprintf(fpOut,"{");
//  for (int j=1; j<=nDim-1; ++j) {
//    fprintf(fpOut,FORMAT",",element[j-1]);
//  }
//  if (nDim>0) {
//    fprintf(fpOut,FORMAT"}\n",element[nDim-1]);
//  } else {
//    fprintf(fpOut,"  }\n");
//  }
//}
