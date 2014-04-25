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
// Here is the start of ``example1.r.cpp''

/*
example1.dat:

"Example 1: mDim = 3, nBLOCK = 1, {2}"
   3  =  mDIM
   1  =  nBOLCK
   2  = bLOCKsTRUCT
{48, -8, 20}
{ {-11,  0},{ 0,  23} } 
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

#include "rsdpa_class.h"
#include "rsdpa_lib.h"

#include <cstdio>
using namespace std;

void PrintMatrix(int nRow, int nCol, double* element, FILE* fpOut);
void PrintVector(int nDim, double* element, FILE* fpOut);

int main ()
{
  rSdpaLib Problem1;

  FILE* Display = stdout;
  // FILE* Display = NULL; // When you may not want to display
  
  Problem1.setDisplay(Display);

  Problem1.setDefaultParameter();
  Problem1.setParameterMaxIteration(50);
  Problem1.setParameterEpsilonStar (1.0e-7);
  Problem1.setParameterLambdaStar  (1.0e+2);
  Problem1.setParameterOmegaStar   (2.0);
  Problem1.setParameterLowerBound  (-1.0e+5);
  Problem1.setParameterUpperBound  (1.0e+5);
  Problem1.setParameterBetaStar    (0.1);
  Problem1.setParameterBetaBar     (0.2);
  Problem1.setParameterGammaStar   (0.9);
  Problem1.setParameterEpsilonDash (1.0e-7);

  int           m  = 3;
  int      nBlock  = 1;
  int* blockStruct = new int[nBlock];
  blockStruct[0]   = 2;
  bool initialPoint = true;

  Problem1.initialize1(m,nBlock,blockStruct,initialPoint);

  Problem1.countUpperTriangle(0,1,2);
  Problem1.countUpperTriangle(1,1,2);
  Problem1.countUpperTriangle(2,1,1);
  Problem1.countUpperTriangle(3,1,2);

  Problem1.initialize2();

  // c = {48, -8, 20}
  Problem1.inputCVec(1,48);
  Problem1.inputCVec(2,-8);
  Problem1.inputCVec(3,20);

  // F_0 = { { -11, 0 }, { 0 , 23} }
  Problem1.inputElement(0,1,1,1,-11);
  Problem1.inputElement(0,1,2,2, 23);

  // F_1 = { { 10, 4 }, { 4 , 0} }
  Problem1.inputElement(1,1,1,1, 10); 
  Problem1.inputElement(1,1,1,2,  4);

  // F_2 = { {  0, 0 }, { 0 , -8} }
  Problem1.inputElement(2,1,2,2, -8);

  // F_3 = { {  0, -8 }, { -8 , -2} }
  Problem1.inputElement(3,1,1,2, -8);
  Problem1.inputElement(3,1,2,2, -2);

  // x0 = { 0, -4, 0}
  Problem1.inputInitXVec(2,-4);
  
  // X0 = { {  11, 0 }, { 0 , 9} }
  Problem1.inputInitXMat(1,1,1,    11);
  Problem1.inputInitXMat(1,2,2,     9);
  
  // Y0 = { { 5.9, -1.375 }, { -1.375 , 1} }
  Problem1.inputInitYMat(1,1,1,   5.9);
  Problem1.inputInitYMat(1,1,2,-1.375);
  Problem1.inputInitYMat(1,2,2,     1);


  Problem1.dumpData("example1.r.produce.dat-s");
  Problem1.dumpInit("example1.r.produce.ini-s");
  
  Problem1.solve();

  // delete information which are not relevant to solutions
  Problem1.delete1();

  if (Display) {
    fprintf(Display, "\nStop iteration = %d\n",
	    Problem1.getIteration());
    fprintf(Display, "time           = %10.6e\n",
	    Problem1.getTime());	
    fprintf(Display, "objValPrimal   = %10.6e\n",
	    Problem1.getPrimalObj());
    fprintf(Display, "objValDual     = %10.6e\n",
	    Problem1.getDualObj());
    fprintf(Display, "p. feas. error = %10.6e\n",
	    Problem1.getPrimalError());
    fprintf(Display, "d. feas. error = %10.6e\n",
	    Problem1.getDualError());	

    rSdpaLib::phaseType phaseType = Problem1.getPhaseValue();
    // If you use this value, you will compare such as
    // 'phaseType == rSparseLib::pdOPT'
  
    char phaseString[20];
    // if length is 20, it is enough.
    Problem1.stringPhaseValue(phaseString);
    fprintf(Display, "phase          = %s\n\n",phaseString);
  }
  
  double* element;
  if (Display) {
    fprintf(Display, "xVec = \n");
    // Problem1.printResultXVec(Display);
    element = Problem1.getResultXVec();
    PrintVector(m,element,Display);

    fprintf(Display, "xMat = \n");
    // Problem1.printResultXMat(Display);
    fprintf(Display, "{\n");
    for (int l=0; l<nBlock; ++l) {
      element = Problem1.getResultXMat(l);
      int nRow = blockStruct[l];
      int nCol = nRow;
      if (nRow<0) {
        nCol = -nCol;
        nRow = 1;
      }
      PrintMatrix(nRow,nCol,element,Display);
    }
    fprintf(Display, "}\n");
    
    
    fprintf(Display, "yMat = \n");
    // Problem1.printResultYMat(Display);
    fprintf(Display, "{\n");
    for (int l=0; l<nBlock; ++l) {
      element = Problem1.getResultYMat(l);
      int nRow = blockStruct[l];
      int nCol = nRow;
      if (nRow<0) {
        nCol = -nCol;
        nRow = 1;
      }
      PrintMatrix(nRow,nCol,element,Display);
    }
    fprintf(Display, "}\n");
  } // end of 'if (Display)'

  // You have to delete 'blockStruct' by yourself.
  delete[] blockStruct;
  // delete information abount solution
  Problem1.delete2();
  return EXIT_SUCCESS;
}

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
