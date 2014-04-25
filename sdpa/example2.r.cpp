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

#include <cstdio>
#include <cstdlib>
#include "rsdpa_class.h"
#include "rsdpa_lib.h"
using namespace std;

void PrintMatrix(int nRow, int nCol, double* element, FILE* fpOut);
void PrintVector(int nDim, double* element, FILE* fpOut);

extern "C" int sdpa_EX2r ()
{
  rSdpaLib Problem1;

  Problem1.setDisplay(stdout);

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
  

  int           m  = 5;
  int      nBlock  = 3;
  int* blockStruct = new int[nBlock];
  blockStruct[0]   = 2;
  blockStruct[1]   = 3;
  blockStruct[2]   =-2;
  bool initialPoint = false;

  Problem1.initialize1(m,nBlock,blockStruct,initialPoint);

  for (int k=0; k<=m; ++k) {
    Problem1.countUpperTriangle(k,1,3);
    Problem1.countUpperTriangle(k,2,6);
    
    // You do not have to assign the number of
    // diagonal parts.
    // The next line is not be needed.
    Problem1.countUpperTriangle(k,3,2);
  }

  Problem1.initialize2();

//      Input c

//      {1.1, -10, 6.6 , 19 , 4.1}
  Problem1.inputCVec(1,1.1);
  Problem1.inputCVec(2,-10);
  Problem1.inputCVec(3,6.6);
  Problem1.inputCVec(4, 19);
  Problem1.inputCVec(5,4.1);


//      Input F_0

//      1st block
//      { { -1.4, -3.2 },
//        { -3.2,-28   }   }
  Problem1.inputElement( 0, 1, 1, 1,-1.4);
  Problem1.inputElement( 0, 1, 1, 2,-3.2);
  Problem1.inputElement( 0, 1, 2, 2, -28);

//      2nd block
//      { { 15,  -12,    2.1 },
//        {-12,   16,   -3.8 },
//        {  2.1, -3.8, 15   }   }
  Problem1.inputElement( 0, 2, 1, 1,   15);
  Problem1.inputElement( 0, 2, 1, 2,  -12);
  Problem1.inputElement( 0, 2, 1, 3,  2.1);
  Problem1.inputElement( 0, 2, 2, 2,   16);
  Problem1.inputElement( 0, 2, 2, 3, -3.8);
  Problem1.inputElement( 0, 2, 3, 3,   15);

//      3rd block
//      {  1.8, -4.0 }
  Problem1.inputElement( 0, 3, 1, 1,  1.8);
  Problem1.inputElement( 0, 3, 2, 2, -4.0);


//      Input F_1

//      1st block
//      { {  0.5,  5.2 },
//        {  5.2, -5.3 }   }
  Problem1.inputElement( 1, 1, 1, 1, 0.5);
  Problem1.inputElement( 1, 1, 1, 2, 5.2);
  Problem1.inputElement( 1, 1, 2, 2, -5.3);

  //      2nd block
  //      { {  7.8, -2.4,  6.0 },
  //        { -2.4,  4.2,  6.5 },
  //        {  6.0,  6.5,  2.1 }   }
  Problem1.inputElement( 1, 2, 1, 1, 7.8);
  Problem1.inputElement( 1, 2, 1, 2, -2.4);
  Problem1.inputElement( 1, 2, 1, 3, 6.0);
  Problem1.inputElement( 1, 2, 2, 2, 4.2);
  Problem1.inputElement( 1, 2, 2, 3, 6.5);
  Problem1.inputElement( 1, 2, 3, 3, 2.1);

  //      3rd block
  //      { -4.5, -3.5 }
  Problem1.inputElement( 1, 3, 1, 1, -4.5);
  Problem1.inputElement( 1, 3, 2, 2, -3.5);

  //      Input F_2

//      1st block
//      { { 1.7,  7.0 },
//        { 7.0, -9.3 }   }
  Problem1.inputElement( 2, 1, 1, 1, 1.7);
  Problem1.inputElement( 2, 1, 1, 2, 7.0);
  Problem1.inputElement( 2, 1, 2, 2, -9.3);

  //      2nd block
  //      { {-1.9, -0.9, -1.3 },
  //        {-0.9, -0.8, -2.1 },
  //        {-1.3, -2.1,  4.0 }   }
  Problem1.inputElement( 2, 2, 1, 1, -1.9);
  Problem1.inputElement( 2, 2, 1, 2, -0.9);
  Problem1.inputElement( 2, 2, 1, 3, -1.3);
  Problem1.inputElement( 2, 2, 2, 2, -0.8);
  Problem1.inputElement( 2, 2, 2, 3, -2.1);
  Problem1.inputElement( 2, 2, 3, 3, 4.0);

  //      3rd blcok
  //      {-0.2, -3.7 }
  Problem1.inputElement( 2, 3, 1, 1, -0.2);
  Problem1.inputElement( 2, 3, 2, 2, -3.7);

  //      Input F_3

//      1st block
//      { { 6.3, -7.5 },
//        {-7.5, -3.3 }   }
  Problem1.inputElement( 3, 1, 1, 1, 6.3);
  Problem1.inputElement( 3, 1, 1, 2, -7.5);
  Problem1.inputElement( 3, 1, 2, 2, -3.3);

  //      2nd block
  //      { { 0.2,  8.8,  5.4 },
  //        { 8.8,  3.4, -0.4 },
  //        { 5.4, -0.4,  7.5 }   }
  Problem1.inputElement( 3, 2, 1, 1, 0.2);
  Problem1.inputElement( 3, 2, 1, 2, 8.8);
  Problem1.inputElement( 3, 2, 1, 3, 5.4);
  Problem1.inputElement( 3, 2, 2, 2, 3.4);
  Problem1.inputElement( 3, 2, 2, 3, -0.4);
  Problem1.inputElement( 3, 2, 3, 3, 7.5);

  //      3rd blcok
  //      {-3.3, -4.0 }
  Problem1.inputElement( 3, 3, 1, 1, -3.3);
  Problem1.inputElement( 3, 3, 2, 2, -4.0);

  //      Input F_4

//      1st block
//      { { -2.4, -2.5 },
//        { -2.5, -2.9 }   }
  Problem1.inputElement( 4, 1, 1, 1, -2.4);
  Problem1.inputElement( 4, 1, 1, 2, -2.5);
  Problem1.inputElement( 4, 1, 2, 2, -2.9);

  //      2nd block
  //      { {  3.4, -3.2, -4.5 },
  //        { -3.2,  3.0, -4.8 },
  //        { -4.5, -4.8,  3.6 }   }
  Problem1.inputElement( 4, 2, 1, 1, 3.4);
  Problem1.inputElement( 4, 2, 1, 2, -3.2);
  Problem1.inputElement( 4, 2, 1, 3, -4.5);
  Problem1.inputElement( 4, 2, 2, 2, 3.0);
  Problem1.inputElement( 4, 2, 2, 3, -4.8);
  Problem1.inputElement( 4, 2, 3, 3, 3.6);

  //      3rd block
  //      {  4.8 , 9.7 }
  Problem1.inputElement( 4, 3, 1, 1, 4.8);
  Problem1.inputElement( 4, 3, 2, 2, 9.7);

  //      Input F_5

//      1st block
//      { { -6.5, -5.4 },
//        { -5.4, -6.6 }   }
  Problem1.inputElement( 5, 1, 1, 1, -6.5);
  Problem1.inputElement( 5, 1, 1, 2, -5.4);
  Problem1.inputElement( 5, 1, 2, 2, -6.6);

  //      2nd block
  //      { {  6.7, -7.2, -3.6 },
  //        { -7.2,  7.3, -3.0 },
  //        { -3.6, -3.0, -1.4 }   }
  Problem1.inputElement( 5, 2, 1, 1, 6.7);
  Problem1.inputElement( 5, 2, 1, 2, -7.2);
  Problem1.inputElement( 5, 2, 1, 3, -3.6);
  Problem1.inputElement( 5, 2, 2, 2, 7.3);
  Problem1.inputElement( 5, 2, 2, 3, -3.0);
  Problem1.inputElement( 5, 2, 3, 3, -1.4);

  //      3rd block
  //      {  6.1, -1.5 }
  Problem1.inputElement( 5, 3, 1, 1, 6.1);
  Problem1.inputElement( 5, 3, 2, 2, -1.5);

  Problem1.dumpData("dump2.dat-s");
  Problem1.dumpInit("Init2.ini-s");
  
  Problem1.solve();

  // Delete information which are not relevant to solutions.
  // If you have enough memory, you don't have to call
  // delete1();.
  Problem1.delete1();

  fprintf(stdout, "\nStop iteration = %d\n",
	  Problem1.getIteration());
  fprintf(stdout, "objValPrimal   = %10.6e\n",
	  Problem1.getPrimalObj());
  fprintf(stdout, "objValDual     = %10.6e\n",
	  Problem1.getDualObj());
  fprintf(stdout, "p. feas. error = %10.6e\n",
	  Problem1.getPrimalError());
  fprintf(stdout, "d. feas. error = %10.6e\n",
	  Problem1.getDualError());	
  fprintf(stdout, "mu             = %10.6e\n",
	  Problem1.getMu());	
  fprintf(stdout, "time           = %10.6e\n",
	  Problem1.getTime());	
  #if 1
  rSdpaLib::phaseType phaseType = Problem1.getPhaseValue();
  // If you use this value,
  // you will compare such as
  // 'phaseType == rSparseLib::pdOPT'
  #endif
  
  char phaseString[20];
  // Length 20 is sutisfied.
  Problem1.stringPhaseValue(phaseString);
  fprintf(stdout, "phase          = %s\n\n",
	  phaseString);

  double* element;
  fprintf(stdout, "xVec = \n");
  // Problem1.printResultYVec(stdout);
  element = Problem1.getResultXVec();
  PrintVector(m,element,stdout);

  fprintf(stdout, "xMat = \n");
  // Problem1.printResultXMat(stdout);
  fprintf(stdout, "{\n");
  for (int l=0; l<nBlock; ++l) {
    element = Problem1.getResultXMat(l);
    int nRow = blockStruct[l];
    int nCol = nRow;
    if (nRow<0) {
      nCol = -nCol;
      nRow = 1;
    }
    PrintMatrix(nRow,nCol,element,stdout);
  }
  fprintf(stdout, "}\n");
    
  
  fprintf(stdout, "YMat = \n");
  // Problem1.printResultYMat(stdout);
  fprintf(stdout, "{\n");
  for (int l=0; l<nBlock; ++l) {
    element = Problem1.getResultYMat(l);
    int nRow = blockStruct[l];
    int nCol = nRow;
    if (nRow<0) {
      nCol = -nCol;
      nRow = 1;
    }
    PrintMatrix(nRow,nCol,element,stdout);
  }
  fprintf(stdout, "}\n");

  Problem1.printTime(stdout);

  // You may have to delete 'blockStruct' by yourself.
  delete[] blockStruct;
  // Delete all informaion about Problem1.
  Problem1.delete2();
  return EXIT_SUCCESS;
}

#define FORMAT "%+8.16e"

//void PrintMatrix(int nRow, int nCol, double* element, FILE* fpOut)
//{
//  fprintf(fpOut,"{");
//  for (int i=0; i<nRow-1; ++i) {
//    if (i==0) {
//      fprintf(fpOut," ");
//    } else {
//      fprintf(fpOut,"  ");
//    }
//    fprintf(fpOut,"{");
//    for (int j=0; j<nCol-1; ++j) {
//      fprintf(fpOut, FORMAT",",element[i+nRow*j]);
//    }
//    fprintf(fpOut,FORMAT" },\n",element[i+nRow*(nCol-1)]);
//  }
//  if (nRow>1) {
//    fprintf(fpOut,"  {");
//  }
//  for (int j=0; j<nCol-1; ++j) {
//    fprintf(fpOut,FORMAT",",element[(nRow-1)+nRow*j]);
//  }
//  fprintf(fpOut,FORMAT" }",element[(nRow-1)+nRow*(nCol-1)]);
//  if (nRow>1) {
//    fprintf(fpOut,"   }\n");
//  } else {
//    fprintf(fpOut,"\n");
//  }
//}
//
//void PrintVector(int nDim, double* element, FILE* fpOut)
//{
//  fprintf(fpOut,"{");
//  for (int j=0; j<nDim-1; ++j) {
//    fprintf(fpOut,FORMAT",",element[j]);
//  }
//  if (nDim>0) {
//    fprintf(fpOut,FORMAT"}\n",element[nDim-1]);
//  } else {
//    fprintf(fpOut,"  }\n");
//  }
//}
