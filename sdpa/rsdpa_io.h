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
/*-------------------------------------------------
  rsdpa_io.h
-------------------------------------------------*/

#ifndef __rsdpa_io_h__
#define __rsdpa_io_h__

#include "rsdpa_parts.h"

#define lengthOfString 256

class rIO
{
public:
  static void read(FILE* fpData, FILE* fpout, int& m, char* str);
  static void read(FILE* fpData, int& nBlock);
  static void read(FILE* fpData, int nBlock, int* blockStruct);
  static void read(FILE* fpData, rVector& b);
  static void read(FILE* fpData, rBlockDenseMatrix& xMat,
		   rVector& yVec,
		   rBlockDenseMatrix& zMat, int nBlock,
		   int* blockStruct, int inputSparse);
  static void read(FILE* fpData, rBlockSparseMatrix& C,
		   rBlockSparseMatrix* A, int m, int nBlock,
		   int* blockStruct);
  static void read(FILE* fpData, int m, int nBlock,
		   int* blockStruct, int* CNonZeroCount,
		   int* ANonZeroElement,bool isDataSparse);
  static void read(FILE* fpData, rBlockSparseMatrix& C,
		   rBlockSparseMatrix* A,int m, int nBlock,
		   int* blockStruct, long position,
		   bool isDataSparse);

  static void printHeader(FILE* fpout, FILE* Display);

  static void printOneIteration(int pIteration,
				rAverageComplementarity& mu,
				rRatioInitResCurrentRes& theta,
				rSolveInfo& solveInfo,
				rStepLength& alpha,
				rDirectionParameter& beta,
				rResiduals& currentRes,
				FILE* fpout,
				FILE* Display);
  static void printLastInfo(int pIteration,
			    rAverageComplementarity& mu,
			    rRatioInitResCurrentRes& theta,
			    rSolveInfo& solveInfo,
			    rStepLength& alpha,
			    rDirectionParameter& beta,
			    rResiduals& currentRes,
			    rPhase & phase,
			    rSolutions& currentPt,
			    double cputime,
			    int nDim,
			    rVector& b,
			    rBlockSparseMatrix& C,
			    rBlockSparseMatrix* A,
			    rComputeTime& com,
			    rParameter& param,
			    FILE* fpout,
			    FILE* Display,
			    bool printTime = true);
};
#endif // __rsdpa_io_h__
