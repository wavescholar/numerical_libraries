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
  rsdpa_struct.h
-------------------------------------------------*/

#ifndef __rsdpa_struct_h__
#define __rsdpa_struct_h__

#include "rsdpa_include.h"

class rVector
{
public:
  int nDim;
  double* ele;

  rVector();
  rVector(int nDim, double value = 0.0);
  ~rVector();

  void initialize(int nDim, double value = 0.0);
  void initialize(double value);
  void setZero();
  
  void display(FILE* fpout = stdout);
  void display(FILE* fpout,double scalar);
  bool copyFrom(rVector& other);
};

class rBlockVector
{
public:
  int  nBlock;
  int* blockStruct;

  rVector* ele;
  
  rBlockVector();
  rBlockVector(int nBlock, int* blockStruct, double value = 0.0);
  ~rBlockVector();
  
  void initialize(int nBlock, int* blockStruct, double value = 0.0);
  void initialize(double value);
  void setZero();
  
  void display(FILE* fpout = stdout);
  bool copyFrom(rBlockVector& other);
};

class rSparseMatrix
{
public:
  int nRow, nCol;

  enum rSpMat_Sp_De_Di { SPARSE, DENSE ,DIAGONAL};
  rSpMat_Sp_De_Di Sp_De_Di;
  // flag of Sparse or Dense or Diagonal
  
  int NonZeroNumber;
  // for memory
  int NonZeroCount;
  // currentry stored
  int NonZeroEffect;
  // use for calculation of F1,F2,F3 

  // for Dense
  double* de_ele;

  // for Sparse
  int*    row_index;
  int*    column_index;
  double* sp_ele;

  // for Diagonal
  double* di_ele;

  rSparseMatrix();
  rSparseMatrix(int nRow,int nCol, rSpMat_Sp_De_Di Sp_De_Di,
		int NonZeroNumber);
  ~rSparseMatrix();

  void initialize(int nRow,int nCol, rSpMat_Sp_De_Di Sp_De_Di,
		  int NonZeroNumber);

  void display(FILE* fpout = stdout);
  bool copyFrom(rSparseMatrix& other);

  void changeToDense(bool forceChange = false);
  void setZero();
  void setIdentity(double scalar = 1.0);

  bool sortSparseIndex(int&i, int& j);
};

class rDenseMatrix
{
public:
  int nRow, nCol;

  enum rDeMat_De_Di { DENSE ,DIAGONAL};
  rDeMat_De_Di De_Di;
  // flag of Dense or Diagonal
  
  // for Dense
  double* de_ele;

  // for Diagonal
  double* di_ele;

  rDenseMatrix();
  rDenseMatrix(int nRow,int nCol, rDeMat_De_Di De_Di);
  ~rDenseMatrix();

  void initialize(int nRow,int nCol, rDeMat_De_Di De_Di);
  
  void display(FILE* fpout = stdout);
  bool copyFrom(rDenseMatrix& other);
  bool copyFrom(rSparseMatrix& other);

  void setZero();
  void setIdentity(double scalar = 1.0);
};

class rBlockSparseMatrix
{
public:
  int  nBlock;
  int* blockStruct;

  rSparseMatrix* ele;
  
  rBlockSparseMatrix();
  rBlockSparseMatrix(int nBlock,int* blockStruct);
  ~rBlockSparseMatrix();

  void initialize(int nBlock,int* blockStruct);
  
  void display(FILE* fpout = stdout);
  bool copyFrom(rBlockSparseMatrix& other);
  
  void setZero();
  void setIdentity(double scalar = 1.0);
  void changeToDense(bool forceChange=false);
  bool sortSparseIndex(int&l , int& i, int& j);
};

class rBlockDenseMatrix
{
public:
  int  nBlock;
  int* blockStruct;

  rDenseMatrix* ele;

  rBlockDenseMatrix();
  rBlockDenseMatrix(int nBlock,int* blockStruct);
  ~rBlockDenseMatrix();

  void initialize(int nBlock,int* blockStruct);
  
  void display(FILE* fpout = stdout);
  bool copyFrom(rBlockSparseMatrix& other);
  bool copyFrom(rBlockDenseMatrix& other);
  
  void setZero();
  void setIdentity(double scalar = 1.0);
};

#endif // __rsdpa_struct_h__
