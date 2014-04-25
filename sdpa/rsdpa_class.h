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
/*-----------------------------------------
   rsdpa_class.h                            
 -----------------------------------------*/
 #ifndef __rsdpa_class_h__                  
 #define __rsdpa_class_h__                  

/*--------------------------------------------------
  rsdpa_right.h
--------------------------------------------------*/

#ifndef __rsdpa_right_h__
#define __rsdpa_right_h__

/*---------------------------------------------------
  the first character is r, that means 'rosemary'.
---------------------------------------------------*/

static const char rsdpa_right[] =
  "'rsdpa' has been written by Makoto Yamashita since 2002.03.28";


#endif // __rsdpa_right_h__
/*--------------------------------------------------

  rsdpa_tool.h

--------------------------------------------------*/



#ifndef __rsdpa_tool_h__

#define __rsdpa_tool_h__





#include <iostream>

#include <time.h>

#include <sys/types.h>

#include <sys/timeb.h>

// #include <sys/time.h>

#include <string>



#if 1

#define rMessage(message) cout << message << " :: line " << __LINE__  << " in " << __FILE__ << endl;

#else

#define rMessage(message)

#endif



#define rError(message) cout << message << " :: line " << __LINE__ << " in " << __FILE__ << endl;  exit(false);



#if 0

#define rNewCheck() rMessage("new invoked");

#else

#define rNewCheck() ;

#endif



#define REVERSE_PRIMAL_DUAL 1





// These are constant. Do NOT change

extern int IZERO   ; // =  0;

extern int IONE    ; // =  1;

extern int IMONE   ; // = -1;

extern double DZERO; // =  0.0;

extern double DONE ; // =  1.0;

extern double DMONE; // = -1.0;



struct rrealtime {

  time_t ltime;

  _timeb tstruct;

};



class rTime

{

public:

  static double rGetUserTime();

  static void rSetTimeVal(rrealtime& targetVal);

  static double rGetRealTime(rrealtime& start,

			     rrealtime& end);

};



#if 1 // count time with process time

#define rTimeStart(START__)  static clock_t START__; START__ = clock();

#define rTimeEnd(END__)   static clock_t END__;   END__ = clock();

#define rTimeCal(START__,END__) (((double)(END__ - START__)) / CLOCKS_PER_SEC);



#else // count time with real time

#define rTimeStart(START__) \

   static rrealtime START__; rTime::rSetTimeVal(START__)

#define rTimeEnd(END__) \

   static rrealtime END__; rTime::rSetTimeVal(END__)

#define rTimeCal(START__,END__) rTime::rGetRealTime(START__,END__)

#endif



#endif // __rsdpa_tool_h__

/*--------------------------------------------------
  rsdpa_include.h
--------------------------------------------------*/

#ifndef __rsdpa_include_h__
#define __rsdpa_include_h__


// if you use ATLAS, you need to set 0
// otherwise (for example, BLAS in clapack.tgz), set 1
//           and edit Makefile to change LAPACK_LIB

#define NON_ATLAS_SDPA 1

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
	#include "mkl_lapack.h"
#include "mkl_blas.h"
extern "C" {
#include "f2c.h"

#if NON_ATLAS_SDPA//bbcrevisit edited in port to Intel BLAS
//#include "blaswrap.h"
#endif
//#include "fblaswr.h"
//#include "cblas.h"

//#include "clapack.h"
};





using namespace std;

#define _SUCCESS true
#define FAILURE false


#endif // __rsdpa_include_h__
/*-------------------------------------------------
  rsdpa_struct.h
-------------------------------------------------*/

#ifndef __rsdpa_struct_h__
#define __rsdpa_struct_h__


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
/*-----------------------------------------------
  rsdpa_dpotrf.cpp
  modification of ATL_dpotrfL
  int rATL_dpotrfL(int N, double *A,int lda)
-----------------------------------------------*/

#ifndef __rsdpa_dpotrf_h__
#define __rsdpa_dpotrf_h__

#ifdef __cplusplus
extern "C" int rATL_dpotrfL(int N, double *A,int lda);
#else
extern int rATL_dpotrfL(int N, double *A,int lda);
#endif

#endif // __rsdpa_dpotrf_h__


/*-------------------------------------------------
  rsdpa_algebra.h
-------------------------------------------------*/

#ifndef __rsdpa_algebra_h__
#define __rsdpa_algebra_h__


class rAl
{
public:
  static double getMinEigenValue(rDenseMatrix& aMat,
				 rVector& eigenVec,
				 rVector& workVec);
  static double getMinEigenValue(rBlockDenseMatrix& aMat,
				 rBlockVector& eigenVec,
				 rBlockVector& workVec);

  static bool getInnerProduct(double& ret,
			      rVector& aVec, rVector& bVec);
  static bool getInnerProduct(double& ret,
			      rBlockVector& aVec,
			      rBlockVector& bVec);
  static bool getInnerProduct(double& ret,
			      rDenseMatrix& aMat,
			      rDenseMatrix& bMat);
  static bool getInnerProduct(double& ret,
			      rSparseMatrix& aMat,
			      rDenseMatrix&  bMat);
  static bool getInnerProduct(double& ret,
			      rBlockDenseMatrix& aMat,
			      rBlockDenseMatrix& bMat);
  static bool getInnerProduct(double& ret,
			      rBlockSparseMatrix& aMat,
			      rBlockDenseMatrix&  bMat);

  static bool getCholesky(rDenseMatrix& retMat, rDenseMatrix& aMat);

  static bool getInvLowTriangularMatrix(rDenseMatrix& retMat,
					rDenseMatrix& aMat);
  static bool getCholeskyAndInv(rBlockDenseMatrix& choleskyMat,
				rBlockDenseMatrix& inverseMat,
				rBlockDenseMatrix& aMat);

  static bool getSymmetrize(rDenseMatrix& aMat);
  static bool getSymmetrize(rBlockDenseMatrix& aMat);

  static bool getTranspose(rDenseMatrix& retMat,
			   rDenseMatrix& aMat);
  static bool getTranspose(rBlockDenseMatrix& retMat,
			   rBlockDenseMatrix& aMat);
  
  static int rdpotf2_(char*uplo, int *n, double *a, int *lda, int *info);
  static int rdpotrf_(char *uplo, int *n, double *a, int *lda, int *info);
  static bool choleskyFactorWithAdjust(rDenseMatrix& aMat);
  
  static bool solveSystems(rVector& xVec,
			   rDenseMatrix& aMat, rVector& bVec);
  // solve aMat * xVec = bVec
  // aMat must be Cholesky Factorized.

  static bool multiply(rDenseMatrix& retMat,
		       rDenseMatrix& aMat, rDenseMatrix& bMat,
		       double* scalar = NULL);
  static bool multiply(rDenseMatrix& retMat,
		       rSparseMatrix& aMat, rDenseMatrix& bMat,
		       double* scalar = NULL);
  static bool multiply(rDenseMatrix& retMat,
		       rDenseMatrix& aMat, rSparseMatrix& bMat,
		       double* scalar = NULL);
  static bool multiply(rDenseMatrix& retMat,
		       rDenseMatrix& aMat, double* scalar = NULL);
  static bool multiply(rBlockDenseMatrix& retMat,
		       rBlockDenseMatrix& aMat,
		       double* scalar = NULL);
  static bool multiply(rVector& retVec,
		       rVector& aVec, double* scalar = NULL);
  static bool multiply(rBlockVector& retVec,
		       rBlockVector& aVec,
		       double* scalar = NULL);
  static bool multiply(rVector& retVec,
		       rDenseMatrix& aMat, rVector& bVec,
		       double* scalar = NULL);
  static bool multiply(rBlockVector& retVec,
		       rBlockDenseMatrix& aMat,
		       rBlockVector& bVec,
		       double* scalar = NULL);
  static bool multiply(rBlockDenseMatrix& retMat,
		       rBlockDenseMatrix& aMat,
		       rBlockDenseMatrix& bMat,
		       double* scalar = NULL);
  static bool multiply(rBlockDenseMatrix& retMat,
		       rBlockSparseMatrix& aMat,
		       rBlockDenseMatrix& bMat,
		       double* scalar = NULL);
  static bool multiply(rBlockDenseMatrix& retMat,
		       rBlockDenseMatrix& aMat,
		       rBlockSparseMatrix& bMat,
		       double* scalar = NULL);
  // ret = aMat**T * bMat
  static bool tran_multiply(rDenseMatrix& retMat,
			    rDenseMatrix& aMat, rDenseMatrix& bMat,
			    double* scalar = NULL);
  static bool tran_multiply(rBlockDenseMatrix& retMat,
			    rBlockDenseMatrix& aMat,
			    rBlockDenseMatrix& bMat,
			    double* scalar = NULL);
  // ret = aMat * bMat**T
  static bool multiply_tran(rDenseMatrix& retMat,
			    rDenseMatrix& aMat, rDenseMatrix& bMat,
			    double* scalar = NULL);
  static bool multiply_tran(rBlockDenseMatrix& retMat,
			    rBlockDenseMatrix& aMat,
			    rBlockDenseMatrix& bMat,
			    double* scalar = NULL);
  // ret = a + (*scalar)*b
  static bool plus(rVector& retVec, rVector& aVec,
		   rVector& bVec, double* scalar = NULL);
  static bool plus(rDenseMatrix& retMat,
		   rDenseMatrix& aMat, rDenseMatrix& bMat,
		   double* scalar = NULL);
  static bool plus(rDenseMatrix& retMat,
		   rSparseMatrix& aMat, rDenseMatrix& bMat,
		   double* scalar = NULL);
  static bool plus(rDenseMatrix& retMat,
		   rDenseMatrix& aMat, rSparseMatrix& bMat,
		   double* scalar = NULL);
  
  static bool plus(rBlockVector& retVec,
		   rBlockVector& aVec,
		   rBlockVector& bVec, double* scalar = NULL);
  static bool plus(rBlockDenseMatrix& retMat,
		   rBlockDenseMatrix& aMat,
		   rBlockDenseMatrix& bMat,
		   double* scalar = NULL);
  static bool plus(rBlockDenseMatrix& retMat,
		   rBlockSparseMatrix& aMat,
		   rBlockDenseMatrix& bMat,
		   double* scalar = NULL);
  static bool plus(rBlockDenseMatrix& retMat,
		   rBlockDenseMatrix& aMat,
		   rBlockSparseMatrix& bMat,
		   double* scalar = NULL);

  // ret = a '*' (*scalar)
  static bool let(rVector& retVec, const char eq,
		  rVector& aVec, const char op,
		  double* scalar = NULL);

  // ret = a '*' (*scalar)
  static bool let(rBlockVector& retVec, const char eq,
		  rBlockVector& aVec, const char op,
		  double* scalar = NULL);

  // ret = a '*' (*scalar)
  static bool let(rDenseMatrix& retMat, const char eq,
		  rDenseMatrix& aMat, const char op,
		  double* scalar = NULL);

  // ret = a '*' (*scalar)
  static bool let(rBlockDenseMatrix& retMat, const char eq,
		  rBlockDenseMatrix& aMat, const char op,
		  double* scalar = NULL);

  // ret = a '+' '-' b*(*scalar)
  static bool let(rVector& retVec, const char eq,
		  rVector& aVec, const char op,
		  rVector& bVec, double* scalar = NULL);

  // ret = a '+' '-' '*' 't' 'T' b*(*scalar)
  static bool let(rDenseMatrix& retMat, const char eq,
		  rDenseMatrix& aMat, const char op,
		  rDenseMatrix& bMat, double* scalar = NULL);

  // ret = a '+' '-' '*' b*(*scalar)
  static bool let(rDenseMatrix& retMat, const char eq,
		  rSparseMatrix& aMat, const char op,
		  rDenseMatrix& bMat, double* scalar = NULL);

  // ret = a '+' '-' '*' b*(*scalar)
  static bool let(rDenseMatrix& retMat, const char eq,
		  rDenseMatrix& aMat, const char op,
		  rSparseMatrix& bMat, double* scalar = NULL);

  // ret = a '+' '-' '*' 't' 'T' b*(*scalar)
  static bool let(rBlockDenseMatrix& retMat, const char eq,
		  rBlockDenseMatrix& aMat, const char op,
		  rBlockDenseMatrix& bMat, double* scalar = NULL);

  // ret = a '+' '-' '*' b*(*scalar)
  static bool let(rBlockDenseMatrix& retMat, const char eq,
		  rBlockSparseMatrix& aMat, const char op,
		  rBlockDenseMatrix& bMat, double* scalar = NULL);

  // ret = a '+' '-' '*' b*(*scalar)
  static bool let(rBlockDenseMatrix& retMat, const char eq,
		  rBlockDenseMatrix& aMat, const char op,
		  rBlockSparseMatrix& bMat, double* scalar = NULL);


  // ret = aMat '*' '/' bVec
  static bool let(rVector& rVec, const char eq,
		  rDenseMatrix& aMat, const char op,
		  rVector& bVec);

  // ret = inner_product(a,b) // op = '.'
  static bool let(double& ret, const char eq,
		  rVector& aVec, const char op,
		  rVector& bVec);
  
  // ret = inner_product(a,b) // op = '.'
  static bool let(double& ret, const char eq,
		  rDenseMatrix& aMat, const char op,
		  rDenseMatrix& bMat);
  
  // ret = inner_product(a,b) // op = '.'
  static bool let(double& ret, const char eq,
		  rDenseMatrix& aMat, const char op,
		  rSparseMatrix& bMat);
  
  // ret = inner_product(a,b) // op = '.'
  static bool let(double& ret, const char eq,
		  rSparseMatrix& aMat, const char op,
		  rDenseMatrix& bMat);

  // ret = inner_product(a,b) // op = '.'
  static bool let(double& ret, const char eq,
		  rBlockVector& aVec, const char op,
		  rBlockVector& bVec);
  
  // ret = inner_product(a,b) // op = '.'
  static bool let(double& ret, const char eq,
		  rBlockDenseMatrix& aMat, const char op,
		  rBlockDenseMatrix& bMat);

  // ret = inner_product(a,b) // op = '.'
  static bool let(double& ret, const char eq,
		  rBlockSparseMatrix& aMat, const char op,
		  rBlockDenseMatrix& bMat);

  // ret = inner_product(a,b) // op = '.'
  static bool let(double& ret, const char eq,
		  rBlockDenseMatrix& aMat, const char op,
		  rBlockSparseMatrix& bMat);

};

#endif // __rsdpa_algebra_h__
/*-------------------------------------------------
  rsdpa_parts.h
-------------------------------------------------*/

#ifndef __rsdpa_parts_h__
#define __rsdpa_parts_h__


class rComputeTime;
class rParameter;
class rSolutions;
class rNewton;
class rResiduals;
class rStepLength;
class rDirectionParameter;
class rSwitch;

class rRatioInitResCurrentRes;
class rSolveInfo;
class rPhase;
class rAverageComplementarity;
class rLanczos;

class rComputeTime
{
public:
  double Predictor;
  double Corrector;
  double StepPredictor;
  double StepCorrector;
  double xMatTime;
  double zMatTime;
  double invzMatTime;
  double xMatzMatTime;
  double EigxMatTime;
  double EigzMatTime;
  double EigxMatzMatTime;
  double makerMat;
  double makebMat;
  double B_DIAG;
  double B_F1;
  double B_F2;
  double B_F3;
  double B_PRE;
  double makegVecMul;
  double makegVec;
  double choleskybMat;
  double solve;
  double sumDz;
  double makedX;
  double symmetriseDx;
  double makedXdZ;
  double updateRes;
  double MainLoop;
  double FileRead;
  double FileCheck;
  double FileChange;
  double TotalTime;
  rComputeTime();
  ~rComputeTime();
  void display(FILE* fpout=stdout);
};

class rParameter
{
public:
  enum parameterType {PARAMETER_DEFAULT,
		      PARAMETER_AGGRESSIVE,
		      PARAMETER_STABLE};
  int    maxIteration;
  double epsilonStar;
  double lambdaStar;
  double omegaStar;
  double lowerBound;
  double upperBound;
  double betaStar;
  double betaBar;
  double gammaStar;
  double epsilonDash;
  rParameter();
  rParameter(FILE* parameterFile);
  ~rParameter();
  void setDefaultParameter(parameterType type
			   = PARAMETER_DEFAULT);
  void readFile(FILE* parameterFile);
  void display(FILE* fpout=stdout);
};

class rSolutions
{
public:
  rBlockDenseMatrix xMat;
  rBlockDenseMatrix zMat;
  rVector           yVec;
  rBlockDenseMatrix xMatzMat;

  rBlockDenseMatrix choleskyX;
  rBlockDenseMatrix invCholeskyX;
  rBlockDenseMatrix choleskyZ;
  rBlockDenseMatrix invCholeskyZ;

  rBlockVector xzEigenValues;
  double xzMinEigenValue;

  rBlockDenseMatrix workMat;
  rBlockDenseMatrix workMat2;
  rBlockVector      workVec;

  rSolutions();
  rSolutions(int m,int nBlock, int* blockStruct, double lambdaStar,
	     rComputeTime& com);
  ~rSolutions();
  void copyFrom(rSolutions& other);

  void initialize(int m, int nBlock, int* blockStruct,
		  double lambda,rComputeTime& com);

  // if we set initial point,
  // call initializeZero before we set initial point,
  // and call initializeResetup after we set initial point
  void initializeZero(int m, int nBlock,
		      int* blockStruct,
		      rComputeTime& com);
  void initializeResetup(int m, int nBlock,
			 int* blockStruct,
			 rComputeTime& com);
  
  bool update(rStepLength& alpha, rNewton& newton,
	      rComputeTime& com);
  bool update_last(rComputeTime& com);
  void display(FILE* fpout=stdout);
};

class rNewton
{
public:
  rDenseMatrix bMat; // the coefficent of Schur complement
  rVector      gVec; // the right hand side of Schur complement

  rBlockDenseMatrix DxMat;
  rVector           DyVec;
  rBlockDenseMatrix DzMat;
  
  rBlockDenseMatrix rMat;
  
  int* upNonZeroCount;
  enum FormulaType {F1,F2,F3};
  FormulaType* useFormula;

  rBlockDenseMatrix invzMat;

  // temporary matrix
  rBlockDenseMatrix fMat;
  rBlockDenseMatrix gMat;
  rBlockDenseMatrix DxMatDzMat;
  

  rNewton();
  rNewton(int m,int nBlock,int* blockStruct);
  ~rNewton();
  
  void initialize(int m,int nBlock,int* blockStruct);

  void computeFormula(int m, rBlockSparseMatrix* A,
		      double DenseRatio,double Kappa);

  void calF1(double& ret, rDenseMatrix& G,
	     rSparseMatrix& Aj);
  void calF2(double& ret, rDenseMatrix& F, rDenseMatrix& G,
	     rDenseMatrix& X, rSparseMatrix& Aj, bool& hasF2Gcal);
  void calF3(double& ret,
	     rDenseMatrix& F, rDenseMatrix& G,
	     rDenseMatrix& X, rDenseMatrix& invZ,
	     rSparseMatrix& Ai, rSparseMatrix& Aj);

  // B_{i,j} = (X A_i Z^{-1}) \bullet A_j
  void compute_bMat(int m, rBlockSparseMatrix* A,
		    // rBlockSparseMatrix& C,
		    rBlockDenseMatrix& xMat,
		    rBlockDenseMatrix& invzMat,
		    rComputeTime& com);

  enum WHICH_DIRECTION {PREDICTOR, CORRECTOR};
  void compute_rMat(WHICH_DIRECTION direction,
		    rAverageComplementarity& mu,
		    rDirectionParameter& beta,
		    rSolutions& cuurentPt);
  bool Mehrotra(WHICH_DIRECTION direction,
		int m,
		rBlockSparseMatrix* A,
		rBlockSparseMatrix& C,
		rAverageComplementarity& mu,
		rDirectionParameter& beta,
		rSwitch& reduction,
		rPhase& phase,
		rSolutions& currentPt,
		rResiduals& currentRes,
		rComputeTime& com);

  void display(FILE* fpout=stdout);
};

class rResiduals
{
public:
  rVector           primalVec;
  rBlockDenseMatrix dualMat;
  double            normPrimalVec;
  double            normDualMat;
  double            centerNorm;

  rResiduals();
  rResiduals(int m, int nBlock, int* blockStruct,
	     rVector& b, rBlockSparseMatrix& C,
	     rBlockSparseMatrix* A,rSolutions& currentPt);
  ~rResiduals();

  void initialize(int m,int nBlock, int* blockStruct,
		  rVector& b, rBlockSparseMatrix& C,
		  rBlockSparseMatrix* A,rSolutions& currentPt);
  void copyFrom(rResiduals& other);
  
  double computeMaxNorm(rVector& primalVec);
  double computeMaxNorm(rBlockDenseMatrix& dualMat);

  void update(int m, int nBlock, int* blockStruct,
	      rVector& b, rBlockSparseMatrix& C,
	      rBlockSparseMatrix* A,
	      rResiduals& initResidual,
	      rRatioInitResCurrentRes& theta,
	      rSolutions& currentPt,
	      rPhase& phase,
	      rAverageComplementarity& mu,
	      rComputeTime& com);
  void compute(int m, int nBlock, int* blockStruct,
	       rVector& b, rBlockSparseMatrix& C,
	       rBlockSparseMatrix* A, rSolutions& currentPt,
	       rAverageComplementarity& mu);
  void display(FILE* fpout = stdout);

};

class rStepLength
{
public:
  rBlockDenseMatrix workMat1;
  rBlockDenseMatrix workMat2;
  rBlockVector xInvDxEigenValues;
  rBlockVector zInvDzEigenValues;
  rBlockVector workVec;

  double primal;
  double dual;
  rStepLength();
  rStepLength(double alphaP, double alphaD, int nBlock,
	      int* blockStruct);
  ~rStepLength();
  void initialize(double alphaP, double alphaD, int nBlock,
		  int* blockStruct);
  
  static double minBlockVector(rBlockVector& aVec);

  void MehrotraPredictor(rVector& b, rBlockSparseMatrix& C,
			 rBlockSparseMatrix* A,
			 rSolutions& currentPt, rPhase& phase,
			 rNewton& newton, rLanczos& lanczos,
			 rComputeTime& com);
  void MehrotraCorrector(int nDim, rVector& b, rBlockSparseMatrix& C,
			 rBlockSparseMatrix* A,
			 rSolutions& currentPt, rPhase& phase,
			 rSwitch& reduction, rNewton& newton,
			 rAverageComplementarity& mu,
			 rRatioInitResCurrentRes& theta,
			 rLanczos& lanczos,
			 rParameter& param, rComputeTime& com);
  void display(FILE* fpout = stdout);
};

class rDirectionParameter
{
public:
  double value;
  rDirectionParameter(double betaStar=0.0);
  ~rDirectionParameter();
  void initialize(double betaStar=0.0);
  
  void MehrotraPredictor(rPhase& phase, rSwitch& reduction,
			 rParameter& param);
  void MehrotraCorrector(int nDim, rPhase& phase, rStepLength& alpha,
			 rSolutions& currentPt, rNewton& newton,
			 rAverageComplementarity& mu,
			 rParameter& param);
  void display(FILE* fpout = stdout);
};

class rSwitch
{
public:
  enum rSwitchType {ON,OFF};
  rSwitchType switchType;

  rSwitch(rSwitchType switchType=ON);
  ~rSwitch();
  void initialize(rSwitchType switchType=ON);
  
  void MehrotraPredictor(rPhase& phase);
  void display(FILE* fpout = stdout);

};

class rAverageComplementarity
{
public:
  double initial;
  double current;
  rAverageComplementarity(double lambdaStar = 0.0);
  ~rAverageComplementarity();
  void initialize(double lambdaStar = 0.0);
  void initialize(int nDim, rSolutions& initPt);
  void update(int nDim, rSolutions& currentPt);
  void display(FILE* fpout = stdout);
};

class rRatioInitResCurrentRes
{
public:
  double primal;
  double dual;
  
  rRatioInitResCurrentRes();
  rRatioInitResCurrentRes(rParameter& param, rResiduals& initRes);
  ~rRatioInitResCurrentRes();

  void initialize(rParameter& param, rResiduals& initRes);

  void update(rSwitch& reduction, rStepLength& alpha);
  void update_exact(rResiduals& initRes, rResiduals& currentRes);
  void display(FILE* fpout = stdout);
};

class rSolveInfo
{
public:
  enum phaseType { noINFO,pFEAS,dFEAS,pdFEAS,pdINF,pFEAS_dINF,
		   pINF_dFEAS,pdOPT,pUNBD,dUNBD};

  double rho;
  double etaPrimal;
  double etaDual;
  double objValPrimal;
  double objValDual;

  rSolveInfo();
  rSolveInfo(int nDim, rVector& b, rBlockSparseMatrix& C,
	     rBlockSparseMatrix* A, rSolutions& initPt,
	     double mu0, double omegaStar);
  ~rSolveInfo();

  void initialize(int nDim, rVector& b, rBlockSparseMatrix& C,
		  rBlockSparseMatrix* A, rSolutions& initPt,
		  double mu0, double omegaStar);

  void update(int nDim, rVector& b, rBlockSparseMatrix& C,
	      rSolutions& initPt, rSolutions& currentPt,
	      rResiduals& currentRes,
	      rAverageComplementarity& mu,
	      rRatioInitResCurrentRes& theta,
	      rParameter& param);
  void display(FILE* fpout = stdout);
};

class rPhase
{
public:
  int nDim;
  rSolveInfo::phaseType value;

  rPhase();
  rPhase(rResiduals& initRes, rSolveInfo& solveInfo,
	 rParameter& param, int nDim);
  ~rPhase();

  bool initialize(rResiduals& initRes,
		  rSolveInfo& solveInfo,
		  rParameter& param, int nDim);
  bool updateCheck(rResiduals& currentRes,
		   rSolveInfo& solveInfo,
		   rParameter& param);
  void reverse();
  void display(FILE* fpout = stdout);
};

class rLanczos
{
public:
  // temporary matricies
  rBlockDenseMatrix Q;
  rBlockDenseMatrix A;
  rBlockVector      out;
  rBlockVector      b;
  rBlockVector      diag;
  rBlockVector      bet;
  rBlockVector      r;
  rBlockVector      q;
  rBlockVector      qold;
  rBlockVector      w;
  rBlockVector      tmp;
  rBlockVector      diagVec;
  rBlockVector      diagVec2;
  rBlockVector      workVec;

  rLanczos();
  rLanczos(int nBlock, int* blockStruct);
  void initialize(int nBlock, int* blockStruct);
  ~rLanczos();

  // calculate the minimum eigen value of lMat*xMat*(lMat^T)
  // by Lanczos methods.
  // lMat is lower triangular¡¢xMat is symmetric
  double getMinEigen(rBlockDenseMatrix& lMat,
		     rBlockDenseMatrix& xMat);
  // NonDiagonal
  double getMinEigen(rDenseMatrix& lMat, rDenseMatrix& xMat,
		     rDenseMatrix& Q,
		     rVector& out, rVector& b,  rVector& r,
		     rVector& q, rVector& qold, 
		     rVector& w, rVector& tmp,
		     rVector& diagVec, rVector& diagVec2,
		     rVector& workVec);
  // Diagonal
  double getMinEigen(rDenseMatrix& lMat, rDenseMatrix& xMat); 
};


#endif // __rsdpa_parts_h__
/*-------------------------------------------------
  rsdpa_io.h
-------------------------------------------------*/

#ifndef __rsdpa_io_h__
#define __rsdpa_io_h__


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
#endif // __rsdpa_class_h__                

