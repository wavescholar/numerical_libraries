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
  rsdpa_parts.h
-------------------------------------------------*/

#ifndef __rsdpa_parts_h__
#define __rsdpa_parts_h__

#include "rsdpa_algebra.h"

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
