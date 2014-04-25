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
  rsdpa_parts.cpp
-------------------------------------------------*/

#include "rsdpa_parts.h"

rComputeTime::rComputeTime()
{
  Predictor = 0.0;
  Corrector = 0.0;
  
  StepPredictor = 0.0;
  StepCorrector = 0.0;
  
  xMatTime = 0.0;
  zMatTime = 0.0;
  xMatzMatTime = 0.0;
  invzMatTime = 0.0;
  
  EigxMatTime = 0.0;
  EigzMatTime = 0.0;
  EigxMatzMatTime = 0.0;

  makebMat = 0.0;
  B_DIAG   = 0.0;
  B_F1     = 0.0;
  B_F2     = 0.0;
  B_F3     = 0.0;
  B_PRE    = 0.0;

  makegVecMul = 0.0;
  makegVec = 0.0;
  makerMat = 0.0;
  choleskybMat = 0.0;
  
  solve    = 0.0;
  sumDz    = 0.0;
  makedX    = 0.0;
  symmetriseDx = 0.0;
  makedXdZ = 0.0;
  updateRes = 0.0;

  MainLoop  = 0.0;
  FileRead  = 0.0;
  FileCheck = 0.0;
  FileChange= 0.0;
  TotalTime = 0.0;
}

rComputeTime::~rComputeTime()
{
  // Nothing needs.
}


void rComputeTime::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }
  fprintf(fpout,"\n");
  #if 0
  if (TotalTime <= 0.0) {
    return;
  }
  #endif
  fprintf(fpout, "                         Time(sec) ");
  fprintf(fpout," Ratio(%% : MainLoop) \n");
  fprintf(fpout, " Predictor time  =       %f,  %f\n",
	  Predictor, Predictor/MainLoop*100.0);
  fprintf(fpout, " Corrector time  =       %f,  %f\n",
	  Corrector, Corrector/MainLoop*100.0);
  fprintf(fpout, " Make bMat time  =       %f,  %f\n",
	  makebMat, makebMat/MainLoop*100.0);
  fprintf(fpout, " Make bDia time  =       %f,  %f\n",
	  B_DIAG,B_DIAG/MainLoop*100.0);
  fprintf(fpout, " Make bF1  time  =       %f,  %f\n",
	  B_F1,B_F1/MainLoop*100.0);
  fprintf(fpout, " Make bF2  time  =       %f,  %f\n",
	  B_F2,B_F2/MainLoop*100.0);
  fprintf(fpout, " Make bF3  time  =       %f,  %f\n",
	  B_F3,B_F3/MainLoop*100.0);
  fprintf(fpout, " Make bPRE time  =       %f,  %f\n",
	  B_PRE,B_PRE/MainLoop*100.0);
  fprintf(fpout, " Make rMat time  =       %f,  %f\n",
	  makerMat, makerMat/MainLoop*100.0);
  fprintf(fpout, " Make gVec Mul   =       %f,  %f\n",
	  makegVecMul, makegVecMul/MainLoop*100.0);
  fprintf(fpout, " Make gVec time  =       %f,  %f\n",
	  makegVec, makegVec/MainLoop*100.0);
  fprintf(fpout, " Cholesky bMat   =       %f,  %f\n",
	  choleskybMat, choleskybMat/MainLoop*100.0);
  fprintf(fpout, " Ste Pre time    =       %f,  %f\n",
	  StepPredictor, StepPredictor/MainLoop*100.0);
  fprintf(fpout, " Ste Cor time    =       %f,  %f\n",
	  StepCorrector, StepCorrector/MainLoop*100.0);
  fprintf(fpout, " solve           =       %f,  %f\n",
	  solve, solve/MainLoop*100.0);
  fprintf(fpout, " sumDz           =       %f,  %f\n",
  	  sumDz, sumDz/MainLoop*100.0);
  fprintf(fpout, " makedX          =       %f,  %f\n",
  	  makedX, makedX/MainLoop*100.0);
  fprintf(fpout, " symmetriseDx    =       %f,  %f\n",
  	  symmetriseDx, symmetriseDx/MainLoop*100.0);
  fprintf(fpout, " makedXdZ        =       %f,  %f\n",
	  makedXdZ, makedXdZ/MainLoop*100.0);
  fprintf(fpout, " xMatTime        =       %f,  %f\n",
	  xMatTime, xMatTime/MainLoop*100.0);
  fprintf(fpout, " zMatTime        =       %f,  %f\n",
	  zMatTime, zMatTime/MainLoop*100.0);
  fprintf(fpout, " invzMatTime     =       %f,  %f\n",
  	  invzMatTime, invzMatTime/MainLoop*100.0);
  fprintf(fpout, " xMatzMatTime    =       %f,  %f\n",
	  xMatzMatTime, xMatzMatTime/MainLoop*100.0);
  fprintf(fpout, " EigxMatTime     =       %f,  %f\n",
	  EigxMatTime, EigxMatTime/MainLoop*100.0);
  fprintf(fpout, " EigzMatTime     =       %f,  %f\n",
	  EigzMatTime, EigzMatTime/MainLoop*100.0);
  fprintf(fpout, " EigxMatzMatTime =       %f,  %f\n",
	  EigxMatzMatTime, EigxMatzMatTime/MainLoop*100.0);
  fprintf(fpout, " updateRes       =       %f,  %f\n",
	  updateRes, updateRes/MainLoop*100.0);
  double total_eigen = EigxMatTime + EigzMatTime + EigxMatzMatTime;
  fprintf(fpout, " EigTime         =       %f,  %f\n",
	  total_eigen, total_eigen/MainLoop*100.0);
  double sub_total_bMat = MainLoop - makebMat;
  fprintf(fpout, " sub_total_bMat  =       %f,  %f\n",
	  sub_total_bMat, sub_total_bMat/MainLoop*100.0);
  fprintf(fpout, " Main Loop       =       %f,  %f\n",
	  MainLoop, MainLoop/MainLoop*100.0);
  fprintf(fpout, " File Check      =       %f,  %f\n",
	  FileCheck, FileCheck/MainLoop*100.0);
  fprintf(fpout, " File Change     =       %f,  %f\n",
	  FileChange, FileChange/MainLoop*100.0);
  fprintf(fpout, " File Read       =       %f,  %f\n",
	  FileRead, FileRead/MainLoop*100.0);
  fprintf(fpout, " Total           =       %f,  %f\n",
	  TotalTime, TotalTime/MainLoop*100.0);
  fprintf(fpout, "\n");

  return;
}

//-------------------------------------------------------------

rParameter::rParameter()
{
  // setDefaultParameter();
}
rParameter::rParameter(FILE* parameterFile)
{
  readFile(parameterFile);
}

rParameter::~rParameter()
{
  // Nothings needs.
}

void rParameter::setDefaultParameter(rParameter::parameterType type)
{
  if (type == PARAMETER_STABLE) {
    maxIteration =  1000;
    epsilonStar  =  1.0e-7;
    lambdaStar   =  1.0e+2;
    omegaStar    =  2.0;
    lowerBound   = -1.0e+5;
    upperBound   =  1.0e+5;
    betaStar     =  0.2;
    betaBar      =  0.4;
    gammaStar    =  0.5;
    epsilonDash  =  1.0e-7;
  }
  else if (type == PARAMETER_AGGRESSIVE) {
    maxIteration =  100;
    epsilonStar  =  1.0e-7;
    lambdaStar   =  1.0e+2;
    omegaStar    =  2.0;
    lowerBound   = -1.0e+5;
    upperBound   =  1.0e+5;
    betaStar     =  0.01;
    betaBar      =  0.02;
    gammaStar    =  0.98;
    epsilonDash  =  1.0e-7;
  }
  else {
    maxIteration =  100;
    epsilonStar  =  1.0e-7;
    lambdaStar   =  1.0e+2;
    omegaStar    =  2.0;
    lowerBound   = -1.0e+5;
    upperBound   =  1.0e+5;
    betaStar     =  0.1;
    betaBar      =  0.2;
    gammaStar    =  0.9;
    epsilonDash  =  1.0e-7;
  }    
}
void rParameter::readFile(FILE* parameterFile)
{
  fscanf(parameterFile,"%d%*[^\n]",&maxIteration);
  fscanf(parameterFile,"%lf%*[^\n]",&epsilonStar);
  fscanf(parameterFile,"%lf%*[^\n]",&lambdaStar);
  fscanf(parameterFile,"%lf%*[^\n]",&omegaStar);
  fscanf(parameterFile,"%lf%*[^\n]",&lowerBound);
  fscanf(parameterFile,"%lf%*[^\n]",&upperBound);
  fscanf(parameterFile,"%lf%*[^\n]",&betaStar);
  fscanf(parameterFile,"%lf%*[^\n]",&betaBar);
  fscanf(parameterFile,"%lf%*[^\n]",&gammaStar);
  fscanf(parameterFile,"%lf%*[^\n]",&epsilonDash);
}

void rParameter::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }

  fprintf(fpout, "maxIteration =    %d\n",maxIteration);
  fprintf(fpout, "epsilonStar  = %8.3e\n",epsilonStar );
  fprintf(fpout, "lambdaStar   = %8.3e\n",lambdaStar  );
  fprintf(fpout, "omegaStar    = %8.3e\n",omegaStar   );
  fprintf(fpout, "lowerBound   = %8.3e\n",lowerBound  );
  fprintf(fpout, "upperBound   = %8.3e\n",upperBound  );
  fprintf(fpout, "betaStar     = %8.3e\n",betaStar    );
  fprintf(fpout, "betaBar      = %8.3e\n",betaBar     );
  fprintf(fpout, "gammaStar    = %8.3e\n",gammaStar   );
  fprintf(fpout, "epsilonDash  = %8.3e\n",epsilonDash );
  return;
}

//----------------------------------------------------------

rSolutions::rSolutions()
{
  // Nothings needs.
}

rSolutions::~rSolutions()
{
  xMat.~rBlockDenseMatrix();
  zMat.~rBlockDenseMatrix();
  yVec.~rVector();
  xMatzMat.~rBlockDenseMatrix();

  choleskyX.~rBlockDenseMatrix();
  invCholeskyX.~rBlockDenseMatrix();
  choleskyZ.~rBlockDenseMatrix();
  invCholeskyZ.~rBlockDenseMatrix();

  xzEigenValues.~rBlockVector();

  workMat.~rBlockDenseMatrix();
  workMat2.~rBlockDenseMatrix();
}

rSolutions::rSolutions(int m, int nBlock, int* blockStruct,
		       double lambda,rComputeTime& com)
{
  initialize(m,nBlock,blockStruct,lambda,com);
}

void rSolutions::initialize(int m, int nBlock, int* blockStruct,
		double lambda,rComputeTime& com)
{
  xMat.initialize(nBlock,blockStruct);
  xMat.setIdentity(lambda);
  zMat.initialize(nBlock,blockStruct);
  zMat.setIdentity(lambda);
  yVec.initialize(m);
  yVec.setZero();
  xMatzMat.initialize(nBlock,blockStruct);
  xMatzMat.setIdentity(lambda*lambda);

  choleskyX.initialize(nBlock,blockStruct);
  choleskyX.setIdentity(sqrt(lambda));
  choleskyZ.initialize(nBlock,blockStruct);
  choleskyZ.setIdentity(sqrt(lambda));
  invCholeskyX.initialize(nBlock,blockStruct);
  invCholeskyX.setIdentity(1.0/sqrt(lambda));
  invCholeskyZ.initialize(nBlock,blockStruct);
  invCholeskyZ.setIdentity(1.0/sqrt(lambda));

  xzEigenValues.initialize(nBlock,blockStruct);
  
  workMat.initialize(nBlock,blockStruct);
  workMat2.initialize(nBlock,blockStruct);
  
  int* workStruct;
  workStruct = NULL;
  rNewCheck();
  workStruct = new int[nBlock];
  if (workStruct==NULL) {
    rError("rSolutions :: memory exhauseted ");
  }
  for (int l=0; l<nBlock; ++l) {
    if (blockStruct[l]>0) {
      workStruct[l] = 3*blockStruct[l]-1;
    } else {
      workStruct[l] = - 3*blockStruct[l]-1;
    }
    // rMessage("workStruct = " << workStruct[l]);
  }
  workVec.initialize(nBlock,workStruct);
  delete[] workStruct;
  workStruct = NULL;

  // -----------------------------------

#if 0
  // These works have already done in initialization step.
  rTimeStart(START1);
  rAl::getCholeskyAndInv(choleskyX,invCholeskyX,xMat);
  rTimeEnd(END1);
  com.xMatTime += rTimeCal(START1,END1);
  rTimeStart(START2);
  rAl::getCholeskyAndInv(choleskyZ,invCholeskyZ,zMat);
  rTimeEnd(END2);
  com.zMatTime += rTimeCal(START2,END2);

  rTimeStart(START3);
  #if 0
  rAl::let(workMat2,'=',choleskyX,'t',zMat);
  rAl::let(workMat,'=',workMat2,'*',choleskyX);

  // rMessage("eigen test");
  
  xzMinEigenValue =
    rAl::getMinEigenValue(workMat,xzEigenValues,workVec);
  if (xzMinEigenValue < 0.0) {
    rMessage(" xzMinEigenValue < 0.0");
  }
  #endif

  rAl::let(xMatzMat,'=',xMat,'*',zMat);
  rTimeEnd(END3);
  com.xMatzMatTime += rTimeCal(START3,END3);
#endif
}
void rSolutions::initializeZero(int m, int nBlock,
				int* blockStruct,
				rComputeTime& com)
{
  // if we set initial point,
  // we malloc only the space of xMat, yVec, zMat
  xMat.initialize(nBlock,blockStruct);
  xMat.setZero();
  zMat.initialize(nBlock,blockStruct);
  zMat.setZero();
  yVec.initialize(m);
  yVec.setZero();
}

void rSolutions::initializeResetup(int m, int nBlock,
				   int* blockStruct,
				   rComputeTime& com)
{
  // Do the rest work , which are made by initalizeZero
  xMatzMat.initialize(nBlock,blockStruct);
  choleskyX.initialize(nBlock,blockStruct);
  // choleskyX.setIdentity(sqrt(lambda));
  choleskyZ.initialize(nBlock,blockStruct);
  // choleskyZ.setIdentity(sqrt(lambda));
  invCholeskyX.initialize(nBlock,blockStruct);
  // invCholeskyX.setIdentity(1.0/sqrt(lambda));
  invCholeskyZ.initialize(nBlock,blockStruct);
  // invCholeskyZ.setIdentity(1.0/sqrt(lambda));

  xzEigenValues.initialize(nBlock,blockStruct);
  
  workMat.initialize(nBlock,blockStruct);
  workMat2.initialize(nBlock,blockStruct);
  
  int* workStruct;
  workStruct = NULL;
  rNewCheck();
  workStruct = new int[nBlock];
  if (workStruct==NULL) {
    rError("rSolutions :: memory exhauseted ");
  }
  for (int l=0; l<nBlock; ++l) {
    if (blockStruct[l]>0) {
      workStruct[l] = 3*blockStruct[l]-1;
    } else {
      workStruct[l] = - 3*blockStruct[l]-1;
    }
    // rMessage("workStruct = " << workStruct[l]);
  }
  workVec.initialize(nBlock,workStruct);
  delete[] workStruct;
  workStruct = NULL;

  // -----------------------------------
  
  rTimeStart(START1);
  rAl::getCholeskyAndInv(choleskyX,invCholeskyX,xMat);
  rTimeEnd(END1);
  com.xMatTime += rTimeCal(START1,END1);
  rTimeStart(START2);
  rAl::getCholeskyAndInv(choleskyZ,invCholeskyZ,zMat);
  rTimeEnd(END2);
  com.zMatTime += rTimeCal(START2,END2);

  rTimeStart(START3);
  #if 0
  rAl::let(workMat2,'=',choleskyX,'t',zMat);
  rAl::let(workMat,'=',workMat2,'*',choleskyX);

  // rMessage("eigen test");
  
  xzMinEigenValue =
    rAl::getMinEigenValue(workMat,xzEigenValues,workVec);
  if (xzMinEigenValue < 0.0) {
    rMessage(" xzMinEigenValue < 0.0");
  }
  #endif

  rAl::let(xMatzMat,'=',xMat,'*',zMat);
  rTimeEnd(END3);
  com.xMatzMatTime += rTimeCal(START3,END3);
}

void rSolutions::copyFrom(rSolutions& other)
{
  if (this == &other) {
    return;
  }
  xMat.copyFrom(other.xMat);
  yVec.copyFrom(other.yVec);
  zMat.copyFrom(other.zMat);
  xMatzMat.copyFrom(other.xMatzMat);

  choleskyX.copyFrom(other.choleskyX);
  invCholeskyX.copyFrom(other.invCholeskyX);
  choleskyZ.copyFrom(other.choleskyZ);
  invCholeskyZ.copyFrom(other.invCholeskyZ);

  xzEigenValues.copyFrom(other.xzEigenValues);
  workMat.copyFrom(other.workMat);
  workMat2.copyFrom(other.workMat2);
  workVec.copyFrom(other.workVec);
}


bool rSolutions::update(rStepLength& alpha, rNewton& newton,
			rComputeTime& com)
{

  bool total_judge = _SUCCESS;

  rTimeStart(START1_1);
  rAl::let(xMat,'=',xMat,'+',newton.DxMat,&alpha.primal);
  rTimeEnd(END1_1);
  com.xMatTime += rTimeCal(START1_1,END1_1);
  rAl::let(yVec,'=',yVec,'+',newton.DyVec,&alpha.dual);
  rTimeStart(START1_2);
  rAl::let(zMat,'=',zMat,'+',newton.DzMat,&alpha.dual);
  rTimeEnd(END1_2);
  com.zMatTime += rTimeCal(START1_2,END1_2);

  const double cannot_move = 1.0e-4;
  if (alpha.primal < cannot_move && alpha.dual < cannot_move) {
    rMessage("Step length is too small. ");
    return FAILURE;
  }

  rTimeStart(START1_3);
  if (rAl::getCholeskyAndInv(choleskyX,invCholeskyX,xMat) == false) {
    total_judge = FAILURE;
  }
  rTimeEnd(END1_3);
  com.xMatTime += rTimeCal(START1_3,END1_3);
  // rMessage(" xMat cholesky :: " << rTimeCal(START1_3,END1_3));

  rTimeStart(START1_4); 
  if (rAl::getCholeskyAndInv(choleskyZ,invCholeskyZ,zMat) == false) {
    total_judge = FAILURE;
  }
  rTimeEnd(END1_4);
  // rMessage(" zMat cholesky :: " << rTimeCal(START1_4,END1_4));
  com.zMatTime += rTimeCal(START1_4,END1_4);

  rTimeStart(START1);
  rAl::let(xMatzMat,'=',xMat,'*',zMat);
  rTimeEnd(END1);
  com.xMatzMatTime += rTimeCal(START1,END1);

  #if 0
  rTimeStart(START2);
  rAl::let(workMat2,'=',choleskyX,'t',zMat);
  rAl::let(workMat,'=',workMat2,'*',choleskyX);

  xzMinEigenValue =
    rAl::getMinEigenValue(workMat,xzEigenValues,workVec);
  if (xzMinEigenValue < 0.0) {
    rMessage(" xzMinEigenValue < 0.0");
    total_judge = FAILURE;
  }
  rTimeEnd(END2);
  com.EigxMatzMatTime += rTimeCal(START2,END2);
  #else
  xzMinEigenValue = 1.0;
  #endif
  
  return total_judge;
}

bool rSolutions::update_last(rComputeTime& com)
{
  rTimeStart(START1);
  bool total_judge = _SUCCESS;

#if 1
  // cut calculation of eigenvalues of xMatzMat
  xzMinEigenValue = 0.0;
#else  
  rAl::let(workMat2,'=',choleskyX,'t',zMat);
  rAl::let(workMat,'=',workMat2,'*',choleskyX);

  // rMessage("eigen test");
  xzMinEigenValue =
    rAl::getMinEigenValue(workMat,xzEigenValues,workVec);

  rAl::let(xMatzMat,'=',xMat,'*',zMat);

  if (xzMinEigenValue < 0.0) {
    rMessage(" xzMinEigenValue < 0.0");
  }
#endif
  rTimeEnd(END1);
  com.EigxMatzMatTime += rTimeCal(START1,END1);
  
  return total_judge;
}

void rSolutions::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }

  fprintf(fpout,"xMat = \n");
  xMat.display(fpout);
  fprintf(fpout,"yVec = \n");
  yVec.display(fpout);
  fprintf(fpout,"zMat = \n");
  zMat.display(fpout);
}

// ------------------------------------------------------------

rNewton::rNewton()
{
  upNonZeroCount = NULL;
  useFormula     = NULL;
}

rNewton::rNewton(int m,int nBlock,int* blockStruct)
{
  initialize(m,nBlock,blockStruct);
}

rNewton::~rNewton()
{
  bMat.~rDenseMatrix();
  gVec.~rVector();
  DxMat.~rBlockDenseMatrix();
  DyVec.~rVector();
  DzMat.~rBlockDenseMatrix();
  rMat.~rBlockDenseMatrix();

  if (upNonZeroCount!=NULL) {
    delete[] upNonZeroCount;
  }
  upNonZeroCount = NULL;
  if (useFormula!=NULL) {
    delete[] useFormula;
  }
  useFormula = NULL;

  invzMat.~rBlockDenseMatrix();
  fMat.~rBlockDenseMatrix();
  gMat.~rBlockDenseMatrix();
  DxMatDzMat.~rBlockDenseMatrix();
}

void rNewton::initialize(int m, int nBlock, int* blockStruct)
{
  bMat.initialize(m,m,rDenseMatrix::DENSE);
  gVec.initialize(m);

  DxMat.initialize(nBlock,blockStruct);
  DyVec.initialize(m);
  DzMat.initialize(nBlock,blockStruct);
  rMat.initialize(nBlock,blockStruct);

  invzMat.initialize(nBlock,blockStruct);
  fMat.initialize(nBlock,blockStruct);
  gMat.initialize(nBlock,blockStruct);
  DxMatDzMat.initialize(nBlock,blockStruct);

  rNewCheck();
  upNonZeroCount = new int[m*nBlock];
  if (upNonZeroCount == NULL) {
    rError("rNewton:: memory exhausted ");
  }
  rNewCheck();
  useFormula = new FormulaType[m*nBlock];
  if (useFormula == NULL) {
    rError("rNewton:: memory exhausted ");
  }
}

void rNewton::computeFormula(int m, rBlockSparseMatrix* A,
			     double DenseRatio, double Kappa)
{
  // We have no chance to use DenseRatio
  if (upNonZeroCount == NULL || useFormula == NULL) {
    rError("rNewton:: failed initialization");
  }

  #if 0
  for (int k=0; k<m; ++k) {
    for (int l=0; l<A[0].nBlock; ++l) {
      rMessage("A[" << k << "].ele[" << l << "] ="
	       << A[k].ele[l].NonZeroEffect);
    }
  }
  #endif

  // Count sum of number of elements
  // that each number of elements are less than own.

  int nBlock = A[0].nBlock;
  for (int l=0; l<nBlock; ++l) {
    if (A[0].blockStruct[l] < 0) {
      // in Diagonal case, we don't have to calculate.
      continue;
    }
    for (int k=0; k<m; ++k) {
      int up = A[k].ele[l].NonZeroEffect;
      // rMessage("up = " << up);
      for (int k2=0; k2<m; ++k2) {
	if (A[k2].ele[l].NonZeroEffect < A[k].ele[l].NonZeroEffect) {
	  up += A[k2].ele[l].NonZeroEffect;
	}
	#if 1
	else if (A[k2].ele[l].NonZeroEffect ==
		   A[k ].ele[l].NonZeroEffect
		   && k2<k ) {
	 up += A[k2].ele[l].NonZeroEffect;
	}
	#endif
      }
      upNonZeroCount[k*nBlock + l] = up;
      // rMessage("up = " << up);
    }
  }


  // Determine which formula
  for (int l=0; l<nBlock; ++l) {
    if (A[0].blockStruct[l] < 0) {
      // in Diagonal case, we don't have to calculate.
      continue;
    }
    int countf1,countf2,countf3;
    countf1 = countf2 = countf3 = 0;
    for (int k=0; k<m; ++k) {
      double f1,f2,f3;
      #if 0
      int n       = A[k].ele[l].nRow;
      int up      = upNonZeroCount[k*nBlock + l];
      int nonzero = A[k].ele[l].NonZeroEffect;
      #else
      double n       = A[k].ele[l].nRow;
      double up      = upNonZeroCount[k*nBlock + l];
      double nonzero = A[k].ele[l].NonZeroEffect;
      #endif
      f1 = Kappa*n*nonzero + n*n*n + Kappa*up;
      f2 = Kappa*n*nonzero + Kappa*(n+1)*up;
      #if 0
      f3 = Kappa*(2*Kappa*nonzero+1)*up/Kappa;
      #else
      f3 = Kappa*(2*Kappa*nonzero+1)*up;
      #endif
      // rMessage("up = " << up << " nonzero = " << nonzero);
      // rMessage("f1=" << f1 << " f2=" << f2 << " f3=" << f3);
      // printf("%d %d %lf %lf %lf %lf\n",k,l,nonzero,f1,f2,f3);
      if (A[k].ele[l].Sp_De_Di == rSparseMatrix::DENSE) {
	// if DENSE, we use only F1 or F2,
	// that is we don't use F3
	if (f1<f2) {
	  useFormula[k*nBlock+l] = F1;
	  countf1++;
	} else {
	  useFormula[k*nBlock+l] = F2;
	  countf2++;
	}
      } else {
	// this case is SPARSE
	if (f1<f2 && f1<f3) {
	  // rMessage("line " << k << " is F1");
	  useFormula[k*nBlock+l] = F1;
	  countf1++;
	} else if (f2<f3) {
	  // rMessage("line " << k << " is F2");
	  useFormula[k*nBlock+l] = F2;
	  countf2++;
	} else {
	  // rMessage("line " << k << " is F3");
	  useFormula[k*nBlock+l] = F3;
	  countf3++;
	}
      }
    }
    // rMessage("Kappa = " << Kappa);
    #if 0
    rMessage("count f1 = " << countf1
	     << ":: count f2 = " << countf2
	     << ":: count f3 = " << countf3);
    #endif
  } // end of 'for (int l)'
  return;
}

void rNewton::calF1(double& ret, rDenseMatrix& G,
		    rSparseMatrix& Aj)
{
  rAl::let(ret,'=',Aj,'.',G);
}
void rNewton::calF2(double& ret,
		    rDenseMatrix& F, rDenseMatrix& G,
		    rDenseMatrix& X, rSparseMatrix& Aj,
		    bool& hasF2Gcal)
{
  int alpha,beta;
  double value1,value2;

  int n    = Aj.nRow;
  int index = 0;
  // rMessage(" using F2 ");
  switch (Aj.Sp_De_Di) {
  case rSparseMatrix::SPARSE:
    // rMessage("F2::SPARSE  " << Aj.NonZeroCount);
    ret = 0.0;
    for (index = 0; index < Aj.NonZeroCount; ++index) {
      alpha  = Aj.row_index[index];
      beta   = Aj.column_index[index];
      value1 = Aj.sp_ele[index];

      // value2 = ddot_(&n, &X.de_ele[alpha+n*0], &n,
      //	     &F.de_ele[0+n*beta], &IONE);
      value2 = ddot(&n, &X.de_ele[alpha], &n,
		     &F.de_ele[n*beta], &IONE);
      ret += value1*value2;
      if (alpha!=beta) {
	//value2 = ddot_(&n, &X.de_ele[beta+n*0], &n,
	//       &F.de_ele[0+n*alpha], &IONE);
	value2 = ddot(&n, &X.de_ele[beta], &n,
		       &F.de_ele[n*alpha], &IONE);
	ret += value1*value2;
      }
    }
    break;
  case rSparseMatrix::DENSE:
    // G is temporary matrix
    // rMessage("F2::DENSE");
    if (hasF2Gcal == false) {
      // rMessage(" using F2 changing to F1");
      rAl::let(G,'=',X,'*',F);
      hasF2Gcal = true;
    }
    rAl::let(ret,'=',Aj,'.',G);
    break;
  case rSparseMatrix::DIAGONAL:
    // DIAGONAL is something wrong
    rMessage("F2::DIAGONAL");
    rMessage("calF2:: miss condition");
    break;
  } // end of switch
}

void rNewton::calF3(double& ret,
		    rDenseMatrix& F, rDenseMatrix& G,
		    rDenseMatrix& X, rDenseMatrix& invZ,
		    rSparseMatrix& Ai, rSparseMatrix& Aj)
{
  #if 0
  if (Ai.Sp_De_Di == rSparseMatrix::DIAGONAL
      || Aj.Sp_De_Di == rSparseMatrix::DIAGONAL) {
    rMessage("calF3:: miss condition");
  } // end of switch
  #endif

  #if 0
  if (Ai.Sp_De_Di == rSparseMatrix::DENSE) {
    rMessage("change from F3 to F2");
    rAl::let(F,'=',Ai,'*',invZ);
    calF2(ret,F,G,X,Aj);
    return;
  }
  #endif 

  #if 0
  // Ai must have more elements than Aj.
  // We must not have this followring case.
  if (Ai.Sp_De_Di == rSparseMatrix::SPARSE
      && Aj.Sp_De_Di == rSparseMatrix::DENSE) {
    rError("calF3:: miss condition");
  }
  #endif

  // Ai and Aj are SPARSE
  ret = 0.0;
  int index2=0;
  double sum;
  // rMessage("Aj.NonZeroCount = " << Aj.NonZeroCount);
  for (int index1=0; index1<Aj.NonZeroCount; ++index1) {
    int alpha = Aj.row_index[index1];
    int beta  = Aj.column_index[index1];
    double value1 = Aj.sp_ele[index1];
    sum = 0.0;
    #if 1
    for (index2=0; index2<Ai.NonZeroCount; ++index2) {
      int gamma = Ai.row_index[index2];
      int delta  = Ai.column_index[index2];
      double value2 = Ai.sp_ele[index2];
      double plu = value2*invZ.de_ele[delta+invZ.nCol*beta]
        * X.de_ele[alpha+X.nCol*gamma];
      sum += plu;
      if (gamma!=delta) {
        double plu2 = value2*invZ.de_ele[gamma+invZ.nCol*beta]
          * X.de_ele[alpha+X.nCol*delta];
        sum += plu2;
      }
    }
    ret += value1*sum;
    #else
    int shou = Ai.NonZeroCount / 4;
    int amari = Ai.NonZeroCount % 4;
    for (int index2=0; index2<amari; ++index2) {
      int gamma = Ai.row_index[index2];
      int delta  = Ai.column_index[index2];
      double value2 = Ai.sp_ele[index2];
      double plu = value2*invZ.de_ele[delta+invZ.nCol*beta]
        * X.de_ele[alpha+X.nCol*gamma];
      sum += plu;
      if (gamma!=delta) {
        double plu2 = value2*invZ.de_ele[gamma+invZ.nCol*beta]
          * X.de_ele[alpha+X.nCol*delta];
        sum += plu2;
      }
    }
    for (int index=amari,count=0; count<shou;
         ++count, index+=4) {
      int gamma_1 = Ai.row_index[index];
      int delta_1  = Ai.column_index[index];
      double value_1 = Ai.sp_ele[index];
      double plu_1 = value_1*invZ.de_ele[delta_1+invZ.nCol*beta]
        * X.de_ele[alpha+X.nCol*gamma_1];
      double plu__1 = 0.0;
      if (gamma_1!=delta_1) {
        plu__1 = value_1*invZ.de_ele[gamma_1+invZ.nCol*beta]
          * X.de_ele[alpha+X.nCol*delta_1];
      }
      int gamma_2 = Ai.row_index[index+1];
      int delta_2  = Ai.column_index[index+1];
      double value_2 = Ai.sp_ele[index+1];
      double plu_2 = value_2*invZ.de_ele[delta_2+invZ.nCol*beta]
        * X.de_ele[alpha+X.nCol*gamma_2];
      double plu__2 = 0.0;
      if (gamma_2!=delta_2) {
        plu__2 = value_2*invZ.de_ele[gamma_2+invZ.nCol*beta]
          * X.de_ele[alpha+X.nCol*delta_2];
      }
      int gamma_3 = Ai.row_index[index+2];
      int delta_3  = Ai.column_index[index+2];
      double value_3 = Ai.sp_ele[index+2];
      double plu_3 = value_3*invZ.de_ele[delta_3+invZ.nCol*beta]
        * X.de_ele[alpha+X.nCol*gamma_3];
      double plu__3 = 0.0;
      if (gamma_3!=delta_3) {
        plu__3 = value_3*invZ.de_ele[gamma_3+invZ.nCol*beta]
          * X.de_ele[alpha+X.nCol*delta_3];
      }
      int gamma_4 = Ai.row_index[index+3];
      int delta_4  = Ai.column_index[index+3];
      double value_4 = Ai.sp_ele[index+3];
      double plu_4 = value_4*invZ.de_ele[delta_4+invZ.nCol*beta]
        * X.de_ele[alpha+X.nCol*gamma_4];
      double plu__4 = 0.0;
      if (gamma_4!=delta_4) {
        plu__4 = value_4*invZ.de_ele[gamma_4+invZ.nCol*beta]
          * X.de_ele[alpha+X.nCol*delta_4];
      }
      sum += (plu_1+plu__1 + plu_2+plu__2
              + plu_3+plu__3 + plu_4+plu__4 );
    }
    ret += value1*sum;
    #endif
    if (alpha==beta) {
      continue;
    }
    sum = 0.0;
    #if 1
    for (index2=0; index2<Ai.NonZeroCount; ++index2) {
      int gamma = Ai.row_index[index2];
      int delta  = Ai.column_index[index2];
      double value2 = Ai.sp_ele[index2];
      double plu = value2*invZ.de_ele[delta+invZ.nCol*alpha]
        * X.de_ele[beta+X.nCol*gamma];
      sum += plu;
      if (gamma!=delta) {
        double plu2 = value2*invZ.de_ele[gamma+invZ.nCol*alpha]
          * X.de_ele[beta+X.nCol*delta];
        sum += plu2;
      }
    }
    ret += value1*sum;
    #else
    shou = Ai.NonZeroCount / 4;
    amari = Ai.NonZeroCount % 4;
    for (int index2=0; index2<amari; ++index2) {
      int gamma = Ai.row_index[index2];
      int delta  = Ai.column_index[index2];
      double value2 = Ai.sp_ele[index2];
      double plu = value2*invZ.de_ele[delta+invZ.nCol*alpha]
        * X.de_ele[beta+X.nCol*gamma];
      sum += plu;
      if (gamma!=delta) {
        double plu2 = value2*invZ.de_ele[gamma+invZ.nCol*alpha]
          * X.de_ele[beta+X.nCol*delta];
        sum += plu2;
      }
    }
    for (int index=amari,count=0; count<shou;
         ++count, index+=4) {
      int gamma_1 = Ai.row_index[index];
      int delta_1  = Ai.column_index[index];
      double value_1 = Ai.sp_ele[index];
      double plu_1 = value_1*invZ.de_ele[delta_1+invZ.nCol*alpha]
        * X.de_ele[beta+X.nCol*gamma_1];
      double plu__1 = 0.0;
      if (gamma_1!=delta_1) {
        plu__1 = value_1*invZ.de_ele[gamma_1+invZ.nCol*alpha]
         * X.de_ele[beta+X.nCol*delta_1];
      }
      int gamma_2 = Ai.row_index[index+1];
      int delta_2  = Ai.column_index[index+1];
      double value_2 = Ai.sp_ele[index+1];
      double plu_2 = value_2*invZ.de_ele[delta_2+invZ.nCol*alpha]
        * X.de_ele[beta+X.nCol*gamma_2];
      double plu__2 = 0.0;
      if (gamma_2!=delta_2) {
        plu__2 = value_2*invZ.de_ele[gamma_2+invZ.nCol*alpha]
          * X.de_ele[beta+X.nCol*delta_2];
      }
      int gamma_3 = Ai.row_index[index+2];
      int delta_3  = Ai.column_index[index+2];
      double value_3 = Ai.sp_ele[index+2];
      double plu_3 = value_3*invZ.de_ele[delta_3+invZ.nCol*alpha]
        * X.de_ele[beta+X.nCol*gamma_3];
      double plu__3 = 0.0;
      if (gamma_3!=delta_3) {
        plu__3 = value_3*invZ.de_ele[gamma_3+invZ.nCol*alpha]
          * X.de_ele[beta+X.nCol*delta_3];
      }
      int gamma_4 = Ai.row_index[index+3];
      int delta_4  = Ai.column_index[index+3];
      double value_4 = Ai.sp_ele[index+3];
      double plu_4 = value_4*invZ.de_ele[delta_4+invZ.nCol*alpha]
        * X.de_ele[beta+X.nCol*gamma_4];
      double plu__4 = 0.0;
      if (gamma_4!=delta_4) {
        plu__4 = value_4*invZ.de_ele[gamma_4+invZ.nCol*alpha]
          * X.de_ele[beta+X.nCol*delta_4];
      }
      sum += (plu_1+plu__1 + plu_2+plu__2
              + plu_3+plu__3 + plu_4+plu__4 );
    }
    ret += value1*sum;
    #endif
  } // end of 'for (index1)'
  return;
}

void rNewton::compute_bMat(int m, rBlockSparseMatrix* A,
			   // rBlockSparseMatrix& C,
			   rBlockDenseMatrix& xMat,
			   rBlockDenseMatrix& invzMat,
			   rComputeTime& com)
{
  int nBlock = A[0].nBlock;
  int* blockStruct = A[0].blockStruct;
  
  bMat.setZero();
  for (int l=0; l<nBlock; ++l) {
    if (blockStruct[l]<0) {
      rTimeStart(B_DIAG_START1);
      // case Diagonal
      for (int i=0; i<m; ++i) {
	rAl::let(fMat.ele[l],'=',A[i].ele[l],'*',invzMat.ele[l]);
	rAl::let(gMat.ele[l],'=',xMat.ele[l],'*',fMat.ele[l]);
	for (int j=0; j<=i; ++j) {
	  double value;
	  rAl::let(value,'=',gMat.ele[l],'.',A[j].ele[l]);
	  if (i!=j) {
	    bMat.de_ele[i+bMat.nCol*j] += value;
	    bMat.de_ele[j+bMat.nCol*i] += value;
	  } else {
	    bMat.de_ele[i+bMat.nCol*i] += value;
	  }
	}
      }
      rTimeEnd(B_DIAG_END1);
      com.B_DIAG += rTimeCal(B_DIAG_START1,B_DIAG_END1);
    } else {
      // case not Diagonal
      for (int i=0; i<m; ++i) {
        const int A_i_ele_l_NonZeroEffect = A[i].ele[l].NonZeroEffect;
        if (A_i_ele_l_NonZeroEffect == 0) {
  		continue;
        }
	FormulaType formula = useFormula[i*nBlock + l];
	// ---------------------------------------------------
	// formula = F3; // this is force change
	// ---------------------------------------------------
	rTimeStart(B_NDIAG_START1);
	rTimeStart(B_NDIAG_START2);

	bool hasF2Gcal = false;
	if (formula==F1) {
	  rAl::let(fMat.ele[l],'=',A[i].ele[l],'*',invzMat.ele[l]);
	  rAl::let(gMat.ele[l],'=',xMat.ele[l],'*',fMat.ele[l]);
	} else if (formula==F2) {
	  rAl::let(fMat.ele[l],'=',A[i].ele[l],'*',invzMat.ele[l]);
	  hasF2Gcal = false;
	  // rAl::let(gMat.ele[l],'=',xMat.ele[l],'*',fMat.ele[l]);
	}
	rTimeEnd(B_NDIAG_END2);
	com.B_PRE += rTimeCal(B_NDIAG_START2,B_NDIAG_END2);
	for (int j=0; j<m; ++j) {
          const int A_j_ele_l_NonZeroEffect = A[j].ele[l].NonZeroEffect;
          if (A_j_ele_l_NonZeroEffect == 0) {
            continue;
          }
	  // Select the formula A[i] or the formula A[j].
	  // Use formula that has more NonZeroEffects than others.
	  // We must calculate i==j.
          if (A_i_ele_l_NonZeroEffect < A_j_ele_l_NonZeroEffect
              || ( (A_i_ele_l_NonZeroEffect == A_j_ele_l_NonZeroEffect)
                   && i<j) ) {
            continue;
          }
#if 0
	  rMessage("A[" << i << "].ele[" <<l <<
		   "].NonZeroEffect =" <<
		   A[i].ele[l].NonZeroEffect );
	  rMessage("A[" << j << "].ele[" <<l <<
		   "].NonZeroEffect =" <<
		   A[j].ele[l].NonZeroEffect );
#endif

	  double value;
	  switch (formula) {
	  case F1:
	    // rMessage("calF1");
	    calF1(value,gMat.ele[l],A[j].ele[l]);
	    break;
	  case F2:
	    // rMessage("calF2 ");
	    calF2(value,fMat.ele[l],gMat.ele[l],xMat.ele[l],
		  A[j].ele[l],hasF2Gcal);
	    // calF1(value2,gMat.ele[l],A[j].ele[l]);
	    // rMessage("calF2:  " << (value-value2));
	    break;
	  case F3:
	    // rMessage("calF3");
	    calF3(value,fMat.ele[l],gMat.ele[l],
		  xMat.ele[l],invzMat.ele[l],A[i].ele[l],
		  A[j].ele[l]);
	    break;
	  } // end of switch
	  if (i!=j) {
	    bMat.de_ele[i+bMat.nCol*j] += value;
	    bMat.de_ele[j+bMat.nCol*i] += value;
	  } else {
	    bMat.de_ele[i+bMat.nCol*i] += value;
	  }
	} // end of 'for (int j)'
	rTimeEnd(B_NDIAG_END1);
	double t = rTimeCal(B_NDIAG_START1,B_NDIAG_END1);
	switch (formula) {
	case F1: com.B_F1 += t; break;
	case F2: com.B_F2 += t; break;
	case F3: com.B_F3 += t; break;
	}
      } // end of 'for (int i)'
    } // end of non Diagonal
  } // end of 'for (int l)'
}

void rNewton::compute_rMat(rNewton::WHICH_DIRECTION direction,
			   rAverageComplementarity& mu,
			   rDirectionParameter& beta,
			   rSolutions& currentPt)
{

  //     CORRECTOR ::  R = -XZ -dXdZ + mu I
  // not CORRECTOR :: R = -XZ + mu I
  rAl::let(rMat,'=',currentPt.xMatzMat,'*',&DMONE);
  if (direction == CORRECTOR) {
    rAl::let(DxMatDzMat,'=',DxMat,'*',DzMat);
    rAl::let(rMat,'=',rMat,'+',DxMatDzMat,&DMONE);
  }
  // rMessage("rMat = ");
  // rMat.display();
  double target_mu = beta.value*mu.current;
  for (int l=0; l<rMat.nBlock; ++l) {
    int size = rMat.blockStruct[l];
    if (size < 0) {
      // case Diagonal 
      size = -size;
      double* target = rMat.ele[l].di_ele;
      #if 0
      for (int j=0; j<size; ++j) {
	target[j] += target_mu;
      }
      #else
      int shou = size / 4;
      int amari = size % 4;
      for (int j=0; j<amari; ++j) {
	target[j] += target_mu;
      }
	  int count;
      for (int j=amari,count=0; count<shou; ++count, j+=4) {
	target[j] += target_mu;
	target[j+1] += target_mu;
	target[j+2] += target_mu;
	target[j+3] += target_mu;
      }
      #endif
    } else { // case non Diagonal
      double* target = rMat.ele[l].de_ele;
      #if 0
      for (int j=0; j<size; ++j) {
	target[j+j*size] += target_mu;
      }
      #else
      int shou = size / 4;
      int amari = size % 4;
      for (int j=0; j<amari; ++j) {
	target[j+j*size] += target_mu;
      }
	  int count;
      for (int j=amari,count=0; count<shou; ++count, j+=4) {
	target[j+j*size] += target_mu;
	target[(j+1)+(j+1)*size] += target_mu;
	target[(j+2)+(j+2)*size] += target_mu;
	target[(j+3)+(j+3)*size] += target_mu;
      }
      #endif
    }
  }
}

bool rNewton::Mehrotra(rNewton::WHICH_DIRECTION direction,
		       int m,
		       rBlockSparseMatrix* A,
		       rBlockSparseMatrix& C,
		       rAverageComplementarity& mu,
		       rDirectionParameter& beta,
		       rSwitch& reduction,
		       rPhase& phase,
		       rSolutions& currentPt,
		       rResiduals& currentRes,
		       rComputeTime& com)
{
  if (direction==PREDICTOR) {
    // invzMat = Z^{-1}
    rTimeStart(S1);
    #if 0
    rAl::let(invzMat,'=',currentPt.invCholeskyZ,'t',
             currentPt.invCholeskyZ);
    #else
    // bool total_judge = _SUCCESS;
    for (int l=0; l<currentPt.invCholeskyZ.nBlock; ++l) {
      if (currentPt.invCholeskyZ.ele[l].De_Di
	  == rDenseMatrix::DENSE) {
	invzMat.ele[l].copyFrom(currentPt.invCholeskyZ.ele[l]);
        dtrmm("Left","Lower","Transpose","NonUnitDiag",
               &currentPt.invCholeskyZ.ele[l].nRow,
               &currentPt.invCholeskyZ.ele[l].nCol,
	       &DONE,
               currentPt.invCholeskyZ.ele[l].de_ele,
               &currentPt.invCholeskyZ.ele[l].nRow,
	       invzMat.ele[l].de_ele,
               &invzMat.ele[l].nRow);
      }
      else {
        for (int j=0;j<invzMat.ele[l].nRow; ++j) {
          invzMat.ele[l].di_ele[j] =
            currentPt.invCholeskyZ.ele[l].di_ele[j]
            *currentPt.invCholeskyZ.ele[l].di_ele[j];
        }
      }
    }
    #endif
    rTimeEnd(E1);
    com.invzMatTime += rTimeCal(S1,E1);
    // rMessage("invzMatTime = " << rTimeCal(S1,E1));
  }
  rTimeStart(START1);
  // rMessage("mu = " << mu.current);
  // rMessage("beta = " << beta.value);
  compute_rMat(direction,mu,beta,currentPt);
  // rMessage("rMat = ");
  // rMat.display();
  rTimeEnd(END1);
  com.makerMat += rTimeCal(START1,END1);

  rTimeStart(START2);
  rTimeStart(START_GVEC_MUL);
  if (phase.value == rSolveInfo:: pFEAS
      || phase.value == rSolveInfo::noINFO) {
    // currentPt is infeasilbe, that is the residual
    // dualMat is not 0.
    // fMat, gMat is temporary
    // gMat = R-XD
    rAl::let(fMat,'=',currentPt.xMat,'*',currentRes.dualMat);
    rAl::let(gMat,'=',rMat,'+',fMat,&DMONE);
  } else {
    // dualMat == 0
    gMat.copyFrom(rMat);
  }
  
  // rMessage("currentRes.dualMat = ");
  // currentRes.dualMat.display();
  // rMessage("K = R-XD = ");
  // gMat.display();
  rAl::let(fMat,'=',gMat,'*',invzMat,&DMONE);
  // fMat = (-1)*(R-XD)*Z^{-1}
  // rMessage("- (R-XD)*Z^{-1} =  ");
  // fMat.display();
  // rMessage("currentRes.primalVec =  ");
  // currentRes.primalVec.display();
  // rMessage("currentRes.dualMat =  ");
  // currentRes.dualMat.display();
  rTimeEnd(END_GVEC_MUL);
  com.makegVecMul += rTimeCal(START_GVEC_MUL,END_GVEC_MUL);
    
  for (int k=0; k<m; ++k) {
    rAl::let(gVec.ele[k],'=',fMat,'.',A[k]);
  }
  // rMessage("gVec =  ");
  // gVec.display();
  // rMessage("currentRes.primalVec = ");
  // currentRes.primalVec.display();

  #if 0
  if (phase.value == rSolveInfo:: dFEAS
      || phase.value == rSolveInfo::noINFO) {
  #endif
    rAl::let(gVec,'=',gVec,'+',currentRes.primalVec);
  #if 0
  }
  #endif

  // rMessage("gVec =  ");
  // gVec.display();
  
  rTimeEnd(END2);
  com.makegVec += rTimeCal(START2,END2);

  if (direction == PREDICTOR) {
    rTimeStart(START3);
    compute_bMat(m,A,currentPt.xMat,invzMat,com);
    // rMessage("Cholesky bMat");
    // rMessage("bMat =  ");
    // bMat.display();
    rTimeEnd(END3);
    com.makebMat += rTimeCal(START3,END3);
    rTimeStart(START3_2);
    bool ret = rAl::choleskyFactorWithAdjust(bMat);
    if (ret == FAILURE) {
      return FAILURE;
    }
    rTimeEnd(END3_2);
    com.choleskybMat += rTimeCal(START3_2,END3_2);
  }
  // bMat is already cholesky factorized.
  rTimeStart(START4);
  rAl::let(DyVec,'=',bMat,'/',gVec);
  rTimeEnd(END4);
  com.solve += rTimeCal(START4,END4);
  // rMessage("DyVec =  ");
  // DyVec.display();

  rTimeStart(START5);
  if (phase.value == rSolveInfo:: pFEAS
      || phase.value == rSolveInfo::noINFO) {
    DzMat.copyFrom(currentRes.dualMat);
  } else {
    DzMat.setZero();
  }
  rTimeStart(START_SUMDZ);
  for (int k=0; k<m; ++k) {
    rAl::let(DzMat,'=',DzMat,'-',A[k],&DyVec.ele[k]);
  }
  rTimeEnd(END_SUMDZ);

  // fMat,gMat are temporary.
  rTimeStart(START_DX);
  rAl::let(fMat,'=',currentPt.xMat,'*',DzMat);
  rAl::let(gMat,'=',rMat,'+',fMat,&DMONE);
  rAl::let(DxMat,'=',gMat,'*',invzMat);
  rTimeEnd(END_DX);
  rTimeStart(START_SYMM);
  rAl::getSymmetrize(DxMat);
  rTimeEnd(END_SYMM);
  // rMessage("DxMat =  ");
  // DxMat.display();
  rTimeEnd(END5);
  com.makedXdZ += rTimeCal(START5,END5);
  com.sumDz += rTimeCal(START_SUMDZ,END_SUMDZ);
  com.makedX += rTimeCal(START_DX,END_DX);
  com.symmetriseDx += rTimeCal(START_SYMM,END_SYMM);
  return true;
}

void rNewton::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }

  fprintf(fpout,"rNewton.DxMat = \n");
  DxMat.display(fpout);
  fprintf(fpout,"rNewton.DyVec = \n");
  DyVec.display(fpout);
  fprintf(fpout,"rNewton.DzMat = \n");
  DzMat.display(fpout);
}

//---------------------------------------------------------

rResiduals::rResiduals()
{
  normPrimalVec = 0.0;
  normDualMat   = 0.0;
  centerNorm    = 0.0;
}

rResiduals::rResiduals(int m, int nBlock, int* blockStruct,
		       rVector& b, rBlockSparseMatrix& C,
		       rBlockSparseMatrix* A,rSolutions& currentPt)
{
  initialize(m,nBlock,blockStruct,b,C,A,currentPt);
}

rResiduals::~rResiduals()
{
  primalVec.~rVector();
  dualMat.~rBlockDenseMatrix();
}

void rResiduals::copyFrom(rResiduals& other)
{
  if (this==&other) {
    return;
  }
  primalVec.copyFrom(other.primalVec);
  dualMat.copyFrom(other.dualMat);
  normPrimalVec = other.normPrimalVec;
  normDualMat   = other.normDualMat;
  centerNorm    = other.centerNorm;
}

void rResiduals::initialize(int m, int nBlock, int* blockStruct,
			    rVector& b, rBlockSparseMatrix& C,
			    rBlockSparseMatrix* A,
			    rSolutions& currentPt)
{
  primalVec.initialize(m);
  dualMat.initialize(nBlock,blockStruct);

  // p[k] = b[k] - A[k].X;
  for (int k=0; k<m; ++k) {
    double ip;
    rAl::let(ip,'=',A[k],'.',currentPt.xMat);
    primalVec.ele[k] = b.ele[k] - ip;
  }
  
  // D = C - Z - \sum A[k]y[k]
  rAl::let(dualMat,'=',C,'+',currentPt.zMat,&DMONE);
  
  for (int k=0; k<m; ++k) {
    rAl::let(dualMat,'=',dualMat,'-',A[k],&currentPt.yVec.ele[k]);
  }

  // rMessage("primal residual =");
  // primalVec.display();
  // rMessage("dual residual =");
  // dualMat.display();
  
  normPrimalVec = computeMaxNorm(primalVec);
  normDualMat   = computeMaxNorm(dualMat);
}

double rResiduals::computeMaxNorm(rVector& primalVec)
{
  double ret = 0.0;
  #if 1
  for (int k=0; k<primalVec.nDim; ++k) {
    double tmp = fabs(primalVec.ele[k]);
    if (tmp > ret) {
      ret = tmp;
    }
  }
  #else
  int index = idamax_(&primalVec.nDim,primalVec.ele,&IONE);
  ret = fabs(primalVec.ele[index]);
  #endif
  return ret;
}

double rResiduals::computeMaxNorm(rBlockDenseMatrix& dualMat)
{
  int nBlock = dualMat.nBlock;
  int* blockStruct = dualMat.blockStruct;

  double ret = 0.0;
  double tmp;
  for (int l=0; l<nBlock; ++l) {
    int size = blockStruct[l];
    if (size < 0) {
      // case Diagonal
      #if 1 // ????????? bag in idamax_ ?
      double* target = dualMat.ele[l].di_ele;
      for (int j=0; j<-size; ++j) {
	tmp = fabs(target[j]);
	if (tmp > ret) {
	  ret = tmp;
	}
      }
      #else
      int length = -size;
      int index = idamax_(&length,dualMat.ele[l].di_ele,&IONE);
      tmp = fabs(dualMat.ele[l].di_ele[index]);
      if (tmp > ret) {
	ret = tmp;
      }
      #endif
    } else {
      // case non Diagonal
      #if 1 // ????????????? bag in idamax_ ?
      double* target = dualMat.ele[l].de_ele;
      for (int j=0; j<size*size; ++j) {
	tmp = fabs(target[j]);
	if (tmp > ret) {
	  ret = tmp;
	}
      }
      #else
      int length = size*size;
      int index = idamax_(&length,dualMat.ele[l].de_ele,&IONE);
      tmp = fabs(dualMat.ele[l].de_ele[index]);
      if (tmp > ret) {
	ret = tmp;
      }
      #endif
      // rMessage("ret = " << ret);
    } // end of 'if (size < 0)'
  }
  return ret;
}

void rResiduals::update(int m, int nBlock, int* blockStruct,
			rVector& b, rBlockSparseMatrix& C,
			rBlockSparseMatrix* A,
			rResiduals& initResidual,
			rRatioInitResCurrentRes& theta,
			rSolutions& currentPt,
			rPhase& phase,
			rAverageComplementarity& mu,
			rComputeTime& com)
{
  rTimeStart(UPDATE_START);
  #if 0
  rAl::let(primalVec,'=',initResidual.primalVec,
	   '*',&theta.primal);
  normPrimalVec = initResidual.normPrimalVec * theta.primal;
  #else
  for (int k=0; k<m; ++k) {
    double ip;
    rAl::let(ip,'=',A[k],'.',currentPt.xMat);
    primalVec.ele[k] = b.ele[k] - ip;
  }
  normPrimalVec = computeMaxNorm(primalVec);
  #endif
  
  #if 0
  rAl::let(dualMat,'=',initResidual.dualMat,
	   '*',&theta.dual);
  normDualMat   = initResidual.normDualMat   * theta.dual;
  #else
  rAl::let(dualMat,'=',C,'+',currentPt.zMat,&DMONE);
  for (int k=0; k<m; ++k) {
    rAl::let(dualMat,'=',dualMat,'-',A[k],&currentPt.yVec.ele[k]);
  }
  normDualMat   = computeMaxNorm(dualMat);
  #endif
  rTimeStart(UPDATE_END);
  com.updateRes += rTimeCal(UPDATE_START,UPDATE_END);

  // rMessage("mu.current =" << mu.current);
  // rMessage("currentPt.xzMinEigenValue ="
  // << currentPt.xzMinEigenValue);
  centerNorm = 1.0 - currentPt.xzMinEigenValue/mu.current;
}

void rResiduals::compute(int m, int nBlock, int* blockStruct,
			 rVector& b, rBlockSparseMatrix& C,
			 rBlockSparseMatrix* A,
			 rSolutions& currentPt,
			 rAverageComplementarity& mu)
{
  initialize(m,nBlock,blockStruct,b,C,A,currentPt);
  // normPrimalVec = computeMaxNorm(primalVec);
  // normDualMat   = computeMaxNorm(dualMat);
  // centerNorm = 1.0 - currentPt.xzMinEigenValue/mu.current;
}

void rResiduals::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }
  fprintf(fpout," currentRes.primalVec = \n");
  primalVec.display(fpout);
  fprintf(fpout," currentRes.dualMat = \n");
  dualMat.display(fpout);

  fprintf(fpout," currentRes.normPrimalVec = %8.3e\n",
	  normPrimalVec);
  fprintf(fpout," currentRes.normDualMat = %8.3e\n",
	  normDualMat);
}

//-----------------------------------------------------------

rStepLength::rStepLength()
{
  primal = 0.0;
  dual   = 0.0;

}

rStepLength::rStepLength(double alphaP, double alphaD, int nBlock,
			 int* blockStruct)
{
  initialize(alphaP,alphaD,nBlock,blockStruct);
}

rStepLength::~rStepLength()
{
  workMat1.~rBlockDenseMatrix();
  workMat2.~rBlockDenseMatrix();
  xInvDxEigenValues.~rBlockVector();
  zInvDzEigenValues.~rBlockVector();
  workVec.~rBlockVector();
}

void rStepLength::initialize(double alphaP, double alphaD,
			     int nBlock, int* blockStruct)
{
  primal = alphaP;
  dual   = alphaD;
  workMat1.initialize(nBlock,blockStruct);
  workMat2.initialize(nBlock,blockStruct);
  xInvDxEigenValues.initialize(nBlock,blockStruct);
  zInvDzEigenValues.initialize(nBlock,blockStruct);

  int* workBlockStruct = NULL;
  rNewCheck();
  workBlockStruct = new int[nBlock];
  if (workBlockStruct==NULL) {
    rError("rStepLength:: memory exhausted");
  }
  for (int l=0; l<nBlock; ++l) {
    if (blockStruct[l]>0) {
      workBlockStruct[l] = 3*blockStruct[l]-1;
    } else {
      workBlockStruct[l] = -3*blockStruct[l]-1;
    }
  }
  workVec.initialize(nBlock,workBlockStruct);
  delete[] workBlockStruct;
  workBlockStruct = NULL;
}

double rStepLength::minBlockVector(rBlockVector& aVec)
{
  int nBlock = aVec.nBlock;
  double ret = aVec.ele[0].ele[0];
  double tmp;
  int size = aVec.ele[0].nDim;
  for (int j=1; j<size; ++j) {
    tmp = aVec.ele[0].ele[j];
    if (tmp < ret) {
      ret = tmp;
    }
  }
  for (int k=1; k<nBlock; ++k) {
    size = aVec.ele[k].nDim;
    for (int j=0; j<size; ++j) {
      tmp = aVec.ele[k].ele[j];
      if (tmp < ret) {
	ret = tmp;
      }
    }
  }
  return ret;
}

void rStepLength::MehrotraPredictor(rVector& b,
				    rBlockSparseMatrix& C,
				    rBlockSparseMatrix* A,
				    rSolutions& currentPt,
				    rPhase& phase,
				    rNewton& newton,
				    rLanczos& lanczos,
				    rComputeTime& com)
{
  #if 1
  primal = dual = 0.9;
  #else
  double alphaBD = 100.0;
  double xi      = 3.0;

  // calculate  eigenvalues of X^{-1} dX
  rTimeStart(START1);
  // rMessage("invCholeskyX=");
  // currentPt.invCholeskyX.display();
  // rMessage("Dx=");
  // newton.DxMat.display();
  double minxInvDxEigenValue;
  #define ALL_EIGEN 0
  #if ALL_EIGEN
  rAl::let(workMat2,'=',newton.DxMat,'T',currentPt.invCholeskyX);
  rAl::let(workMat1,'=',currentPt.invCholeskyX,'*',workMat2);
  // rMessage("LDxLT=");
  // workMat1.display();
  // rMessage("eigen test");
  rAl::getMinEigenValue(workMat1,xInvDxEigenValues,workVec);
  minxInvDxEigenValue = minBlockVector(xInvDxEigenValues);
  // rMessage("minxInvDxEigenValue = " << minxInvDxEigenValue);
  #else
  minxInvDxEigenValue
    = lanczos.getMinEigen(currentPt.invCholeskyX,
			  newton.DxMat);
  // rMessage("minxInvDxEigenValue = " << minxInvDxEigenValue);
  #endif
  if (-minxInvDxEigenValue > 1.0 /alphaBD) {
    primal = - 1.0/minxInvDxEigenValue;
    // the limint of primal steplength
  } else {
    primal = alphaBD;
  }
  rTimeEnd(END1);
  com.EigxMatTime += rTimeCal(START1,END1);

  
  // calculate  eigenvalues of Z^{-1} dZ
  rTimeStart(START2);
  // rMessage("invCholeskyZ=");
  // currentPt.invCholeskyZ.display();
  // rMessage("Dz=");
  // newton.DzMat.display();

  double minzInvDzEigenValue;
  #if ALL_EIGEN
  rAl::let(workMat2,'=',newton.DzMat,'T',currentPt.invCholeskyZ);
  rAl::let(workMat1,'=',currentPt.invCholeskyZ,'*',workMat2);
  // rMessage("LDzLT=");
  // workMat1.display();
  // rMessage("eigen test");
  rAl::getMinEigenValue(workMat1,zInvDzEigenValues,workVec);
  minzInvDzEigenValue = minBlockVector(zInvDzEigenValues);
  // rMessage("minzInvDzEigenValue = " << minzInvDzEigenValue);
  #else
  minzInvDzEigenValue
    = lanczos.getMinEigen(currentPt.invCholeskyZ,
			  newton.DzMat);
  // rMessage("minzInvDzEigenValue = " << minzInvDzEigenValue);
  #endif
  if (-minzInvDzEigenValue > 1.0 /alphaBD) {
    dual = - 1.0/minzInvDzEigenValue;
    // the limint of dual steplength
  } else {
    dual = alphaBD;
  }
  rTimeEnd(END2);
  com.EigzMatTime += rTimeCal(START2,END2);
  
  #endif

  if (phase.value==rSolveInfo::noINFO
      || phase.value==rSolveInfo::dFEAS) {
    // primal is infeasible
    if (primal>1.0) {
      primal = 1.0;
    }
  } else {
    // when primal is feasible,
    // check stepP1 is effective or not.
    double incPrimalObj;
    rAl::let(incPrimalObj,'=',C,'.',newton.DxMat);
    if (incPrimalObj>0.0) {
      if (primal>dual) {
	primal = dual;
      }
      if (primal>1.0) {
	primal = 1.0;
      }
    }
  }
  if (phase.value==rSolveInfo::noINFO
      || phase.value==rSolveInfo::pFEAS) {
    // dual is infeasible
    if (dual>1.0) {
      dual = 1.0;
    }
  } else {
    // when dual is feasible
    // check stepD1 is effective or not.
    double incDualObj;
    rAl::let(incDualObj,'=',b,'.',newton.DyVec);
    if(incDualObj<0.0) {
      if (dual>primal) {
	dual = primal;
      }
      if (dual>1.0) {
	dual = 1.0;
      }
    }
  }
}

void rStepLength::MehrotraCorrector(int nDim, rVector& b,
				    rBlockSparseMatrix& C,
				    rBlockSparseMatrix* A,
				    rSolutions& currentPt,
				    rPhase& phase,
				    rSwitch& reduction,
				    rNewton& newton,
				    rAverageComplementarity& mu,
				    rRatioInitResCurrentRes& theta,
				    rLanczos& lanczos,
				    rParameter& param,
				    rComputeTime& com)
{
  double alphaBD = 100.0;
  double xi      = 3.0;

  // calculate  eigenvalues of X^{-1} dX
  rTimeStart(START1);
  // rMessage("invCholeskyX=");
  // currentPt.invCholeskyX.display();
  // rMessage("Dx=");
  // newton.DxMat.display();
  double minxInvDxEigenValue;
  #define ALL_EIGEN 0
  #if ALL_EIGEN
  rAl::let(workMat2,'=',newton.DxMat,'T',currentPt.invCholeskyX);
  rAl::let(workMat1,'=',currentPt.invCholeskyX,'*',workMat2);
  // rMessage("LDxLT=");
  // workMat1.display();
  // rMessage("eigen test");
  rAl::getMinEigenValue(workMat1,xInvDxEigenValues,workVec);
  minxInvDxEigenValue = minBlockVector(xInvDxEigenValues);
  // rMessage("minxInvDxEigenValue = " << minxInvDxEigenValue);
  #else
  minxInvDxEigenValue
    = lanczos.getMinEigen(currentPt.invCholeskyX,
			  newton.DxMat);
  // rMessage("minxInvDxEigenValue = " << minxInvDxEigenValue);
  #endif
  if (-minxInvDxEigenValue > 1.0 /alphaBD) {
    primal = - 1.0/minxInvDxEigenValue;
    // the limint of primal steplength
  } else {
    primal = alphaBD;
  }
  rTimeEnd(END1);
  com.EigxMatTime += rTimeCal(START1,END1);

  
  // calculate  eigenvalues of Z^{-1} dZ
  rTimeStart(START2);
  // rMessage("invCholeskyZ=");
  // currentPt.invCholeskyZ.display();
  // rMessage("Dz=");
  // newton.DzMat.display();

  double minzInvDzEigenValue;
  #if ALL_EIGEN
  rAl::let(workMat2,'=',newton.DzMat,'T',currentPt.invCholeskyZ);
  rAl::let(workMat1,'=',currentPt.invCholeskyZ,'*',workMat2);
  // rMessage("LDzLT=");
  // workMat1.display();
  // rMessage("eigen test");
  rAl::getMinEigenValue(workMat1,zInvDzEigenValues,workVec);
  minzInvDzEigenValue = minBlockVector(zInvDzEigenValues);
  // rMessage("minzInvDzEigenValue = " << minzInvDzEigenValue);
  #else
  minzInvDzEigenValue
    = lanczos.getMinEigen(currentPt.invCholeskyZ,
			  newton.DzMat);
  // rMessage("minzInvDzEigenValue = " << minzInvDzEigenValue);
  #endif
  if (-minzInvDzEigenValue > 1.0 /alphaBD) {
    dual = - 1.0/minzInvDzEigenValue;
    // the limint of dual steplength
  } else {
    dual = alphaBD;
  }
  rTimeEnd(END2);
  com.EigzMatTime += rTimeCal(START2,END2);

  // adjust steplength with param.gammaStar
  // param.gammaStar = 0.5;
  primal = param.gammaStar * primal;
  dual   = param.gammaStar * dual;
  // rMessage("primal = " << primal);
  // rMessage("dual = " << dual);
  // rMessage("phase = ");
  // phase.display();
  
  if (phase.value==rSolveInfo::noINFO
      || phase.value==rSolveInfo::dFEAS) {
    // primal is infeasible.
    if (primal>1.0) {
      primal = 1.0;
    }
  } else {
    double incPrimalObj;
    rAl::let(incPrimalObj,'=',C,'.',newton.DxMat);
    if(incPrimalObj>0.0) {
      // when primal is feasible
      // check stepD1 is effective or not.
      if (primal>dual) {
	primal = dual;
      }
      if (primal>1.0) {
	primal = 1.0;
      }
    }
  }
  if (phase.value==rSolveInfo::noINFO
      || phase.value==rSolveInfo::pFEAS) {
    // dual is infeasible
    if (dual>1.0) {
      dual = 1.0;
    }
  } else {
    // when dual is feasible
    // check stepD1 is effective or not.
    double incDualObj;
    rAl::let(incDualObj,'=',b,'.',newton.DyVec);
    if(incDualObj<0.0) {
      if (dual>primal) {
	// change because noneffective
	dual = primal;
      }
      if (dual>1.0) {
	dual = 1.0;
      }
    }
  }
#if 1
  // attain feasibility before mu reduction
  if (reduction.switchType==rSwitch::ON
      && (phase.value == rSolveInfo::noINFO
	  || phase.value == rSolveInfo::pFEAS
	  || phase.value == rSolveInfo::dFEAS) ) {
    double xMatvMat;
    rAl::let(xMatvMat,'=',currentPt.xMat,'.',newton.DzMat);
    double uMatzMat;
    rAl::let(uMatzMat,'=',newton.DxMat,'.',currentPt.zMat);
    double uMatvMat;
    rAl::let(uMatvMat,'=',newton.DxMat,'.',newton.DzMat);

    double thetaMax = max((1.0-primal)*theta.primal,
			  (1.0-dual  )*theta.dual);
    double muNew = mu.current
      + (primal*uMatzMat + dual*xMatvMat
	 + primal*dual*uMatvMat) / nDim;
    double alphaMax;
    while (thetaMax*mu.initial > xi*muNew) {
      alphaMax = 0.95 * max(primal,dual);
      primal = min(primal,alphaMax);
      dual   = min(dual  ,alphaMax);
      thetaMax = max((1.0-primal)*theta.primal,
		     (1.0-dual  )*theta.dual);
      muNew = mu.current + (primal*uMatzMat + dual*xMatvMat
			    + primal*dual*uMatvMat) / nDim;
      // if "too short step", then break down the algorithm.
      if (primal < 1.0e-6 && dual < 1.0e-6) {
	break;
      }
    }
  }
#endif
}

void rStepLength::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }

  fprintf(fpout,"alpha.primal = %8.3e\n",primal);
  fprintf(fpout,"alpha.dual   = %8.3e\n",dual);
}

//-------------------------------------------------
rDirectionParameter::rDirectionParameter(double betaStar)
{
  initialize(betaStar);
}

rDirectionParameter::~rDirectionParameter()
{
  // Nothing needs.
}

void rDirectionParameter::initialize(double betaStar)
{
  value = betaStar;
}

void rDirectionParameter::MehrotraPredictor(rPhase& phase,
					    rSwitch& reduction,
					    rParameter& param)
{
  const double nu = 2.0;
  if (phase.value == rSolveInfo::pdFEAS) {
    value = 0.0;
  } else {
    value = param.betaBar;
    if (reduction.switchType==rSwitch::OFF) {
      value = nu;
    }
  }
}

void rDirectionParameter::
MehrotraCorrector(int nDim,rPhase& phase,rStepLength& alpha,
		  rSolutions& currentPt,rNewton& newton,
		  rAverageComplementarity& mu,rParameter& param)
{
  double xMatvMat;
  rAl::let(xMatvMat,'=',currentPt.xMat,'.',newton.DzMat);
  double uMatzMat;
  rAl::let(uMatzMat,'=',newton.DxMat,'.',currentPt.zMat);
  double uMatvMat;
  rAl::let(uMatvMat,'=',newton.DxMat,'.',newton.DzMat);

  double muTarget = mu.current
    + (alpha.primal*uMatzMat + alpha.dual*xMatvMat
       + alpha.primal*alpha.dual*uMatvMat) / nDim;
  // rMessage("muTarget : " << muTarget);
  value = muTarget/mu.current;
  // rMessage("muValue : " << value);
  if (value < 1.0) {
    value = value*value;
  }
  if (phase.value==rSolveInfo::pdFEAS) {
    // rMessage("MehrotraCorrector : pdFEAS" << value);
    if (value < param.betaStar) {
      value = param.betaStar;
    }
    if (value > 1.0) {
      value = 1.0;
    }
  } else {
    if (value < param.betaBar) {
      value = param.betaBar;
    }
  }
  // rMessage("MehrotraCorrector : " << value);
}

void rDirectionParameter::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }
  fprintf(fpout,"beta.value = %8.3e\n",value);
}

//---------------------------------------------------

rSwitch::rSwitch(rSwitchType switchType)
{
  initialize(switchType);
}

rSwitch::~rSwitch()
{
  // Nothing needs.
}

void rSwitch::initialize(rSwitchType switchType)
{
  this->switchType = switchType;
}

void rSwitch::MehrotraPredictor(rPhase& phase)
{
  if (phase.value==rSolveInfo::noINFO
      || phase.value==rSolveInfo::pFEAS
      || phase.value==rSolveInfo::dFEAS) {
    // At least one of primal or dual is infeasible.
    switchType = ON;
  } else {
    switchType = OFF;
  }
}

void rSwitch::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }

  if (switchType == ON) {
    fprintf(fpout,"reduction.switchType == ON\n");
  } else {
    fprintf(fpout,"reduction.switchType == OFF\n");
  }
}

// ----------------------------------------

rAverageComplementarity::rAverageComplementarity(double lambdaStar)
{
  initialize(lambdaStar);
}

rAverageComplementarity::~rAverageComplementarity()
{
  // Nothing needs.
}

void rAverageComplementarity::initialize(double lambdaStar)
{
  initial = lambdaStar*lambdaStar;
  current = initial;
  // rMessage("initial average = " << initial);
}

void rAverageComplementarity::initialize(int nDim,
					 rSolutions& initPt)
{
  rAl::let(initial,'=',initPt.xMat,'.',initPt.zMat);
  initial /= nDim;
  current = initial;
}

void rAverageComplementarity::update(int nDim,
				     rSolutions& currentPt)
{
  rAl::let(current,'=',currentPt.xMat,'.',currentPt.zMat);
  current /= nDim;
}

void rAverageComplementarity::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }

  fprintf(fpout,"mu0 = %8.3e\n",initial);
  fprintf(fpout,"mu  = %8.3e\n",current);
}

//--------------------------------------------------


rRatioInitResCurrentRes::rRatioInitResCurrentRes()
{
  primal = 0.0;
  dual   = 0.0;
}

rRatioInitResCurrentRes::~rRatioInitResCurrentRes()
{
  // Nothing needs.
}

rRatioInitResCurrentRes::rRatioInitResCurrentRes(rParameter& param,
						 rResiduals& initRes)
{
  initialize(param,initRes);
}

void rRatioInitResCurrentRes::initialize(rParameter& param,
					 rResiduals& initRes)
{
  double accuracy = param.epsilonDash;
  if (initRes.normPrimalVec < accuracy) {
    primal = 0.0;
  } else {
    primal = 1.0;
  }
  if (initRes.normDualMat < accuracy) {
    dual = 0.0;
  } else {
    dual = 1.0;
  }
}

void rRatioInitResCurrentRes::update(rSwitch& reduction,
				     rStepLength& alpha)
{
  if (reduction.switchType==rSwitch::ON) {
    // At least one of primal or dual is infeasible
    primal = fabs((1.0-alpha.primal)*primal);
    dual   = fabs((1.0-alpha.dual  )*dual  );
  }
}

void rRatioInitResCurrentRes::update_exact(rResiduals& initRes,
					   rResiduals& currentRes)
{
  if (initRes.normPrimalVec > 1.0e-10) {
    primal = currentRes.normPrimalVec / initRes.normPrimalVec;
  }
  else {
    primal = 0.0;
  }
  if (initRes.normDualMat > 1.0e-10) {
    dual = currentRes.normDualMat / initRes.normDualMat;
  }
  else {
    dual = 0.0;
  }
}

void rRatioInitResCurrentRes::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }

  fprintf(fpout,"theta.primal = %8.3e\n",primal);
  fprintf(fpout,"theta.dual   = %8.3e\n",dual);
}

//---------------------------------------------------

rSolveInfo::rSolveInfo()
{
  rho = 0.0;
  etaPrimal = 0.0;
  etaDual   = 0.0;
  objValPrimal = 0.0;
  objValDual   = 0.0;
}

rSolveInfo::rSolveInfo(int nDim, rVector& b, rBlockSparseMatrix& C,
		       rBlockSparseMatrix* A, rSolutions& initPt,
		       double mu0, double omegaStar)
{
  initialize(nDim,b,C,A,initPt,mu0,omegaStar);
}

rSolveInfo::~rSolveInfo()
{
  // Nothing needs.
}

void rSolveInfo::initialize(int nDim, rVector& b,
			    rBlockSparseMatrix& C,
			    rBlockSparseMatrix* A,
			    rSolutions& initPt,
			    double mu0, double omegaStar)
{
  rho = 1.0;
  etaPrimal = omegaStar * nDim * mu0;
  etaDual   = omegaStar * nDim * mu0;
  rAl::let(objValPrimal,'=',C,'.',initPt.xMat);
  rAl::let(objValDual  ,'=',b,'.',initPt.yVec);
}

void rSolveInfo::update(int nDim, rVector& b, rBlockSparseMatrix& C,
			rSolutions& initPt, rSolutions& currentPt,
			rResiduals& currentRes,
			rAverageComplementarity& mu,
			rRatioInitResCurrentRes& theta,
			rParameter& param)
{
  rAl::let(objValPrimal,'=',C,'.',currentPt.xMat);
  rAl::let(objValDual  ,'=',b,'.',currentPt.yVec);
  double primal = theta.primal;
  double dual   = theta.dual;
  double omega  = param.omegaStar;
  rho = 0.0;
  double x0z0     = nDim*mu.initial;
  double xMatzMat = nDim*mu.current;
  double x0zMat   = 0.0;
  double xMatz0   = 0.0;
  rAl::let(x0zMat,'=',initPt.xMat,'.',currentPt.zMat);
  rAl::let(xMatz0,'=',currentPt.xMat,'.',initPt.zMat);

  double accuracy = param.epsilonDash;

  if (currentRes.normPrimalVec <= accuracy) {
    // rMessage("primal accuracy");
    if (xMatz0 < etaPrimal) {
      etaPrimal = xMatz0;
    }
  }
  if (currentRes.normDualMat <= accuracy) {
    // rMessage("dual accuracy");
    if (x0zMat < etaDual) {
      etaDual = x0zMat;
    }
  }

  // primal is infeasible and dual is feasible
  if (currentRes.normPrimalVec > accuracy
      && currentRes.normDualMat <= accuracy) {
    rho = primal*x0zMat
      / ((primal+(1.0-primal)*omega)*etaDual + xMatzMat);
  }

  // primal is feasible and dual is infeasible
  if (currentRes.normPrimalVec <= accuracy
      && currentRes.normDualMat > accuracy) {
    rho = dual*xMatz0
      / ((dual+(1.0-dual)*omega)* etaPrimal + xMatzMat);
  }
  
  // primal and dual are infeasible
  if (currentRes.normPrimalVec > accuracy
      && currentRes.normDualMat > accuracy) {
    rho = (dual*xMatz0+primal*x0zMat)
      / ((primal*dual
	  + omega *(primal*(1.0-dual) + (1.0-primal)*dual))* x0z0
	 + xMatzMat);
  }
  // rMessage("eta Primal = " << etaPrimal);
  // rMessage("eta Dual = " << etaDual);
}

void rSolveInfo::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }

  fprintf(fpout,"rSolveInfo.rho          = %8.3e\n",rho);
  fprintf(fpout,"rSolveInfo.etaPrimal    = %8.3e\n",etaPrimal);
  fprintf(fpout,"rSolveInfo.etaDual      = %8.3e\n",etaDual);
  fprintf(fpout,"rSolveInfo.objValPrimal = %8.3e\n",objValPrimal);
  fprintf(fpout,"rSolveInfo.objValDual   = %8.3e\n",objValDual);
}

// ----------------------------------------------------

rPhase::rPhase()
{
  nDim = 0;
  value = rSolveInfo::noINFO;
}

rPhase::rPhase(rResiduals& initRes, rSolveInfo& solveInfo,
	       rParameter& param, int nDim)
{
  initialize(initRes,solveInfo,param,nDim);
}

rPhase::~rPhase()
{
  // Nothing needs.
}

bool rPhase::initialize(rResiduals& initRes,
			rSolveInfo& solveInfo,
			rParameter& param, int nDim)
{
  this->nDim = nDim;
  return updateCheck(initRes,solveInfo,param);
}

bool rPhase::updateCheck(rResiduals& currentRes,
			  rSolveInfo& solveInfo,
			  rParameter& param)
{
  const double NONZERO = 1.0e-6;
  double accuracy = param.epsilonDash;
  value = rSolveInfo::noINFO;

  if (currentRes.normPrimalVec <= accuracy) {
    if (currentRes.normDualMat <= accuracy) {
      value = rSolveInfo::pdFEAS;
    } else {
      value = rSolveInfo::pFEAS;
    }
  }
  if (value==rSolveInfo::noINFO
      && currentRes.normDualMat <= accuracy) {
    value = rSolveInfo::dFEAS;
  }
  if (value==rSolveInfo::pdFEAS) {
    double mean = (fabs(solveInfo.objValPrimal)+
		   fabs(solveInfo.objValDual)) / 2.0;
    double PDgap = fabs(solveInfo.objValPrimal
			- solveInfo.objValDual);

    double dominator;
    if (mean < 1.0) {
      dominator = 1.0;
    } else {
      dominator = mean;
    }
    #if 0
    rMessage("PDgap = " << PDgap);
    rMessage("dominator = " << dominator);
    rMessage("PDgap/dominator = " << PDgap/dominator);
    #endif
    if (PDgap/dominator <= param.epsilonStar) {
      value = rSolveInfo::pdOPT;
      return false;
    }
  }
  if (value == rSolveInfo::noINFO
      && solveInfo.rho > 1.0+NONZERO) {
    value = rSolveInfo::pdINF;
    return false;
  }
  if (value == rSolveInfo::pFEAS) {
    #if REVERSE_PRIMAL_DUAL
    if (solveInfo.objValPrimal<=-param.upperBound) {
      value = rSolveInfo::pUNBD;
      return false;
    }
    #else
    if (solveInfo.objValPrimal<=param.lowerBound) {
      value = rSolveInfo::pUNBD;
      return false;
    }
    #endif
    if (solveInfo.rho > 1.0+NONZERO) {
      value = rSolveInfo::pFEAS_dINF;
      return false;
    }
  }

  if (value == rSolveInfo::dFEAS) {
    #if REVERSE_PRIMAL_DUAL
    if (solveInfo.objValDual>=-param.lowerBound) {
      value = rSolveInfo::dUNBD;
      return false;
    }
    #else
    if (solveInfo.objValDual>=param.upperBound) {
      value = rSolveInfo::dUNBD;
      return false;
    }
    #endif
    if (solveInfo.rho > 1.0+NONZERO) {
      value = rSolveInfo::pINF_dFEAS;
      return false;
    }
  }
  #if 0
  rMessage("phase =");
  display();
  #endif
  return true;
}

void rPhase::reverse()
{
  #if REVERSE_PRIMAL_DUAL
  switch (value) {
  case rSolveInfo::noINFO    :                               ; break;
  case rSolveInfo::pFEAS     : value = rSolveInfo::dFEAS     ; break;
  case rSolveInfo::dFEAS     : value = rSolveInfo::pFEAS     ; break;
  case rSolveInfo::pdFEAS    :                               ; break;
  case rSolveInfo::pdINF     :                               ; break;
  case rSolveInfo::pFEAS_dINF: value = rSolveInfo::pINF_dFEAS; break;
  case rSolveInfo::pINF_dFEAS: value = rSolveInfo::pFEAS_dINF; break;
  case rSolveInfo::pdOPT     :                               ; break;
  case rSolveInfo::pUNBD     : value = rSolveInfo::dUNBD     ; break;
  case rSolveInfo::dUNBD     : value = rSolveInfo::pUNBD     ; break;
  default: break;
  }
  #else
  // do nothing
  #endif
}


void rPhase::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }
  char* str;
  switch (value) {
  case rSolveInfo::noINFO    : str = "noINFO    "; break;
  case rSolveInfo::pFEAS     : str = "pFEAS     "; break;
  case rSolveInfo::dFEAS     : str = "dFEAS     "; break;
  case rSolveInfo::pdFEAS    : str = "pdFEAS    "; break;
  case rSolveInfo::pdINF     : str = "pdINF     "; break;
  case rSolveInfo::pFEAS_dINF: str = "pFEAS_dINF"; break;
  case rSolveInfo::pINF_dFEAS: str = "pINF_dFEAS"; break;
  case rSolveInfo::pdOPT     : str = "pdOPT     "; break;
  case rSolveInfo::pUNBD     : str = "pUNBD     "; break;
  case rSolveInfo::dUNBD     : str = "dUNBD     "; break;
  default:
    str = "phase error";
    rMessage("rPhase:: phase error");
    break;
  }
  fprintf(fpout,"phase.value = %s\n",str);
}

//---------------------------------------------------------

rLanczos::rLanczos()
{
  // Nothing needs.
}

rLanczos::rLanczos(int nBlock,int* blockStruct)
{
  initialize(nBlock,blockStruct);
}

void rLanczos::initialize(int nBlock,int* blockStruct)
{
  Q.initialize(nBlock,blockStruct);
  A.initialize(nBlock,blockStruct);
  out.initialize(nBlock,blockStruct);
  b.initialize(nBlock,blockStruct);
  r.initialize(nBlock,blockStruct);
  q.initialize(nBlock,blockStruct);
  qold.initialize(nBlock,blockStruct);
  w.initialize(nBlock,blockStruct);
  tmp.initialize(nBlock,blockStruct);
  diagVec.initialize(nBlock,blockStruct);
  diagVec2.initialize(nBlock,blockStruct);

  int* workStruct = NULL;
  rNewCheck();
  workStruct = new int[nBlock];
  if (workStruct == NULL) {
    rMessage("rLanczos :: memory exhausted");
  }
  for (int l=0; l<nBlock; ++l) {
    workStruct[l] = max(1,2*blockStruct[l]-2);
  }
  workVec.initialize(nBlock,workStruct);
  delete[] workStruct;
  workStruct = NULL;
}

rLanczos::~rLanczos()
{
  Q.~rBlockDenseMatrix();
  A.~rBlockDenseMatrix();
  out.~rBlockVector();
  b.~rBlockVector();
  r.~rBlockVector();
  q.~rBlockVector();
  qold.~rBlockVector();
  w.~rBlockVector();
  tmp.~rBlockVector();
  diagVec.~rBlockVector();
  diagVec2.~rBlockVector();
  workVec.~rBlockVector();
}

double rLanczos::getMinEigen(rBlockDenseMatrix& lMat,
			     rBlockDenseMatrix& xMat)
{
  const int nBlock = xMat.nBlock;
  const int* blockStruct = xMat.blockStruct;
  double min;
  if (blockStruct[0] < 0) {
    min = getMinEigen(lMat.ele[0],xMat.ele[0]);
  } else {
    min = getMinEigen(lMat.ele[0],xMat.ele[0],Q.ele[0],
		      out.ele[0],b.ele[0],r.ele[0],q.ele[0],
		      qold.ele[0],w.ele[0],tmp.ele[0],
		      diagVec.ele[0],diagVec2.ele[0],
		      workVec.ele[0]);
  }
  double value;
  for (int l=1; l<nBlock; ++l) {
    if (blockStruct[l] < 0) {
      value = getMinEigen(lMat.ele[l],xMat.ele[l]);
    } else {
      value = getMinEigen(lMat.ele[l],xMat.ele[l],Q.ele[l],
			  out.ele[l],b.ele[l],r.ele[l],q.ele[l],
			  qold.ele[l],w.ele[l],tmp.ele[l],
			  diagVec.ele[l],diagVec2.ele[l],
			  workVec.ele[l]);
    }
    if (value < min) {
      min = value;
    }
  } // end of 'for (int l)'

  return min;
}

double rLanczos::getMinEigen(rDenseMatrix& lMat,
			     rDenseMatrix& xMat,
			     rDenseMatrix& Q,
			     rVector& out, rVector& b, rVector& r,
			     rVector& q, rVector& qold,
			     rVector& w, rVector& tmp,
			     rVector& diagVec, rVector& diagVec2,
			     rVector& workVec)
{
  double alpha,beta,value;
  double min = 1.0e+51, min_old = 1.0e+52;
  double error = 1.0e+10;

  int nDim = xMat.nRow;
  int k = 0, kk = 0;
  
  diagVec.initialize(1.0e+50);
  diagVec2.setZero();
  q.setZero();
  r.initialize(1.0);
  beta = sqrt((double)nDim);  // norm of "r"

  while (k<nDim && k<sqrt((double)nDim)+10
	 && beta > 1.0e-16
	 && (fabs(min-min_old) > (1.0e-5)*fabs(min)+(1.0e-8)
	     // && (fabs(min-min_old) > (1.0e-3)*fabs(min)+(1.0e-6)
	     || fabs(error*beta) > (1.0e-2)*fabs(min)+(1.0e-4) )) {
    // rMessage("k = " << k);
    qold.copyFrom(q);
    value = 1.0/beta;
    rAl::let(q,'=',r,'*',&value);

    // w = (lMat^T)*q
    w.copyFrom(q);
    dtrmv("Lower","Transpose","NotUnit",&nDim,
	   lMat.de_ele,&nDim,w.ele,&IONE);
    rAl::let(tmp,'=',xMat,'*',w);
    w.copyFrom(tmp);
    dtrmv("Lower","NoTranspose","NotUnit",&nDim,
	   lMat.de_ele,&nDim,w.ele,&IONE);
    // w = lMat*xMat*(lMat^T)*q
    // rMessage("w = ");
    // w.display();
    rAl::let(alpha,'=',q,'.',w);
    diagVec.ele[k] = alpha;
    rAl::let(r,'=',w,'-',q,&alpha);
    rAl::let(r,'=',r,'-',qold,&beta);
    // rMessage("r = ");
    // r.display();

    if ( kk>=sqrt((double)k) || k==nDim-1 || k>sqrt((double)nDim+9) ) {
      kk = 0;
      out.copyFrom(diagVec);
      b.copyFrom(diagVec2);
      out.ele[nDim-1] = diagVec.ele[k];
      b.ele[nDim-1]   = 0.0;
      
      // rMessage("out = ");
      // out.display();
      // rMessage("b = ");
      // b.display();

      int info;
      int kp1 = k+1;
      dsteqr_("I_withEigenvalues",&kp1,out.ele,b.ele,
	      Q.de_ele, &Q.nRow, workVec.ele, &info);
      if (info < 0) {
	rError(" rLanczos :: bad argument " << -info
	       << " Q.nRow = " << Q.nRow
	       << ": nDim = " << nDim
	       << ": kp1 = " << kp1);
      } else if (info > 0) {
	rMessage(" rLanczos :: cannot converge " << info);
	break;
      }
      
      // rMessage("out = ");
      // out.display();
      // rMessage("Q = ");
      // Q.display();
      
      min_old = min;
      #if 0
      min = 1.0e+50;
      error = 1.0e+10;
      for (int i=0; i<k+1; ++i) {
	if (min>out.ele[i]){
	  min = out.ele[i];
	  error = Q.de_ele[k+Q.nCol*i];
	}
      }
      #else
      // out have eigen values with ascending order.
      min = out.ele[0];
      error = Q.de_ele[k];
      #endif

    } // end of 'if ( kk>=sqrt(k) ...)'
    // printf("\n");

    rAl::let(value,'=',r,'.',r);
    beta = sqrt(value);
    diagVec2.ele[k] = beta;
    ++k;
    ++kk;
  } // end of while
  // rMessage("k = " << k);
  return min - fabs(error*beta);
}

double rLanczos::getMinEigen(rDenseMatrix& lMat,
			     rDenseMatrix& xMat)
{
  // lMat, xMat is Diagonal
  const int nDim = xMat.nRow;
  const double* l_ele = lMat.di_ele;
  const double* x_ele = xMat.di_ele;

  double min = l_ele[0]*x_ele[0]*l_ele[0];
  double value;
  #if 0
  for (int j=1; j<nDim; ++j) {
    value = l_ele[j]*x_ele[j]*l_ele[j];
    if (value < min) {
      min = value;
    }
  }
  #else
  int shou = nDim / 4;
  int amari = nDim % 4;
  for (int j=1; j<amari; ++j) {
    value = l_ele[j]*x_ele[j]*l_ele[j];
    if (value < min) {
      min = value;
    }
  }
  int count;
  for (int j=amari,count=0; count < shou; ++count ,j+=4) {
    double value1 = l_ele[j]*x_ele[j]*l_ele[j];
    double value2 = l_ele[j+1]*x_ele[j+1]*l_ele[j+1];
    double value3 = l_ele[j+2]*x_ele[j+2]*l_ele[j+2];
    double value4 = l_ele[j+3]*x_ele[j+3]*l_ele[j+3];
    double tmp1,tmp2;
    if (value1 < value2) {
      tmp1 = value1;
    }
    else {
      tmp1 = value2;
    }
    if (value3 < value4) {
      tmp2 = value3;
    }
    else {
      tmp2 = value4;
    }
    if (tmp1 < tmp2) {
      if (tmp1<min) {
	min = tmp1;
      }
    }
    else if (tmp2<min) {
      min = tmp2;
    }
  }
  #endif
  return min;
}
