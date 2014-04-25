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
  rsdpa_lib.cpp
-------------------------------------------------*/

#include "rsdpa_io.h"
#include "rsdpa_lib.h"

#define KAPPA 2.2

rSdpaLib::rSdpaLib()
{
  DisplayInformation = NULL;
  m = 0;
  nBlock = 0;
  blockStruct = NULL;
  nonZeroNumber = 0;

  nDim = 0;
  A = NULL;
  hasSolved = false;
  hasDelete1 = false;
  setDefaultParameter();

  // for compability SDPA
  ParameterFile = NULL;
  OutputFile    = NULL;
  InitialFile   = NULL;
  InputFile     = NULL;
  CheckMatrix   = false;
  mDIM         = 0;
  nBLOCK       = 0;
  bLOCKsTRUCT  = NULL;
  InitialPoint = false;
}

rSdpaLib::~rSdpaLib()
{
  // delete information which are not relevant to solutions
  if (hasDelete1 == false) {
    delete1();
  }

  // delete information which are relevant to solutions.
  mu.~rAverageComplementarity();
  phase.~rPhase();
  solveInfo.~rSolveInfo();
  currentPt.~rSolutions();
  currentRes.~rResiduals();
  pARAM.~rParameter();
  com.~rComputeTime();
  if (blockStruct) {
    delete[] blockStruct;
    blockStruct = NULL;
  }
}

void rSdpaLib::delete1()
{
  // delete information which are not relevant to solutions
  b.~rVector();
  C.~rBlockSparseMatrix();
  if (A) {
    for (int k=0; k<m; ++k) {
      A[k].~rBlockSparseMatrix();
    }
    delete[] A;
    A = NULL;
  }
  
  initPt.~rSolutions();
  initRes.~rResiduals();

  newton.~rNewton();
  alpha.~rStepLength();
  beta.~rDirectionParameter();
  reduction.~rSwitch();
  lanczos.~rLanczos();
  theta.~rRatioInitResCurrentRes();
  hasDelete1 = true;
}

void rSdpaLib::delete2()
{
  // delete all from memory,
  // that is the same work as destructor
  this->~rSdpaLib();
}

void rSdpaLib::setDisplay(FILE* Display)
{
  this->DisplayInformation = Display;
}

void rSdpaLib::setDefaultParameter(parameterType type)
{
  if (type == ::PARAMETER_AGGRESSIVE) {
    pARAM.setDefaultParameter(rParameter::PARAMETER_AGGRESSIVE);
  }
  if (type == ::PARAMETER_STABLE) {
    pARAM.setDefaultParameter(rParameter::PARAMETER_STABLE);
  }
  else {
    pARAM.setDefaultParameter(rParameter::PARAMETER_DEFAULT);
  }
}

void rSdpaLib::setDefaultParameter(rParameterType type)
{
  if (type == rPARAMETER_AGGRESSIVE) {
    pARAM.setDefaultParameter(rParameter::PARAMETER_AGGRESSIVE);
  }
  if (type == rPARAMETER_STABLE) {
    pARAM.setDefaultParameter(rParameter::PARAMETER_STABLE);
  }
  else {
    pARAM.setDefaultParameter(rParameter::PARAMETER_DEFAULT);
  }
}

void rSdpaLib::setParameterMaxIteration(int maxIteration)
{
  pARAM.maxIteration = maxIteration;
}
void rSdpaLib::setParameterEpsilonStar(double epsilonStar)
{
  pARAM.epsilonStar = epsilonStar;
}
void rSdpaLib::setParameterLambdaStar(double lambdaStar)
{
  pARAM.lambdaStar = lambdaStar;
}
void rSdpaLib::setParameterOmegaStar(double omegaStar)
{
  pARAM.omegaStar = omegaStar;
}
void rSdpaLib::setParameterLowerBound(double lowerBound)
{
  pARAM.lowerBound = lowerBound;
}
void rSdpaLib::setParameterUpperBound(double upperBound)
{
  pARAM.upperBound = upperBound;
}
void rSdpaLib::setParameterBetaStar(double betaStar)
{
  pARAM.betaStar = betaStar;
}
void rSdpaLib::setParameterBetaBar(double betaBar)
{
  pARAM.betaBar = betaBar;
}
void rSdpaLib::setParameterGammaStar(double gammaStar)
{
  pARAM.gammaStar = gammaStar;
}
void rSdpaLib::setParameterEpsilonDash(double epsilonDash)
{
  pARAM.epsilonDash = epsilonDash;
}

void rSdpaLib::initialize1(int m, int nBlock, int* blockStruct,
			   bool initialPoint)
{
  time_t ltime;
  time( &ltime );
  char date_string[256];
  sprintf(date_string,"SDPA library start at %s", ctime(&ltime));
 if (DisplayInformation) {
    fprintf(DisplayInformation,date_string);
    fprintf(DisplayInformation,"let me see your ...\n");
    #if 0
    fprintf(DisplayInformation,"   mu      thetaP  thetaD  objP"
	    "      objD "
            "     alphaP  alphaD  beta \n");
    #endif
  }
  
  this->m = m;
  this->nBlock = nBlock;
  if (this->blockStruct) {
    delete[] this->blockStruct;
    this->blockStruct = NULL;
  }
  this->blockStruct = new int[nBlock];
  if (this->blockStruct==NULL) {
    rError("rSdpaLib::initialize1 memory exhausted");
  }
  nDim = 0;
  for (int l=0; l<nBlock; ++l) {
    this->blockStruct[l] = blockStruct[l];
    if (blockStruct[l] < 0 ){
      nDim -= blockStruct[l];
    } else {
      nDim += blockStruct[l];
    }
  }
  if (nonZeroNumber) {
    delete[] nonZeroNumber;
    nonZeroNumber = NULL;
  }
  nonZeroNumber = new int[(m+1)*nBlock];
  if (nonZeroNumber==NULL) {
    rError("rSdpaLib::initialize1 memory exhausted");
  }
  for (int j=0; j<(m+1)*nBlock; ++j) {
    nonZeroNumber[j] = 0;
  }

  if (initialPoint) {
    // initialize the initial solution with 0
    initPt.initializeZero(m,nBlock,blockStruct,com);
    InitialPoint= true;
  } else {
    initPt.initialize(m,nBlock,blockStruct,pARAM.lambdaStar,com);
    InitialPoint = false;
  }
  b.initialize(m);
  
}

void rSdpaLib::countUpperTriangle(int k, int l, int nonzero)
{
  k--;
  if (k<-1 || k>=m) {
    #if REVERSE_PRIMAL_DUAL
    rMessage("Over index of F:: " << k+1);
    rError("Length of F :: " << m);
    #else
    rMessage("Over index of A or C:: " << k+1);
    rError("Length of A :: " << m);
    #endif
  }
  l--;
  if (l<0 || l>=nBlock) {
    #if REVERSE_PRIMAL_DUAL
    rMessage("Over nBlock of F:: " << l+1);
    rError("nBlock of F :: " << nBlock);
    #else
    rMessage("Over nBlock of A or C:: " << l+1);
    rError("nBlock of A :: " << nBlock);
    #endif
  }
  nonZeroNumber[(k+1)*nBlock+l] = nonzero;
}

void rSdpaLib::initialize2()
{
  C.initialize(nBlock,blockStruct);
  for (int l=0; l<nBlock; ++l) {
    int size = blockStruct[l];
    if (size > 0) {
      C.ele[l].initialize(size,size,rSparseMatrix::SPARSE,
			  nonZeroNumber[l]);
      // rMessage("C.ele[l].NonZeroNumber ="
      //       << C.ele[l].NonZeroNumber);
      C.ele[l].NonZeroCount = 0;
    } else {
      C.ele[l].initialize(-size,-size,rSparseMatrix::DIAGONAL,
			  -size);
    }
  }
  if (A) {
    for (int k=0; k<m; ++k) {
      A[k].~rBlockSparseMatrix();
    }
    delete[] A;
  }
  A = NULL;
  A = new rBlockSparseMatrix[m];
  if (A==NULL) {
    rError("rSdpaLib::initialize2 memory exhausted");
  }
  for (int k=0; k<m; ++k) {
    A[k].initialize(nBlock,blockStruct);
    for (int l=0; l<nBlock; ++l) {
      int size = blockStruct[l];
      if (size > 0) {
	A[k].ele[l].initialize(size,size,rSparseMatrix::SPARSE,
			       nonZeroNumber[(k+1)*nBlock+l]);
	// rMessage("A["<<k<< "].ele["<<l<<"].NonZeroNumber ="
	//      << A[k].ele[l].NonZeroNumber);
	A[k].ele[l].NonZeroCount  = 0;
	A[k].ele[l].NonZeroEffect = 0;
      } else {
	A[k].ele[l].initialize(-size,-size,rSparseMatrix::DIAGONAL,
			       -size);
      }
    }
  }
  delete[] nonZeroNumber;
  nonZeroNumber = NULL;
}

#if REVERSE_PRIMAL_DUAL
void rSdpaLib::inputCVec(int k, double value)
{
  k--;
  if (k<0 || k>=m) {
    rMessage("Over index of c:: " << k+1);
    rError("Length of c:: " << b.nDim);
  }
  b.ele[k] = value;
}

#else //REVERSE_PRIMAL_DUAL

void rSdpaLib::inputBVec(int k, double value)
{
  k--;
  if (k<0 || k>=m) {
    rMessage("Over index of b:: " << k+1);
    rError("Length of b:: " << b.nDim);
  }
  b.ele[k] = value;
}

void rSdpaLib::inputCMat(int l, int i, int j, double value)
{
  l--;
  if (l<0 || l>=nBlock) {
    rMessage("Over nBlock of C:: " << l+1);
    rError("nBlock of C :: " << nBlock);
  }
  i--;
  j--;
  int size = blockStruct[l];
  rSparseMatrix& target = C.ele[l];
  int count = target.NonZeroCount;
  if (target.Sp_De_Di == rSparseMatrix::SPARSE
      && count >= target.NonZeroNumber) {
    rError("C.ele[" << (l+1) << "]"
	   << " is too much element which is assigned");
  }
  if (size>0) {
    target.row_index[count] = i;
    target.column_index[count] = j;
    target.sp_ele[count] = value;
    target.NonZeroCount++;
    if (i==j) {
      target.NonZeroEffect++;
    } else {
      target.NonZeroEffect += 2;
    }
  } else {
    target.di_ele[j] = value;
  }
}
#endif // REVERSE_PRIMAL_DUAL

void rSdpaLib::inputElement(int k,int l, int i, int j, double value)
{
  k--;
  if (k<-1 || k>=m) {
    #if REVERSE_PRIMAL_DUAL
    rMessage("Over index of F:: " << k+1);
    rError("Length of F :: " << m);
    #else
    rMessage("Over index of A or C:: " << k+1);
    rError("Length of A :: " << m);
    #endif
  }
  l--;
  if (l<0 || l>=nBlock) {
    #if REVERSE_PRIMAL_DUAL
    rMessage("Over nBlock of F:: " << l+1);
    rError("nBlock of F :: " << nBlock);
    #else
    rMessage("Over nBlock of A or C:: " << l+1);
    rError("nBlock of A :: " << nBlock);
    #endif
  }
  i--;
  j--;
  int size = blockStruct[l];
  if (k==-1) {
    rSparseMatrix& target = C.ele[l];
    int count = target.NonZeroCount;
    if (target.Sp_De_Di == rSparseMatrix::SPARSE
	&& count >= target.NonZeroNumber) {
      #if REVERSE_PRIMAL_DUAL
      rError("C.ele[" << (l+1) << "]"
	     << " is too much element which is assigned");
      #else
      rError("F[0].ele[" << (l+1) << "]"
	     << " is too much element which is assigned");
      #endif
    }
    if (size>0) {
      target.row_index[count] = i;
      target.column_index[count] = j;
      #if REVERSE_PRIMAL_DUAL
      target.sp_ele[count] = -value;
      #else
      target.sp_ele[count] = value;
      #endif
      target.NonZeroCount++;
      if (i==j) {
	target.NonZeroEffect++;
      } else {
	target.NonZeroEffect += 2;
      }
    } else {
      #if REVERSE_PRIMAL_DUAL
      target.di_ele[j] = -value;
      #else
      target.di_ele[j] = value;
      #endif
    }
  } else { // k>=0
    rSparseMatrix& target = A[k].ele[l];
    int count = target.NonZeroCount;
    if (target.Sp_De_Di == rSparseMatrix::SPARSE
	&& count >= target.NonZeroNumber) {
      #if REVERSE_PRIMAL_DUAL
      rError("F["<<(k+1)<<"].ele[" << (l+1) << "]"
	     << " is too much element which is assigned");
      #else
      rError("A["<<(k+1)<<"].ele[" << (l+1) << "]"
	     << " is too much element which is assigned");
      #endif
    }
    if (size>0) {
      target.row_index[count] = i;
      target.column_index[count] = j;
      target.sp_ele[count] = value;
      target.NonZeroCount++;
      if (i==j) {
	target.NonZeroEffect++;
      } else {
	target.NonZeroEffect += 2;
      }
    } else {
      target.di_ele[j] = value;
    }
  }
}

#if REVERSE_PRIMAL_DUAL


void rSdpaLib::inputInitXVec(int k, double value)
{
  k--;
  if (k<0 || k>=m) {
    rMessage("Over index of x:: " << k+1);
    rError("Length of x:: " << initPt.yVec.nDim);
  }
  initPt.yVec.ele[k] = -value;
}

void rSdpaLib::inputInitXMat(int l, int i, int j, double value)
{
  l--;
  if (l<0 || l>=nBlock) {
    rMessage("Over nBlock of X:: " << l+1);
    rError("nBlock of X :: " << initPt.zMat.nBlock);
  }
  i--;
  j--;
  int size = blockStruct[l];
  rDenseMatrix& target = initPt.zMat.ele[l];
  if (size>0) {
    if (i!=j) {
      target.de_ele[i+target.nRow*j] =value;
      target.de_ele[j+target.nRow*i] =value;
    } else {
      target.de_ele[i+target.nRow*i] =value;
    }
  } else {
    target.di_ele[j] = value;
  }
}

void rSdpaLib::inputInitYMat(int l, int i, int j, double value)
{
  l--;
  if (l<0 || l>=nBlock) {
    rMessage("Over nBlock of Y:: " << l+1);
    rError("nBlock of Y :: " << initPt.xMat.nBlock);
  }
  i--;
  j--;
  int size = blockStruct[l];
  rDenseMatrix& target = initPt.xMat.ele[l];
  if (size>0) {
    if (i!=j) {
      target.de_ele[i+target.nRow*j] =value;
      target.de_ele[j+target.nRow*i] =value;
    } else {
      target.de_ele[i+target.nRow*i] =value;
    }
  } else {
    target.di_ele[j] = value;
  }
}

#else //REVERSE_PRIMAL_DUAL

void rSdpaLib::inputInitXMat(int l, int i, int j, double value)
{
  l--;
  if (l<0 || l>=nBlock) {
    rMessage("Over nBlock of X:: " << l+1);
    rError("nBlock of X :: " << initPt.xMat.nBlock);
  }
  i--;
  j--;
  int size = blockStruct[l];
  rDenseMatrix& target = initPt.xMat.ele[l];
  if (size>0) {
    if (i!=j) {
      target.de_ele[i+target.nRow*j] =value;
      target.de_ele[j+target.nRow*i] =value;
    } else {
      target.de_ele[i+target.nRow*i] =value;
    }
  } else {
    target.di_ele[j] = value;
  }
}

void rSdpaLib::inputInitYVec(int k, double value)
{
  k--;
  if (k<0 || k>=m) {
    rMessage("Over index of y:: " << k+1);
    rError("Length of y:: " << initPt.yVec.nDim);
  }
  initPt.yVec.ele[k] = value;
}

void rSdpaLib::inputInitZMat(int l, int i, int j, double value)
{
  l--;
  if (l<0 || l>=nBlock) {
    rMessage("Over nBlock of Z:: " << l+1);
    rError("nBlock of Z :: " << initPt.zMat.nBlock);
  }
  i--;
  j--;
  int size = blockStruct[l];
  rDenseMatrix& target = initPt.zMat.ele[l];
  if (size>0) {
    if (i!=j) {
      target.de_ele[i+target.nRow*j] =value;
      target.de_ele[j+target.nRow*i] =value;
    } else {
      target.de_ele[i+target.nRow*i] =value;
    }
  } else {
    target.di_ele[j] = value;
  }
}

#endif // REVERSE_PRIMAL_DUAL

bool rSdpaLib::dumpData(const char* filename)
{
  ofstream output;
  output.open(filename);
  if (output.fail()) {
    rError("Cannot Open " << filename);
  }
  output << m << endl;
  output << nBlock << endl;
  for (int l = 0; l<nBlock ; ++l) {
    output << blockStruct[l] << " " ;
  }
  output << '\n';
  output << "{ ";
  for (int k=0; k<m; ++k) {
    output << b.ele[k] << " ";
  }
  output << "}";
  output << '\n';
  for (int l=0; l<nBlock; ++l) {
	  int index = 0;
    switch (C.ele[l].Sp_De_Di) {
    case rSparseMatrix::SPARSE:
      for (index = 0; index < C.ele[l].NonZeroCount; ++index) {
	int i = C.ele[l].row_index[index];
	int j = C.ele[l].column_index[index];
	double value = C.ele[l].sp_ele[index];
	if (value!=0.0) {
	  output << "0 " << l+1 << " "
		 << i+1 << " " << j+1 << " "
		 << -value << '\n';
	}
      }
      break;
    case rSparseMatrix::DENSE:
      // a Matrix must not be DENSE,
      // since it was read as SPARSE or DIAGONAL
      break;
    case rSparseMatrix::DIAGONAL:
      for (index = 0; index < C.ele[l].nRow; ++index) {
	double value = C.ele[l].di_ele[index];
	if (value!=0.0) {
	  output << "0 " << l+1 << " "
		 << index+1 << " " << index+1 << " "
		 << -value << '\n';
	}
      }
	break;
    } // end of switch
  }// end of 'for (int l)'

  for (int k=0; k<m; ++k) {
    for (int l=0; l<nBlock; ++l) {
		int index = 0;
      switch (A[k].ele[l].Sp_De_Di) {
      case rSparseMatrix::SPARSE:
	for (index = 0; index < A[k].ele[l].NonZeroCount;
	     ++index) {
	  int i = A[k].ele[l].row_index[index];
	  int j = A[k].ele[l].column_index[index];
	  double value = A[k].ele[l].sp_ele[index];
	  if (value!=0.0) {
	    output << k+1 << " "  << l+1 << " "
		   << i+1 << " " << j+1 << " "
		   << value << '\n';
	  }
	}
	break;
      case rSparseMatrix::DENSE:
	// a Matrix must not be DENSE,
	// since it was read as SPARSE or DIAGONAL
	break;
      case rSparseMatrix::DIAGONAL:
	for (index = 0; index < A[k].ele[l].nRow; ++index) {
	  double value = A[k].ele[l].di_ele[index];
	  if (value!=0.0) {
	    output << k+1 << " " << l+1 << " "
		   << index+1 << " " << index+1 << " "
		   << value << '\n';
	  }
	}
	break;
      } // end of switch
    } // end of 'for (int l)'
  } // end of 'for (int k)'
  output.close();
  return true;
}

bool rSdpaLib::dumpInit(const char* filename)
{
  ofstream output;
  output.open(filename);
  if (output.fail()) {
    rError("Cannot Open " << filename);
  }
  output << "{ ";
  for (int k = 0; k<m ; ++k) {
    output << -initPt.yVec.ele[k] << " ";
  }
  output << "}";
  output << '\n';
  for (int l=0; l<nBlock; ++l) {
    int size;
	int i = 0;
    switch (initPt.zMat.ele[l].De_Di) {
    case rDenseMatrix::DENSE:
      size = initPt.zMat.ele[l].nRow;
      for (i=0; i<size; ++i) {
	for (int j=0; j<size; ++j) {
	  double value = initPt.zMat.ele[l].de_ele[i+size*j];
	  if (i<=j && value!=0.0) {
	    output << "1 " << l+1 << " "
		   << i+1 << " " << j+1 << " "
		   << value << '\n';
	  }
	}
      }
      break;
    case rDenseMatrix::DIAGONAL:
      for (int index = 0; index < initPt.zMat.ele[l].nRow; ++index) {
	double value = initPt.zMat.ele[l].di_ele[index];
	if (value!=0.0) {
	  output << "1 " << l+1 << " "
		 << index+1 << " " << index+1 << " "
		 << value << '\n';
	}
      }
	break;
    } // end of switch
  }// end of 'for (int l)'
  for (int l=0; l<nBlock; ++l) {
    int size;
	int i = 0;
    switch (initPt.xMat.ele[l].De_Di) {
    case rDenseMatrix::DENSE:
      size = initPt.xMat.ele[l].nRow;
      for (i=0; i<size; ++i) {
	for (int j=0; j<size; ++j) {
	  double value = initPt.xMat.ele[l].de_ele[i+size*j];
	  if (i<=j && value!=0.0) {
	    output << "2 " << l+1 << " "
		   << i+1 << " " << j+1 << " "
		   << value << '\n';
	  }
	}
      }
      break;
    case rDenseMatrix::DIAGONAL:
      for (int index = 0; index < initPt.xMat.ele[l].nRow; ++index) {
	double value = initPt.xMat.ele[l].di_ele[index];
	if (value!=0.0) {
	  output << "2 " << l+1 << " "
		 << index+1 << " " << index+1 << " "
		 << value << '\n';
	}
      }
	break;
    } // end of switch
  }// end of 'for (int l)'
  output.close();
  return true;
}

bool rSdpaLib::checkData(int& k, int& l, int& i, int& j)
{
  rTimeStart(FILE_CHECK_START1);
  if (C.sortSparseIndex(l,i,j)==FAILURE) {
    k=0;
    l++;
    i++;
    j++;
    return false;
  }
  for (k=0; k<m; ++k) {
    if (A[k].sortSparseIndex(l,i,j)==FAILURE) {
      k++;
      l++;
      i++;
      j++;
      return false;
    }
  }
  rTimeEnd(FILE_CHECK_END1);
  com.FileCheck += rTimeCal(FILE_CHECK_START1,
			    FILE_CHECK_END1);
  return true;
}

#if REVERSE_PRIMAL_DUAL

double* rSdpaLib::getResultXVec()
{
  return currentPt.yVec.ele;
}

double* rSdpaLib::getResultXMat(int l)
{
  if (l<0 || l>=nBlock) {
    rMessage("Over nBlock of X:: " << l);
    rError("nBlock of X :: " << nBlock);
  }
  if (blockStruct[l] > 0) {
    return currentPt.zMat.ele[l].de_ele;
  } else {
    return currentPt.zMat.ele[l].di_ele;
  }
}

double* rSdpaLib::getResultYMat(int l)
{
  if (l<0 || l>=nBlock) {
    rMessage("Over nBlock of Y:: " << l);
    rError("nBlock of Y :: " << nBlock);
  }
  if (blockStruct[l] > 0) {
    return currentPt.xMat.ele[l].de_ele;
  } else {
    return currentPt.xMat.ele[l].di_ele;
  }
}

void rSdpaLib::printResultXVec(FILE* fpOut)
{
  currentPt.yVec.display(fpOut);
}
void rSdpaLib::printResultXMat(FILE* fpOut)
{
  currentPt.zMat.display(fpOut);
}
void rSdpaLib::printResultYMat(FILE* fpOut)
{
  currentPt.xMat.display(fpOut);
}

#else

double* rSdpaLib::getResultXMat(int l)
{
  if (l<0 || l>=nBlock) {
    rMessage("Over nBlock of X:: " << l);
    rError("nBlock of X :: " << nBlock);
  }
  if (blockStruct[l] > 0) {
    return currentPt.xMat.ele[l].de_ele;
  } else {
    return currentPt.xMat.ele[l].di_ele;
  }
}

double* rSdpaLib::getResultYVec()
{
  return currentPt.yVec.ele;
}

double* rSdpaLib::getResultZMat(int l)
{
  if (l<0 || l>=nBlock) {
    rMessage("Over nBlock of Z:: " << l);
    rError("nBlock of Z :: " << nBlock);
  }
  if (blockStruct[l] > 0) {
    return currentPt.zMat.ele[l].de_ele;
  } else {
    return currentPt.zMat.ele[l].di_ele;
  }
}

void rSdpaLib::printResultXMat(FILE* fpOut)
{
  currentPt.xMat.display(fpOut);
}
void rSdpaLib::printResultYVec(FILE* fpOut)
{
  currentPt.yVec.display(fpOut);
}
void rSdpaLib::printResultZMat(FILE* fpOut)
{
  currentPt.zMat.display(fpOut);
}

#endif // REVERSE_PRIMAL_DUAL

double rSdpaLib::getPrimalObj()
{
  #if REVERSE_PRIMAL_DUAL
  return -solveInfo.objValDual;
  #else
  return solveInfo.objValPrimal;
  #endif
}

double rSdpaLib::getDualObj()
{
  #if REVERSE_PRIMAL_DUAL
  return -solveInfo.objValPrimal;
  #else
  return solveInfo.objValDual;
  #endif
}

double rSdpaLib::getPrimalError()
{
  #if REVERSE_PRIMAL_DUAL
  return currentRes.normDualMat;
  #else
  return currentRes.normPrimalVec;
  #endif
}

double rSdpaLib::getDualError()
{
  #if REVERSE_PRIMAL_DUAL
  return currentRes.normPrimalVec;
  #else
  return currentRes.normDualMat;
  #endif
}

double rSdpaLib::getDigits()
{
  double mean = (fabs(solveInfo.objValPrimal)
		 + fabs(solveInfo.objValDual)) / 2.0;
  double PDgap = fabs(solveInfo.objValPrimal
		      - solveInfo.objValDual);
  double digits = 1000; // 1000 means infinity in this case
  digits = -log10(fabs(PDgap/mean));
  return digits;
}

double rSdpaLib::getMu()
{
  return mu.current;
}

double rSdpaLib::getTime()
{
  return com.MainLoop;
}

int rSdpaLib::getIteration()
{
  return iteration;
}

rSdpaLib::phaseType rSdpaLib::getPhaseValue()
{
  switch (phase.value) {
  case rSolveInfo::noINFO    : return noINFO;
  case rSolveInfo::pFEAS     : return pFEAS;
  case rSolveInfo::dFEAS     : return dFEAS;
  case rSolveInfo::pdFEAS    : return pdFEAS;
  case rSolveInfo::pdINF     : return pdINF;
  case rSolveInfo::pFEAS_dINF: return pFEAS_dINF;
  case rSolveInfo::pINF_dFEAS: return pINF_dFEAS;
  case rSolveInfo::pdOPT     : return pdOPT;
  case rSolveInfo::pUNBD     : return pUNBD;
  case rSolveInfo::dUNBD     : return dUNBD;
  }
  rError("Convert Error");
}

void rSdpaLib::stringPhaseValue(char* str)
{
  // str needs stirng with at least 15 length.
  char* copy = "";
  switch (phase.value) {
  case rSolveInfo::noINFO    : copy = "noINFO    "; break;
  case rSolveInfo::pFEAS     : copy = "pFEAS     "; break;
  case rSolveInfo::dFEAS     : copy = "dFEAS     "; break;
  case rSolveInfo::pdFEAS    : copy = "pdFEAS    "; break;
  case rSolveInfo::pdINF     : copy = "pdINF     "; break;
  case rSolveInfo::pFEAS_dINF: copy = "pFEAS_dINF"; break;
  case rSolveInfo::pINF_dFEAS: copy = "pINF_dFEAS"; break;
  case rSolveInfo::pdOPT     : copy = "pdOPT     "; break;
  case rSolveInfo::pUNBD     : copy = "pUNBD     "; break;
  case rSolveInfo::dUNBD     : copy = "dUNBD     "; break;
  default:   rError("Convert Error"); break;
  }
  strcpy(str,copy);
}

void rSdpaLib::printTime(FILE* fpOut)
{
  com.display(fpOut);
}

void rSdpaLib::solve()
{
  if (CheckMatrix) {
    rTimeStart(FILE_CHECK_START1);
      int i=0,j=0,k=0,l=0;
    // this method needs long time
    checkData(k,l,i,j);
    if (i>0 || j>0) {
      rError("checkData stops");
      cout << "constraint " << k <<":"
	   << "block " << l <<":"
	   << "row " << i <<":"
	   << "column " << j << endl;
    }
    rTimeEnd(FILE_CHECK_END1);
    com.FileCheck += rTimeCal(FILE_CHECK_START1,
			      FILE_CHECK_END1);
  }
#if 1
  rTimeStart(FILE_CHANGE_START1);
  // if possible, change C and A to Dense type matrices.
  C.changeToDense();
  for (int k=0; k<m; ++k) {
    A[k].changeToDense();
  }
  rTimeEnd(FILE_CHANGE_END1);
  com.FileChange += rTimeCal(FILE_CHANGE_START1,
			     FILE_CHANGE_END1);
#endif

  if (InitialPoint) {
    initPt.initializeResetup(m,nBlock,blockStruct,com);
    currentPt.copyFrom(initPt);
    // rMessage("initialize? ");
  } else {
    currentPt.initialize(m,nBlock,blockStruct,pARAM.lambdaStar,com);
  }
  // rMessage("initPt = ");
  // initPt.display();
  // rMessage("currentPt = ");
  // currentPt.display();
  initRes.initialize(m, nBlock, blockStruct, b, C, A, currentPt);
  currentRes.copyFrom(initRes);
  // rMessage("initial currentRes = ");
  // currentRes.display();

  newton.initialize(m, nBlock, blockStruct);
  newton.computeFormula(m,A,0.0,KAPPA);

  alpha.initialize(1.0,1.0,nBlock, blockStruct);
  beta.initialize(pARAM.betaStar);
  reduction.initialize(rSwitch::ON);
  mu.initialize(pARAM.lambdaStar);
  lanczos.initialize(nBlock,blockStruct);

  if (InitialPoint) {
    mu.initialize(nDim,initPt);
  }

  theta.initialize(pARAM,initRes);
  solveInfo.initialize(nDim, b, C, A, initPt, mu.initial,
		       pARAM.omegaStar);
  phase.initialize(initRes, solveInfo, pARAM, nDim);

  rIO::printHeader(OutputFile,DisplayInformation);
  
  int pIteration = 0;

  // -----------------------------------------------------
  // Here is MAINLOOP
  // -----------------------------------------------------

  rTimeStart(MAIN_LOOP_START1);

  // explisit maxIteration
  // pARAM.maxIteration = 100;
  while (phase.updateCheck(currentRes, solveInfo, pARAM)
	 && pIteration < pARAM.maxIteration) {
    // rMessage(" turn hajimari " << pIteration );

    #if 0
    if (alpha.primal<1.0e-5 && alpha.dual<1.0e-5) {
      break;
    }
    #endif

    // Mehrotra's Predictor
    rTimeStart(MEHROTRA_PREDICTOR_START1);

    // calculate variables of Mehrotra
    reduction.MehrotraPredictor(phase);
    beta.MehrotraPredictor(phase, reduction, pARAM);
    // rMessage("reduction = ");
    // reduction.display();
    // rMessage("phase = ");
    // phase.display();
    // rMessage("beta.predictor.value = " << beta.value);
    
    // rMessage("xMat = ");
    // currentPt.xMat.display();
    // rMessage("zMat = ");
    // currentPt.zMat.display();
    // rMessage(" mu = " << mu.current);

    bool isSuccessCholesky;
    isSuccessCholesky = newton.Mehrotra(rNewton::PREDICTOR,
					m, A, C, mu,
					beta, reduction,
					phase, currentPt,
					currentRes, com);
    if (isSuccessCholesky == false) {
      break;
    }
      
    // rMessage("newton predictor = ");
    // newton.display();
    // rMessage("newton Dy predictor = ");
    // newton.DyVec.display();
    // newton.bMat.display();
    rTimeEnd(MEHROTRA_PREDICTOR_END1);
    com.Predictor += rTimeCal(MEHROTRA_PREDICTOR_START1,
			      MEHROTRA_PREDICTOR_END1);
    
    rTimeStart(STEP_PRE_START1);
    alpha.MehrotraPredictor(b,C,A, currentPt, phase,
			    newton, lanczos, com);
    // rMessage("xMat = ");
    // currentPt.xMat.display();
    // rMessage("zMat = ");
    // currentPt.zMat.display();
    // rMessage("alpha predictor = ");
    // alpha.display();
    // phase.display();
    // rMessage("newton predictor = ");
    // newton.display();
    // rMessage("currentPt = ");
    // currentPt.display();
    rTimeStart(STEP_PRE_END1);
    com.StepPredictor += rTimeCal(STEP_PRE_START1,STEP_PRE_END1);

    // rMessage("alphaStar = " << pARAM.alphaStar);
    // Mehrotra's Corrector
    // rMessage(" Corrector ");
    rTimeStart(CORRECTOR_START1);
    beta.MehrotraCorrector(nDim,phase,alpha,currentPt,
			   newton,mu,pARAM);
    // rMessage("beta corrector = " << beta.value);
    newton.Mehrotra(rNewton::CORRECTOR,m,A,C,mu,beta,reduction,
		    phase, currentPt, currentRes, com);
    // rMessage("currentPt = ");
    // currentPt.display();
    // rMessage("newton corrector = ");
    // newton.display();
    // rMessage("newton Dy corrector = ");
    // newton.DyVec.display();

    rTimeEnd(CORRECTOR_END1);
    com.Corrector += rTimeCal(CORRECTOR_START1,
			      CORRECTOR_END1);
      
    rTimeStart(CORRECTOR_STEP_START1);
    alpha.MehrotraCorrector(nDim, b, C, A, currentPt, phase,
			    reduction, newton, mu, theta,
			    lanczos, pARAM, com);
    // rMessage("alpha corrector = ");
    // alpha.display();
    rTimeEnd(CORRECTOR_STEP_END1);
    com.StepCorrector += rTimeCal(CORRECTOR_STEP_START1,
				  CORRECTOR_STEP_END1);
    // the end of Corrector
    
    rIO::printOneIteration(pIteration, mu, theta, solveInfo,
			   alpha, beta, currentRes, OutputFile,
			   DisplayInformation);

    if (currentPt.update(alpha,newton,com)==false) {
      // if step length is too short,
      // the algorithm ends.
      // rMessage("cannot move");
      pIteration++;
      break;
    }

    // rMessage("currentPt = ");
    // currentPt.display();
    // rMessage("newton = ");
    // newton.display();

    // rMessage("updated");
    theta.update(reduction,alpha);
    mu.update(nDim,currentPt);
    currentRes.update(m,nBlock,blockStruct,b,C,A,
		      initRes, theta, currentPt, phase, mu,com);
    
    theta.update_exact(initRes,currentRes);
    solveInfo.update(nDim, b, C, initPt, currentPt,
		     currentRes, mu, theta, pARAM);
    pIteration++;

  } // end of MAIN_LOOP

  rTimeEnd(MAIN_LOOP_END1);

  com.MainLoop = rTimeCal(MAIN_LOOP_START1,
			  MAIN_LOOP_END1);
  currentPt.update_last(com);
  currentRes.compute(m,nBlock,blockStruct,b,C,A,currentPt,mu);
  rIO::printLastInfo(pIteration, mu, theta, solveInfo, alpha, beta,
		     currentRes, phase, currentPt, com.TotalTime,
		     nDim,b,C,A,com,pARAM, OutputFile,
		     DisplayInformation,false);

  #if REVERSE_PRIMAL_DUAL
  // reverse the sign of y and phase
  rAl::let(currentPt.yVec,'=',currentPt.yVec,'*',&DMONE);
  phase.reverse();
  #endif
  
  iteration = pIteration;

  // for compability SDPA
  PrimalObj    = getPrimalObj();
  DualObj      = getDualObj();
  PrimalError  = getPrimalError();
  DualError    = getDualError();
  Iteration    = getIteration();
  switch (phase.value) {
  case rSolveInfo::noINFO    : Value = ::noINFO;      break;
  case rSolveInfo::pFEAS     : Value = ::pFEAS;       break;
  case rSolveInfo::dFEAS     : Value = ::dFEAS;       break;
  case rSolveInfo::pdFEAS    : Value = ::pdFEAS;      break;
  case rSolveInfo::pdINF     : Value = ::pdINF;       break;
  case rSolveInfo::pFEAS_dINF: Value = ::pFEAS_dINF;  break;
  case rSolveInfo::pINF_dFEAS: Value = ::pINF_dFEAS;  break;
  case rSolveInfo::pdOPT     : Value = ::pdOPT;       break;
  case rSolveInfo::pUNBD     : Value = ::pUNBD;       break;
  case rSolveInfo::dUNBD     : Value = ::dUNBD;       break;
  }
}

void rSdpaLib::initializeFromFile()
{
  rTimeStart(FILE_READ_START1);
  bool isDataSparse = false;
  int len = strlen(InputFileName);
  if (InputFileName[len-1] == 's'
      && InputFileName[len-2] == '-') {
    isDataSparse = true;
  }

  // initialize b,C,A
  char titleAndComment[1024];
  rIO::read(InputFile,OutputFile,m,titleAndComment);
  if (OutputFile) {
    fprintf(OutputFile,"data      is %s\n",InputFileName);
    fprintf(OutputFile,"parameter is %s\n",ParameterFileName);
    fprintf(OutputFile,"initial   is %s\n",InitialFileName);
    fprintf(OutputFile,"out       is %s\n",OutputFileName);
  }
  mDIM = m;
  rIO::read(InputFile,nBlock);
  blockStruct = NULL;
  blockStruct = new int[nBlock];
  if (blockStruct==NULL) {
    rError("Memory exhausted about blockStruct");
  }
  nBLOCK = nBlock;
  bLOCKsTRUCT = blockStruct;
  rIO::read(InputFile,nBlock,blockStruct);
  nDim = 0;
  for (int l=0; l<nBlock; ++l) {
    nDim += abs(blockStruct[l]);
  }
  
  b.initialize(m);
  rIO::read(InputFile,b);
  A = new rBlockSparseMatrix[m];
  if (A==NULL) {
    rError("Memory exhausted about blockStruct");
  }

  long position = ftell(InputFile);
  // C,A must be accessed "twice".

  // count numbers of elements of C and A
  int* CNonZeroCount = NULL;
  CNonZeroCount = new int[nBlock];
  if (CNonZeroCount==NULL) {
    rError("Memory exhausted about blockStruct");
  }
  int* ANonZeroCount = NULL;
  ANonZeroCount = new int[nBlock*m];
  if (ANonZeroCount==NULL) {
    rError("Memory exhausted about blockStruct");
  }
  // initialize C and A
  rIO::read(InputFile,m,nBlock,blockStruct,
	    CNonZeroCount,ANonZeroCount,isDataSparse);
  // rMessage(" C and A count over");
  C.initialize(nBlock,blockStruct);
  for (int l=0; l<nBlock; ++l) {
    int size = blockStruct[l];
    if (size > 0) {
      C.ele[l].initialize(size,size,rSparseMatrix::SPARSE,
			  CNonZeroCount[l]);
    } else {
      C.ele[l].initialize(-size,-size,rSparseMatrix::DIAGONAL,
			  -size);
    }
  }
  for (int k=0; k<m; ++k) {
    A[k].initialize(nBlock,blockStruct);
    for (int l=0; l<nBlock; ++l) {
      int size = blockStruct[l];
      if (size > 0) {
	A[k].ele[l].initialize(size,size,rSparseMatrix::SPARSE,
			       ANonZeroCount[k*nBlock+l]);
      } else {
	A[k].ele[l].initialize(-size,-size,rSparseMatrix::DIAGONAL,
			       -size);
      }
    }
  }
  delete[] CNonZeroCount;
  CNonZeroCount = NULL;
  delete[] ANonZeroCount;
  ANonZeroCount = NULL;

  rIO::read(InputFile, C, A, m, nBlock, blockStruct,
	    position, isDataSparse);

  if (InitialFile != NULL && InitialPoint == true) {
    // isInitPoint = true;
    bool isInitSparse = false;
    int len = strlen(InitialFileName);
    if (InitialFileName[len-1] == 's'
	&& InitialFileName[len-2] == '-') {
      isInitSparse = true;
    }
    initPt.initializeZero(m,nBlock,blockStruct,com);
    rIO::read(InitialFile,initPt.xMat,initPt.yVec,initPt.zMat, nBlock,
	      blockStruct, isInitSparse);
    initPt.initializeResetup(m,nBlock,blockStruct,com);
  } else {
    initPt.initialize(m,nBlock,blockStruct,pARAM.lambdaStar,com);
  }
 
  rTimeEnd(FILE_READ_END1);
  com.FileRead += rTimeCal(FILE_READ_START1,
			   FILE_READ_END1);

}



void rSdpaLib::Delete()
{
  delete2();
  #if 0
  // bLOCKsTRUCT cannot release successfully
  if (bLOCKsTRUCT) {
    delete[] bLOCKsTRUCT;
  }
  bLOCKsTRUCT = NULL;
  #endif
}

// for compability SDPA
extern bool SDPA_initialize(rSdpaLib& SDP)
{
  if (SDP.ParameterFile != NULL) {
    SDP.pARAM.readFile(SDP.ParameterFile);
  }
  if (SDP.InputFile != NULL) {
    SDP.initializeFromFile();
  }
  return true;
}

extern bool SDPA_initialize2(rSdpaLib& SDP)
{
  if (SDP.InputFile == NULL) {
    SDP.initialize1(SDP.mDIM,SDP.nBLOCK,SDP.bLOCKsTRUCT,
		    SDP.InitialPoint);
  }
  return true;
}

bool SDPA_Input_cVECT(rSdpaLib& SDP, int i, double value)
{
  SDP.inputCVec(i,value);
  return true;
}

bool SDPA_CountUpperTriangle(rSdpaLib& SDP, int i, int j,
				int number)
{
  SDP.countUpperTriangle(i,j,number);
  return true;
}

bool SDPA_Make_sfMAT(rSdpaLib& SDP)
{
  SDP.initialize2();
  return true;
}

bool SDPA_InputElement(rSdpaLib& SDP,
			  int i, int j, int k, int ell,
			  double value)
{
  SDP.inputElement(i,j,k,ell,value);
  return true;
}

bool SDPA_Input_IniXMat(rSdpaLib& SDP, int j, int k, int ell,
			   double value)
{
  SDP.inputInitXMat(j,k,ell,value);
  return true;
}

bool SDPA_Input_InixVec(rSdpaLib& SDP, int i, double value)
{
  SDP.inputInitXVec(i,value); 
  return true;
}

bool SDPA_Input_IniYMat(rSdpaLib& SDP, int j, int k, int ell,
			   double value)
{
  SDP.inputInitYMat(j,k,ell,value);
  return true;
}

bool SDPA_Check_sfMAT(rSdpaLib& SDP)
{
  int k=0,l=0,i=0,j=0;
  SDP.checkData(k,l,i,j);
  if (i>0 || j>0) {
    return false;
  }
  return true;
}

bool SDPA_Solve(rSdpaLib& SDP)
{
  SDP.solve();
  return true;
}

bool SDPA_Copy_Current_To_Ini(rSdpaLib& SDP)
{
  SDP.initPt.copyFrom(SDP.currentPt);
  #if REVERSE_PRIMAL_DUAL
  rAl::let(SDP.initPt.yVec,'=',SDP.initPt.yVec,'*',&DMONE);
  #endif
  return true;
}


