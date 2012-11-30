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
  rsdpa_lib.h
-------------------------------------------------*/

#ifndef __rsdpa_lib_h__
#define __rsdpa_lib_h__

enum phaseType {noINFO,pFEAS,dFEAS,pdFEAS,pdINF,pFEAS_dINF,
		pINF_dFEAS,pdOPT,pUNBD,dUNBD};
enum parameterType {PARAMETER_DEFAULT,
		    PARAMETER_AGGRESSIVE,
		    PARAMETER_STABLE};
enum which_method {KSH, NT, AHO};

class rSdpaLib
{
public:
  enum phaseType {noINFO,pFEAS,dFEAS,pdFEAS,pdINF,pFEAS_dINF,
		  pINF_dFEAS,pdOPT,pUNBD,dUNBD};
  enum rParameterType {rPARAMETER_DEFAULT,
		       rPARAMETER_AGGRESSIVE,
		       rPARAMETER_STABLE};
  rSdpaLib();
  ~rSdpaLib();
  void delete1();
  void delete2();

  void setDisplay(FILE* Display=stdout);

  void setDefaultParameter(parameterType type
			   = PARAMETER_DEFAULT);
  void setDefaultParameter(rParameterType type);
  void setParameterMaxIteration(int maxIteration);
  void setParameterEpsilonStar (double epsilonStar);
  void setParameterLambdaStar  (double lambdaStar);
  void setParameterOmegaStar   (double omegaStar);
  void setParameterLowerBound  (double lowerBound);
  void setParameterUpperBound  (double upperBound);
  void setParameterBetaStar    (double betaStar);
  void setParameterBetaBar     (double betaBar);
  void setParameterGammaStar   (double gammaStar);
  void setParameterEpsilonDash (double epsilonDash);
  
  void initialize1(int m, int nBlock,int* blockStruct,
		   bool initialPoint=false);
  void countUpperTriangle(int k, int l, int nonzero);
  void initialize2();

#if REVERSE_PRIMAL_DUAL
  void inputCVec(int k, double value);
#else
  void inputBVec(int k, double value);
  void inputCMat(int l, int i, int j, double value);
#endif
  void inputElement(int k, int l, int i, int j, double value);

#if REVERSE_PRIMAL_DUAL
  void inputInitXVec(int k,double value); 
  void inputInitXMat(int l,int i,int j, double value);
  void inputInitYMat(int l,int i,int j, double value);
#else
  void inputInitXMat(int l,int i,int j, double value);
  void inputInitYVec(int k,double value);
  void inputInitZMat(int l,int i,int j, double value);
#endif
  bool checkData(int& k,int& l,int& i, int& j);
  bool dumpData(const char* filename);
  bool dumpInit(const char* filename);
  
  void solve();
  
#if REVERSE_PRIMAL_DUAL
  double* getResultXVec();
  double* getResultXMat(int l);
  double* getResultYMat(int l);
  void printResultXVec(FILE* fpOut = stdout);
  void printResultXMat(FILE* fpOut = stdout);
  void printResultYMat(FILE* fpOut = stdout);
#else
  double* getResultXMat(int l);
  double* getResultYVec();
  double* getResultZMat(int l);
  void printResultXMat(FILE* fpOut = stdout);
  void printResultYVec(FILE* fpOut = stdout);
  void printResultZMat(FILE* fpOut = stdout);
#endif
  double getPrimalObj();
  double getDualObj();
  double getPrimalError();
  double getDualError();
  double getDigits();
  int    getIteration();
  double getMu();
  phaseType getPhaseValue();
  void   stringPhaseValue(char* str);
  double getTime();
  
  void printTime(FILE* fpOut=stdout);

  rParameter pARAM;

  /*----------------------------------------*/
  // for compatibility of SDPA
  double PrimalObj;
  double DualObj;
  double PrimalError;
  double DualError;
  int    Iteration;
  FILE* DisplayInformation;

  int mDIM;
  int nBLOCK;
  int* bLOCKsTRUCT;
  bool InitialPoint;
  int Method;
  bool CheckMatrix;

  char ParameterFileName[1024];
  char OutputFileName[1024];
  char InitialFileName[1024];
  char InputFileName[1024];

  FILE* ParameterFile;
  FILE* OutputFile;
  FILE* InitialFile;
  FILE* InputFile;
  
  void initializeFromFile();

  ::phaseType Value;
  void Delete();
  /*----------------------------------------*/
  
 
private:

  int m;
  int nBlock;
  int* blockStruct;

  int nDim;
  int* nonZeroNumber;
  rComputeTime com;

  rVector b;
  rBlockSparseMatrix C;
  rBlockSparseMatrix* A;

  rSolutions initPt;
  rSolutions currentPt;
  // bool isInitPoint;

  rResiduals initRes;
  rResiduals currentRes;

  rNewton newton;
  rStepLength alpha;
  rDirectionParameter beta;
  rSwitch reduction;
  rAverageComplementarity mu;
  rLanczos lanczos;

  rRatioInitResCurrentRes theta;
  rSolveInfo solveInfo;
  rPhase phase;
  bool hasSolved;
  bool hasDelete1;
  int iteration;

  // for compability SDPA
  friend bool SDPA_Copy_Current_To_Ini(rSdpaLib& SDP);
};

bool SDPA_initialize(rSdpaLib& SDP);
bool SDPA_initialize2(rSdpaLib& SDP);
bool SDPA_Input_cVECT(rSdpaLib& SDP, int i, double value);
bool SDPA_CountUpperTriangle(rSdpaLib& SDP, int i, int j,
				int number);
bool SDPA_Make_sfMAT(rSdpaLib& SDP);
bool SDPA_InputElement(rSdpaLib& SDP,
			  int i, int j, int k, int ell,
			  double value);
bool SDPA_Input_IniXMat(rSdpaLib& SDP, int j, int k, int ell,
			   double value);
bool SDPA_Input_InixVec(rSdpaLib& SDP, int i, double value);
bool SDPA_Input_IniYMat(rSdpaLib& SDP, int j, int k, int ell,
			   double value);
bool SDPA_Check_sfMAT(rSdpaLib& SDP);
bool SDPA_Solve(rSdpaLib& SDP);
bool SDPA_Copy_Current_To_Ini(rSdpaLib& SDP);

typedef rSdpaLib SDPA;
#endif // __rsdpa_lib_h__
