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
  rsdpa_main.cpp
-------------------------------------------------*/

#include "rsdpa_io.h"
#include <time.h>
#define LengthOfBuffer 1024
static double KAPPA = 2.2;

bool pinpal(char* dataFile, char* initFile, char* outFile,
	    char* paraFile, bool isInitFile, bool isInitSparse,
	    bool isDataSparse, bool isParameter,
	    rParameter::parameterType parameterType,
	    FILE* Display)
{

  rTimeStart(TOTAL_TIME_START1);
  rTimeStart(FILE_READ_START1);
  rComputeTime com;
  
  FILE* fpData      = NULL;
  FILE* fpOut       = NULL;

  if ((fpOut=fopen(outFile,"w"))==NULL) {
    rError("Cannot open out file " << outFile);
  }
  rParameter param;
  param.setDefaultParameter(parameterType);
  if (isParameter) {
    FILE* fpParameter = NULL;
    if ((fpParameter=fopen(paraFile,"r"))==NULL) {
      fprintf(Display,"Cannot open parameter file %s \n",
	      paraFile);
      exit(0);
    } else {
      param.readFile(fpParameter);
      fclose(fpParameter);
    }
  }
  // param.display(Display);

  if ((fpData=fopen(dataFile,"r"))==NULL) {
    rError("Cannot open data file " << dataFile);
  }
  char titleAndComment[LengthOfBuffer];
  int m;
  time_t ltime;
  time( &ltime );
  fprintf(fpOut,"SDPA start at %s",ctime(&ltime));
  rIO::read(fpData,fpOut,m,titleAndComment);
  fprintf(fpOut,"data      is %s\n",dataFile);
  if (paraFile) {
    fprintf(fpOut,"parameter is %s\n",paraFile);
  }
  if (initFile) {
    fprintf(fpOut,"initial   is %s\n",initFile);
  }
  fprintf(fpOut,"out       is %s\n",outFile);

  int nBlock;
  rIO::read(fpData,nBlock);
  int* blockStruct = NULL;
  blockStruct = new int[nBlock];
  if (blockStruct==NULL) {
    rError("Memory exhausted about blockStruct");
  }
  rIO::read(fpData,nBlock,blockStruct);
  int nDim = 0;
  for (int l=0; l<nBlock; ++l) {
    nDim += abs(blockStruct[l]);
  }
  
  // rMessage("b has not been read yet , m = " << m);
  rVector b(m);
  rIO::read(fpData,b);
  // rMessage("b has been read");
  
  rBlockSparseMatrix C;
  rBlockSparseMatrix* A = NULL;
  A = new rBlockSparseMatrix[m];
  if (A==NULL) {
    rError("Memory exhausted about blockStruct");
  }

  long position = ftell(fpData);
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
  rIO::read(fpData,m,nBlock,blockStruct,
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
    
  // rMessage(" C and A initialize over");
  rIO::read(fpData, C, A, m, nBlock, blockStruct, position, isDataSparse);
  // rMessage(" C and A have been read");
  fclose(fpData);

#if 0
  fprintf(Display,"C = \n");
  C.display(Display);
  for (int k=0; k<m; ++k) {
    fprintf(Display,"A[%d] = \n",k);
    A[k].display(Display);
  }
#endif

#if 0
  // write  C and A in SDPA sparse data format to file
  ofstream output;
  output.open("dumped.rsdpa.dat-s");
  if (output.fail()) {
    rError("Cannot Open dumped.rsdpa.dat-s");
  }
  output << m << endl;
  output << nBlock << endl;
  for (l = 0; l<nBlock ; ++l) {
    output << blockStruct[l] << " " ;
  }
  output << endl;
  for (k=0; k<m; ++k) {
	output << b.ele[k] << " ";
  }
  output << endl;
  int index=0;
  for (l=0; l<nBlock; ++l) {
    switch (C.ele[l].Sp_De_Di) {
    case rSparseMatrix::SPARSE:
      for (index = 0; index < C.ele[l].NonZeroCount; ++index) {
	int i = C.ele[l].row_index[index];
	int j = C.ele[l].column_index[index];
	double value = C.ele[l].sp_ele[index];
	if (value!=0.0) {
	  output << "0 " << l+1 << " "
		 << i+1 << " " << j+1 << " "
		 << -value << endl;
	}
      }
      break;
    case rSparseMatrix::DENSE:
      break;
    case rSparseMatrix::DIAGONAL:
      for (int index = 0; index < C.ele[l].nRow; ++index) {
	double value = C.ele[l].di_ele[index];
	if (value!=0.0) {
	  output << "0 " << l+1 << " "
		 << index+1 << " " << index+1 << " "
		 << -value << endl;
	}
      }
	break;
    } // end of switch
  }// end of 'for (int l)'

  for (k=0; k<m; ++k) {
    for (int l=0; l<nBlock; ++l) {
      switch (A[k].ele[l].Sp_De_Di) {
      case rSparseMatrix::SPARSE:
	for (index = 0; index < A[k].ele[l].NonZeroCount; ++index) {
	  int i = A[k].ele[l].row_index[index];
	  int j = A[k].ele[l].column_index[index];
	  double value = A[k].ele[l].sp_ele[index];
	  if (value!=0.0) {
	    output << k+1 << " "  << l+1 << " "
		   << i+1 << " " << j+1 << " "
		   << value << endl;
	  }
	}
	break;
      case rSparseMatrix::DENSE:
	break;
      case rSparseMatrix::DIAGONAL:
	for (int index = 0; index < A[k].ele[l].nRow; ++index) {
	  double value = A[k].ele[l].di_ele[index];
	  if (value!=0.0) {
	    output << k+1 << " " << l+1 << " "
		   << index+1 << " " << index+1 << " "
		   << value << endl;
	  }
	}
	break;
      } // end of switch
    } // end of 'for (int l)'
  } // end of 'for (int k)'
  output.close();
#endif
  
#if 0
  rTimeStart(FILE_CHECK_START1);
  // check whether C,A are symmetric or not.
  int lin,iin,jin;
  if (C.sortSparseIndex(lin,iin,jin)==FAILURE) {
    fprintf(Display,"C is not symmetric, block %d,"
	    "(%d,%d) ", lin+1,iin+1,jin+1);
    exit(0);
  }
  for (int k=0; k<m; ++k) {
    if (A[k].sortSparseIndex(lin,iin,jin)==FAILURE) {
      fprintf(Display,"A[%d] is not symmetric, block %d,"
	      "(%d,%d) ", k+1,lin+1,iin+1,jin+1);
      exit(0);
    }
  }
  rTimeEnd(FILE_CHECK_END1);
  com.FileCheck += rTimeCal(FILE_CHECK_START1,
			    FILE_CHECK_END1);
#endif
  
#if 1
  rTimeStart(FILE_CHANGE_START1);
  // if possible , change C and A to Dense
  C.changeToDense();
  for (int k=0; k<m; ++k) {
    A[k].changeToDense();
  }
  rTimeEnd(FILE_CHANGE_END1);
  com.FileChange += rTimeCal(FILE_CHANGE_START1,
			     FILE_CHANGE_END1);
#endif

  // rMessage("C = ");
  // C.display(Display);
  // for (int k=0; k<m; ++k) {
  //   rMessage("A["<<k<<"] = ");
  //   A[k].display(Display);
  //   }
  
  // the end of initialization of C and A

  // set initial solutions.
  rSolutions initPt;
  rSolutions currentPt;

  if (isInitFile) {
    initPt.initializeZero(m,nBlock,blockStruct,com);
    FILE* fpInit = NULL;
    if ((fpInit=fopen(initFile,"r"))==NULL) {
      rError("Cannot open init file " << initFile);
    }
    rIO::read(fpInit,initPt.xMat,initPt.yVec,initPt.zMat, nBlock,
	      blockStruct, isInitSparse);
    fclose(fpInit);
    initPt.initializeResetup(m,nBlock,blockStruct,com);
    currentPt.copyFrom(initPt);
  } else {
    initPt.initialize(m,nBlock,blockStruct,param.lambdaStar,com);
    currentPt.initialize(m,nBlock,blockStruct,param.lambdaStar,com);
  }
  // rMessage("initial pt = ");
  // initPt.display(Display);
  // rMessage("current pt = ");
  // currentPt.display(Display);
  
  
  rTimeEnd(FILE_READ_END1);
  com.FileRead += rTimeCal(FILE_READ_START1,
			   FILE_READ_END1);

  // -------------------------------------------------------------
  // the end of file read
  // -------------------------------------------------------------
  
  rResiduals initRes(m, nBlock, blockStruct, b, C, A, currentPt);
  rResiduals currentRes;
  currentRes.copyFrom(initRes);
  // rMessage("initial currentRes = ");
  // currentRes.display(Display);

  rNewton newton(m, nBlock, blockStruct);
  newton.computeFormula(m,A,0.0,KAPPA);

  rStepLength alpha(1.0,1.0,nBlock, blockStruct);
  rDirectionParameter beta(param.betaStar);
  rSwitch reduction(rSwitch::ON);
  rAverageComplementarity mu(param.lambdaStar);
  rLanczos lanczos(nBlock,blockStruct);
  
  // rMessage("init mu");
  // mu.display();
  if (isInitFile) {
    mu.initialize(nDim, initPt);
  }
  rRatioInitResCurrentRes theta(param, initRes);
  rSolveInfo solveInfo(nDim, b, C, A, initPt, mu.initial,
		       param.omegaStar);
  rPhase phase(initRes, solveInfo, param, nDim);

  int pIteration = 0;
  rIO::printHeader(fpOut, Display);
  // -----------------------------------------------------
  // Here is MAINLOOP
  // -----------------------------------------------------

  rTimeStart(MAIN_LOOP_START1);

  // explicit maxIteration
  // param.maxIteration = 2;
  while (phase.updateCheck(currentRes, solveInfo, param)
	 && pIteration < param.maxIteration) {
    // rMessage(" turn hajimari " << pIteration );
    // Mehrotra's Predictor
    rTimeStart(MEHROTRA_PREDICTOR_START1);

    // set variable of Mehrotra
    reduction.MehrotraPredictor(phase);
    beta.MehrotraPredictor(phase, reduction, param);
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
    alpha.MehrotraPredictor(b,C,A, currentPt, phase, newton,
			    lanczos, com);
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

    // rMessage("alphaStar = " << param.alphaStar);
    // Mehrotra's Corrector
    // rMessage(" Corrector ");
    rTimeStart(CORRECTOR_START1);
    beta.MehrotraCorrector(nDim,phase,alpha,currentPt,
			   newton,mu,param);
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
			    lanczos, param, com);
    // rMessage("alpha corrector = ");
    // alpha.display();
    rTimeEnd(CORRECTOR_STEP_END1);
    com.StepCorrector += rTimeCal(CORRECTOR_STEP_START1,
				  CORRECTOR_STEP_END1);
    // the end of Corrector
    
    rIO::printOneIteration(pIteration, mu, theta, solveInfo,
			   alpha, beta, currentRes, fpOut, Display);

    if (currentPt.update(alpha,newton,com)==false) {
      // if step length is too short,
      // we finish algorithm
      rMessage("cannot move");
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
		     currentRes, mu, theta, param);
    pIteration++;

  } // end of MAIN_LOOP

  rTimeEnd(MAIN_LOOP_END1);

  com.MainLoop = rTimeCal(MAIN_LOOP_START1,
			  MAIN_LOOP_END1);
  currentPt.update_last(com);
  currentRes.compute(m,nBlock,blockStruct,b,C,A,currentPt,mu);
  
  rTimeEnd(TOTAL_TIME_END1);
  
  com.TotalTime = rTimeCal(TOTAL_TIME_START1,
			   TOTAL_TIME_END1);

  #if REVERSE_PRIMAL_DUAL
  phase.reverse();
  #endif
#if 1
  rIO::printLastInfo(pIteration, mu, theta, solveInfo, alpha, beta,
		     currentRes, phase, currentPt, com.TotalTime,
		     nDim, b, C, A, com, param, fpOut, Display);
#endif
  // com.display(fpOut);

  if (blockStruct) {
    delete[] blockStruct;
    blockStruct = NULL;
  }

  C.~rBlockSparseMatrix();
  for (int k=0; k<m; ++k) {
    A[k].~rBlockSparseMatrix();
  }
  delete[] A;
  A = NULL;

  
  fprintf(Display,   "  main loop time = %.6f\n",com.MainLoop);
  fprintf(fpOut,   "    main loop time = %.6f\n",com.MainLoop);
  fprintf(Display,   "      total time = %.6f\n",com.TotalTime);
  fprintf(fpOut,   "        total time = %.6f\n",com.TotalTime);
  #if 0
  fprintf(Display,   "file  check time = %.6f\n",com.FileCheck);
  fprintf(fpOut,   "  file  check time = %.6f\n",com.FileCheck);
  fprintf(Display,   "file change time = %.6f\n",com.FileChange);
  fprintf(fpOut,   "  file change time = %.6f\n",com.FileChange);
  #endif
  fprintf(Display,   "file   read time = %.6f\n",com.FileRead);
  fprintf(fpOut,   "  file   read time = %.6f\n",com.FileRead);
  fclose(fpOut);
  
  return true;
}

static void message(char* argv0)
{
  cout << endl;
  cout << "*** Please assign data file and output file.***" << endl;
  cout << endl;
  cout << "---- option type 1 ------------" << endl;
  cout << argv0 <<" DataFile OutputFile [InitialPtFile]"
    " [-pt parameters]"<< endl;
  cout << "parameters = 0 default, 1 aggressive,"
    " 2 stable" << endl;
  cout << "example1-1: " << argv0
       << " example1.dat example1.result" << endl;
  cout << "example1-2: " << argv0
       << " example1.dat-s example1.result" << endl;
  cout << "example1-3: " << argv0
       << " example1.dat example1.result example1.ini" << endl;
  cout << "example1-4: " << argv0
       << " example1.dat example1.result -pt 2" << endl;

  cout << endl;
  cout << "---- option type 2 ------------" << endl;
  cout << argv0 << " [option filename]+ " << endl;
  cout << "  -dd : data dense :: -ds : data sparse     " << endl;
  cout << "  -id : init dense :: -is : init sparse     " << endl;
  cout << "  -o  : output     :: -p  : parameter       " << endl;
  cout << "  -pt : parameters , 0 default, 1 aggressive" << endl;
  cout << "                     2 stable               " << endl;
  // cout << "  -k  : Kappa(RealValue)" << endl;
  cout << "example2-1: " << argv0
       << " -o example1.result -dd example1.dat" << endl;
  cout << "example2-2: " << argv0
       << " -ds example1.dat-s -o example2.result "
       << "-p param.sdpa" << endl;
  cout << "example2-3: " << argv0
       << " -ds example1.dat-s -o example3.result "
       << "-pt 2" << endl;
  exit(1);
}
  
int main(int argc, char** argv)
{
  FILE* Display = stdout;
  setbuf(Display,NULL);

  time_t ltime;
  time( &ltime );
  cout << "SDPA start at    " << ctime(&ltime);
  // << "... (built at "<< __DATE__ << " " <<__TIME__ ")" << endl;
  // cout << "let me see your ..." << endl;

  bool isInitFile   = false;
  bool isInitSparse = false;
  bool isOutFile    = false;
  bool isDataSparse = false;
  bool isParameter  = false;
  
  char* dataFile = NULL;
  char* initFile = NULL;
  char* outFile  = NULL;
  char* paraFile = NULL;

  rParameter::parameterType parameterType =
    rParameter::PARAMETER_DEFAULT;
  
  if (argc == 1) {
    message(argv[0]);
  }
  if (argv[1][0] == '-') {
    // rsdpa argument
    
    for (int index = 0; index < argc; ++index) {
      char* target = argv[index];
      if (strcmp(target,"-dd")==0 && index+1 < argc) {
	dataFile = argv[index+1];
	index++;
	continue;
      }
      if (strcmp(target,"-ds")==0 && index+1 < argc) {
	dataFile = argv[index+1];
	index++;
	isDataSparse = true;
	continue;
      }
      if (strcmp(target,"-id")==0 && index+1 < argc) {
	initFile = argv[index+1];
	index++;
	isInitFile = true;
	continue;
      }
      if (strcmp(target,"-is")==0 && index+1 < argc) {
	initFile = argv[index+1];
	index++;
	isInitFile   = true;
	isInitSparse = true;
	continue;
      }
      if (strcmp(target,"-o")==0 && index+1 < argc) {
	outFile = argv[index+1];
	index++;
	isOutFile = true;
	continue;
      }
      if (strcmp(target,"-p")==0 && index+1 < argc) {
	paraFile = argv[index+1];
	index++;
	isParameter = true;
	continue;
      }
      if (strcmp(target,"-k")==0 && index+1 < argc) {
	KAPPA = atof(argv[index+1]);
	rMessage("Kappa = " << KAPPA);
	index++;
	continue;
      }
      if (strcmp(target,"-pt")==0 && index+1 < argc) {
	int tmp = atoi(argv[index+1]);
	switch (tmp) {
	case 0:
	  parameterType = rParameter::PARAMETER_DEFAULT;
	  break;
	case 1:
	  parameterType = rParameter::PARAMETER_AGGRESSIVE;
	  break;
	case 2:
	  parameterType = rParameter::PARAMETER_STABLE;
	  break;
	default:
	  parameterType = rParameter::PARAMETER_DEFAULT;
	}
	index++;
	paraFile = NULL;
	isParameter = false;
	continue;
      }
    }
  }
  else { // SDPA argument
    dataFile = argv[1];
    int len = strlen(dataFile);
    if (dataFile[len-1] == 's'
	&& dataFile[len-2] == '-') {
      isDataSparse = true;
    }
	
    outFile = argv[2];

    paraFile = "./param.sdpa";
    isParameter = true;

    for (int index=3; index<argc; ++index) {
      if (strcmp(argv[index],"-pt")==0 && index+1 < argc) {
	int tmp = atoi(argv[index+1]);
	switch (tmp) {
	case 0:
	  parameterType = rParameter::PARAMETER_DEFAULT;
	  break;
	case 1:
	  parameterType = rParameter::PARAMETER_AGGRESSIVE;
	  break;
	case 2:
	  parameterType = rParameter::PARAMETER_STABLE;
	  break;
	default:
	  parameterType = rParameter::PARAMETER_DEFAULT;
	}
	index++;
	paraFile = NULL;
	isParameter = false;
      } // end of "-pt"
      else {
	initFile = argv[index];
	isInitFile = true;
	int len = strlen(initFile);
	if (initFile[len-1] == 's'
	    && initFile[len-2] == '-') {
	  isInitSparse = true;
	}
      }
    } // end of 'for'
    
  }
  
  if (dataFile == NULL || outFile == NULL) {
    message(argv[0]);
  }

  cout << "data      is " << dataFile;
  if (isDataSparse) {
    cout << " : sparse" << endl;
  } else {
    cout << " : dense" << endl;
  }
  if (paraFile) {
    cout << "parameter is " << paraFile << endl;
  }
  if (outFile) {
    cout << "out       is " << outFile <<endl;
  }
  if (initFile) {
    cout << "initial   is " << initFile;
  }
  if (isInitFile) {
    if (isInitSparse) {
      cout << " : sparse" << endl;
    } else {
      cout << " : dense" << endl;
    }
  } else {
    cout << endl;
  }
  if (paraFile == NULL) {
    if (parameterType == rParameter::PARAMETER_DEFAULT) {
      cout << "set       is DEFAULT" << endl;
    }
    else if (parameterType == rParameter::PARAMETER_AGGRESSIVE) {
      cout << "set       is AGGRESSIVE" << endl;
    }
    else if (parameterType == rParameter::PARAMETER_STABLE) {
      cout << "set       is STABLE" << endl;
    }
  }
  pinpal(dataFile, initFile, outFile, paraFile, isInitFile,
  	 isInitSparse, isDataSparse, isParameter,
	 parameterType, Display);
  return 0;
}

