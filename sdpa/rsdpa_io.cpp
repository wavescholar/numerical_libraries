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
  rsdpa_io.cpp
-------------------------------------------------*/

#define DIMACS_PRINT 0

#include "rsdpa_io.h"

void rIO::read(FILE* fpData, FILE* fpout, int& m, char* str)
{
  while (true) {
    fgets(str,lengthOfString,fpData);
    if (str[0]=='*' || str[0]=='"') {
      fprintf(fpout,"%s",str);
    } else {
      sscanf(str,"%d",&m);
      break;
    }
  }
}

void rIO::read(FILE* fpData, int& nBlock)
{
  fscanf(fpData,"%d",&nBlock);
}

void rIO::read(FILE* fpData, int nBlock, int* blockStruct)
{
  for (int l=0; l<nBlock; ++l) {
    fscanf(fpData,"%*[^0-9+-]%d",&blockStruct[l]);
    if (blockStruct[l] == 1) {
      // matrices 1x1 are thougt as diagonal matrices.
      blockStruct[l] = -1;
    }
  }
}

void rIO::read(FILE* fpData, rVector& b)
{
  for (int k=0; k<b.nDim; ++k) {
    fscanf(fpData,"%*[^0-9+-]%lf",&b.ele[k]);
  }
}

void rIO::read(FILE* fpData, rBlockDenseMatrix& xMat,
	       rVector& yVec,
	       rBlockDenseMatrix& zMat, int nBlock,
	       int* blockStruct, int inputSparse)
{
  // read initial point 

  // yVec is opposite sign
  for (int k=0; k<yVec.nDim; ++k) {
    double tmp;
    fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
    yVec.ele[k] = -tmp;
    // rMessage("yVec.ele[" << k << "] = " << tmp);
  }
  
  if (inputSparse) {
    // sparse case , zMat , xMat in this order
    int i,j,l,target;
    double value;
    while (true) {
      if (fscanf(fpData,"%*[^0-9+-]%d",&target)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&l)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&i)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&j)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%lf",&value)<=0) {
	break;
      }
      #if 0
      rMessage("target = " << target
	       << ": l " << l
	       << ": i " << i
	       << ": j " << j
	       << ": value " <<value);
      #endif
      if (blockStruct[l-1]>=0) {
	if (target==1) {
	  int nCol = zMat.ele[l-1].nCol;
	  zMat.ele[l-1].de_ele[(i-1)+nCol*(j-1)] = value;
	  zMat.ele[l-1].de_ele[(j-1)+nCol*(i-1)] = value;
	} else {
	  int nCol = xMat.ele[l-1].nCol;
	  xMat.ele[l-1].de_ele[(i-1)+nCol*(j-1)] = value;
	  xMat.ele[l-1].de_ele[(j-1)+nCol*(i-1)] = value;
	}
      } else {
	if (target==1) {
	  zMat.ele[l-1].di_ele[i-1] = value;
	} else {
	  xMat.ele[l-1].di_ele[i-1] = value;
	}
      }
    } // end of 'while (true)'
  } else {
    // dense case , zMat , xMat in this order
    for (int l=0; l<nBlock; ++l) {
      int size = blockStruct[l];
      double value;
      if (size>0) {
	for (int j=0; j<size*size; ++j) {
	  fscanf(fpData,"%*[^0-9+-]%lf",&value);
	  zMat.ele[l].de_ele[j] = value;
	}
      } else {
	for (int j=0; j<-size; ++j) {
	  fscanf(fpData,"%*[^0-9+-]%lf",&value);
	  zMat.ele[l].di_ele[j] = value;
	}
      }
    }
    for (int l=0; l<nBlock; ++l) {
      int size = blockStruct[l];
      double value;
      if (size>0) {
	for (int j=0; j<size*size; ++j) {
	  fscanf(fpData,"%*[^0-9+-]%lf",&value);
	  xMat.ele[l].de_ele[j] = value;
	}
      } else {
	for (int j=0; j<-size; ++j) {
	  fscanf(fpData,"%*[^0-9+-]%lf",&value);
	  xMat.ele[l].di_ele[j] = value;
	}
      }
    }

  } // end of 'if (inputSparse)'
}

void rIO::read(FILE* fpData, rBlockSparseMatrix& C,
	       rBlockSparseMatrix* A, int m, int nBlock,
	       int* blockStruct)
{
  // read C,A
  // only case  C,A are dense
  double value;
  for (int l=0; l<nBlock; ++l) {
    int size = blockStruct[l];
    if (size>=0) {
      for (int j=0; j<size*size;++j) {
	fscanf(fpData,"%*[^0-9+-]%lf",&value);
	// C must be opposite sign
	C.ele[l].de_ele[j] = -value;
      }
    } else {
      for (int j=0; j<-size;++j) {
	fscanf(fpData,"%*[^0-9+-]%lf",&value);
	// C must be opposite sign
	C.ele[l].di_ele[j] = -value;
      }
    }
  }
  for (int k=0; k<m; ++k) {
    for (int l=0; l<nBlock; ++l) {
      int size = blockStruct[l];
      if (size>=0) {
	for (int j=0; j<size*size;++j) {
	  fscanf(fpData,"%*[^0-9+-]%lf",&value);
	  A[k].ele[l].de_ele[j] = value;
	}
      } else {
	for (int j=0; j<-size;++j) {
	  fscanf(fpData,"%*[^0-9+-]%lf",&value);
	  A[k].ele[l].di_ele[j] = value;
	}
      }
    }
  }
}

void rIO::read(FILE* fpData, int m, int nBlock,
	       int* blockStruct, int* CNonZeroCount,
	       int* ANonZeroCount,bool isDataSparse)
{
  // only count the numbers of C,A[k]
  for (int l=0; l<nBlock; ++l) {
    CNonZeroCount[l] = 0;
  }
  for (int k=0; k<m; ++k) {
    for (int l=0; l<nBlock; ++l) {
      ANonZeroCount[k*nBlock + l] = 0;
    }
  }
  if (isDataSparse) {
    int i,j,k,l;
    double value;
    while (true) {
      if (fscanf(fpData,"%*[^0-9+-]%d",&k)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&l)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&i)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&j)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%lf",&value)<=0) {
	break;
      }
      if (k==0) {
	CNonZeroCount[l-1]++;
      } else {
	ANonZeroCount[(k-1)*nBlock+(l-1)]++;
      }
      #if 0
      rMessage("k=" << k << ": l=" << l
	       << ": i=" << i << ": j=" << j
	       << ": value=" << value);
      #endif
    }// end of 'while (true)'
  } else { // isDataSparse == false
    // k==0
    for (int l=0; l<nBlock; ++l) {
      int size = blockStruct[l];
      if (size > 0) {
	for (int i=0; i<size; ++i) {
	  for (int j=0; j<size; ++j) {
	    double tmp;
	    fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	    if (i<=j && tmp!=0.0) {
	      CNonZeroCount[l]++;
	    }
	  }
	}
      } else {
	for (int j=0; j<-size; ++j) {
	  double tmp;
	  fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	  CNonZeroCount[l]++;
	}
      }
    }
    for (int k=0; k<m; ++k) {
      for (int l=0; l<nBlock; ++l) {
	int size = blockStruct[l];
	if (size > 0) {
	  for (int i=0; i<size; ++i) {
	    for (int j=0; j<size; ++j) {
	      double tmp;
	      fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	      if (i<=j && tmp!=0.0) {
		ANonZeroCount[k*nBlock+l]++;
	      }
	    }
	  }
	} else {
	  for (int j=0; j<-size; ++j) {
	    double tmp;
	    fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	    ANonZeroCount[k*nBlock+l]++;
	  }
	}
      }
    }
  } // end of 'if (isDataSparse)'
}

void rIO::read(FILE* fpData, rBlockSparseMatrix& C,
	       rBlockSparseMatrix* A,int m, int nBlock,
	       int* blockStruct, long position,bool isDataSparse)
{
  // in Sparse, read C,A[k]

  // seed the positon of C in the fpData
  fseek(fpData, position, 0);

  for (int l=0; l<nBlock; ++l) {
    if (C.ele[l].Sp_De_Di == rSparseMatrix:: SPARSE) {
      C.ele[l].NonZeroCount  = 0;
      C.ele[l].NonZeroEffect = 0;
    }
  }
  for (int k=0; k<m; ++k) {
    for (int l=0; l<nBlock; ++l) {
      if (A[k].ele[l].Sp_De_Di == rSparseMatrix:: SPARSE) {
	A[k].ele[l].NonZeroCount  = 0;
	A[k].ele[l].NonZeroEffect = 0;
      }
    }
  }

  if (isDataSparse) {
    int i,j,k,l;
    double value;
    while (true) {
      if (fscanf(fpData,"%*[^0-9+-]%d",&k)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&l)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&i)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&j)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%lf",&value)<=0) {
	break;
      }
#if 0
      rMessage("input k:" << k <<
	       " l:" << l <<
	       " i:" << i <<
	       " j:" << j);
#endif     
      if (k==0) {
	rSparseMatrix& target = C.ele[l-1];
	int count = target.NonZeroCount;
	if (target.Sp_De_Di == rSparseMatrix::SPARSE
	    && count >= target.NonZeroNumber) {
	  rError("C.ele[" << l << "]"
		 << " is too much element which is assigned");
	
	}
	if (C.blockStruct[l-1]>0) {
	  target.row_index[count] = i-1;
	  target.column_index[count] = j-1;
	  // C must be opposite sign.
	  target.sp_ele[count] = -value;
	  target.NonZeroCount++;
	  if (i==j) {
	    target.NonZeroEffect++;
	  } else {
	    target.NonZeroEffect += 2;
	  }
	} else {
	  target.di_ele[j-1] = -value;
	}
      } else {
	rSparseMatrix& target = A[k-1].ele[l-1];
	int count = target.NonZeroCount;
	if (target.Sp_De_Di == rSparseMatrix::SPARSE
	    && count >= target.NonZeroNumber) {
	  rError("A[" << k << "].ele[" << l << "]"
		 << " is too much element which is assigned");
	}
	if (A[k-1].blockStruct[l-1]>0) {
	  target.row_index[count] = i-1;
	  target.column_index[count] = j-1;
	  target.sp_ele[count] = value;
	  target.NonZeroCount++;
	  if (i==j) {
	    target.NonZeroEffect++;
	  } else {
	    target.NonZeroEffect += 2;
	  }
	} else {
	  target.di_ele[j-1] = value;
	}
      }
    }// end of 'while (true)'
  } else { // isDataSparse == false
    // k==0
    for (int l=0; l<nBlock; ++l) {
      int size = blockStruct[l];
      if (size > 0) {
	for (int i=0; i<size; ++i) {
	  for (int j=0; j<size; ++j) {
	    double tmp = 0.0;
	    fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	    if (i<=j && tmp!=0.0) {
	      rSparseMatrix& target = C.ele[l];
	      int count = target.NonZeroCount;
	      target.row_index[count] = i;
	      target.column_index[count] = j;
	      // C must be opposite sign.
	      target.sp_ele[count] = -tmp;
	      target.NonZeroCount++;
	      if (i==j) {
		target.NonZeroEffect++;
	      } else {
		target.NonZeroEffect += 2;
	      }
	    }
	  }
	}
      } else { // size < 0
	for (int j=0; j<-size; ++j) {
	  double tmp;
	  fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	  // C must be opposite sign.
	  C.ele[l].di_ele[j] = -tmp;
	}
      }
    } // end of 'for (int l)'

    for (int k=0; k<m; ++k) {
      for (int l=0; l<nBlock; ++l) {
	int size = blockStruct[l];
	if (size > 0) {
	  for (int i=0; i<size; ++i) {
	    for (int j=0; j<size; ++j) {
	      double tmp;
	      fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	      if (i<=j && tmp!=0.0) {
		rSparseMatrix& target = A[k].ele[l];
		int count = target.NonZeroCount;
		target.row_index[count] = i;
		target.column_index[count] = j;
		target.sp_ele[count] = tmp;
		target.NonZeroCount++;
		if (i==j) {
		  target.NonZeroEffect++;
		} else {
		  target.NonZeroEffect += 2;
		}
	      }
	    }
	  }
	} else { // size < 0
	  for (int j=0; j<-size; ++j) {
	    double tmp;
	    fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	    A[k].ele[l].di_ele[j] = tmp;
	  }
	}
      } // end of 'for (int l)'
	
    } // end of 'for (int k)'
  } // end of 'if (isDataSparse)'
}

void rIO::printHeader(FILE* fpout, FILE* Display)
{
  if (fpout) {
    fprintf(fpout,"   mu       thetaP   thetaD   objP       objD "

	    "      alphaP   alphaD   beta \n");

  }
  if (Display) {
    fprintf(Display,"   mu       thetaP   thetaD   objP       objD "
	    "      alphaP   alphaD   beta \n");
  }
}

void rIO::printOneIteration(int pIteration,
			    rAverageComplementarity& mu,
			    rRatioInitResCurrentRes& theta,
			    rSolveInfo& solveInfo,
			    rStepLength& alpha,
			    rDirectionParameter& beta,
			    rResiduals& currentRes,
			    FILE* fpout,
			    FILE* Display)
{
  #if REVERSE_PRIMAL_DUAL
  if (Display) {
    fprintf(Display,"%2d %4.1e %4.1e %4.1e %+7.2e %+7.2e"
	    " %4.1e %4.1e %4.2e\n", pIteration, mu.current,
	    theta.dual, theta.primal,
	    -solveInfo.objValDual,-solveInfo.objValPrimal,
	    alpha.dual, alpha.primal, beta.value);
  }
  if (fpout) {
    fprintf(fpout,"%2d %4.1e %4.1e %4.1e %+7.2e %+7.2e"
	    " %4.1e %4.1e %4.2e\n", pIteration, mu.current,
	    theta.dual, theta.primal,
	    -solveInfo.objValDual,-solveInfo.objValPrimal,
	    alpha.dual, alpha.primal, beta.value);
  }
  #else
  if (Display) {
    fprintf(Display,"%2d %4.1e %4.1e %4.1e %+7.2e %+7.2e"
	    " %4.1e %4.1e %4.2e\n", pIteration, mu.current,
	    theta.primal, theta.dual,
	    solveInfo.objValPrimal, solveInfo.objValDual,
	    alpha.primal, alpha.dual, beta.value);
  }
  if (fpout) {
    fprintf(fpout,"%2d %4.1e %4.1e %4.1e %+7.2e %+7.2e"
	    " %4.1e %4.1e %4.2e\n", pIteration, mu.current,
	    theta.primal, theta.dual,
	    solveInfo.objValPrimal, solveInfo.objValDual,
	    alpha.primal, alpha.dual, beta.value);
  }
  #endif
}

void rIO::printLastInfo(int pIteration,
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
			bool printTime)
{
  printOneIteration(pIteration,mu,theta,solveInfo,alpha,
		    beta,currentRes, fpout, Display);

  double mean = (fabs(solveInfo.objValPrimal)
		 + fabs(solveInfo.objValDual)) / 2.0;
  double PDgap = fabs(solveInfo.objValPrimal
		      - solveInfo.objValDual);
  // double dominator;
  double relgap;
  if (mean < 1.0) {
    relgap = PDgap;
  } else {
    relgap = PDgap/mean;
  }

  double gap    = mu.current*nDim;
  double digits = 1000; // 1000 means infinity in this case
  digits = -log10(fabs(PDgap/mean));

  #if DIMACS_PRINT
  double b1 = 0.0;
  for (int k=0; k<b.nDim; ++k) {
    double tmp = fabs(b.ele[k]);
    if (b1 < tmp) {
      b1 = tmp;
    }
  }
  double c1 = 0.0;
  for (int l=0; l<C.nBlock; ++l) {
    rSparseMatrix& Cl = C.ele[l];
    if (Cl.Sp_De_Di == rSparseMatrix::SPARSE) {
      for (int i=0; i<Cl.NonZeroCount; ++i) {
	double tmp  = fabs(Cl.sp_ele[i]);
	if (c1 < tmp) {
	  c1 = tmp;
	}
      }
    }
    else if (Cl.Sp_De_Di == rSparseMatrix::DENSE) {
      for (int i=0; i<Cl.nRow*Cl.nCol; ++i) {
	double tmp  = fabs(Cl.de_ele[i]);
	if (c1 < tmp) {
	  c1 = tmp;
	}
      }
    }
    else if (Cl.Sp_De_Di == rSparseMatrix::DIAGONAL) {
      for (int i=0; i<Cl.nCol; ++i) {
	double tmp  = fabs(Cl.di_ele[i]);
	if (c1 < tmp) {
	  c1 = tmp;
	}
      }
    }
  }
  double p_norm = dnrm2_(&currentRes.primalVec.nDim,
			 currentRes.primalVec.ele, &IONE);
  double d_norm = 0.0;
  for (int l=0; l<currentRes.dualMat.nBlock; ++l) {
    int size = C.blockStruct[l];
    double tmp;
    if (size < 0) {
      size = -size;
      tmp = dnrm2_(&size,currentRes.dualMat.ele[l].di_ele, &IONE);
      d_norm += tmp;
    }
    else {
      size = size * size;
      tmp = dnrm2_(&size,currentRes.dualMat.ele[l].de_ele, &IONE);
      d_norm += tmp;
    }
  }

  currentPt.workMat.copyFrom(currentPt.xMat);
  double x_min =  rAl::getMinEigenValue(currentPt.workMat,
					currentPt.xzEigenValues,
					currentPt.workVec);
  currentPt.workMat.copyFrom(currentPt.zMat);
  double z_min =  rAl::getMinEigenValue(currentPt.workMat,
					currentPt.xzEigenValues,
					currentPt.workVec);
  double ctx = solveInfo.objValPrimal;
  double bty = solveInfo.objValDual;
  double xtz = 0.0;
  rAl::let(xtz,'=',currentPt.xMat,'.',currentPt.zMat);

  double err1 = p_norm / (1+b1);
  double err2 = max( 0.0, - x_min / (1+b1));
  double err3 = d_norm / (1+c1);
  double err4 = max( 0.0, - z_min / (1+c1));
  double err5 = (ctx - bty) / (1 + fabs(ctx) + fabs(bty));
  double err6 = xtz / (1 + fabs(ctx) + fabs(bty));

    
  #endif

  
  
  
  if (Display) {
    fprintf(Display, "\n");
    phase.display(Display);
    fprintf(Display, "   Iteration = %d\n",  pIteration);
    fprintf(Display, "          mu = %4.16e\n",  mu.current);
    fprintf(Display, "relative gap = %4.16e\n", relgap);
    fprintf(Display, "         gap = %4.16e\n",  gap);
    fprintf(Display, "      digits = %4.16e\n",  digits);
    #if REVERSE_PRIMAL_DUAL
    fprintf(Display, "objValPrimal = %10.16e\n",
	    -solveInfo.objValDual);
    fprintf(Display, "objValDual   = %10.16e\n",
	    -solveInfo.objValPrimal);
    fprintf(Display, "p.feas.error = %10.16e\n",
	    currentRes.normDualMat);
    fprintf(Display, "d.feas.error = %10.16e\n",
	    currentRes.normPrimalVec);
    #else
    fprintf(Display, "objValPrimal = %10.16e\n",
	    solveInfo.objValPrimal);
    fprintf(Display, "objValDual   = %10.16e\n",
	    solveInfo.objValDual);
    fprintf(Display, "p.feas.error = %10.16e\n",
	    currentRes.normPrimalVec);
    fprintf(Display, "d.feas.error = %10.16e\n",
	    currentRes.normDualMat);
    #endif
    if (printTime == true) {
      fprintf(Display, "total time   = %.3f\n",cputime);
    }
    #if DIMACS_PRINT
    fprintf(Display, "\n");
    fprintf(Display, "* DIMACS_ERROS * \n");
    fprintf(Display, "err1 = %4.16e  [%40s]\n",
	    err1, "||Ax-b|| / (1+||b||_1) ");
    fprintf(Display, "err2 = %4.16e  [%40s]\n",
	    err2, "max(0, -lambda(x) / (1+||b||_1))");
    fprintf(Display, "err3 = %4.16e  [%40s]\n",
	    err3, "||A^Ty + z - c || / (1+||c||_1) ");
    fprintf(Display, "err4 = %4.16e  [%40s]\n",
	    err4, "max(0, -lambda(z) / (1+||c||_1))");
    fprintf(Display, "err5 = %4.16e  [%40s]\n",
	    err5,"(<c,x> - by) / (1 + |<c,x>| + |by|)");
    fprintf(Display, "err6 = %4.16e  [%40s]\n",
	    err6,"<x,z> / (1 + |<c,x>| + |by|)");
    fprintf(Display, "\n");
    #endif
  }
  if (fpout) {
    fprintf(fpout, "\n");
    phase.display(fpout);
    fprintf(fpout, "   Iteration = %d\n",  pIteration);
    fprintf(fpout, "          mu = %4.16e\n",  mu.current);
    fprintf(fpout, "relative gap = %4.16e\n", relgap);
    fprintf(fpout, "         gap = %4.16e\n",  gap);
    fprintf(fpout, "      digits = %4.16e\n",  digits);
    #if REVERSE_PRIMAL_DUAL
    fprintf(fpout, "objValPrimal = %10.10e\n",
	    -solveInfo.objValDual);
    fprintf(fpout, "objValDual   = %10.10e\n",
	    -solveInfo.objValPrimal);
    fprintf(fpout, "p.feas.error = %10.10e\n",
	    currentRes.normDualMat);
    fprintf(fpout, "d.feas.error = %10.10e\n",
	    currentRes.normPrimalVec);
    #else
    fprintf(fpout, "objValPrimal = %10.10e\n",
	    solveInfo.objValPrimal);
    fprintf(fpout, "objValDual   = %10.10e\n",
	    solveInfo.objValDual);
    fprintf(fpout, "p.feas.error = %10.10e\n",
	    currentRes.normPrimalVec);
    fprintf(fpout, "d.feas.error = %10.10e\n",
	    currentRes.normDualMat);
    #endif
    fprintf(fpout, "total time   = %.3f\n",cputime);
    #if DIMACS_PRINT
    fprintf(fpout, "\n");
    fprintf(fpout, "* DIMACS_ERROS * \n");
    fprintf(fpout, "err1 = %4.16e  [%40s]\n",
	    err1, "||Ax-b|| / (1+||b||_1) ");
    fprintf(fpout, "err2 = %4.16e  [%40s]\n",
	    err2, "max(0, -lambda(x) / (1+||b||_1))");
    fprintf(fpout, "err3 = %4.16e  [%40s]\n",
	    err3, "||A^Ty + z - c || / (1+||c||_1) ");
    fprintf(fpout, "err4 = %4.16e  [%40s]\n",
	    err4, "max(0, -lambda(z) / (1+||c||_1))");
    fprintf(fpout, "err5 = %4.16e  [%40s]\n",
	    err5,"(<c,x> - by) / (1 + |<c,x>| + |by|)");
    fprintf(fpout, "err6 = %4.16e  [%40s]\n",
	    err6,"<x,z> / (1 + |<c,x>| + |by|)");
    fprintf(fpout, "\n");
    #endif

    fprintf(fpout, "\n\nParameters are\n");
    param.display(fpout);
    com.display(fpout);

    #if 1
    #if REVERSE_PRIMAL_DUAL
    fprintf(fpout,"xVec = \n");
    currentPt.yVec.display(fpout,-1.0);
    fprintf(fpout,"xMat = \n");
    currentPt.zMat.display(fpout);
    fprintf(fpout,"yMat = \n");
    currentPt.xMat.display(fpout);
    #else
    currentPt.display(fpout);
    #endif
    #endif
  }
}

