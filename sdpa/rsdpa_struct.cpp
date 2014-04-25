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
/*----------------------------------------
  rsdpa_struct.cpp
----------------------------------------*/

#include "rsdpa_struct.h"
// printing presicion of such as vector 
#define P_FORMAT "%+8.3e"

#if NON_ATLAS_SDPA
#define catlas_dset(dset_length,dset_value,dset_pointer,dset_step) \
{for (int dset_i=0,dset_index = 0; dset_i<dset_length; ++dset_i) { \
  dset_pointer[dset_index] = dset_value; \
  dset_index += dset_step; \
}}
#endif

rVector::rVector()
{
  nDim = 0;
  ele  = NULL;
}

rVector::rVector(int nDim, double value)
{
  ele  = NULL;
  initialize(nDim,value);
}

rVector::~rVector()
{
  if (ele) {
    delete[] ele;
  }
  ele = NULL;
}

void rVector::initialize(int nDim,double value)
{
  // rMessage("rVector initialize");
  if (ele && this->nDim!=nDim) {
    if (ele) {
      delete[] ele;
      ele = NULL;
    }
    if (nDim<=0) {
      rError("rVector:: nDim is nonpositive");
    }
  }
  this->nDim = nDim;
  if (ele==NULL) {
    ele = NULL;
    rNewCheck();
    ele = new double[nDim];
    if (ele==NULL) {
      rError("rVector:: memory exhausted");
    }
  }
  catlas_dset(nDim,value,ele,IONE);
}

void rVector::initialize(double value)
{
  if (nDim<=0) {
    rError("rVector:: nDim is nonpositive");
  }
  if (ele==NULL) {
    rNewCheck();
    ele = new double[nDim];
    if (ele==NULL) {
      rError("rVector:: memory exhausted");
    }
  }
  catlas_dset(nDim,value,ele,IONE);
}

void rVector::setZero()
{
  initialize(0.0);
}

void rVector::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }
  fprintf(fpout,"{");
  for (int j=0; j<nDim-1; ++j) {
    fprintf(fpout,P_FORMAT",",ele[j]);
  }
  if (nDim>0) {
    fprintf(fpout,P_FORMAT"}\n",ele[nDim-1]);
  } else {
    fprintf(fpout,"  }\n");
  }
}

void rVector::display(FILE* fpout,double scalar)
{
  if (fpout == NULL) {
    return;
  }
  fprintf(fpout,"{");
  for (int j=0; j<nDim-1; ++j) {
    fprintf(fpout,P_FORMAT",",ele[j]*scalar);
  }
  if (nDim>0) {
    fprintf(fpout,P_FORMAT"}\n",ele[nDim-1]*scalar);
  } else {
    fprintf(fpout,"  }\n");
  }
}

bool rVector::copyFrom(rVector& other)
{
  if (this == &other) {
    return _SUCCESS;
  }
  if (nDim != other.nDim && ele!=NULL) {
    delete[] ele;
    ele = NULL;
  }
  nDim = other.nDim;
  if (nDim<=0) {
    rError("rVector:: nDim is nonpositive");
  }
  if (ele==NULL) {
    rNewCheck();
    ele = new double[nDim];
    if (ele==NULL) {
      rError("rVector:: memory exhausted");
    }
  }
  //bbc
  dcopy(&nDim,other.ele,&IONE,ele,&IONE);
  return _SUCCESS;
}

rBlockVector::rBlockVector()
{
  nBlock = 0;
  blockStruct = NULL;
  ele = NULL;
}

rBlockVector::rBlockVector(int nBlock, int* blockStruct,
			   double value)
{
  initialize(nBlock,blockStruct,value);
}

rBlockVector::~rBlockVector()
{
  if (ele && blockStruct && nBlock>=0) {
    for (int l=0; l<nBlock; ++l) {
      ele[l].~rVector();
    }
    delete[] ele;
    ele = NULL;

    delete[] blockStruct;
    blockStruct = NULL;
  }
}

void rBlockVector::initialize(int nBlock, int* blockStruct,
			      double value)
{
  // rMessage("rBlockVector initialize");
  this->nBlock = nBlock;
  if (nBlock<=0) {
    rError("rBlockVector:: nBlock is nonpositive");
  }
  this->blockStruct = NULL;
  rNewCheck();
  this->blockStruct = new int[nBlock];
  if (this->blockStruct==NULL) {
    rError("rBlockVector:: memory exhausted");
  }
  for (int l=0; l<nBlock; ++l) {
    this->blockStruct[l] = blockStruct[l];
  }

  ele = NULL;
  rNewCheck();
  ele = new rVector[nBlock];
  if (ele==NULL) {
    rError("rBlockVector:: memory exhausted");
  }
  for (int l=0; l<nBlock; ++l) {
    int size = blockStruct[l];
    if (size<0) {
      size = -size;
    }
    ele[l].initialize(size,value);
  }
}

void rBlockVector::initialize(double value)
{
  if (nBlock>0 && blockStruct && ele) {
    for (int l=0; l<nBlock; ++l) {
      ele[l].initialize(value);
    }
  }
}

void rBlockVector::setZero()
{
  if (nBlock>0 && blockStruct && ele) {
    for (int l=0; l<nBlock; ++l) {
      ele[l].setZero();
    }
  }
}

void rBlockVector::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }
  fprintf(fpout,"{ ");
  if (nBlock>0 && blockStruct && ele) {
    for (int l=0; l<nBlock; ++l) {
      ele[l].display(fpout);
    }
  }
  fprintf(fpout,"} \n");
}

bool rBlockVector::copyFrom(rBlockVector& other)
{
  if (this == &other) {
    return _SUCCESS;
  }
  
  if (other.nBlock<=0) {
    rError("rBlockVector:: nBlock is nonpositive");
  }
  if (nBlock!=other.nBlock && blockStruct) {
    delete[] blockStruct;
    blockStruct = NULL;
    delete[] ele;
    ele = NULL;
  }
  if (blockStruct==NULL) {
    nBlock = other.nBlock;
    rNewCheck();
    blockStruct = new int[nBlock];
    if (blockStruct==NULL) {
      rError("rBlockVector:: memory exhausted");
    }
    for (int l=0; l<nBlock; ++l) {
      blockStruct[l] = other.blockStruct[l];
    }
  }
  if (ele==NULL) {
    rNewCheck();
    ele = new rVector[nBlock];
    if (ele==NULL) {
      rError("rBlockVector:: memory exhausted");
    }
  }
  for (int l=0; l<nBlock; ++l) {
    ele[l].copyFrom(other.ele[l]);
  }
  return _SUCCESS;
}

rSparseMatrix::rSparseMatrix()
{
  nRow = 0;
  nCol = 0;
  Sp_De_Di = SPARSE;

  NonZeroNumber = 0;
  
  de_ele = NULL;

  row_index     = NULL;
  column_index  = NULL;
  sp_ele        = NULL;
  NonZeroCount  = 0;
  NonZeroEffect = 0;

  di_ele = NULL;
}

rSparseMatrix::rSparseMatrix(int nRow, int nCol,
			     rSparseMatrix::rSpMat_Sp_De_Di Sp_De_Di,
			     int NonZeroNumber)
{
  initialize(nRow, nCol, Sp_De_Di, NonZeroNumber);
}

rSparseMatrix::~rSparseMatrix()
{
  if (de_ele) {
    delete[] de_ele;
    de_ele = NULL;
  }
  if (row_index) {
    delete[] row_index;
    row_index = NULL;
  }
  if (column_index) {
    delete[] column_index;
    column_index = NULL;
  }
  if (sp_ele) {
    delete[] sp_ele;
    sp_ele = NULL;
  }
  if (di_ele) {
    delete[] di_ele;
    di_ele = NULL;
  }
}

void rSparseMatrix::
initialize(int nRow, int nCol,
	   rSparseMatrix::rSpMat_Sp_De_Di Sp_De_Di,
	   int NonZeroNumber)
{
  // rMessage("rSparseMatrix initialize");

  rSparseMatrix();
  if (nRow<=0 || nCol<=0) {
    rError("rSparseMatrix:: Dimensions are nonpositive");
  }
  this->nRow          = nRow;
  this->nCol          = nCol;
  this->Sp_De_Di      = Sp_De_Di;

  int length;
	row_index = NULL;
	column_index = NULL;
	sp_ele = NULL;
  switch(Sp_De_Di) {
  case SPARSE:
    this->NonZeroNumber  = NonZeroNumber;
    this->NonZeroCount   = 0;
    this->NonZeroEffect  = 0;
    if (NonZeroNumber > 0) {
      rNewCheck();
      row_index    = new int[NonZeroNumber];
      rNewCheck();
      column_index = new int[NonZeroNumber];
      rNewCheck();
      sp_ele       = new double[NonZeroNumber];
      if (row_index==NULL || column_index==NULL
	  || sp_ele==NULL) {
	rError("rSparseMatrix:: memory exhausted");
      }
    }
    break;
  case DENSE:
    this->NonZeroNumber = nRow*nCol;
    this->NonZeroCount  = nRow*nCol;
    this->NonZeroEffect = nRow*nCol;
    rNewCheck();
    de_ele = new double[NonZeroNumber];
    if (de_ele==NULL) {
      rError("rSparseMatrix:: memory exhausted");
    }
    length = nRow*nCol;
    catlas_dset(length,DZERO,de_ele,IONE);
    // all elements are 0.
    break;
  case DIAGONAL:
    if (nRow!=nCol) {
      rError("rSparseMatrix:: Diagonal must be Square matrix");
    }
    this->NonZeroNumber = nCol;
    this->NonZeroCount  = nCol;
    this->NonZeroEffect = nCol;
    if (di_ele==NULL) {
      rNewCheck();
      di_ele = new double[NonZeroNumber];
      if (di_ele==NULL) {
	rError("rSparseMatrix:: memory exhausted");
      }
    }
    catlas_dset(nCol,DZERO,di_ele,IONE);
    // all elements are 0.

    break;
  }
}

void rSparseMatrix::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }
  int index=0,j=0,i=0;
  switch(Sp_De_Di) {
  case SPARSE:
    fprintf(fpout,"{");
    for (index=0; index<NonZeroCount; ++index) {
      int i        = row_index[index];
      int j        = column_index[index];
      double value = sp_ele[index];
      fprintf(fpout,"val[%d,%d] = "P_FORMAT"\n", i,j,value);
    }
    fprintf(fpout,"}\n");
    break;
  case DENSE:
    fprintf(fpout,"{\n");
    for (i=0; i<nRow-1; ++i) {
      if (i==0) {
	fprintf(fpout," ");
      } else {
	fprintf(fpout,"  ");
      }
      fprintf(fpout,"{");
      for (int j=0; j<nCol-1; ++j) {
	fprintf(fpout, P_FORMAT",",de_ele[i+nCol*j]);
      }
      fprintf(fpout,P_FORMAT" },\n",de_ele[i+nCol*(nCol-1)]);
    }
    if (nRow>1) {
      fprintf(fpout,"  {");
    }
    for (j=0; j<nCol-1; ++j) {
      fprintf(fpout,P_FORMAT",",de_ele[(nRow-1)+nCol*j]);
    }
    fprintf(fpout,P_FORMAT" }",de_ele[(nRow-1)+nCol*(nCol-1)]);
    if (nRow>1) {
      fprintf(fpout,"   }\n");
    } else {
      fprintf(fpout,"\n");
    }
    break;
  case DIAGONAL:
    fprintf(fpout,"{");
    for (j=0; j<nCol-1; ++j) {
      fprintf(fpout, P_FORMAT",",di_ele[j]);
    }
    if (nCol>0) {
      fprintf(fpout, P_FORMAT",",di_ele[nCol-1]);
    }
    fprintf(fpout,"}\n");
    break;
  }
}

bool rSparseMatrix::copyFrom(rSparseMatrix& other)
{
  if (Sp_De_Di != other.Sp_De_Di || nRow != other.nRow
      || nCol != other.nCol) {
    this->~rSparseMatrix();
    initialize(other.nRow,other.nCol,other.Sp_De_Di,
	       NonZeroNumber);
    NonZeroCount  = other.NonZeroCount;
    NonZeroEffect = other.NonZeroEffect;
    int length,index=0;
    switch(Sp_De_Di) {
    case SPARSE:
      for (index = 0; index<NonZeroCount;++index) {
	row_index[index]    = other.row_index[index];
	column_index[index] = other.column_index[index];
	sp_ele[index]       = other.sp_ele[index];
      }
      break;
    case DENSE:
      length = nRow*nCol;
      dcopy(&length,other.de_ele,&IONE,de_ele,&IONE);
      break;
    case DIAGONAL:
      dcopy(&nCol,other.di_ele,&IONE,di_ele,&IONE);
      break;
    }
  } else { // Sp_De_Di == other.Sp_De_Di
           // && nRow == other.nRow && nCol == other.nCol
    NonZeroCount  = other.NonZeroCount;
    NonZeroEffect = other.NonZeroEffect;
    int length,index=0;
    switch(Sp_De_Di) {
    case SPARSE:
      if (NonZeroNumber!=other.NonZeroNumber) {
	delete[] row_index;
	delete[] column_index;
	delete[] sp_ele;
	row_index = column_index = NULL;
	sp_ele = NULL;
	rNewCheck();
	row_index    = new int[NonZeroNumber];
	rNewCheck();
	column_index = new int[NonZeroNumber];
	rNewCheck();
	sp_ele       = new double[NonZeroNumber];
	if (row_index==NULL || column_index==NULL
	    || sp_ele==NULL) {
	  rError("rSparseMatrix:: memory exhausted");
	}
      }
      for (index = 0; index<NonZeroCount;++index) {
	row_index[index]    = other.row_index[index];
	column_index[index] = other.column_index[index];
	sp_ele[index]       = other.sp_ele[index];
      }
      break;
    case DENSE:
      length = nRow*nCol;
      dcopy(&length,other.de_ele,&IONE,de_ele,&IONE);
      break;
    case DIAGONAL:
      dcopy(&nCol,other.di_ele,&IONE,di_ele,&IONE);
      break;
    } // end of switch
  } // end of else
  return _SUCCESS;
}

void rSparseMatrix::changeToDense(bool forceChange)
{
  if (Sp_De_Di!=SPARSE) {
    return;
  }
  // if (false)
  // rMessage(" NonZeroCount " << NonZeroCount);
  // rMessage(" nRow*nCol*0.2 " << nRow*nCol*0.2);
  if (forceChange == false && NonZeroCount < (nRow*nCol) * 0.20) {
    // if the number of elements are less than 20 percent,
    // we don't change to Dense.
    return;
  }
  // rMessage("change");
  Sp_De_Di = DENSE;
  de_ele = NULL;
  int length = nRow*nCol;
  rNewCheck();
  de_ele = new double[length];
  if (de_ele==NULL) {
    rError("rSparseMatrix:: memory exhausted");
  }
  catlas_dset(length,DZERO,de_ele,IONE);
  // all elements are set 0.
  for (int index=0; index<NonZeroCount; ++index) {
    int        i = row_index[index];
    int        j = column_index[index];
    double value = sp_ele[index];
    if (i==j) {
      de_ele[i+nCol*j] = value;
    } else {
      de_ele[i+nCol*j] = de_ele[j+nCol*i] = value;
    }
  }
  NonZeroCount = NonZeroNumber = NonZeroEffect = length;

  if (row_index != NULL) {
		 delete[] row_index;
  }
  if (column_index != NULL) {
		 delete[] column_index;
  }
  if (sp_ele != NULL) {
		 delete[] sp_ele;
  }
  row_index = NULL;
  column_index = NULL;
  sp_ele = NULL;
}

void rSparseMatrix::setZero()
{
  int length;
  switch(Sp_De_Di) {
  case SPARSE:
    NonZeroCount  = 0;
    NonZeroEffect = 0;
    // No element is stored.
    break;
  case DENSE:
    length = nRow*nCol;
    catlas_dset(length,DZERO,de_ele,IONE);
    break;
  case DIAGONAL:
    catlas_dset(nCol,DZERO,di_ele,IONE);
    break;
  }
}

void rSparseMatrix::setIdentity(double scalar)
{
  if (nRow != nCol) {
    rError("rSparseMatrix:: Identity matrix must be square matrix");
  }
  int length,step,index=0;
  switch(Sp_De_Di) {
  case SPARSE:
    if (nCol > NonZeroNumber) {
      rError("rSparseMatrix:: cannot store over NonZeroNumber");
      // the number of Diagonal elements equals nCol.
    }
    NonZeroCount  = nCol;
    NonZeroEffect = nCol;
    for (index=0; index< NonZeroCount; ++index) {
      row_index[index]    = index;
      column_index[index] = index;
      sp_ele[index]       = scalar;
    }
    break;
  case DENSE:
    length = nRow*nCol;
    catlas_dset(length,DZERO,de_ele,IONE);
    step = nCol+1;
    catlas_dset(nCol,scalar,de_ele,step);
    // only diagonal elements are the value of scalar.
    break;
  case DIAGONAL:
    catlas_dset(nCol,scalar,di_ele,IONE);
    break;
  }
}
    
bool rSparseMatrix::sortSparseIndex(int& i, int& j)
{
  // if this matrix is not symmetric,
  // return the index(i,j) whose values are not symmetric.
  i = -1;
  j = -1;
  const double tolerance = 1.0e-8;
  int i1=0;
  switch(Sp_De_Di) {
  case SPARSE:
    // Make matrix as Upper Triangluar
    for (i1=0; i1<NonZeroCount; ++i1) {
      int tmpi = row_index[i1];
      int tmpj = column_index[i1];
      if (tmpi>tmpj) {
	row_index   [i1] = tmpj;
	column_index[i1] = tmpi;
      }
    }
    // simple sort
    for (i1=0; i1<NonZeroCount; ++i1) {
      for (int i2=0; i2<i1; ++i2) {
	int index1 = row_index[i1]+nCol*column_index[i1];
	int index2 = row_index[i2]+nCol*column_index[i2];
	if (index1<index2) {
	  int         tmpi = row_index   [i2];
	  int         tmpj = column_index[i2];
	  double      tmpv = sp_ele      [i2];
	  row_index   [i2] = row_index   [i1];
	  column_index[i2] = column_index[i1];
	  sp_ele      [i2] = sp_ele      [i1];
	  row_index   [i1] = tmpi;
	  column_index[i1] = tmpj;
	  sp_ele      [i1] = tmpv;
	}
      }
    }
    // the process for the same index
    for (i1=0; i1<NonZeroCount-1; ++i1) {
      int index1 = row_index[i1  ]+nCol*column_index[i1  ];
      int index2 = row_index[i1+1]+nCol*column_index[i1+1];
      if (index1 == index2) {
	if (fabs(sp_ele[index1] - sp_ele[index2]) > tolerance) {
	  // Here must not be symmetric
	  if (i<0 || j<0) {
	    i = row_index   [i1];
	    j = column_index[i1];
	  }
	}

	// remove redudunt
	for (int i2 = i1+1; i2<NonZeroCount-2;++i2) {
	  row_index   [i2] = row_index   [i2+1];
	  column_index[i2] = column_index[i2+1];
	  sp_ele      [i2] = sp_ele      [i2+1];
	}
	NonZeroCount--;
	if (i==j) {
	  NonZeroEffect--;
	} else {
	  NonZeroEffect -= 2;
	}
      } // end of 'if (index1==index2)'
    }
    break;
  case DENSE:
    if (nRow!=nCol) {
      return FAILURE;
    }
    for (j=1; j<nCol; ++j) {
      for (i=0; i<j; ++i) {
	if (fabs(de_ele[i+nCol*j]-de_ele[j+nCol*i]) > tolerance) {
	  return FAILURE;
	}
      }
    }
    break;
  case DIAGONAL:
    // Nothing needs.
    break;
  }
  return _SUCCESS;
}


rDenseMatrix::rDenseMatrix()
{
  nRow = 0;
  nCol = 0;
  De_Di = DENSE;

  de_ele = NULL;
  di_ele = NULL;
}

rDenseMatrix::rDenseMatrix(int nRow, int nCol,
			   rDenseMatrix::rDeMat_De_Di Sp_De_Di)
{
  initialize(nRow, nCol, Sp_De_Di);
}

rDenseMatrix::~rDenseMatrix()
{
  if (de_ele) {
    delete[] de_ele;
    de_ele = NULL;
  }
  if (di_ele) {
    delete[] di_ele;
    di_ele = NULL;
  }
}

void rDenseMatrix::
initialize(int nRow, int nCol,
	   rDenseMatrix::rDeMat_De_Di De_Di)
{
  // rMessage("rDenseMatrix::initialize");

  rDenseMatrix();
  if (nRow<=0 || nCol<=0) {
    rError("rDenseMatrix:: Dimensions are nonpositive");
  }
  int old_length = this->nRow*this->nCol;
  if (this->De_Di==DIAGONAL) {
    old_length = this->nRow;
  }
  this->nRow  = nRow;
  this->nCol  = nCol;
  this->De_Di = De_Di;

  int length;
  switch(De_Di) {
  case DENSE:
    length = nRow*nCol;
    if (de_ele && old_length!=length) {
      delete[] de_ele;
      de_ele = NULL;
    }
    if (de_ele==NULL) {
      rNewCheck();
      de_ele = new double[length];
      if (de_ele==NULL) {
	rError("rDenseMatrix:: memory exhausted");
      }
    }
    catlas_dset(length,DZERO,de_ele,IONE);
    break;
  case DIAGONAL:
    if (nRow!=nCol) {
      rError("rDenseMatrix:: Diagonal must be Square matrix");
    }
    if (di_ele && old_length!=nRow) {
      delete[] di_ele;
      di_ele = NULL;
    }
    if (di_ele==NULL) {
      rNewCheck();
      di_ele = new double[nCol];
      if (di_ele==NULL) {
	rError("rDenseMatrix:: memory exhausted");
      }
    }
    catlas_dset(nCol,DZERO,di_ele,IONE);
    break;
  }
}

void rDenseMatrix::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }
  int i=0,j=0;
  switch(De_Di) {
  case DENSE:
    fprintf(fpout,"{");
    for (i=0; i<nRow-1; ++i) {
      if (i==0) {
	fprintf(fpout," ");
      } else {
	fprintf(fpout,"  ");
      }
      fprintf(fpout,"{");
      for (int j=0; j<nCol-1; ++j) {
	fprintf(fpout, P_FORMAT",",de_ele[i+nCol*j]);
      }
      fprintf(fpout,P_FORMAT" },\n",de_ele[i+nCol*(nCol-1)]);
    }
    if (nRow>1) {
      fprintf(fpout,"  {");
    }
    for (j=0; j<nCol-1; ++j) {
      fprintf(fpout,P_FORMAT",",de_ele[(nRow-1)+nCol*j]);
    }
    fprintf(fpout,P_FORMAT" }",de_ele[(nRow-1)+nCol*(nCol-1)]);
    if (nRow>1) {
      fprintf(fpout,"   }\n");
    } else {
      fprintf(fpout,"\n");
    }
    break;
  case DIAGONAL:
    fprintf(fpout,"{");
    for (int j=0; j<nCol-1; ++j) {
      fprintf(fpout, P_FORMAT",",di_ele[j]);
    }
    if (nCol>0) {
      fprintf(fpout, P_FORMAT"}\n",di_ele[nCol-1]);
    }
    break;
  }
}

bool rDenseMatrix::copyFrom(rSparseMatrix& other)
{
  int length,index=0;
  switch(other.Sp_De_Di) {
  case rSparseMatrix::SPARSE:
    De_Di = DENSE;
    if (de_ele) {
      delete[] de_ele;
    }
    de_ele = NULL;
    nRow = other.nRow;
    nCol = other.nCol;
    rNewCheck();
    de_ele = new double[nRow*nCol];
    if (de_ele==NULL) {
      rError("rDenseMatrix:: memory exhausted");
    }
    length = nRow*nCol;
    catlas_dset(length,DZERO,de_ele,IONE);
    for (index = 0; index<other.NonZeroCount; ++index) {
      int i = other.row_index[index];
      int j = other.column_index[index];
      double value = other.sp_ele[index];
      de_ele[i+nCol*j] = de_ele[j+nCol*i] = value;
    }
    break;
  case rSparseMatrix::DENSE:
    De_Di = DENSE;
    if (de_ele && (other.nRow!=nRow || other.nCol!=nCol)) {
      delete[] de_ele;
      de_ele = NULL;
    }
    nRow = other.nRow;
    nCol = other.nCol;
    rNewCheck();
    de_ele = new double[nRow*nCol];
    if (de_ele==NULL) {
      rError("rDenseMatrix:: memory exhausted");
    }
    length = nRow*nCol;
    dcopy(&length,other.de_ele,&IONE,de_ele,&IONE);
    break;
  case rSparseMatrix::DIAGONAL:
    De_Di = DIAGONAL;
    if (di_ele && (other.nRow!=nRow || other.nCol!=nCol)) {
      delete[] di_ele;
      di_ele = NULL;
    }
    nRow = other.nRow;
    nCol = other.nCol;
    if (di_ele==NULL) {
      rNewCheck();
      di_ele = new double[nCol];
      if (di_ele==NULL) {
	rError("rDenseMatrix:: memory exhausted");
      }
    }
    dcopy(&nCol,other.di_ele,&IONE,di_ele,&IONE);
    break;
  }
  return _SUCCESS;
}

bool rDenseMatrix::copyFrom(rDenseMatrix& other)
{
  if (this == &other) {
    return _SUCCESS;
  }
  int length;
  switch(other.De_Di) {
  case DENSE:
    De_Di = DENSE;
    if (de_ele && (other.nRow!=nRow || other.nCol!=nCol)) {
      delete[] de_ele;
      de_ele = NULL;
    }
    nRow = other.nRow;
    nCol = other.nCol;
    if (de_ele==NULL) {
      rNewCheck();
      de_ele = new double[nRow*nCol];
      if (de_ele==NULL) {
	rError("rDenseMatrix:: memory exhausted");
      }
    }
    length = nRow*nCol;
    dcopy(&length,other.de_ele,&IONE,de_ele,&IONE);
    break;
  case DIAGONAL:
    De_Di = DIAGONAL;
    if (di_ele && (other.nRow!=nRow || other.nCol!=nCol)) {
      delete[] di_ele;
      di_ele = NULL;
    }	
    nRow = other.nRow;
    nCol = other.nCol;
    if (di_ele==NULL) {
      rNewCheck();
      di_ele = new double[nCol];
      if (di_ele==NULL) {
	rError("rDenseMatrix:: memory exhausted");
      }
    }
    dcopy(&nCol,other.di_ele,&IONE,di_ele,&IONE);
    break;
  }
  return _SUCCESS;
}


void rDenseMatrix::setZero()
{
  int length;
  switch(De_Di) {
  case DENSE:
    length = nRow*nCol;
    catlas_dset(length,DZERO,de_ele,IONE);
    break;
  case DIAGONAL:
    catlas_dset(nCol,DZERO,di_ele,IONE);
    break;
  }
}

void rDenseMatrix::setIdentity(double scalar)
{
  if (nRow != nCol) {
    rError("rSparseMatrix:: Identity matrix must be square matrix");
  }
  int length,step;
  switch(De_Di) {
  case DENSE:
    length = nRow*nCol;
    catlas_dset(length,DZERO,de_ele,IONE);
    step = nCol+1;
    catlas_dset(nCol,scalar,de_ele,step);
    break;
  case DIAGONAL:
    catlas_dset(nCol,scalar,di_ele,IONE);
    break;
  }
}
    
rBlockSparseMatrix::rBlockSparseMatrix()
{
  nBlock = 0;
  blockStruct = NULL;
  ele = NULL;
}

rBlockSparseMatrix::rBlockSparseMatrix(int nBlock,
				       int* blockStruct)
{
  initialize(nBlock,blockStruct);
}

rBlockSparseMatrix::~rBlockSparseMatrix()
{
  if (ele && blockStruct && nBlock>=0) {
    for (int l=0; l<nBlock; ++l) {
      ele[l].~rSparseMatrix();
    }
    delete[] ele;
    ele = NULL;
    
    delete[] blockStruct;
    blockStruct = NULL;
  }
}

void rBlockSparseMatrix::initialize(int nBlock,
				    int* blockStruct)
{
  // rMessage("rBlockSparseMatrix::initialize");

  this->nBlock = nBlock;
  if (nBlock<=0) {
    rError("rBlockSparseMatrix:: nBlock is nonpositive");
  }
  this->blockStruct = NULL;
  rNewCheck();
  this->blockStruct = new int[nBlock];
  if (this->blockStruct==NULL) {
    rError("rBlockSparseMatrix:: memory exhausted");
  }
  for (int l=0; l<nBlock; ++l) {
    this->blockStruct[l] = blockStruct[l];
  }

  ele = NULL;
  rNewCheck();
  ele = new rSparseMatrix[nBlock];
  if (ele==NULL) {
    rError("rBlockSparseMatrix:: memory exhausted");
  }
  // ATTENSION
  // You have to initialize after you count NonZeroNumber.
  
}

void rBlockSparseMatrix::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }
  fprintf(fpout,"{\n");
  if (nBlock>0 && blockStruct && ele) {
    for (int l=0; l<nBlock; ++l) {
      ele[l].display(fpout);
    }
  }
  fprintf(fpout,"} \n");
}

bool rBlockSparseMatrix::copyFrom(rBlockSparseMatrix& other)
{
  if (this == &other) {
    return _SUCCESS;
  }
  if (other.nBlock<=0) {
    rError("rBlockSparseMatix:: nBlock is nonpositive");
  }
  if (nBlock!=other.nBlock && blockStruct) {
    delete[] blockStruct;
    blockStruct = NULL;
    delete[] ele;
    ele = NULL;
  }
  nBlock = other.nBlock;
  if (blockStruct==NULL) {
    rNewCheck();
    blockStruct = new int[nBlock];
    if (blockStruct==NULL) {
      rError("rBlockSparseMatrix:: memory exhausted");
    }
    for (int l=0; l<nBlock; ++l) {
      blockStruct[l] = other.blockStruct[l];
    }
  }
  if (ele==NULL) {
    rNewCheck();
    ele = new rSparseMatrix[nBlock];
    if (ele==NULL) {
      rError("rBlockSparseMatrix:: memory exhausted");
    }
  }
  bool total_judge = _SUCCESS;
  for (int l=0; l<nBlock; ++l) {
    total_judge = ele[l].copyFrom(other.ele[l]);
  }
  if (total_judge==FAILURE) {
    rError("rBlockSparseMatrix:: copy miss");
  }
  return total_judge;
}

void rBlockSparseMatrix::setZero()
{
  if (nBlock>0 && blockStruct && ele) {
    for (int l=0; l<nBlock; ++l) {
      ele[l].setZero();
    }
  }
}

void rBlockSparseMatrix::setIdentity(double scalar)
{
  if (nBlock>0 && blockStruct && ele) {
    for (int l=0; l<nBlock; ++l) {
      ele[l].setIdentity(scalar);
    }
  }
}

void rBlockSparseMatrix::changeToDense(bool forceChange)
{
  if (nBlock>0 && blockStruct && ele) {
    for (int l=0; l<nBlock; ++l) {
      ele[l].changeToDense(forceChange);
    }
  }
}

bool rBlockSparseMatrix::sortSparseIndex(int& l, int& i, int& j)
{
  bool total_judge = _SUCCESS;
  l = -1;
  int i_in,j_in; 
  if (nBlock>0 && blockStruct && ele) {
    for (int l_in=0; l_in<nBlock; ++l_in) {
      total_judge = ele[l_in].sortSparseIndex(i_in,j_in);
      if (total_judge==FAILURE && l<0) {
	l = l_in;
	i = i_in;
	j = j_in;
      }
    }
  }
  return total_judge;
}

rBlockDenseMatrix::rBlockDenseMatrix()
{
  nBlock = 0;
  blockStruct = NULL;
  ele = NULL;
}

rBlockDenseMatrix::rBlockDenseMatrix(int nBlock,
				       int* blockStruct)
{
  initialize(nBlock,blockStruct);
}

rBlockDenseMatrix::~rBlockDenseMatrix()
{
  if (ele && blockStruct && nBlock>=0) {
    for (int l=0; l<nBlock; ++l) {
      ele[l].~rDenseMatrix();
    }
    delete[] ele;
    ele = NULL;
    
    delete[] blockStruct;
    blockStruct = NULL;
  }
}

void rBlockDenseMatrix::initialize(int nBlock,
				    int* blockStruct)
{
  // rMessage("rBlockDenseMatrix::initialize");
  if (this->blockStruct && this->nBlock!=nBlock) {
    delete[] this->blockStruct;
    this->blockStruct = NULL;
    if (ele!=NULL) {
      delete[] ele;
    }
    ele = NULL;
  }
  this->nBlock = nBlock;
  if (nBlock<=0) {
    rError("rBlockDenseMatrix:: nBlock is nonpositive");
  }
  if (this->blockStruct==NULL) {
    rNewCheck();
    this->blockStruct = new int[nBlock];
    if (this->blockStruct==NULL) {
      rError("rBlockDenseMatrix:: memory exhausted");
    }
  }
  for (int l=0; l<nBlock; ++l) {
    this->blockStruct[l] = blockStruct[l];
  }

  if (ele==NULL) {
    rNewCheck();
    ele = new rDenseMatrix[nBlock];
    if (ele==NULL) {
      rError("rBlockDenseMatrix:: memory exhausted");
    }
  }

  for (int l=0; l<nBlock; ++l) {
    int size = blockStruct[l];
    if (size>0) {
      ele[l].initialize(size,size,rDenseMatrix::DENSE);
    } else {
      size = -size;
      ele[l].initialize(size,size,rDenseMatrix::DIAGONAL);
    }
  }
}

void rBlockDenseMatrix::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }
  fprintf(fpout,"{\n");
  if (nBlock>0 && blockStruct && ele) {
    for (int l=0; l<nBlock; ++l) {
      ele[l].display(fpout);
    }
  }
  fprintf(fpout,"} \n");
}

bool rBlockDenseMatrix::copyFrom(rBlockSparseMatrix& other)
{
  if (other.nBlock<=0) {
    rError("rBlockDenseMatix:: nBlock is nonpositive");
  }
  if (nBlock!=other.nBlock && blockStruct) {
    delete[] blockStruct;
    blockStruct = NULL;
    delete[] ele;
    ele = NULL;
  }
  nBlock = other.nBlock;
  rNewCheck();
  blockStruct = new int[nBlock];
  if (blockStruct==NULL) {
    rError("rBlockDenseMatrix:: memory exhausted");
  }
  rNewCheck();
  ele = new rDenseMatrix[nBlock];
  if (ele==NULL) {
    rError("rBlockDenseMatrix:: memory exhausted");
  }
  for (int l=0; l<nBlock; ++l) {
    blockStruct[l] = other.blockStruct[l];
  }
  bool total_judge = _SUCCESS;
  for (int l=0; l<nBlock; ++l) {
    total_judge = ele[l].copyFrom(other.ele[l]);
  }
  if (total_judge==FAILURE) {
    rError("rBlockSparseMatrix:: copy miss");
  }
  return total_judge;
}

bool rBlockDenseMatrix::copyFrom(rBlockDenseMatrix& other)
{
  if (this == &other) {
    return _SUCCESS;
  }
  if (other.nBlock<=0) {
    rError("rBlockDenseMatix:: nBlock is nonpositive");
  }
  if (nBlock!=other.nBlock && blockStruct) {
    delete[] blockStruct;
    blockStruct = NULL;
    delete[] ele;
    ele = NULL;
  }
  if (blockStruct == NULL) {
    nBlock = other.nBlock;
    rNewCheck();
    blockStruct = new int[nBlock];
    if (blockStruct==NULL) {
      rError("rBlockDenseMatrix:: memory exhausted");
    }
    for (int l=0; l<nBlock; ++l) {
      blockStruct[l] = other.blockStruct[l];
    }
  }
  if (ele == NULL) {
    rNewCheck();
    ele = new rDenseMatrix[nBlock];
    if (ele==NULL) {
      rError("rBlockDenseMatrix:: memory exhausted");
    }
  }
  bool total_judge = _SUCCESS;
  for (int l=0; l<nBlock; ++l) {
    total_judge = ele[l].copyFrom(other.ele[l]);
  }
  if (total_judge==FAILURE) {
    rError("rBlockSparseMatrix:: copy miss");
  }
  return total_judge;
}

void rBlockDenseMatrix::setZero()
{
  if (nBlock>0 && blockStruct && ele) {
    for (int l=0; l<nBlock; ++l) {
      ele[l].setZero();
    }
  }
}

void rBlockDenseMatrix::setIdentity(double scalar)
{
  if (nBlock>0 && blockStruct && ele) {
    for (int l=0; l<nBlock; ++l) {
      ele[l].setIdentity(scalar);
    }
  }
}

