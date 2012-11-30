
#ifndef _testregressunit_h
#define _testregressunit_h

#include "ap.h"
#include "ialglib.h"

#include "descriptivestatistics.h"
#include "gammafunc.h"
#include "normaldistr.h"
#include "igammaf.h"
#include "hblas.h"
#include "reflections.h"
#include "creflections.h"
#include "sblas.h"
#include "ablasf.h"
#include "ablas.h"
#include "ortfac.h"
#include "blas.h"
#include "rotations.h"
#include "bdsvd.h"
#include "svd.h"
#include "linreg.h"


bool testlinregression(bool silent);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testregressunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testregressunit_test();


#endif

