
#ifndef _testminlmunit_h
#define _testminlmunit_h

#include "ap.h"
#include "ialglib.h"

#include "blas.h"
#include "reflections.h"
#include "creflections.h"
#include "hqrnd.h"
#include "matgen.h"
#include "ablasf.h"
#include "ablas.h"
#include "trfac.h"
#include "trlinsolve.h"
#include "safesolve.h"
#include "rcond.h"
#include "matinv.h"
#include "hblas.h"
#include "sblas.h"
#include "ortfac.h"
#include "rotations.h"
#include "bdsvd.h"
#include "svd.h"
#include "xblas.h"
#include "densesolver.h"
#include "linmin.h"
#include "minlbfgs.h"
#include "minlm.h"


bool testminlm(bool silent);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testminlmunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testminlmunit_test();


#endif

