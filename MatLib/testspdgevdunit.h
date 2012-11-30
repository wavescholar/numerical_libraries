
#ifndef _testspdgevdunit_h
#define _testspdgevdunit_h

#include "ap.h"
#include "ialglib.h"

#include "reflections.h"
#include "creflections.h"
#include "hqrnd.h"
#include "matgen.h"
#include "ablasf.h"
#include "ablas.h"
#include "trfac.h"
#include "sblas.h"
#include "blas.h"
#include "trlinsolve.h"
#include "safesolve.h"
#include "rcond.h"
#include "matinv.h"
#include "hblas.h"
#include "ortfac.h"
#include "rotations.h"
#include "hsschur.h"
#include "evd.h"
#include "spdgevd.h"


/*************************************************************************
Testing bidiagonal SVD decomposition subroutine
*************************************************************************/
bool testspdgevd(bool silent);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testspdgevdunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testspdgevdunit_test();


#endif

