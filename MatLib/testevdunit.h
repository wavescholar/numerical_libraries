
#ifndef _testevdunit_h
#define _testevdunit_h

#include "ap.h"
#include "ialglib.h"

#include "hblas.h"
#include "reflections.h"
#include "creflections.h"
#include "sblas.h"
#include "ablasf.h"
#include "ablas.h"
#include "ortfac.h"
#include "blas.h"
#include "rotations.h"
#include "hsschur.h"
#include "evd.h"


/*************************************************************************
Testing symmetric EVD subroutine
*************************************************************************/
bool testevd(bool silent);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testevdunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testevdunit_test();


#endif

