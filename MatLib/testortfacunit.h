
#ifndef _testortfacunit_h
#define _testortfacunit_h

#include "ap.h"
#include "ialglib.h"

#include "hblas.h"
#include "reflections.h"
#include "creflections.h"
#include "sblas.h"
#include "ablasf.h"
#include "ablas.h"
#include "ortfac.h"


/*************************************************************************
Main unittest subroutine
*************************************************************************/
bool testortfac(bool silent);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testortfacunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testortfacunit_test();


#endif

