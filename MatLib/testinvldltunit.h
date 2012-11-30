
#ifndef _testinvldltunit_h
#define _testinvldltunit_h

#include "ap.h"
#include "ialglib.h"

#include "sblas.h"
#include "ldlt.h"
#include "sinverse.h"


bool testinvldlt(bool silent);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testinvldltunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testinvldltunit_test();


#endif

