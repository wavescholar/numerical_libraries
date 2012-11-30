
#ifndef _testxblasunit_h
#define _testxblasunit_h

#include "ap.h"
#include "ialglib.h"

#include "xblas.h"


bool testxblas(bool silent);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testxblasunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testxblasunit_test();


#endif

