
#ifndef _testblasunit_h
#define _testblasunit_h

#include "ap.h"
#include "ialglib.h"

#include "blas.h"


bool testblas(bool silent);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testblasunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testblasunit_test();


#endif

