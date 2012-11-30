
#ifndef _testsblasunit_h
#define _testsblasunit_h

#include "ap.h"
#include "ialglib.h"

#include "sblas.h"


bool testsblas(bool silent);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testsblasunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testsblasunit_test();


#endif

