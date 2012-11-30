
#ifndef _testsafesolveunit_h
#define _testsafesolveunit_h

#include "ap.h"
#include "ialglib.h"

#include "safesolve.h"


/*************************************************************************
Main unittest subroutine
*************************************************************************/
bool testsafesolve(bool silent);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testsafesolveunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testsafesolveunit_test();


#endif

