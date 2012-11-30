
#ifndef _testhblasunit_h
#define _testhblasunit_h

#include "ap.h"
#include "ialglib.h"

#include "hblas.h"


bool testhblas(bool silent);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testhblasunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testhblasunit_test();


#endif

