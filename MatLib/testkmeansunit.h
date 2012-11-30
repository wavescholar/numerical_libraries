
#ifndef _testkmeansunit_h
#define _testkmeansunit_h

#include "ap.h"
#include "ialglib.h"

#include "blas.h"
#include "kmeans.h"


bool testkmeans(bool silent);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testkmeansunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testkmeansunit_test();


#endif

