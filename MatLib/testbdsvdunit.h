
#ifndef _testbdsvdunit_h
#define _testbdsvdunit_h

#include "ap.h"
#include "ialglib.h"

#include "rotations.h"
#include "bdsvd.h"


/*************************************************************************
Testing bidiagonal SVD decomposition subroutine
*************************************************************************/
bool testbdsvd(bool silent);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testbdsvdunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testbdsvdunit_test();


#endif

