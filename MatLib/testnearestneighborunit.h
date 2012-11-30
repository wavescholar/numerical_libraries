
#ifndef _testnearestneighborunit_h
#define _testnearestneighborunit_h

#include "ap.h"
#include "ialglib.h"

#include "tsort.h"
#include "nearestneighbor.h"


/*************************************************************************
Testing Nearest Neighbor Search
*************************************************************************/
bool testnearestneighbor(bool silent);


/*************************************************************************
Silent unit test
*************************************************************************/
bool testnearestneighborunit_test_silent();


/*************************************************************************
Unit test
*************************************************************************/
bool testnearestneighborunit_test();


#endif

