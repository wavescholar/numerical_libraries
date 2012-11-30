
 
#include <stdio.h>
#include "testlinminunit.h"

bool testlinmin(bool silent)
{
    bool result;
    bool waserrors;

    waserrors = false;
    if( !silent )
    {
        printf("TESTING LINMIN\n");
        if( waserrors )
        {
            printf("TEST FAILED\n");
        }
        else
        {
            printf("TEST PASSED\n");
        }
        printf("\n\n");
    }
    result = !waserrors;
    return result;
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testlinminunit_test_silent()
{
    bool result;

    result = testlinmin(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testlinminunit_test()
{
    bool result;

    result = testlinmin(false);
    return result;
}




