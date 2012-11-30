
 
#include <stdio.h>
#include "testsstunit.h"

static void makeacopy(const ap::real_2d_array& a,
     int m,
     int n,
     ap::real_2d_array& b);

/*************************************************************************
Main unittest subroutine
*************************************************************************/
bool testsst(bool silent)
{
    bool result;
    int maxmn;
    int passcount;
    double threshold;
    ap::real_2d_array aeffective;
    ap::real_2d_array aparam;
    ap::real_1d_array xe;
    ap::real_1d_array b;
    int n;
    int pass;
    int i;
    int j;
    int cnts;
    int cntu;
    int cntt;
    int cntm;
    bool waserrors;
    bool isupper;
    bool istrans;
    bool isunit;
    double v;
    double s;

    waserrors = false;
    maxmn = 15;
    passcount = 15;
    threshold = 1000*ap::machineepsilon;
    
    //
    // Different problems
    //
    for(n = 1; n <= maxmn; n++)
    {
        aeffective.setbounds(0, n-1, 0, n-1);
        aparam.setbounds(0, n-1, 0, n-1);
        xe.setbounds(0, n-1);
        b.setbounds(0, n-1);
        for(pass = 1; pass <= passcount; pass++)
        {
            for(cnts = 0; cnts <= 1; cnts++)
            {
                for(cntu = 0; cntu <= 1; cntu++)
                {
                    for(cntt = 0; cntt <= 1; cntt++)
                    {
                        for(cntm = 0; cntm <= 2; cntm++)
                        {
                            isupper = cnts==0;
                            isunit = cntu==0;
                            istrans = cntt==0;
                            
                            //
                            // Skip meaningless combinations of parameters:
                            // (matrix is singular) AND (matrix is unit diagonal)
                            //
                            if( cntm==2&&isunit )
                            {
                                continue;
                            }
                            
                            //
                            // Clear matrices
                            //
                            for(i = 0; i <= n-1; i++)
                            {
                                for(j = 0; j <= n-1; j++)
                                {
                                    aeffective(i,j) = 0;
                                    aparam(i,j) = 0;
                                }
                            }
                            
                            //
                            // Prepare matrices
                            //
                            if( isupper )
                            {
                                for(i = 0; i <= n-1; i++)
                                {
                                    for(j = i; j <= n-1; j++)
                                    {
                                        aeffective(i,j) = 0.9*(2*ap::randomreal()-1);
                                        aparam(i,j) = aeffective(i,j);
                                    }
                                    aeffective(i,i) = (2*ap::randominteger(2)-1)*(0.8+ap::randomreal());
                                    aparam(i,i) = aeffective(i,i);
                                }
                            }
                            else
                            {
                                for(i = 0; i <= n-1; i++)
                                {
                                    for(j = 0; j <= i; j++)
                                    {
                                        aeffective(i,j) = 0.9*(2*ap::randomreal()-1);
                                        aparam(i,j) = aeffective(i,j);
                                    }
                                    aeffective(i,i) = (2*ap::randominteger(2)-1)*(0.8+ap::randomreal());
                                    aparam(i,i) = aeffective(i,i);
                                }
                            }
                            if( isunit )
                            {
                                for(i = 0; i <= n-1; i++)
                                {
                                    aeffective(i,i) = 1;
                                    aparam(i,i) = 0;
                                }
                            }
                            if( istrans )
                            {
                                if( isupper )
                                {
                                    for(i = 0; i <= n-1; i++)
                                    {
                                        for(j = i+1; j <= n-1; j++)
                                        {
                                            aeffective(j,i) = aeffective(i,j);
                                            aeffective(i,j) = 0;
                                        }
                                    }
                                }
                                else
                                {
                                    for(i = 0; i <= n-1; i++)
                                    {
                                        for(j = i+1; j <= n-1; j++)
                                        {
                                            aeffective(i,j) = aeffective(j,i);
                                            aeffective(j,i) = 0;
                                        }
                                    }
                                }
                            }
                            
                            //
                            // Prepare task, solve, compare
                            //
                            for(i = 0; i <= n-1; i++)
                            {
                                xe(i) = 2*ap::randomreal()-1;
                            }
                            for(i = 0; i <= n-1; i++)
                            {
                                v = ap::vdotproduct(&aeffective(i, 0), 1, &xe(0), 1, ap::vlen(0,n-1));
                                b(i) = v;
                            }
                            rmatrixtrsafesolve(aparam, n, b, s, isupper, istrans, isunit);
                            ap::vmul(&xe(0), 1, ap::vlen(0,n-1), s);
                            ap::vsub(&xe(0), 1, &b(0), 1, ap::vlen(0,n-1));
                            v = ap::vdotproduct(&xe(0), 1, &xe(0), 1, ap::vlen(0,n-1));
                            v = sqrt(v);
                            waserrors = waserrors||ap::fp_greater(v,threshold);
                        }
                    }
                }
            }
        }
    }
    
    //
    // report
    //
    if( !silent )
    {
        printf("TESTING RMatrixTRSafeSolve\n");
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
Copy
*************************************************************************/
static void makeacopy(const ap::real_2d_array& a,
     int m,
     int n,
     ap::real_2d_array& b)
{
    int i;
    int j;

    b.setbounds(0, m-1, 0, n-1);
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            b(i,j) = a(i,j);
        }
    }
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testsstunit_test_silent()
{
    bool result;

    result = testsst(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testsstunit_test()
{
    bool result;

    result = testsst(false);
    return result;
}




