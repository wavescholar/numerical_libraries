
 
#include <stdio.h>
#include "testsblasunit.h"

bool testsblas(bool silent)
{
    bool result;
    ap::real_2d_array a;
    ap::real_2d_array ua;
    ap::real_2d_array la;
    ap::real_1d_array x;
    ap::real_1d_array y1;
    ap::real_1d_array y2;
    ap::real_1d_array y3;
    int n;
    int maxn;
    int i;
    int j;
    int i1;
    int i2;
    int gpass;
    bool waserrors;
    double mverr;
    double threshold;
    double alpha;
    double v;

    mverr = 0;
    waserrors = false;
    maxn = 10;
    threshold = 1000*ap::machineepsilon;
    
    //
    // Test MV
    //
    for(n = 2; n <= maxn; n++)
    {
        a.setbounds(1, n, 1, n);
        ua.setbounds(1, n, 1, n);
        la.setbounds(1, n, 1, n);
        x.setbounds(1, n);
        y1.setbounds(1, n);
        y2.setbounds(1, n);
        y3.setbounds(1, n);
        
        //
        // fill A, UA, LA
        //
        for(i = 1; i <= n; i++)
        {
            a(i,i) = 2*ap::randomreal()-1;
            for(j = i+1; j <= n; j++)
            {
                a(i,j) = 2*ap::randomreal()-1;
                a(j,i) = a(i,j);
            }
        }
        for(i = 1; i <= n; i++)
        {
            for(j = 1; j <= n; j++)
            {
                ua(i,j) = 0;
            }
        }
        for(i = 1; i <= n; i++)
        {
            for(j = i; j <= n; j++)
            {
                ua(i,j) = a(i,j);
            }
        }
        for(i = 1; i <= n; i++)
        {
            for(j = 1; j <= n; j++)
            {
                la(i,j) = 0;
            }
        }
        for(i = 1; i <= n; i++)
        {
            for(j = 1; j <= i; j++)
            {
                la(i,j) = a(i,j);
            }
        }
        
        //
        // test on different I1, I2
        //
        for(i1 = 1; i1 <= n; i1++)
        {
            for(i2 = i1; i2 <= n; i2++)
            {
                
                //
                // Fill X, choose Alpha
                //
                for(i = 1; i <= i2-i1+1; i++)
                {
                    x(i) = 2*ap::randomreal()-1;
                }
                alpha = 2*ap::randomreal()-1;
                
                //
                // calculate A*x, UA*x, LA*x
                //
                for(i = i1; i <= i2; i++)
                {
                    v = ap::vdotproduct(&a(i, i1), 1, &x(1), 1, ap::vlen(i1,i2));
                    y1(i-i1+1) = alpha*v;
                }
                symmetricmatrixvectormultiply(ua, true, i1, i2, x, alpha, y2);
                symmetricmatrixvectormultiply(la, false, i1, i2, x, alpha, y3);
                
                //
                // Calculate error
                //
                ap::vsub(&y2(1), 1, &y1(1), 1, ap::vlen(1,i2-i1+1));
                v = ap::vdotproduct(&y2(1), 1, &y2(1), 1, ap::vlen(1,i2-i1+1));
                mverr = ap::maxreal(mverr, sqrt(v));
                ap::vsub(&y3(1), 1, &y1(1), 1, ap::vlen(1,i2-i1+1));
                v = ap::vdotproduct(&y3(1), 1, &y3(1), 1, ap::vlen(1,i2-i1+1));
                mverr = ap::maxreal(mverr, sqrt(v));
            }
        }
    }
    
    //
    // report
    //
    waserrors = ap::fp_greater(mverr,threshold);
    if( !silent )
    {
        printf("TESTING SYMMETRIC BLAS\n");
        printf("MV error:                                %5.3le\n",
            double(mverr));
        printf("Threshold:                               %5.3le\n",
            double(threshold));
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
bool testsblasunit_test_silent()
{
    bool result;

    result = testsblas(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testsblasunit_test()
{
    bool result;

    result = testsblas(false);
    return result;
}




