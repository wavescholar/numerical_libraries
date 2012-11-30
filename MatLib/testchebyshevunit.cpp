
 
#include <stdio.h>
#include "testchebyshevunit.h"

bool testchebyshev(bool silent)
{
    bool result;
    double err;
    double sumerr;
    double cerr;
    double ferr;
    double threshold;
    double x;
    double v;
    double t;
    int pass;
    int i;
    int j;
    int k;
    int n;
    int maxn;
    ap::real_1d_array c;
    ap::real_1d_array p1;
    ap::real_1d_array p2;
    ap::real_2d_array a;
    bool waserrors;

    err = 0;
    sumerr = 0;
    cerr = 0;
    ferr = 0;
    threshold = 1.0E-9;
    waserrors = false;
    
    //
    // Testing Chebyshev polynomials of the first kind
    //
    err = ap::maxreal(err, fabs(chebyshevcalculate(1, 0, 0.00)-1));
    err = ap::maxreal(err, fabs(chebyshevcalculate(1, 0, 0.33)-1));
    err = ap::maxreal(err, fabs(chebyshevcalculate(1, 0, -0.42)-1));
    x = 0.2;
    err = ap::maxreal(err, fabs(chebyshevcalculate(1, 1, x)-0.2));
    x = 0.4;
    err = ap::maxreal(err, fabs(chebyshevcalculate(1, 1, x)-0.4));
    x = 0.6;
    err = ap::maxreal(err, fabs(chebyshevcalculate(1, 1, x)-0.6));
    x = 0.8;
    err = ap::maxreal(err, fabs(chebyshevcalculate(1, 1, x)-0.8));
    x = 1.0;
    err = ap::maxreal(err, fabs(chebyshevcalculate(1, 1, x)-1.0));
    x = 0.2;
    err = ap::maxreal(err, fabs(chebyshevcalculate(1, 2, x)+0.92));
    x = 0.4;
    err = ap::maxreal(err, fabs(chebyshevcalculate(1, 2, x)+0.68));
    x = 0.6;
    err = ap::maxreal(err, fabs(chebyshevcalculate(1, 2, x)+0.28));
    x = 0.8;
    err = ap::maxreal(err, fabs(chebyshevcalculate(1, 2, x)-0.28));
    x = 1.0;
    err = ap::maxreal(err, fabs(chebyshevcalculate(1, 2, x)-1.00));
    n = 10;
    err = ap::maxreal(err, fabs(chebyshevcalculate(1, n, 0.2)-0.4284556288));
    n = 11;
    err = ap::maxreal(err, fabs(chebyshevcalculate(1, n, 0.2)+0.7996160205));
    n = 12;
    err = ap::maxreal(err, fabs(chebyshevcalculate(1, n, 0.2)+0.7483020370));
    
    //
    // Testing Chebyshev polynomials of the second kind
    //
    n = 0;
    err = ap::maxreal(err, fabs(chebyshevcalculate(2, n, 0.2)-1.0000000000));
    n = 1;
    err = ap::maxreal(err, fabs(chebyshevcalculate(2, n, 0.2)-0.4000000000));
    n = 2;
    err = ap::maxreal(err, fabs(chebyshevcalculate(2, n, 0.2)+0.8400000000));
    n = 3;
    err = ap::maxreal(err, fabs(chebyshevcalculate(2, n, 0.2)+0.7360000000));
    n = 4;
    err = ap::maxreal(err, fabs(chebyshevcalculate(2, n, 0.2)-0.5456000000));
    n = 10;
    err = ap::maxreal(err, fabs(chebyshevcalculate(2, n, 0.2)-0.6128946176));
    n = 11;
    err = ap::maxreal(err, fabs(chebyshevcalculate(2, n, 0.2)+0.6770370970));
    n = 12;
    err = ap::maxreal(err, fabs(chebyshevcalculate(2, n, 0.2)+0.8837094564));
    
    //
    // Testing Clenshaw summation
    //
    maxn = 20;
    c.setbounds(0, maxn);
    for(k = 1; k <= 2; k++)
    {
        for(pass = 1; pass <= 10; pass++)
        {
            x = 2*ap::randomreal()-1;
            v = 0;
            for(n = 0; n <= maxn; n++)
            {
                c(n) = 2*ap::randomreal()-1;
                v = v+chebyshevcalculate(k, n, x)*c(n);
                sumerr = ap::maxreal(sumerr, fabs(v-chebyshevsum(c, k, n, x)));
            }
        }
    }
    
    //
    // Testing coefficients
    //
    chebyshevcoefficients(0, c);
    cerr = ap::maxreal(cerr, fabs(c(0)-1));
    chebyshevcoefficients(1, c);
    cerr = ap::maxreal(cerr, fabs(c(0)-0));
    cerr = ap::maxreal(cerr, fabs(c(1)-1));
    chebyshevcoefficients(2, c);
    cerr = ap::maxreal(cerr, fabs(c(0)+1));
    cerr = ap::maxreal(cerr, fabs(c(1)-0));
    cerr = ap::maxreal(cerr, fabs(c(2)-2));
    chebyshevcoefficients(3, c);
    cerr = ap::maxreal(cerr, fabs(c(0)-0));
    cerr = ap::maxreal(cerr, fabs(c(1)+3));
    cerr = ap::maxreal(cerr, fabs(c(2)-0));
    cerr = ap::maxreal(cerr, fabs(c(3)-4));
    chebyshevcoefficients(4, c);
    cerr = ap::maxreal(cerr, fabs(c(0)-1));
    cerr = ap::maxreal(cerr, fabs(c(1)-0));
    cerr = ap::maxreal(cerr, fabs(c(2)+8));
    cerr = ap::maxreal(cerr, fabs(c(3)-0));
    cerr = ap::maxreal(cerr, fabs(c(4)-8));
    chebyshevcoefficients(9, c);
    cerr = ap::maxreal(cerr, fabs(c(0)-0));
    cerr = ap::maxreal(cerr, fabs(c(1)-9));
    cerr = ap::maxreal(cerr, fabs(c(2)-0));
    cerr = ap::maxreal(cerr, fabs(c(3)+120));
    cerr = ap::maxreal(cerr, fabs(c(4)-0));
    cerr = ap::maxreal(cerr, fabs(c(5)-432));
    cerr = ap::maxreal(cerr, fabs(c(6)-0));
    cerr = ap::maxreal(cerr, fabs(c(7)+576));
    cerr = ap::maxreal(cerr, fabs(c(8)-0));
    cerr = ap::maxreal(cerr, fabs(c(9)-256));
    
    //
    // Testing FromChebyshev
    //
    maxn = 10;
    a.setbounds(0, maxn, 0, maxn);
    for(i = 0; i <= maxn; i++)
    {
        for(j = 0; j <= maxn; j++)
        {
            a(i,j) = 0;
        }
        chebyshevcoefficients(i, c);
        ap::vmove(&a(i, 0), 1, &c(0), 1, ap::vlen(0,i));
    }
    c.setbounds(0, maxn);
    p1.setbounds(0, maxn);
    for(n = 0; n <= maxn; n++)
    {
        for(pass = 1; pass <= 10; pass++)
        {
            for(i = 0; i <= n; i++)
            {
                p1(i) = 0;
            }
            for(i = 0; i <= n; i++)
            {
                c(i) = 2*ap::randomreal()-1;
                v = c(i);
                ap::vadd(&p1(0), 1, &a(i, 0), 1, ap::vlen(0,i), v);
            }
            fromchebyshev(c, n, p2);
            for(i = 0; i <= n; i++)
            {
                ferr = ap::maxreal(ferr, fabs(p1(i)-p2(i)));
            }
        }
    }
    
    //
    // Reporting
    //
    waserrors = ap::fp_greater(err,threshold)||ap::fp_greater(sumerr,threshold)||ap::fp_greater(cerr,threshold)||ap::fp_greater(ferr,threshold);
    if( !silent )
    {
        printf("TESTING CALCULATION OF THE CHEBYSHEV POLYNOMIALS\n");
        printf("Max error against table                   %5.2le\n",
            double(err));
        printf("Summation error                           %5.2le\n",
            double(sumerr));
        printf("Coefficients error                        %5.2le\n",
            double(cerr));
        printf("FrobChebyshev error                       %5.2le\n",
            double(ferr));
        printf("Threshold                                 %5.2le\n",
            double(threshold));
        if( !waserrors )
        {
            printf("TEST PASSED\n");
        }
        else
        {
            printf("TEST FAILED\n");
        }
    }
    result = !waserrors;
    return result;
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testchebyshevunit_test_silent()
{
    bool result;

    result = testchebyshev(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testchebyshevunit_test()
{
    bool result;

    result = testchebyshev(false);
    return result;
}




