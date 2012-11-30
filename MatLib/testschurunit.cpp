
 
#include <stdio.h>
#include "testschurunit.h"

static void fillsparsea(ap::real_2d_array& a, int n, double sparcity);
static void testschurproblem(const ap::real_2d_array& a,
     int n,
     double& materr,
     double& orterr,
     bool& errstruct,
     bool& wfailed);

/*************************************************************************
Testing Schur decomposition subroutine
*************************************************************************/
bool testschur(bool silent)
{
    bool result;
    ap::real_2d_array a;
    int n;
    int maxn;
    int i;
    int j;
    int pass;
    int passcount;
    bool waserrors;
    bool errstruct;
    bool wfailed;
    double materr;
    double orterr;
    double threshold;

    materr = 0;
    orterr = 0;
    errstruct = false;
    wfailed = false;
    waserrors = false;
    maxn = 70;
    passcount = 1;
    threshold = 5*100*ap::machineepsilon;
    a.setbounds(0, maxn-1, 0, maxn-1);
    
    //
    // zero matrix, several cases
    //
    for(i = 0; i <= maxn-1; i++)
    {
        for(j = 0; j <= maxn-1; j++)
        {
            a(i,j) = 0;
        }
    }
    for(n = 1; n <= maxn; n++)
    {
        if( n>30&&n%2==0 )
        {
            continue;
        }
        testschurproblem(a, n, materr, orterr, errstruct, wfailed);
    }
    
    //
    // Dense matrix
    //
    for(pass = 1; pass <= passcount; pass++)
    {
        for(n = 1; n <= maxn; n++)
        {
            if( n>30&&n%2==0 )
            {
                continue;
            }
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    a(i,j) = 2*ap::randomreal()-1;
                }
            }
            testschurproblem(a, n, materr, orterr, errstruct, wfailed);
        }
    }
    
    //
    // Sparse matrices, very sparse matrices, incredible sparse matrices
    //
    for(pass = 1; pass <= 1; pass++)
    {
        for(n = 1; n <= maxn; n++)
        {
            if( n>30&&n%3!=0 )
            {
                continue;
            }
            fillsparsea(a, n, 0.8);
            testschurproblem(a, n, materr, orterr, errstruct, wfailed);
            fillsparsea(a, n, 0.9);
            testschurproblem(a, n, materr, orterr, errstruct, wfailed);
            fillsparsea(a, n, 0.95);
            testschurproblem(a, n, materr, orterr, errstruct, wfailed);
            fillsparsea(a, n, 0.997);
            testschurproblem(a, n, materr, orterr, errstruct, wfailed);
        }
    }
    
    //
    // report
    //
    waserrors = ap::fp_greater(materr,threshold)||ap::fp_greater(orterr,threshold)||errstruct||wfailed;
    if( !silent )
    {
        printf("TESTING SCHUR DECOMPOSITION\n");
        printf("Schur decomposition error:               %5.3le\n",
            double(materr));
        printf("Schur orthogonality error:               %5.3le\n",
            double(orterr));
        printf("T matrix structure:                      ");
        if( !errstruct )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("Always converged:                        ");
        if( !wfailed )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
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


static void fillsparsea(ap::real_2d_array& a, int n, double sparcity)
{
    int i;
    int j;

    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            if( ap::fp_greater_eq(ap::randomreal(),sparcity) )
            {
                a(i,j) = 2*ap::randomreal()-1;
            }
            else
            {
                a(i,j) = 0;
            }
        }
    }
}


static void testschurproblem(const ap::real_2d_array& a,
     int n,
     double& materr,
     double& orterr,
     bool& errstruct,
     bool& wfailed)
{
    ap::real_2d_array s;
    ap::real_2d_array t;
    ap::real_1d_array sr;
    ap::real_1d_array astc;
    ap::real_1d_array sastc;
    int i;
    int j;
    int k;
    double v;
    double locerr;

    sr.setbounds(0, n-1);
    astc.setbounds(0, n-1);
    sastc.setbounds(0, n-1);
    
    //
    // Schur decomposition, convergence test
    //
    t.setbounds(0, n-1, 0, n-1);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            t(i,j) = a(i,j);
        }
    }
    if( !rmatrixschur(t, n, s) )
    {
        wfailed = true;
        return;
    }
    
    //
    // decomposition error
    //
    locerr = 0;
    for(j = 0; j <= n-1; j++)
    {
        ap::vmove(&sr(0), 1, &s(j, 0), 1, ap::vlen(0,n-1));
        for(k = 0; k <= n-1; k++)
        {
            v = ap::vdotproduct(&t(k, 0), 1, &sr(0), 1, ap::vlen(0,n-1));
            astc(k) = v;
        }
        for(k = 0; k <= n-1; k++)
        {
            v = ap::vdotproduct(&s(k, 0), 1, &astc(0), 1, ap::vlen(0,n-1));
            sastc(k) = v;
        }
        for(k = 0; k <= n-1; k++)
        {
            locerr = ap::maxreal(locerr, fabs(sastc(k)-a(k,j)));
        }
    }
    materr = ap::maxreal(materr, locerr);
    
    //
    // orthogonality error
    //
    locerr = 0;
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&s(0, i), s.getstride(), &s(0, j), s.getstride(), ap::vlen(0,n-1));
            if( i!=j )
            {
                locerr = ap::maxreal(locerr, fabs(v));
            }
            else
            {
                locerr = ap::maxreal(locerr, fabs(v-1));
            }
        }
    }
    orterr = ap::maxreal(orterr, locerr);
    
    //
    // T matrix structure
    //
    for(j = 0; j <= n-1; j++)
    {
        for(i = j+2; i <= n-1; i++)
        {
            if( ap::fp_neq(t(i,j),0) )
            {
                errstruct = true;
            }
        }
    }
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testschurunit_test_silent()
{
    bool result;

    result = testschur(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testschurunit_test()
{
    bool result;

    result = testschur(false);
    return result;
}




