
 
#include <stdio.h>
#include "testbdsvdunit.h"

static int failcount;
static int succcount;

static void fillidentity(ap::real_2d_array& a, int n);
static void fillsparsede(ap::real_1d_array& d,
     ap::real_1d_array& e,
     int n,
     double sparcity);
static void getbdsvderror(const ap::real_1d_array& d,
     const ap::real_1d_array& e,
     int n,
     bool isupper,
     const ap::real_2d_array& u,
     const ap::real_2d_array& c,
     const ap::real_1d_array& w,
     const ap::real_2d_array& vt,
     double& materr,
     double& orterr,
     bool& wsorted);
static void testbdsvdproblem(const ap::real_1d_array& d,
     const ap::real_1d_array& e,
     int n,
     double& materr,
     double& orterr,
     bool& wsorted,
     bool& wfailed);

/*************************************************************************
Testing bidiagonal SVD decomposition subroutine
*************************************************************************/
bool testbdsvd(bool silent)
{
    bool result;
    ap::real_1d_array d;
    ap::real_1d_array e;
    ap::real_2d_array mempty;
    int n;
    int maxn;
    int i;
    int j;
    int gpass;
    int pass;
    bool waserrors;
    bool wsorted;
    bool wfailed;
    bool failcase;
    double materr;
    double orterr;
    double threshold;
    double failthreshold;
    double failr;

    failcount = 0;
    succcount = 0;
    materr = 0;
    orterr = 0;
    wsorted = true;
    wfailed = false;
    waserrors = false;
    maxn = 15;
    threshold = 5*100*ap::machineepsilon;
    failthreshold = 1.0E-2;
    d.setbounds(0, maxn-1);
    e.setbounds(0, maxn-2);
    
    //
    // special case: fail matrix
    //
    n = 5;
    d(0) = -8.27448347422711894000e-01;
    d(1) = -8.16705832087160854600e-01;
    d(2) = -2.53974358904729382800e-17;
    d(3) = -1.24626684881972815700e+00;
    d(4) = -4.64744131545637651000e-01;
    e(0) = -3.25785088656270038800e-01;
    e(1) = -1.03732413708914436580e-01;
    e(2) = -9.57365642262031357700e-02;
    e(3) = -2.71564153973817390400e-01;
    failcase = rmatrixbdsvd(d, e, n, true, false, mempty, 0, mempty, 0, mempty, 0);
    
    //
    // special case: zero divide matrix
    // unfixed LAPACK routine should fail on this problem
    //
    n = 7;
    d(0) = -6.96462904751731892700e-01;
    d(1) = 0.00000000000000000000e+00;
    d(2) = -5.73827770385971991400e-01;
    d(3) = -6.62562624399371191700e-01;
    d(4) = 5.82737148001782223600e-01;
    d(5) = 3.84825263580925003300e-01;
    d(6) = 9.84087420830525472200e-01;
    e(0) = -7.30307931760612871800e-02;
    e(1) = -2.30079042939542843800e-01;
    e(2) = -6.87824621739351216300e-01;
    e(3) = -1.77306437707837570600e-02;
    e(4) = 1.78285126526551632000e-15;
    e(5) = -4.89434737751289969400e-02;
    rmatrixbdsvd(d, e, n, true, false, mempty, 0, mempty, 0, mempty, 0);
    
    //
    // zero matrix, several cases
    //
    for(i = 0; i <= maxn-1; i++)
    {
        d(i) = 0;
    }
    for(i = 0; i <= maxn-2; i++)
    {
        e(i) = 0;
    }
    for(n = 1; n <= maxn; n++)
    {
        testbdsvdproblem(d, e, n, materr, orterr, wsorted, wfailed);
    }
    
    //
    // Dense matrix
    //
    for(n = 1; n <= maxn; n++)
    {
        for(pass = 1; pass <= 10; pass++)
        {
            for(i = 0; i <= maxn-1; i++)
            {
                d(i) = 2*ap::randomreal()-1;
            }
            for(i = 0; i <= maxn-2; i++)
            {
                e(i) = 2*ap::randomreal()-1;
            }
            testbdsvdproblem(d, e, n, materr, orterr, wsorted, wfailed);
        }
    }
    
    //
    // Sparse matrices, very sparse matrices, incredible sparse matrices
    //
    for(n = 1; n <= maxn; n++)
    {
        for(pass = 1; pass <= 10; pass++)
        {
            fillsparsede(d, e, n, 0.5);
            testbdsvdproblem(d, e, n, materr, orterr, wsorted, wfailed);
            fillsparsede(d, e, n, 0.8);
            testbdsvdproblem(d, e, n, materr, orterr, wsorted, wfailed);
            fillsparsede(d, e, n, 0.9);
            testbdsvdproblem(d, e, n, materr, orterr, wsorted, wfailed);
            fillsparsede(d, e, n, 0.95);
            testbdsvdproblem(d, e, n, materr, orterr, wsorted, wfailed);
        }
    }
    
    //
    // report
    //
    failr = double(failcount)/double(succcount+failcount);
    waserrors = ap::fp_greater(materr,threshold)||ap::fp_greater(orterr,threshold)||!wsorted||ap::fp_greater(failr,failthreshold);
    if( !silent )
    {
        printf("TESTING BIDIAGONAL SVD DECOMPOSITION\n");
        printf("SVD decomposition error:                 %5.3le\n",
            double(materr));
        printf("SVD orthogonality error:                 %5.3le\n",
            double(orterr));
        printf("Singular values order:                   ");
        if( wsorted )
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
            printf("YES\n");
        }
        else
        {
            printf("NO\n");
            printf("Fail ratio:                              %5.3lf\n",
                double(failr));
        }
        printf("Fail matrix test:                        ");
        if( !failcase )
        {
            printf("AS EXPECTED\n");
        }
        else
        {
            printf("CONVERGED (UNEXPECTED)\n");
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


static void fillidentity(ap::real_2d_array& a, int n)
{
    int i;
    int j;

    a.setbounds(0, n-1, 0, n-1);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            if( i==j )
            {
                a(i,j) = 1;
            }
            else
            {
                a(i,j) = 0;
            }
        }
    }
}


static void fillsparsede(ap::real_1d_array& d,
     ap::real_1d_array& e,
     int n,
     double sparcity)
{
    int i;
    int j;

    d.setbounds(0, n-1);
    e.setbounds(0, ap::maxint(0, n-2));
    for(i = 0; i <= n-1; i++)
    {
        if( ap::fp_greater_eq(ap::randomreal(),sparcity) )
        {
            d(i) = 2*ap::randomreal()-1;
        }
        else
        {
            d(i) = 0;
        }
    }
    for(i = 0; i <= n-2; i++)
    {
        if( ap::fp_greater_eq(ap::randomreal(),sparcity) )
        {
            e(i) = 2*ap::randomreal()-1;
        }
        else
        {
            e(i) = 0;
        }
    }
}


static void getbdsvderror(const ap::real_1d_array& d,
     const ap::real_1d_array& e,
     int n,
     bool isupper,
     const ap::real_2d_array& u,
     const ap::real_2d_array& c,
     const ap::real_1d_array& w,
     const ap::real_2d_array& vt,
     double& materr,
     double& orterr,
     bool& wsorted)
{
    int i;
    int j;
    int k;
    double locerr;
    double sm;

    
    //
    // decomposition error
    //
    locerr = 0;
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            sm = 0;
            for(k = 0; k <= n-1; k++)
            {
                sm = sm+w(k)*u(i,k)*vt(k,j);
            }
            if( isupper )
            {
                if( i==j )
                {
                    locerr = ap::maxreal(locerr, fabs(d(i)-sm));
                }
                else
                {
                    if( i==j-1 )
                    {
                        locerr = ap::maxreal(locerr, fabs(e(i)-sm));
                    }
                    else
                    {
                        locerr = ap::maxreal(locerr, fabs(sm));
                    }
                }
            }
            else
            {
                if( i==j )
                {
                    locerr = ap::maxreal(locerr, fabs(d(i)-sm));
                }
                else
                {
                    if( i-1==j )
                    {
                        locerr = ap::maxreal(locerr, fabs(e(j)-sm));
                    }
                    else
                    {
                        locerr = ap::maxreal(locerr, fabs(sm));
                    }
                }
            }
        }
    }
    materr = ap::maxreal(materr, locerr);
    
    //
    // check for C = U'
    // we consider it as decomposition error
    //
    locerr = 0;
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            locerr = ap::maxreal(locerr, fabs(u(i,j)-c(j,i)));
        }
    }
    materr = ap::maxreal(materr, locerr);
    
    //
    // orthogonality error
    //
    locerr = 0;
    for(i = 0; i <= n-1; i++)
    {
        for(j = i; j <= n-1; j++)
        {
            sm = ap::vdotproduct(&u(0, i), u.getstride(), &u(0, j), u.getstride(), ap::vlen(0,n-1));
            if( i!=j )
            {
                locerr = ap::maxreal(locerr, fabs(sm));
            }
            else
            {
                locerr = ap::maxreal(locerr, fabs(sm-1));
            }
            sm = ap::vdotproduct(&vt(i, 0), 1, &vt(j, 0), 1, ap::vlen(0,n-1));
            if( i!=j )
            {
                locerr = ap::maxreal(locerr, fabs(sm));
            }
            else
            {
                locerr = ap::maxreal(locerr, fabs(sm-1));
            }
        }
    }
    orterr = ap::maxreal(orterr, locerr);
    
    //
    // values order error
    //
    for(i = 1; i <= n-1; i++)
    {
        if( ap::fp_greater(w(i),w(i-1)) )
        {
            wsorted = false;
        }
    }
}


static void testbdsvdproblem(const ap::real_1d_array& d,
     const ap::real_1d_array& e,
     int n,
     double& materr,
     double& orterr,
     bool& wsorted,
     bool& wfailed)
{
    ap::real_2d_array u;
    ap::real_2d_array vt;
    ap::real_2d_array c;
    ap::real_1d_array w;
    int i;
    int j;
    int k;
    double v;
    double mx;

    mx = 0;
    for(i = 0; i <= n-1; i++)
    {
        if( ap::fp_greater(fabs(d(i)),mx) )
        {
            mx = fabs(d(i));
        }
    }
    for(i = 0; i <= n-2; i++)
    {
        if( ap::fp_greater(fabs(e(i)),mx) )
        {
            mx = fabs(e(i));
        }
    }
    if( ap::fp_eq(mx,0) )
    {
        mx = 1;
    }
    
    //
    // Upper BDSVD tests
    //
    w.setbounds(0, n-1);
    fillidentity(u, n);
    fillidentity(vt, n);
    fillidentity(c, n);
    for(i = 0; i <= n-1; i++)
    {
        w(i) = d(i);
    }
    if( !rmatrixbdsvd(w, e, n, true, false, u, n, c, n, vt, n) )
    {
        failcount = failcount+1;
        wfailed = true;
        return;
    }
    getbdsvderror(d, e, n, true, u, c, w, vt, materr, orterr, wsorted);
    fillidentity(u, n);
    fillidentity(vt, n);
    fillidentity(c, n);
    for(i = 0; i <= n-1; i++)
    {
        w(i) = d(i);
    }
    if( !rmatrixbdsvd(w, e, n, true, true, u, n, c, n, vt, n) )
    {
        failcount = failcount+1;
        wfailed = true;
        return;
    }
    getbdsvderror(d, e, n, true, u, c, w, vt, materr, orterr, wsorted);
    
    //
    // Lower BDSVD tests
    //
    w.setbounds(0, n-1);
    fillidentity(u, n);
    fillidentity(vt, n);
    fillidentity(c, n);
    for(i = 0; i <= n-1; i++)
    {
        w(i) = d(i);
    }
    if( !rmatrixbdsvd(w, e, n, false, false, u, n, c, n, vt, n) )
    {
        failcount = failcount+1;
        wfailed = true;
        return;
    }
    getbdsvderror(d, e, n, false, u, c, w, vt, materr, orterr, wsorted);
    fillidentity(u, n);
    fillidentity(vt, n);
    fillidentity(c, n);
    for(i = 0; i <= n-1; i++)
    {
        w(i) = d(i);
    }
    if( !rmatrixbdsvd(w, e, n, false, true, u, n, c, n, vt, n) )
    {
        failcount = failcount+1;
        wfailed = true;
        return;
    }
    getbdsvderror(d, e, n, false, u, c, w, vt, materr, orterr, wsorted);
    
    //
    // update counter
    //
    succcount = succcount+1;
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testbdsvdunit_test_silent()
{
    bool result;

    result = testbdsvd(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testbdsvdunit_test()
{
    bool result;

    result = testbdsvd(false);
    return result;
}




