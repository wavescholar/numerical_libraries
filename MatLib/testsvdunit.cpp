
 
#include <stdio.h>
#include "testsvdunit.h"

static int failcount;
static int succcount;

static void fillsparsea(ap::real_2d_array& a, int m, int n, double sparcity);
static void getsvderror(const ap::real_2d_array& a,
     int m,
     int n,
     const ap::real_2d_array& u,
     const ap::real_1d_array& w,
     const ap::real_2d_array& vt,
     double& materr,
     double& orterr,
     bool& wsorted);
static void testsvdproblem(const ap::real_2d_array& a,
     int m,
     int n,
     double& materr,
     double& orterr,
     double& othererr,
     bool& wsorted,
     bool& wfailed);

/*************************************************************************
Testing SVD decomposition subroutine
*************************************************************************/
bool testsvd(bool silent)
{
    bool result;
    ap::real_2d_array a;
    int m;
    int n;
    int maxmn;
    int i;
    int j;
    int gpass;
    int pass;
    bool waserrors;
    bool wsorted;
    bool wfailed;
    double materr;
    double orterr;
    double othererr;
    double threshold;
    double failthreshold;
    double failr;

    failcount = 0;
    succcount = 0;
    materr = 0;
    orterr = 0;
    othererr = 0;
    wsorted = true;
    wfailed = false;
    waserrors = false;
    maxmn = 30;
    threshold = 5*100*ap::machineepsilon;
    failthreshold = 5.0E-3;
    a.setbounds(0, maxmn-1, 0, maxmn-1);
    
    //
    // TODO: div by zero fail, convergence fail
    //
    for(gpass = 1; gpass <= 1; gpass++)
    {
        
        //
        // zero matrix, several cases
        //
        for(i = 0; i <= maxmn-1; i++)
        {
            for(j = 0; j <= maxmn-1; j++)
            {
                a(i,j) = 0;
            }
        }
        for(i = 1; i <= ap::minint(5, maxmn); i++)
        {
            for(j = 1; j <= ap::minint(5, maxmn); j++)
            {
                testsvdproblem(a, i, j, materr, orterr, othererr, wsorted, wfailed);
            }
        }
        
        //
        // Long dense matrix
        //
        for(i = 0; i <= maxmn-1; i++)
        {
            for(j = 0; j <= ap::minint(5, maxmn)-1; j++)
            {
                a(i,j) = 2*ap::randomreal()-1;
            }
        }
        for(i = 1; i <= maxmn; i++)
        {
            for(j = 1; j <= ap::minint(5, maxmn); j++)
            {
                testsvdproblem(a, i, j, materr, orterr, othererr, wsorted, wfailed);
            }
        }
        for(i = 0; i <= ap::minint(5, maxmn)-1; i++)
        {
            for(j = 0; j <= maxmn-1; j++)
            {
                a(i,j) = 2*ap::randomreal()-1;
            }
        }
        for(i = 1; i <= ap::minint(5, maxmn); i++)
        {
            for(j = 1; j <= maxmn; j++)
            {
                testsvdproblem(a, i, j, materr, orterr, othererr, wsorted, wfailed);
            }
        }
        
        //
        // Dense matrices
        //
        for(m = 1; m <= ap::minint(10, maxmn); m++)
        {
            for(n = 1; n <= ap::minint(10, maxmn); n++)
            {
                for(i = 0; i <= m-1; i++)
                {
                    for(j = 0; j <= n-1; j++)
                    {
                        a(i,j) = 2*ap::randomreal()-1;
                    }
                }
                testsvdproblem(a, m, n, materr, orterr, othererr, wsorted, wfailed);
            }
        }
        
        //
        // Sparse matrices, very sparse matrices, incredible sparse matrices
        //
        for(m = 1; m <= 10; m++)
        {
            for(n = 1; n <= 10; n++)
            {
                for(pass = 1; pass <= 2; pass++)
                {
                    fillsparsea(a, m, n, 0.8);
                    testsvdproblem(a, m, n, materr, orterr, othererr, wsorted, wfailed);
                    fillsparsea(a, m, n, 0.9);
                    testsvdproblem(a, m, n, materr, orterr, othererr, wsorted, wfailed);
                    fillsparsea(a, m, n, 0.95);
                    testsvdproblem(a, m, n, materr, orterr, othererr, wsorted, wfailed);
                }
            }
        }
    }
    
    //
    // report
    //
    failr = double(failcount)/double(succcount+failcount);
    waserrors = ap::fp_greater(materr,threshold)||ap::fp_greater(orterr,threshold)||ap::fp_greater(othererr,threshold)||!wsorted||ap::fp_greater(failr,failthreshold);
    if( !silent )
    {
        printf("TESTING SVD DECOMPOSITION\n");
        printf("SVD decomposition error:                 %5.3le\n",
            double(materr));
        printf("SVD orthogonality error:                 %5.3le\n",
            double(orterr));
        printf("SVD with different parameters error:     %5.3le\n",
            double(othererr));
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


static void fillsparsea(ap::real_2d_array& a, int m, int n, double sparcity)
{
    int i;
    int j;

    for(i = 0; i <= m-1; i++)
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


static void getsvderror(const ap::real_2d_array& a,
     int m,
     int n,
     const ap::real_2d_array& u,
     const ap::real_1d_array& w,
     const ap::real_2d_array& vt,
     double& materr,
     double& orterr,
     bool& wsorted)
{
    int i;
    int j;
    int k;
    int minmn;
    double locerr;
    double sm;

    minmn = ap::minint(m, n);
    
    //
    // decomposition error
    //
    locerr = 0;
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            sm = 0;
            for(k = 0; k <= minmn-1; k++)
            {
                sm = sm+w(k)*u(i,k)*vt(k,j);
            }
            locerr = ap::maxreal(locerr, fabs(a(i,j)-sm));
        }
    }
    materr = ap::maxreal(materr, locerr);
    
    //
    // orthogonality error
    //
    locerr = 0;
    for(i = 0; i <= minmn-1; i++)
    {
        for(j = i; j <= minmn-1; j++)
        {
            sm = ap::vdotproduct(&u(0, i), u.getstride(), &u(0, j), u.getstride(), ap::vlen(0,m-1));
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
    for(i = 1; i <= minmn-1; i++)
    {
        if( ap::fp_greater(w(i),w(i-1)) )
        {
            wsorted = false;
        }
    }
}


static void testsvdproblem(const ap::real_2d_array& a,
     int m,
     int n,
     double& materr,
     double& orterr,
     double& othererr,
     bool& wsorted,
     bool& wfailed)
{
    ap::real_2d_array u;
    ap::real_2d_array vt;
    ap::real_2d_array u2;
    ap::real_2d_array vt2;
    ap::real_1d_array w;
    ap::real_1d_array w2;
    int i;
    int j;
    int k;
    int ujob;
    int vtjob;
    int memjob;
    int ucheck;
    int vtcheck;
    double v;
    double mx;

    
    //
    // Main SVD test
    //
    if( !rmatrixsvd(a, m, n, 2, 2, 2, w, u, vt) )
    {
        failcount = failcount+1;
        wfailed = true;
        return;
    }
    getsvderror(a, m, n, u, w, vt, materr, orterr, wsorted);
    
    //
    // Additional SVD tests
    //
    for(ujob = 0; ujob <= 2; ujob++)
    {
        for(vtjob = 0; vtjob <= 2; vtjob++)
        {
            for(memjob = 0; memjob <= 2; memjob++)
            {
                if( !rmatrixsvd(a, m, n, ujob, vtjob, memjob, w2, u2, vt2) )
                {
                    failcount = failcount+1;
                    wfailed = true;
                    return;
                }
                ucheck = 0;
                if( ujob==1 )
                {
                    ucheck = ap::minint(m, n);
                }
                if( ujob==2 )
                {
                    ucheck = m;
                }
                vtcheck = 0;
                if( vtjob==1 )
                {
                    vtcheck = ap::minint(m, n);
                }
                if( vtjob==2 )
                {
                    vtcheck = n;
                }
                for(i = 0; i <= m-1; i++)
                {
                    for(j = 0; j <= ucheck-1; j++)
                    {
                        othererr = ap::maxreal(othererr, fabs(u(i,j)-u2(i,j)));
                    }
                }
                for(i = 0; i <= vtcheck-1; i++)
                {
                    for(j = 0; j <= n-1; j++)
                    {
                        othererr = ap::maxreal(othererr, fabs(vt(i,j)-vt2(i,j)));
                    }
                }
                for(i = 0; i <= ap::minint(m, n)-1; i++)
                {
                    othererr = ap::maxreal(othererr, fabs(w(i)-w2(i)));
                }
            }
        }
    }
    
    //
    // update counter
    //
    succcount = succcount+1;
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testsvdunit_test_silent()
{
    bool result;

    result = testsvd(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testsvdunit_test()
{
    bool result;

    result = testsvd(false);
    return result;
}




