
 
#include <stdio.h>
#include "testdensesolverunit.h"

static bool rmatrixchecksolutionm(const ap::real_2d_array& xe,
     int n,
     int m,
     double threshold,
     int info,
     const densesolverreport& rep,
     const ap::real_2d_array& xs);
static bool rmatrixchecksolution(const ap::real_2d_array& xe,
     int n,
     double threshold,
     int info,
     const densesolverreport& rep,
     const ap::real_1d_array& xs);
static bool rmatrixchecksingularm(int n,
     int m,
     int info,
     const densesolverreport& rep,
     const ap::real_2d_array& xs);
static bool rmatrixchecksingular(int n,
     int info,
     const densesolverreport& rep,
     const ap::real_1d_array& xs);
static bool cmatrixchecksolutionm(const ap::complex_2d_array& xe,
     int n,
     int m,
     double threshold,
     int info,
     const densesolverreport& rep,
     const ap::complex_2d_array& xs);
static bool cmatrixchecksolution(const ap::complex_2d_array& xe,
     int n,
     double threshold,
     int info,
     const densesolverreport& rep,
     const ap::complex_1d_array& xs);
static bool cmatrixchecksingularm(int n,
     int m,
     int info,
     const densesolverreport& rep,
     const ap::complex_2d_array& xs);
static bool cmatrixchecksingular(int n,
     int info,
     const densesolverreport& rep,
     const ap::complex_1d_array& xs);
static void rmatrixmakeacopy(const ap::real_2d_array& a,
     int m,
     int n,
     ap::real_2d_array& b);
static void cmatrixmakeacopy(const ap::complex_2d_array& a,
     int m,
     int n,
     ap::complex_2d_array& b);
static void rmatrixdrophalf(ap::real_2d_array& a, int n, bool droplower);
static void cmatrixdrophalf(ap::complex_2d_array& a, int n, bool droplower);
static void testrsolver(int maxn,
     int maxm,
     int passcount,
     double threshold,
     bool& rerrors,
     bool& rfserrors);
static void testspdsolver(int maxn,
     int maxm,
     int passcount,
     double threshold,
     bool& spderrors,
     bool& rfserrors);
static void testcsolver(int maxn,
     int maxm,
     int passcount,
     double threshold,
     bool& cerrors,
     bool& rfserrors);
static void testhpdsolver(int maxn,
     int maxm,
     int passcount,
     double threshold,
     bool& hpderrors,
     bool& rfserrors);
static void unset2d(ap::real_2d_array& x);
static void unset1d(ap::real_1d_array& x);
static void cunset2d(ap::complex_2d_array& x);
static void cunset1d(ap::complex_1d_array& x);
static void unsetrep(densesolverreport& r);
static void unsetlsrep(densesolverlsreport& r);

/*************************************************************************
Test
*************************************************************************/
bool testdensesolver(bool silent)
{
    bool result;
    ap::real_2d_array a;
    ap::real_2d_array lua;
    ap::real_2d_array atmp;
    ap::integer_1d_array p;
    ap::real_2d_array xe;
    ap::real_2d_array b;
    ap::real_1d_array bv;
    int i;
    int j;
    int k;
    int n;
    int m;
    int pass;
    int taskkind;
    double mx;
    double v;
    double verr;
    int info;
    densesolverreport rep;
    densesolverlsreport repls;
    ap::real_2d_array x;
    ap::real_1d_array xv;
    ap::real_1d_array y;
    ap::real_1d_array tx;
    int maxn;
    int maxm;
    int passcount;
    double threshold;
    bool rerrors;
    bool cerrors;
    bool spderrors;
    bool hpderrors;
    bool rfserrors;
    bool waserrors;

    maxn = 10;
    maxm = 5;
    passcount = 5;
    threshold = 10000*ap::machineepsilon;
    rfserrors = false;
    rerrors = false;
    cerrors = false;
    spderrors = false;
    hpderrors = false;
    testrsolver(maxn, maxm, passcount, threshold, rerrors, rfserrors);
    testspdsolver(maxn, maxm, passcount, threshold, spderrors, rfserrors);
    testcsolver(maxn, maxm, passcount, threshold, cerrors, rfserrors);
    testhpdsolver(maxn, maxm, passcount, threshold, hpderrors, rfserrors);
    waserrors = rerrors||cerrors||spderrors||hpderrors||rfserrors;
    if( !silent )
    {
        printf("TESTING DENSE SOLVER\n");
        printf("* REAL:                                   ");
        if( rerrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* COMPLEX:                                ");
        if( cerrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* SPD:                                    ");
        if( spderrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* HPD:                                    ");
        if( hpderrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* ITERATIVE IMPROVEMENT:                  ");
        if( rfserrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        if( waserrors )
        {
            printf("TEST FAILED\n");
        }
        else
        {
            printf("TEST PASSED\n");
        }
    }
    result = !waserrors;
    return result;
}


/*************************************************************************
Checks whether solver results are correct solution.
Returns True on success.
*************************************************************************/
static bool rmatrixchecksolutionm(const ap::real_2d_array& xe,
     int n,
     int m,
     double threshold,
     int info,
     const densesolverreport& rep,
     const ap::real_2d_array& xs)
{
    bool result;
    int i;
    int j;

    result = true;
    if( info<=0 )
    {
        result = false;
    }
    else
    {
        result = result&&!(ap::fp_less(rep.r1,100*ap::machineepsilon)||ap::fp_greater(rep.r1,1+1000*ap::machineepsilon));
        result = result&&!(ap::fp_less(rep.rinf,100*ap::machineepsilon)||ap::fp_greater(rep.rinf,1+1000*ap::machineepsilon));
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                result = result&&ap::fp_less_eq(fabs(xe(i,j)-xs(i,j)),threshold);
            }
        }
    }
    return result;
}


/*************************************************************************
Checks whether solver results are correct solution.
Returns True on success.
*************************************************************************/
static bool rmatrixchecksolution(const ap::real_2d_array& xe,
     int n,
     double threshold,
     int info,
     const densesolverreport& rep,
     const ap::real_1d_array& xs)
{
    bool result;
    ap::real_2d_array xsm;

    xsm.setlength(n, 1);
    ap::vmove(&xsm(0, 0), xsm.getstride(), &xs(0), 1, ap::vlen(0,n-1));
    result = rmatrixchecksolutionm(xe, n, 1, threshold, info, rep, xsm);
    return result;
}


/*************************************************************************
Checks whether solver results indicate singular matrix.
Returns True on success.
*************************************************************************/
static bool rmatrixchecksingularm(int n,
     int m,
     int info,
     const densesolverreport& rep,
     const ap::real_2d_array& xs)
{
    bool result;
    int i;
    int j;

    result = true;
    if( info!=-3&&info!=1 )
    {
        result = false;
    }
    else
    {
        result = result&&!(ap::fp_less(rep.r1,0)||ap::fp_greater(rep.r1,1000*ap::machineepsilon));
        result = result&&!(ap::fp_less(rep.rinf,0)||ap::fp_greater(rep.rinf,1000*ap::machineepsilon));
        if( info==-3 )
        {
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= m-1; j++)
                {
                    result = result&&ap::fp_eq(xs(i,j),0);
                }
            }
        }
    }
    return result;
}


/*************************************************************************
Checks whether solver results indicate singular matrix.
Returns True on success.
*************************************************************************/
static bool rmatrixchecksingular(int n,
     int info,
     const densesolverreport& rep,
     const ap::real_1d_array& xs)
{
    bool result;
    ap::real_2d_array xsm;

    xsm.setlength(n, 1);
    ap::vmove(&xsm(0, 0), xsm.getstride(), &xs(0), 1, ap::vlen(0,n-1));
    result = rmatrixchecksingularm(n, 1, info, rep, xsm);
    return result;
}


/*************************************************************************
Checks whether solver results are correct solution.
Returns True on success.
*************************************************************************/
static bool cmatrixchecksolutionm(const ap::complex_2d_array& xe,
     int n,
     int m,
     double threshold,
     int info,
     const densesolverreport& rep,
     const ap::complex_2d_array& xs)
{
    bool result;
    int i;
    int j;

    result = true;
    if( info<=0 )
    {
        result = false;
    }
    else
    {
        result = result&&!(ap::fp_less(rep.r1,100*ap::machineepsilon)||ap::fp_greater(rep.r1,1+1000*ap::machineepsilon));
        result = result&&!(ap::fp_less(rep.rinf,100*ap::machineepsilon)||ap::fp_greater(rep.rinf,1+1000*ap::machineepsilon));
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                result = result&&ap::fp_less_eq(ap::abscomplex(xe(i,j)-xs(i,j)),threshold);
            }
        }
    }
    return result;
}


/*************************************************************************
Checks whether solver results are correct solution.
Returns True on success.
*************************************************************************/
static bool cmatrixchecksolution(const ap::complex_2d_array& xe,
     int n,
     double threshold,
     int info,
     const densesolverreport& rep,
     const ap::complex_1d_array& xs)
{
    bool result;
    ap::complex_2d_array xsm;

    xsm.setlength(n, 1);
    ap::vmove(&xsm(0, 0), xsm.getstride(), &xs(0), 1, "N", ap::vlen(0,n-1));
    result = cmatrixchecksolutionm(xe, n, 1, threshold, info, rep, xsm);
    return result;
}


/*************************************************************************
Checks whether solver results indicate singular matrix.
Returns True on success.
*************************************************************************/
static bool cmatrixchecksingularm(int n,
     int m,
     int info,
     const densesolverreport& rep,
     const ap::complex_2d_array& xs)
{
    bool result;
    int i;
    int j;

    result = true;
    if( info!=-3&&info!=1 )
    {
        result = false;
    }
    else
    {
        result = result&&!(ap::fp_less(rep.r1,0)||ap::fp_greater(rep.r1,1000*ap::machineepsilon));
        result = result&&!(ap::fp_less(rep.rinf,0)||ap::fp_greater(rep.rinf,1000*ap::machineepsilon));
        if( info==-3 )
        {
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= m-1; j++)
                {
                    result = result&&xs(i,j)==0;
                }
            }
        }
    }
    return result;
}


/*************************************************************************
Checks whether solver results indicate singular matrix.
Returns True on success.
*************************************************************************/
static bool cmatrixchecksingular(int n,
     int info,
     const densesolverreport& rep,
     const ap::complex_1d_array& xs)
{
    bool result;
    ap::complex_2d_array xsm;

    xsm.setlength(n, 1);
    ap::vmove(&xsm(0, 0), xsm.getstride(), &xs(0), 1, "N", ap::vlen(0,n-1));
    result = cmatrixchecksingularm(n, 1, info, rep, xsm);
    return result;
}


/*************************************************************************
Copy
*************************************************************************/
static void rmatrixmakeacopy(const ap::real_2d_array& a,
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
Copy
*************************************************************************/
static void cmatrixmakeacopy(const ap::complex_2d_array& a,
     int m,
     int n,
     ap::complex_2d_array& b)
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
Drops upper or lower half of the matrix - fills it by special pattern
which may be used later to ensure that this part wasn't changed
*************************************************************************/
static void rmatrixdrophalf(ap::real_2d_array& a, int n, bool droplower)
{
    int i;
    int j;

    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            if( droplower&&i>j||!droplower&&i<j )
            {
                a(i,j) = 1+2*i+3*j;
            }
        }
    }
}


/*************************************************************************
Drops upper or lower half of the matrix - fills it by special pattern
which may be used later to ensure that this part wasn't changed
*************************************************************************/
static void cmatrixdrophalf(ap::complex_2d_array& a, int n, bool droplower)
{
    int i;
    int j;

    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            if( droplower&&i>j||!droplower&&i<j )
            {
                a(i,j) = 1+2*i+3*j;
            }
        }
    }
}


/*************************************************************************
Real test
*************************************************************************/
static void testrsolver(int maxn,
     int maxm,
     int passcount,
     double threshold,
     bool& rerrors,
     bool& rfserrors)
{
    ap::real_2d_array a;
    ap::real_2d_array lua;
    ap::real_2d_array atmp;
    ap::integer_1d_array p;
    ap::real_2d_array xe;
    ap::real_2d_array b;
    ap::real_1d_array bv;
    int i;
    int j;
    int k;
    int n;
    int m;
    int pass;
    int taskkind;
    double mx;
    double v;
    double verr;
    int info;
    densesolverreport rep;
    densesolverlsreport repls;
    ap::real_2d_array x;
    ap::real_1d_array xv;
    ap::real_1d_array y;
    ap::real_1d_array tx;

    
    //
    // General square matrices:
    // * test general solvers
    // * test least squares solver
    //
    for(pass = 1; pass <= passcount; pass++)
    {
        for(n = 1; n <= maxn; n++)
        {
            for(m = 1; m <= maxm; m++)
            {
                
                //
                // ********************************************************
                // WELL CONDITIONED TASKS
                // ability to find correct solution is tested
                // ********************************************************
                //
                // 1. generate random well conditioned matrix A.
                // 2. generate random solution vector xe
                // 3. generate right part b=A*xe
                // 4. test different methods on original A
                //
                rmatrixrndcond(n, double(1000), a);
                rmatrixmakeacopy(a, n, n, lua);
                rmatrixlu(lua, n, n, p);
                xe.setlength(n, m);
                for(i = 0; i <= n-1; i++)
                {
                    for(j = 0; j <= m-1; j++)
                    {
                        xe(i,j) = 2*ap::randomreal()-1;
                    }
                }
                b.setlength(n, m);
                for(i = 0; i <= n-1; i++)
                {
                    for(j = 0; j <= m-1; j++)
                    {
                        v = ap::vdotproduct(&a(i, 0), 1, &xe(0, j), xe.getstride(), ap::vlen(0,n-1));
                        b(i,j) = v;
                    }
                }
                
                //
                // Test solvers
                //
                info = 0;
                unsetrep(rep);
                unset2d(x);
                rmatrixsolvem(a, n, b, m, ap::fp_greater(ap::randomreal(),0.5), info, rep, x);
                rerrors = rerrors||!rmatrixchecksolutionm(xe, n, m, threshold, info, rep, x);
                info = 0;
                unsetrep(rep);
                unset1d(xv);
                bv.setlength(n);
                ap::vmove(&bv(0), 1, &b(0, 0), b.getstride(), ap::vlen(0,n-1));
                rmatrixsolve(a, n, bv, info, rep, xv);
                rerrors = rerrors||!rmatrixchecksolution(xe, n, threshold, info, rep, xv);
                info = 0;
                unsetrep(rep);
                unset2d(x);
                rmatrixlusolvem(lua, p, n, b, m, info, rep, x);
                rerrors = rerrors||!rmatrixchecksolutionm(xe, n, m, threshold, info, rep, x);
                info = 0;
                unsetrep(rep);
                unset1d(xv);
                bv.setlength(n);
                ap::vmove(&bv(0), 1, &b(0, 0), b.getstride(), ap::vlen(0,n-1));
                rmatrixlusolve(lua, p, n, bv, info, rep, xv);
                rerrors = rerrors||!rmatrixchecksolution(xe, n, threshold, info, rep, xv);
                info = 0;
                unsetrep(rep);
                unset2d(x);
                rmatrixmixedsolvem(a, lua, p, n, b, m, info, rep, x);
                rerrors = rerrors||!rmatrixchecksolutionm(xe, n, m, threshold, info, rep, x);
                info = 0;
                unsetrep(rep);
                unset1d(xv);
                bv.setlength(n);
                ap::vmove(&bv(0), 1, &b(0, 0), b.getstride(), ap::vlen(0,n-1));
                rmatrixmixedsolve(a, lua, p, n, bv, info, rep, xv);
                rerrors = rerrors||!rmatrixchecksolution(xe, n, threshold, info, rep, xv);
                
                //
                // Test DenseSolverRLS():
                // * test on original system A*x = b
                // * test on overdetermined system with the same solution: (A' A')'*x = (b' b')'
                // * test on underdetermined system with the same solution: (A 0 0 0 ) * z = b
                //
                info = 0;
                unsetlsrep(repls);
                unset1d(xv);
                bv.setlength(n);
                ap::vmove(&bv(0), 1, &b(0, 0), b.getstride(), ap::vlen(0,n-1));
                rmatrixsolvels(a, n, n, bv, 0.0, info, repls, xv);
                if( info<=0 )
                {
                    rerrors = true;
                }
                else
                {
                    rerrors = rerrors||ap::fp_less(repls.r2,100*ap::machineepsilon)||ap::fp_greater(repls.r2,1+1000*ap::machineepsilon);
                    rerrors = rerrors||repls.n!=n||repls.k!=0;
                    for(i = 0; i <= n-1; i++)
                    {
                        rerrors = rerrors||ap::fp_greater(fabs(xe(i,0)-xv(i)),threshold);
                    }
                }
                info = 0;
                unsetlsrep(repls);
                unset1d(xv);
                bv.setlength(2*n);
                ap::vmove(&bv(0), 1, &b(0, 0), b.getstride(), ap::vlen(0,n-1));
                ap::vmove(&bv(n), 1, &b(0, 0), b.getstride(), ap::vlen(n,2*n-1));
                atmp.setlength(2*n, n);
                copymatrix(a, 0, n-1, 0, n-1, atmp, 0, n-1, 0, n-1);
                copymatrix(a, 0, n-1, 0, n-1, atmp, n, 2*n-1, 0, n-1);
                rmatrixsolvels(atmp, 2*n, n, bv, 0.0, info, repls, xv);
                if( info<=0 )
                {
                    rerrors = true;
                }
                else
                {
                    rerrors = rerrors||ap::fp_less(repls.r2,100*ap::machineepsilon)||ap::fp_greater(repls.r2,1+1000*ap::machineepsilon);
                    rerrors = rerrors||repls.n!=n||repls.k!=0;
                    for(i = 0; i <= n-1; i++)
                    {
                        rerrors = rerrors||ap::fp_greater(fabs(xe(i,0)-xv(i)),threshold);
                    }
                }
                info = 0;
                unsetlsrep(repls);
                unset1d(xv);
                bv.setlength(n);
                ap::vmove(&bv(0), 1, &b(0, 0), b.getstride(), ap::vlen(0,n-1));
                atmp.setlength(n, 2*n);
                copymatrix(a, 0, n-1, 0, n-1, atmp, 0, n-1, 0, n-1);
                for(i = 0; i <= n-1; i++)
                {
                    for(j = n; j <= 2*n-1; j++)
                    {
                        atmp(i,j) = 0;
                    }
                }
                rmatrixsolvels(atmp, n, 2*n, bv, 0.0, info, repls, xv);
                if( info<=0 )
                {
                    rerrors = true;
                }
                else
                {
                    rerrors = rerrors||ap::fp_neq(repls.r2,0);
                    rerrors = rerrors||repls.n!=2*n||repls.k!=n;
                    for(i = 0; i <= n-1; i++)
                    {
                        rerrors = rerrors||ap::fp_greater(fabs(xe(i,0)-xv(i)),threshold);
                    }
                    for(i = n; i <= 2*n-1; i++)
                    {
                        rerrors = rerrors||ap::fp_greater(fabs(xv(i)),threshold);
                    }
                }
                
                //
                // ********************************************************
                // EXACTLY SINGULAR MATRICES
                // ability to detect singularity is tested
                // ********************************************************
                //
                // 1. generate different types of singular matrices:
                //    * zero
                //    * with zero columns
                //    * with zero rows
                //    * with equal rows/columns
                // 2. generate random solution vector xe
                // 3. generate right part b=A*xe
                // 4. test different methods
                //
                for(taskkind = 0; taskkind <= 4; taskkind++)
                {
                    unset2d(a);
                    if( taskkind==0 )
                    {
                        
                        //
                        // all zeros
                        //
                        a.setlength(n, n);
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = 0; j <= n-1; j++)
                            {
                                a(i,j) = 0;
                            }
                        }
                    }
                    if( taskkind==1 )
                    {
                        
                        //
                        // there is zero column
                        //
                        a.setlength(n, n);
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = 0; j <= n-1; j++)
                            {
                                a(i,j) = 2*ap::randomreal()-1;
                            }
                        }
                        k = ap::randominteger(n);
                        ap::vmul(&a(0, k), a.getstride(), ap::vlen(0,n-1), 0);
                    }
                    if( taskkind==2 )
                    {
                        
                        //
                        // there is zero row
                        //
                        a.setlength(n, n);
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = 0; j <= n-1; j++)
                            {
                                a(i,j) = 2*ap::randomreal()-1;
                            }
                        }
                        k = ap::randominteger(n);
                        ap::vmul(&a(k, 0), 1, ap::vlen(0,n-1), 0);
                    }
                    if( taskkind==3 )
                    {
                        
                        //
                        // equal columns
                        //
                        if( n<2 )
                        {
                            continue;
                        }
                        a.setlength(n, n);
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = 0; j <= n-1; j++)
                            {
                                a(i,j) = 2*ap::randomreal()-1;
                            }
                        }
                        k = 1+ap::randominteger(n-1);
                        ap::vmove(&a(0, 0), a.getstride(), &a(0, k), a.getstride(), ap::vlen(0,n-1));
                    }
                    if( taskkind==4 )
                    {
                        
                        //
                        // equal rows
                        //
                        if( n<2 )
                        {
                            continue;
                        }
                        a.setlength(n, n);
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = 0; j <= n-1; j++)
                            {
                                a(i,j) = 2*ap::randomreal()-1;
                            }
                        }
                        k = 1+ap::randominteger(n-1);
                        ap::vmove(&a(0, 0), 1, &a(k, 0), 1, ap::vlen(0,n-1));
                    }
                    xe.setlength(n, m);
                    for(i = 0; i <= n-1; i++)
                    {
                        for(j = 0; j <= m-1; j++)
                        {
                            xe(i,j) = 2*ap::randomreal()-1;
                        }
                    }
                    b.setlength(n, m);
                    for(i = 0; i <= n-1; i++)
                    {
                        for(j = 0; j <= m-1; j++)
                        {
                            v = ap::vdotproduct(&a(i, 0), 1, &xe(0, j), xe.getstride(), ap::vlen(0,n-1));
                            b(i,j) = v;
                        }
                    }
                    rmatrixmakeacopy(a, n, n, lua);
                    rmatrixlu(lua, n, n, p);
                    
                    //
                    // Test RMatrixSolveM()
                    //
                    info = 0;
                    unsetrep(rep);
                    unset2d(x);
                    rmatrixsolvem(a, n, b, m, ap::fp_greater(ap::randomreal(),0.5), info, rep, x);
                    rerrors = rerrors||!rmatrixchecksingularm(n, m, info, rep, x);
                    
                    //
                    // Test RMatrixSolve()
                    //
                    info = 0;
                    unsetrep(rep);
                    unset2d(x);
                    bv.setlength(n);
                    ap::vmove(&bv(0), 1, &b(0, 0), b.getstride(), ap::vlen(0,n-1));
                    rmatrixsolve(a, n, bv, info, rep, xv);
                    rerrors = rerrors||!rmatrixchecksingular(n, info, rep, xv);
                    
                    //
                    // Test RMatrixLUSolveM()
                    //
                    info = 0;
                    unsetrep(rep);
                    unset2d(x);
                    rmatrixlusolvem(lua, p, n, b, m, info, rep, x);
                    rerrors = rerrors||!rmatrixchecksingularm(n, m, info, rep, x);
                    
                    //
                    // Test RMatrixLUSolve()
                    //
                    info = 0;
                    unsetrep(rep);
                    unset2d(x);
                    bv.setlength(n);
                    ap::vmove(&bv(0), 1, &b(0, 0), b.getstride(), ap::vlen(0,n-1));
                    rmatrixlusolve(lua, p, n, bv, info, rep, xv);
                    rerrors = rerrors||!rmatrixchecksingular(n, info, rep, xv);
                    
                    //
                    // Test RMatrixMixedSolveM()
                    //
                    info = 0;
                    unsetrep(rep);
                    unset2d(x);
                    rmatrixmixedsolvem(a, lua, p, n, b, m, info, rep, x);
                    rerrors = rerrors||!rmatrixchecksingularm(n, m, info, rep, x);
                    
                    //
                    // Test RMatrixMixedSolve()
                    //
                    info = 0;
                    unsetrep(rep);
                    unset2d(x);
                    bv.setlength(n);
                    ap::vmove(&bv(0), 1, &b(0, 0), b.getstride(), ap::vlen(0,n-1));
                    rmatrixmixedsolve(a, lua, p, n, bv, info, rep, xv);
                    rerrors = rerrors||!rmatrixchecksingular(n, info, rep, xv);
                }
            }
        }
    }
    
    //
    // test iterative improvement
    //
    for(pass = 1; pass <= passcount; pass++)
    {
        
        //
        // Test iterative improvement matrices
        //
        // A matrix/right part are constructed such that both matrix
        // and solution components are within (-1,+1). Such matrix/right part
        // have nice properties - system can be solved using iterative
        // improvement with |A*x-b| about several ulps of max(1,|b|).
        //
        n = 100;
        a.setlength(n, n);
        b.setlength(n, 1);
        bv.setlength(n);
        tx.setlength(n);
        xv.setlength(n);
        y.setlength(n);
        for(i = 0; i <= n-1; i++)
        {
            xv(i) = 2*ap::randomreal()-1;
        }
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                a(i,j) = 2*ap::randomreal()-1;
            }
            ap::vmove(&y(0), 1, &a(i, 0), 1, ap::vlen(0,n-1));
            xdot(y, xv, n, tx, v, verr);
            bv(i) = v;
        }
        ap::vmove(&b(0, 0), b.getstride(), &bv(0), 1, ap::vlen(0,n-1));
        
        //
        // Test RMatrixSolveM()
        //
        unset2d(x);
        rmatrixsolvem(a, n, b, 1, true, info, rep, x);
        if( info<=0 )
        {
            rfserrors = true;
        }
        else
        {
            xv.setlength(n);
            ap::vmove(&xv(0), 1, &x(0, 0), x.getstride(), ap::vlen(0,n-1));
            for(i = 0; i <= n-1; i++)
            {
                ap::vmove(&y(0), 1, &a(i, 0), 1, ap::vlen(0,n-1));
                xdot(y, xv, n, tx, v, verr);
                rfserrors = rfserrors||ap::fp_greater(fabs(v-b(i,0)),8*ap::machineepsilon*ap::maxreal(double(1), fabs(b(i,0))));
            }
        }
        
        //
        // Test RMatrixSolve()
        //
        unset1d(xv);
        rmatrixsolve(a, n, bv, info, rep, xv);
        if( info<=0 )
        {
            rfserrors = true;
        }
        else
        {
            for(i = 0; i <= n-1; i++)
            {
                ap::vmove(&y(0), 1, &a(i, 0), 1, ap::vlen(0,n-1));
                xdot(y, xv, n, tx, v, verr);
                rfserrors = rfserrors||ap::fp_greater(fabs(v-bv(i)),8*ap::machineepsilon*ap::maxreal(double(1), fabs(bv(i))));
            }
        }
        
        //
        // Test LS-solver on the same matrix
        //
        rmatrixsolvels(a, n, n, bv, 0.0, info, repls, xv);
        if( info<=0 )
        {
            rfserrors = true;
        }
        else
        {
            for(i = 0; i <= n-1; i++)
            {
                ap::vmove(&y(0), 1, &a(i, 0), 1, ap::vlen(0,n-1));
                xdot(y, xv, n, tx, v, verr);
                rfserrors = rfserrors||ap::fp_greater(fabs(v-bv(i)),8*ap::machineepsilon*ap::maxreal(double(1), fabs(bv(i))));
            }
        }
    }
}


/*************************************************************************
SPD test
*************************************************************************/
static void testspdsolver(int maxn,
     int maxm,
     int passcount,
     double threshold,
     bool& spderrors,
     bool& rfserrors)
{
    ap::real_2d_array a;
    ap::real_2d_array cha;
    ap::real_2d_array atmp;
    ap::integer_1d_array p;
    ap::real_2d_array xe;
    ap::real_2d_array b;
    ap::real_1d_array bv;
    int i;
    int j;
    int k;
    int n;
    int m;
    int pass;
    int taskkind;
    double mx;
    double v;
    double verr;
    bool isupper;
    int info;
    densesolverreport rep;
    densesolverlsreport repls;
    ap::real_2d_array x;
    ap::real_1d_array xv;
    ap::real_1d_array y;
    ap::real_1d_array tx;

    
    //
    // General square matrices:
    // * test general solvers
    // * test least squares solver
    //
    for(pass = 1; pass <= passcount; pass++)
    {
        for(n = 1; n <= maxn; n++)
        {
            for(m = 1; m <= maxm; m++)
            {
                
                //
                // ********************************************************
                // WELL CONDITIONED TASKS
                // ability to find correct solution is tested
                // ********************************************************
                //
                // 1. generate random well conditioned matrix A.
                // 2. generate random solution vector xe
                // 3. generate right part b=A*xe
                // 4. test different methods on original A
                //
                isupper = ap::fp_greater(ap::randomreal(),0.5);
                spdmatrixrndcond(n, double(1000), a);
                rmatrixmakeacopy(a, n, n, cha);
                if( !spdmatrixcholesky(cha, n, isupper) )
                {
                    spderrors = true;
                    return;
                }
                xe.setlength(n, m);
                for(i = 0; i <= n-1; i++)
                {
                    for(j = 0; j <= m-1; j++)
                    {
                        xe(i,j) = 2*ap::randomreal()-1;
                    }
                }
                b.setlength(n, m);
                for(i = 0; i <= n-1; i++)
                {
                    for(j = 0; j <= m-1; j++)
                    {
                        v = ap::vdotproduct(&a(i, 0), 1, &xe(0, j), xe.getstride(), ap::vlen(0,n-1));
                        b(i,j) = v;
                    }
                }
                rmatrixdrophalf(a, n, isupper);
                rmatrixdrophalf(cha, n, isupper);
                
                //
                // Test solvers
                //
                info = 0;
                unsetrep(rep);
                unset2d(x);
                spdmatrixsolvem(a, n, isupper, b, m, info, rep, x);
                spderrors = spderrors||!rmatrixchecksolutionm(xe, n, m, threshold, info, rep, x);
                info = 0;
                unsetrep(rep);
                unset1d(xv);
                bv.setlength(n);
                ap::vmove(&bv(0), 1, &b(0, 0), b.getstride(), ap::vlen(0,n-1));
                spdmatrixsolve(a, n, isupper, bv, info, rep, xv);
                spderrors = spderrors||!rmatrixchecksolution(xe, n, threshold, info, rep, xv);
                info = 0;
                unsetrep(rep);
                unset2d(x);
                spdmatrixcholeskysolvem(cha, n, isupper, b, m, info, rep, x);
                spderrors = spderrors||!rmatrixchecksolutionm(xe, n, m, threshold, info, rep, x);
                info = 0;
                unsetrep(rep);
                unset1d(xv);
                bv.setlength(n);
                ap::vmove(&bv(0), 1, &b(0, 0), b.getstride(), ap::vlen(0,n-1));
                spdmatrixcholeskysolve(cha, n, isupper, bv, info, rep, xv);
                spderrors = spderrors||!rmatrixchecksolution(xe, n, threshold, info, rep, xv);
                
                //
                // ********************************************************
                // EXACTLY SINGULAR MATRICES
                // ability to detect singularity is tested
                // ********************************************************
                //
                // 1. generate different types of singular matrices:
                //    * zero
                //    * with zero columns
                //    * with zero rows
                //    * with equal rows/columns
                // 2. generate random solution vector xe
                // 3. generate right part b=A*xe
                // 4. test different methods
                //
                for(taskkind = 0; taskkind <= 3; taskkind++)
                {
                    unset2d(a);
                    if( taskkind==0 )
                    {
                        
                        //
                        // all zeros
                        //
                        a.setlength(n, n);
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = 0; j <= n-1; j++)
                            {
                                a(i,j) = 0;
                            }
                        }
                    }
                    if( taskkind==1 )
                    {
                        
                        //
                        // there is zero column
                        //
                        a.setlength(n, n);
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = i; j <= n-1; j++)
                            {
                                a(i,j) = 2*ap::randomreal()-1;
                                a(j,i) = a(i,j);
                            }
                        }
                        k = ap::randominteger(n);
                        ap::vmul(&a(0, k), a.getstride(), ap::vlen(0,n-1), 0);
                        ap::vmul(&a(k, 0), 1, ap::vlen(0,n-1), 0);
                    }
                    if( taskkind==2 )
                    {
                        
                        //
                        // there is zero row
                        //
                        a.setlength(n, n);
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = i; j <= n-1; j++)
                            {
                                a(i,j) = 2*ap::randomreal()-1;
                                a(j,i) = a(i,j);
                            }
                        }
                        k = ap::randominteger(n);
                        ap::vmul(&a(k, 0), 1, ap::vlen(0,n-1), 0);
                        ap::vmul(&a(0, k), a.getstride(), ap::vlen(0,n-1), 0);
                    }
                    if( taskkind==3 )
                    {
                        
                        //
                        // equal columns/rows
                        //
                        if( n<2 )
                        {
                            continue;
                        }
                        a.setlength(n, n);
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = i; j <= n-1; j++)
                            {
                                a(i,j) = 2*ap::randomreal()-1;
                                a(j,i) = a(i,j);
                            }
                        }
                        k = 1+ap::randominteger(n-1);
                        ap::vmove(&a(0, 0), a.getstride(), &a(0, k), a.getstride(), ap::vlen(0,n-1));
                        ap::vmove(&a(0, 0), 1, &a(k, 0), 1, ap::vlen(0,n-1));
                    }
                    xe.setlength(n, m);
                    for(i = 0; i <= n-1; i++)
                    {
                        for(j = 0; j <= m-1; j++)
                        {
                            xe(i,j) = 2*ap::randomreal()-1;
                        }
                    }
                    b.setlength(n, m);
                    for(i = 0; i <= n-1; i++)
                    {
                        for(j = 0; j <= m-1; j++)
                        {
                            v = ap::vdotproduct(&a(i, 0), 1, &xe(0, j), xe.getstride(), ap::vlen(0,n-1));
                            b(i,j) = v;
                        }
                    }
                    rmatrixmakeacopy(a, n, n, cha);
                    rmatrixdrophalf(a, n, isupper);
                    rmatrixdrophalf(cha, n, isupper);
                    
                    //
                    // Test SPDMatrixSolveM()
                    //
                    info = 0;
                    unsetrep(rep);
                    unset2d(x);
                    spdmatrixsolvem(a, n, isupper, b, m, info, rep, x);
                    spderrors = spderrors||!rmatrixchecksingularm(n, m, info, rep, x);
                    
                    //
                    // Test SPDMatrixSolve()
                    //
                    info = 0;
                    unsetrep(rep);
                    unset2d(x);
                    bv.setlength(n);
                    ap::vmove(&bv(0), 1, &b(0, 0), b.getstride(), ap::vlen(0,n-1));
                    spdmatrixsolve(a, n, isupper, bv, info, rep, xv);
                    spderrors = spderrors||!rmatrixchecksingular(n, info, rep, xv);
                    
                    //
                    // 'equal columns/rows' are degenerate, but
                    // Cholesky matrix with equal columns/rows IS NOT degenerate,
                    // so it is not used for testing purposes.
                    //
                    if( taskkind!=3 )
                    {
                        
                        //
                        // Test SPDMatrixLUSolveM()
                        //
                        info = 0;
                        unsetrep(rep);
                        unset2d(x);
                        spdmatrixcholeskysolvem(cha, n, isupper, b, m, info, rep, x);
                        spderrors = spderrors||!rmatrixchecksingularm(n, m, info, rep, x);
                        
                        //
                        // Test SPDMatrixLUSolve()
                        //
                        info = 0;
                        unsetrep(rep);
                        unset2d(x);
                        bv.setlength(n);
                        ap::vmove(&bv(0), 1, &b(0, 0), b.getstride(), ap::vlen(0,n-1));
                        spdmatrixcholeskysolve(cha, n, isupper, bv, info, rep, xv);
                        spderrors = spderrors||!rmatrixchecksingular(n, info, rep, xv);
                    }
                }
            }
        }
    }
}


/*************************************************************************
Real test
*************************************************************************/
static void testcsolver(int maxn,
     int maxm,
     int passcount,
     double threshold,
     bool& cerrors,
     bool& rfserrors)
{
    ap::complex_2d_array a;
    ap::complex_2d_array lua;
    ap::complex_2d_array atmp;
    ap::integer_1d_array p;
    ap::complex_2d_array xe;
    ap::complex_2d_array b;
    ap::complex_1d_array bv;
    int i;
    int j;
    int k;
    int n;
    int m;
    int pass;
    int taskkind;
    double mx;
    double verr;
    ap::complex v;
    int info;
    densesolverreport rep;
    densesolverlsreport repls;
    ap::complex_2d_array x;
    ap::complex_1d_array xv;
    ap::complex_1d_array y;
    ap::real_1d_array tx;

    
    //
    // General square matrices:
    // * test general solvers
    // * test least squares solver
    //
    for(pass = 1; pass <= passcount; pass++)
    {
        for(n = 1; n <= maxn; n++)
        {
            for(m = 1; m <= maxm; m++)
            {
                
                //
                // ********************************************************
                // WELL CONDITIONED TASKS
                // ability to find correct solution is tested
                // ********************************************************
                //
                // 1. generate random well conditioned matrix A.
                // 2. generate random solution vector xe
                // 3. generate right part b=A*xe
                // 4. test different methods on original A
                //
                cmatrixrndcond(n, double(1000), a);
                cmatrixmakeacopy(a, n, n, lua);
                cmatrixlu(lua, n, n, p);
                xe.setlength(n, m);
                for(i = 0; i <= n-1; i++)
                {
                    for(j = 0; j <= m-1; j++)
                    {
                        xe(i,j).x = 2*ap::randomreal()-1;
                        xe(i,j).y = 2*ap::randomreal()-1;
                    }
                }
                b.setlength(n, m);
                for(i = 0; i <= n-1; i++)
                {
                    for(j = 0; j <= m-1; j++)
                    {
                        v = ap::vdotproduct(&a(i, 0), 1, "N", &xe(0, j), xe.getstride(), "N", ap::vlen(0,n-1));
                        b(i,j) = v;
                    }
                }
                
                //
                // Test solvers
                //
                info = 0;
                unsetrep(rep);
                cunset2d(x);
                cmatrixsolvem(a, n, b, m, ap::fp_greater(ap::randomreal(),0.5), info, rep, x);
                cerrors = cerrors||!cmatrixchecksolutionm(xe, n, m, threshold, info, rep, x);
                info = 0;
                unsetrep(rep);
                cunset1d(xv);
                bv.setlength(n);
                ap::vmove(&bv(0), 1, &b(0, 0), b.getstride(), "N", ap::vlen(0,n-1));
                cmatrixsolve(a, n, bv, info, rep, xv);
                cerrors = cerrors||!cmatrixchecksolution(xe, n, threshold, info, rep, xv);
                info = 0;
                unsetrep(rep);
                cunset2d(x);
                cmatrixlusolvem(lua, p, n, b, m, info, rep, x);
                cerrors = cerrors||!cmatrixchecksolutionm(xe, n, m, threshold, info, rep, x);
                info = 0;
                unsetrep(rep);
                cunset1d(xv);
                bv.setlength(n);
                ap::vmove(&bv(0), 1, &b(0, 0), b.getstride(), "N", ap::vlen(0,n-1));
                cmatrixlusolve(lua, p, n, bv, info, rep, xv);
                cerrors = cerrors||!cmatrixchecksolution(xe, n, threshold, info, rep, xv);
                info = 0;
                unsetrep(rep);
                cunset2d(x);
                cmatrixmixedsolvem(a, lua, p, n, b, m, info, rep, x);
                cerrors = cerrors||!cmatrixchecksolutionm(xe, n, m, threshold, info, rep, x);
                info = 0;
                unsetrep(rep);
                cunset1d(xv);
                bv.setlength(n);
                ap::vmove(&bv(0), 1, &b(0, 0), b.getstride(), "N", ap::vlen(0,n-1));
                cmatrixmixedsolve(a, lua, p, n, bv, info, rep, xv);
                cerrors = cerrors||!cmatrixchecksolution(xe, n, threshold, info, rep, xv);
                
                //
                // ********************************************************
                // EXACTLY SINGULAR MATRICES
                // ability to detect singularity is tested
                // ********************************************************
                //
                // 1. generate different types of singular matrices:
                //    * zero
                //    * with zero columns
                //    * with zero rows
                //    * with equal rows/columns
                // 2. generate random solution vector xe
                // 3. generate right part b=A*xe
                // 4. test different methods
                //
                for(taskkind = 0; taskkind <= 4; taskkind++)
                {
                    cunset2d(a);
                    if( taskkind==0 )
                    {
                        
                        //
                        // all zeros
                        //
                        a.setlength(n, n);
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = 0; j <= n-1; j++)
                            {
                                a(i,j) = 0;
                            }
                        }
                    }
                    if( taskkind==1 )
                    {
                        
                        //
                        // there is zero column
                        //
                        a.setlength(n, n);
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = 0; j <= n-1; j++)
                            {
                                a(i,j).x = 2*ap::randomreal()-1;
                                a(i,j).y = 2*ap::randomreal()-1;
                            }
                        }
                        k = ap::randominteger(n);
                        ap::vmul(&a(0, k), a.getstride(), ap::vlen(0,n-1), 0);
                    }
                    if( taskkind==2 )
                    {
                        
                        //
                        // there is zero row
                        //
                        a.setlength(n, n);
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = 0; j <= n-1; j++)
                            {
                                a(i,j).x = 2*ap::randomreal()-1;
                                a(i,j).y = 2*ap::randomreal()-1;
                            }
                        }
                        k = ap::randominteger(n);
                        ap::vmul(&a(k, 0), 1, ap::vlen(0,n-1), 0);
                    }
                    if( taskkind==3 )
                    {
                        
                        //
                        // equal columns
                        //
                        if( n<2 )
                        {
                            continue;
                        }
                        a.setlength(n, n);
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = 0; j <= n-1; j++)
                            {
                                a(i,j).x = 2*ap::randomreal()-1;
                                a(i,j).y = 2*ap::randomreal()-1;
                            }
                        }
                        k = 1+ap::randominteger(n-1);
                        ap::vmove(&a(0, 0), a.getstride(), &a(0, k), a.getstride(), "N", ap::vlen(0,n-1));
                    }
                    if( taskkind==4 )
                    {
                        
                        //
                        // equal rows
                        //
                        if( n<2 )
                        {
                            continue;
                        }
                        a.setlength(n, n);
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = 0; j <= n-1; j++)
                            {
                                a(i,j).x = 2*ap::randomreal()-1;
                                a(i,j).y = 2*ap::randomreal()-1;
                            }
                        }
                        k = 1+ap::randominteger(n-1);
                        ap::vmove(&a(0, 0), 1, &a(k, 0), 1, "N", ap::vlen(0,n-1));
                    }
                    xe.setlength(n, m);
                    for(i = 0; i <= n-1; i++)
                    {
                        for(j = 0; j <= m-1; j++)
                        {
                            xe(i,j) = 2*ap::randomreal()-1;
                        }
                    }
                    b.setlength(n, m);
                    for(i = 0; i <= n-1; i++)
                    {
                        for(j = 0; j <= m-1; j++)
                        {
                            v = ap::vdotproduct(&a(i, 0), 1, "N", &xe(0, j), xe.getstride(), "N", ap::vlen(0,n-1));
                            b(i,j) = v;
                        }
                    }
                    cmatrixmakeacopy(a, n, n, lua);
                    cmatrixlu(lua, n, n, p);
                    
                    //
                    // Test CMatrixSolveM()
                    //
                    info = 0;
                    unsetrep(rep);
                    cunset2d(x);
                    cmatrixsolvem(a, n, b, m, ap::fp_greater(ap::randomreal(),0.5), info, rep, x);
                    cerrors = cerrors||!cmatrixchecksingularm(n, m, info, rep, x);
                    
                    //
                    // Test CMatrixSolve()
                    //
                    info = 0;
                    unsetrep(rep);
                    cunset2d(x);
                    bv.setlength(n);
                    ap::vmove(&bv(0), 1, &b(0, 0), b.getstride(), "N", ap::vlen(0,n-1));
                    cmatrixsolve(a, n, bv, info, rep, xv);
                    cerrors = cerrors||!cmatrixchecksingular(n, info, rep, xv);
                    
                    //
                    // Test CMatrixLUSolveM()
                    //
                    info = 0;
                    unsetrep(rep);
                    cunset2d(x);
                    cmatrixlusolvem(lua, p, n, b, m, info, rep, x);
                    cerrors = cerrors||!cmatrixchecksingularm(n, m, info, rep, x);
                    
                    //
                    // Test CMatrixLUSolve()
                    //
                    info = 0;
                    unsetrep(rep);
                    cunset2d(x);
                    bv.setlength(n);
                    ap::vmove(&bv(0), 1, &b(0, 0), b.getstride(), "N", ap::vlen(0,n-1));
                    cmatrixlusolve(lua, p, n, bv, info, rep, xv);
                    cerrors = cerrors||!cmatrixchecksingular(n, info, rep, xv);
                    
                    //
                    // Test CMatrixMixedSolveM()
                    //
                    info = 0;
                    unsetrep(rep);
                    cunset2d(x);
                    cmatrixmixedsolvem(a, lua, p, n, b, m, info, rep, x);
                    cerrors = cerrors||!cmatrixchecksingularm(n, m, info, rep, x);
                    
                    //
                    // Test CMatrixMixedSolve()
                    //
                    info = 0;
                    unsetrep(rep);
                    cunset2d(x);
                    bv.setlength(n);
                    ap::vmove(&bv(0), 1, &b(0, 0), b.getstride(), "N", ap::vlen(0,n-1));
                    cmatrixmixedsolve(a, lua, p, n, bv, info, rep, xv);
                    cerrors = cerrors||!cmatrixchecksingular(n, info, rep, xv);
                }
            }
        }
    }
    
    //
    // test iterative improvement
    //
    for(pass = 1; pass <= passcount; pass++)
    {
        
        //
        // Test iterative improvement matrices
        //
        // A matrix/right part are constructed such that both matrix
        // and solution components magnitudes are within (-1,+1).
        // Such matrix/right part have nice properties - system can
        // be solved using iterative improvement with |A*x-b| about
        // several ulps of max(1,|b|).
        //
        n = 100;
        a.setlength(n, n);
        b.setlength(n, 1);
        bv.setlength(n);
        tx.setlength(2*n);
        xv.setlength(n);
        y.setlength(n);
        for(i = 0; i <= n-1; i++)
        {
            xv(i).x = 2*ap::randomreal()-1;
            xv(i).y = 2*ap::randomreal()-1;
        }
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                a(i,j).x = 2*ap::randomreal()-1;
                a(i,j).y = 2*ap::randomreal()-1;
            }
            ap::vmove(&y(0), 1, &a(i, 0), 1, "N", ap::vlen(0,n-1));
            xcdot(y, xv, n, tx, v, verr);
            bv(i) = v;
        }
        ap::vmove(&b(0, 0), b.getstride(), &bv(0), 1, "N", ap::vlen(0,n-1));
        
        //
        // Test CMatrixSolveM()
        //
        cunset2d(x);
        cmatrixsolvem(a, n, b, 1, true, info, rep, x);
        if( info<=0 )
        {
            rfserrors = true;
        }
        else
        {
            xv.setlength(n);
            ap::vmove(&xv(0), 1, &x(0, 0), x.getstride(), "N", ap::vlen(0,n-1));
            for(i = 0; i <= n-1; i++)
            {
                ap::vmove(&y(0), 1, &a(i, 0), 1, "N", ap::vlen(0,n-1));
                xcdot(y, xv, n, tx, v, verr);
                rfserrors = rfserrors||ap::fp_greater(ap::abscomplex(v-b(i,0)),8*ap::machineepsilon*ap::maxreal(double(1), ap::abscomplex(b(i,0))));
            }
        }
        
        //
        // Test CMatrixSolve()
        //
        cunset1d(xv);
        cmatrixsolve(a, n, bv, info, rep, xv);
        if( info<=0 )
        {
            rfserrors = true;
        }
        else
        {
            for(i = 0; i <= n-1; i++)
            {
                ap::vmove(&y(0), 1, &a(i, 0), 1, "N", ap::vlen(0,n-1));
                xcdot(y, xv, n, tx, v, verr);
                rfserrors = rfserrors||ap::fp_greater(ap::abscomplex(v-bv(i)),8*ap::machineepsilon*ap::maxreal(double(1), ap::abscomplex(bv(i))));
            }
        }
        
        //
        // TODO: Test LS-solver on the same matrix
        //
    }
}


/*************************************************************************
HPD test
*************************************************************************/
static void testhpdsolver(int maxn,
     int maxm,
     int passcount,
     double threshold,
     bool& hpderrors,
     bool& rfserrors)
{
    ap::complex_2d_array a;
    ap::complex_2d_array cha;
    ap::complex_2d_array atmp;
    ap::integer_1d_array p;
    ap::complex_2d_array xe;
    ap::complex_2d_array b;
    ap::complex_1d_array bv;
    int i;
    int j;
    int k;
    int n;
    int m;
    int pass;
    int taskkind;
    double mx;
    ap::complex v;
    bool isupper;
    int info;
    densesolverreport rep;
    densesolverlsreport repls;
    ap::complex_2d_array x;
    ap::complex_1d_array xv;
    ap::complex_1d_array y;
    ap::complex_1d_array tx;

    
    //
    // General square matrices:
    // * test general solvers
    // * test least squares solver
    //
    for(pass = 1; pass <= passcount; pass++)
    {
        for(n = 1; n <= maxn; n++)
        {
            for(m = 1; m <= maxm; m++)
            {
                
                //
                // ********************************************************
                // WELL CONDITIONED TASKS
                // ability to find correct solution is tested
                // ********************************************************
                //
                // 1. generate random well conditioned matrix A.
                // 2. generate random solution vector xe
                // 3. generate right part b=A*xe
                // 4. test different methods on original A
                //
                isupper = ap::fp_greater(ap::randomreal(),0.5);
                hpdmatrixrndcond(n, double(1000), a);
                cmatrixmakeacopy(a, n, n, cha);
                if( !hpdmatrixcholesky(cha, n, isupper) )
                {
                    hpderrors = true;
                    return;
                }
                xe.setlength(n, m);
                for(i = 0; i <= n-1; i++)
                {
                    for(j = 0; j <= m-1; j++)
                    {
                        xe(i,j).x = 2*ap::randomreal()-1;
                        xe(i,j).y = 2*ap::randomreal()-1;
                    }
                }
                b.setlength(n, m);
                for(i = 0; i <= n-1; i++)
                {
                    for(j = 0; j <= m-1; j++)
                    {
                        v = ap::vdotproduct(&a(i, 0), 1, "N", &xe(0, j), xe.getstride(), "N", ap::vlen(0,n-1));
                        b(i,j) = v;
                    }
                }
                cmatrixdrophalf(a, n, isupper);
                cmatrixdrophalf(cha, n, isupper);
                
                //
                // Test solvers
                //
                info = 0;
                unsetrep(rep);
                cunset2d(x);
                hpdmatrixsolvem(a, n, isupper, b, m, info, rep, x);
                hpderrors = hpderrors||!cmatrixchecksolutionm(xe, n, m, threshold, info, rep, x);
                info = 0;
                unsetrep(rep);
                cunset1d(xv);
                bv.setlength(n);
                ap::vmove(&bv(0), 1, &b(0, 0), b.getstride(), "N", ap::vlen(0,n-1));
                hpdmatrixsolve(a, n, isupper, bv, info, rep, xv);
                hpderrors = hpderrors||!cmatrixchecksolution(xe, n, threshold, info, rep, xv);
                info = 0;
                unsetrep(rep);
                cunset2d(x);
                hpdmatrixcholeskysolvem(cha, n, isupper, b, m, info, rep, x);
                hpderrors = hpderrors||!cmatrixchecksolutionm(xe, n, m, threshold, info, rep, x);
                info = 0;
                unsetrep(rep);
                cunset1d(xv);
                bv.setlength(n);
                ap::vmove(&bv(0), 1, &b(0, 0), b.getstride(), "N", ap::vlen(0,n-1));
                hpdmatrixcholeskysolve(cha, n, isupper, bv, info, rep, xv);
                hpderrors = hpderrors||!cmatrixchecksolution(xe, n, threshold, info, rep, xv);
                
                //
                // ********************************************************
                // EXACTLY SINGULAR MATRICES
                // ability to detect singularity is tested
                // ********************************************************
                //
                // 1. generate different types of singular matrices:
                //    * zero
                //    * with zero columns
                //    * with zero rows
                //    * with equal rows/columns
                // 2. generate random solution vector xe
                // 3. generate right part b=A*xe
                // 4. test different methods
                //
                for(taskkind = 0; taskkind <= 3; taskkind++)
                {
                    cunset2d(a);
                    if( taskkind==0 )
                    {
                        
                        //
                        // all zeros
                        //
                        a.setlength(n, n);
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = 0; j <= n-1; j++)
                            {
                                a(i,j) = 0;
                            }
                        }
                    }
                    if( taskkind==1 )
                    {
                        
                        //
                        // there is zero column
                        //
                        a.setlength(n, n);
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = i; j <= n-1; j++)
                            {
                                a(i,j).x = 2*ap::randomreal()-1;
                                a(i,j).y = 2*ap::randomreal()-1;
                                if( i==j )
                                {
                                    a(i,j).y = 0;
                                }
                                a(j,i) = a(i,j);
                            }
                        }
                        k = ap::randominteger(n);
                        ap::vmul(&a(0, k), a.getstride(), ap::vlen(0,n-1), 0);
                        ap::vmul(&a(k, 0), 1, ap::vlen(0,n-1), 0);
                    }
                    if( taskkind==2 )
                    {
                        
                        //
                        // there is zero row
                        //
                        a.setlength(n, n);
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = i; j <= n-1; j++)
                            {
                                a(i,j).x = 2*ap::randomreal()-1;
                                a(i,j).y = 2*ap::randomreal()-1;
                                if( i==j )
                                {
                                    a(i,j).y = 0;
                                }
                                a(j,i) = a(i,j);
                            }
                        }
                        k = ap::randominteger(n);
                        ap::vmul(&a(k, 0), 1, ap::vlen(0,n-1), 0);
                        ap::vmul(&a(0, k), a.getstride(), ap::vlen(0,n-1), 0);
                    }
                    if( taskkind==3 )
                    {
                        
                        //
                        // equal columns/rows
                        //
                        if( n<2 )
                        {
                            continue;
                        }
                        a.setlength(n, n);
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = i; j <= n-1; j++)
                            {
                                a(i,j).x = 2*ap::randomreal()-1;
                                a(i,j).y = 2*ap::randomreal()-1;
                                if( i==j )
                                {
                                    a(i,j).y = 0;
                                }
                                a(j,i) = a(i,j);
                            }
                        }
                        k = 1+ap::randominteger(n-1);
                        ap::vmove(&a(0, 0), a.getstride(), &a(0, k), a.getstride(), "N", ap::vlen(0,n-1));
                        ap::vmove(&a(0, 0), 1, &a(k, 0), 1, "N", ap::vlen(0,n-1));
                    }
                    xe.setlength(n, m);
                    for(i = 0; i <= n-1; i++)
                    {
                        for(j = 0; j <= m-1; j++)
                        {
                            xe(i,j) = 2*ap::randomreal()-1;
                        }
                    }
                    b.setlength(n, m);
                    for(i = 0; i <= n-1; i++)
                    {
                        for(j = 0; j <= m-1; j++)
                        {
                            v = ap::vdotproduct(&a(i, 0), 1, "N", &xe(0, j), xe.getstride(), "N", ap::vlen(0,n-1));
                            b(i,j) = v;
                        }
                    }
                    cmatrixmakeacopy(a, n, n, cha);
                    cmatrixdrophalf(a, n, isupper);
                    cmatrixdrophalf(cha, n, isupper);
                    
                    //
                    // Test SPDMatrixSolveM()
                    //
                    info = 0;
                    unsetrep(rep);
                    cunset2d(x);
                    hpdmatrixsolvem(a, n, isupper, b, m, info, rep, x);
                    hpderrors = hpderrors||!cmatrixchecksingularm(n, m, info, rep, x);
                    
                    //
                    // Test SPDMatrixSolve()
                    //
                    info = 0;
                    unsetrep(rep);
                    cunset2d(x);
                    bv.setlength(n);
                    ap::vmove(&bv(0), 1, &b(0, 0), b.getstride(), "N", ap::vlen(0,n-1));
                    hpdmatrixsolve(a, n, isupper, bv, info, rep, xv);
                    hpderrors = hpderrors||!cmatrixchecksingular(n, info, rep, xv);
                    
                    //
                    // 'equal columns/rows' are degenerate, but
                    // Cholesky matrix with equal columns/rows IS NOT degenerate,
                    // so it is not used for testing purposes.
                    //
                    if( taskkind!=3 )
                    {
                        
                        //
                        // Test SPDMatrixLUSolveM()
                        //
                        info = 0;
                        unsetrep(rep);
                        cunset2d(x);
                        hpdmatrixcholeskysolvem(cha, n, isupper, b, m, info, rep, x);
                        hpderrors = hpderrors||!cmatrixchecksingularm(n, m, info, rep, x);
                        
                        //
                        // Test SPDMatrixLUSolve()
                        //
                        info = 0;
                        unsetrep(rep);
                        cunset2d(x);
                        bv.setlength(n);
                        ap::vmove(&bv(0), 1, &b(0, 0), b.getstride(), "N", ap::vlen(0,n-1));
                        hpdmatrixcholeskysolve(cha, n, isupper, bv, info, rep, xv);
                        hpderrors = hpderrors||!cmatrixchecksingular(n, info, rep, xv);
                    }
                }
            }
        }
    }
}


/*************************************************************************
Unsets real matrix
*************************************************************************/
static void unset2d(ap::real_2d_array& x)
{

    x.setlength(1, 1);
    x(0,0) = 2*ap::randomreal()-1;
}


/*************************************************************************
Unsets real vector
*************************************************************************/
static void unset1d(ap::real_1d_array& x)
{

    x.setlength(1);
    x(0) = 2*ap::randomreal()-1;
}


/*************************************************************************
Unsets real matrix
*************************************************************************/
static void cunset2d(ap::complex_2d_array& x)
{

    x.setlength(1, 1);
    x(0,0) = 2*ap::randomreal()-1;
}


/*************************************************************************
Unsets real vector
*************************************************************************/
static void cunset1d(ap::complex_1d_array& x)
{

    x.setlength(1);
    x(0) = 2*ap::randomreal()-1;
}


/*************************************************************************
Unsets report
*************************************************************************/
static void unsetrep(densesolverreport& r)
{

    r.r1 = -1;
    r.rinf = -1;
}


/*************************************************************************
Unsets report
*************************************************************************/
static void unsetlsrep(densesolverlsreport& r)
{

    r.r2 = -1;
    r.n = -1;
    r.k = -1;
    unset2d(r.cx);
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testdensesolverunit_test_silent()
{
    bool result;

    result = testdensesolver(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testdensesolverunit_test()
{
    bool result;

    result = testdensesolver(false);
    return result;
}




