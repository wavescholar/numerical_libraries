
 
#include <stdio.h>
#include "testortfacunit.h"

static double rmatrixdiff(const ap::real_2d_array& a,
     const ap::real_2d_array& b,
     int m,
     int n);
static void rmatrixmakeacopy(const ap::real_2d_array& a,
     int m,
     int n,
     ap::real_2d_array& b);
static void cmatrixmakeacopy(const ap::complex_2d_array& a,
     int m,
     int n,
     ap::complex_2d_array& b);
static void rmatrixfillsparsea(ap::real_2d_array& a,
     int m,
     int n,
     double sparcity);
static void cmatrixfillsparsea(ap::complex_2d_array& a,
     int m,
     int n,
     double sparcity);
static void internalmatrixmatrixmultiply(const ap::real_2d_array& a,
     int ai1,
     int ai2,
     int aj1,
     int aj2,
     bool transa,
     const ap::real_2d_array& b,
     int bi1,
     int bi2,
     int bj1,
     int bj2,
     bool transb,
     ap::real_2d_array& c,
     int ci1,
     int ci2,
     int cj1,
     int cj2);
static void testrqrproblem(const ap::real_2d_array& a,
     int m,
     int n,
     double threshold,
     bool& qrerrors);
static void testcqrproblem(const ap::complex_2d_array& a,
     int m,
     int n,
     double threshold,
     bool& qrerrors);
static void testrlqproblem(const ap::real_2d_array& a,
     int m,
     int n,
     double threshold,
     bool& lqerrors);
static void testclqproblem(const ap::complex_2d_array& a,
     int m,
     int n,
     double threshold,
     bool& lqerrors);
static void testrbdproblem(const ap::real_2d_array& a,
     int m,
     int n,
     double threshold,
     bool& bderrors);
static void testrhessproblem(const ap::real_2d_array& a,
     int n,
     double threshold,
     bool& hesserrors);
static void testrtdproblem(const ap::real_2d_array& a,
     int n,
     double threshold,
     bool& tderrors);
static void testctdproblem(const ap::complex_2d_array& a,
     int n,
     double threshold,
     bool& tderrors);

/*************************************************************************
Main unittest subroutine
*************************************************************************/
bool testortfac(bool silent)
{
    bool result;
    int maxmn;
    double threshold;
    int passcount;
    int mx;
    ap::real_2d_array ra;
    ap::complex_2d_array ca;
    int m;
    int n;
    int pass;
    int i;
    int j;
    bool rqrerrors;
    bool rlqerrors;
    bool cqrerrors;
    bool clqerrors;
    bool rbderrors;
    bool rhesserrors;
    bool rtderrors;
    bool ctderrors;
    bool waserrors;

    waserrors = false;
    rqrerrors = false;
    rlqerrors = false;
    cqrerrors = false;
    clqerrors = false;
    rbderrors = false;
    rhesserrors = false;
    rtderrors = false;
    ctderrors = false;
    maxmn = 3*ablasblocksize(ra)+1;
    passcount = 1;
    threshold = 5*1000*ap::machineepsilon;
    
    //
    // Different problems
    //
    for(mx = 1; mx <= maxmn; mx++)
    {
        for(pass = 1; pass <= passcount; pass++)
        {
            
            //
            // Rectangular factorizations: QR, LQ, bidiagonal
            // Matrix types: zero, dense, sparse
            //
            n = 1+ap::randominteger(mx);
            m = 1+ap::randominteger(mx);
            if( ap::fp_greater(ap::randomreal(),0.5) )
            {
                n = mx;
            }
            else
            {
                m = mx;
            }
            ra.setlength(m, n);
            ca.setlength(m, n);
            for(i = 0; i <= m-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    ra(i,j) = 0;
                    ca(i,j) = 0;
                }
            }
            testrqrproblem(ra, m, n, threshold, rqrerrors);
            testrlqproblem(ra, m, n, threshold, rlqerrors);
            testcqrproblem(ca, m, n, threshold, cqrerrors);
            testclqproblem(ca, m, n, threshold, clqerrors);
            testrbdproblem(ra, m, n, threshold, rbderrors);
            for(i = 0; i <= m-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    ra(i,j) = 2*ap::randomreal()-1;
                    ca(i,j).x = 2*ap::randomreal()-1;
                    ca(i,j).y = 2*ap::randomreal()-1;
                }
            }
            testrqrproblem(ra, m, n, threshold, rqrerrors);
            testrlqproblem(ra, m, n, threshold, rlqerrors);
            testcqrproblem(ca, m, n, threshold, cqrerrors);
            testclqproblem(ca, m, n, threshold, clqerrors);
            testrbdproblem(ra, m, n, threshold, rbderrors);
            rmatrixfillsparsea(ra, m, n, 0.95);
            cmatrixfillsparsea(ca, m, n, 0.95);
            testrqrproblem(ra, m, n, threshold, rqrerrors);
            testrlqproblem(ra, m, n, threshold, rlqerrors);
            testcqrproblem(ca, m, n, threshold, cqrerrors);
            testclqproblem(ca, m, n, threshold, clqerrors);
            testrbdproblem(ra, m, n, threshold, rbderrors);
            
            //
            // Square factorizations: Hessenberg, tridiagonal
            // Matrix types: zero, dense, sparse
            //
            ra.setlength(mx, mx);
            ca.setlength(mx, mx);
            for(i = 0; i <= mx-1; i++)
            {
                for(j = 0; j <= mx-1; j++)
                {
                    ra(i,j) = 0;
                    ca(i,j) = 0;
                }
            }
            testrhessproblem(ra, mx, threshold, rhesserrors);
            for(i = 0; i <= mx-1; i++)
            {
                for(j = 0; j <= mx-1; j++)
                {
                    ra(i,j) = 2*ap::randomreal()-1;
                    ca(i,j).x = 2*ap::randomreal()-1;
                    ca(i,j).y = 2*ap::randomreal()-1;
                }
            }
            testrhessproblem(ra, mx, threshold, rhesserrors);
            rmatrixfillsparsea(ra, mx, mx, 0.95);
            cmatrixfillsparsea(ca, mx, mx, 0.95);
            testrhessproblem(ra, mx, threshold, rhesserrors);
            
            //
            // Symetric factorizations: tridiagonal
            // Matrix types: zero, dense, sparse
            //
            ra.setlength(mx, mx);
            ca.setlength(mx, mx);
            for(i = 0; i <= mx-1; i++)
            {
                for(j = 0; j <= mx-1; j++)
                {
                    ra(i,j) = 0;
                    ca(i,j) = 0;
                }
            }
            testrtdproblem(ra, mx, threshold, rtderrors);
            testctdproblem(ca, mx, threshold, ctderrors);
            for(i = 0; i <= mx-1; i++)
            {
                for(j = i; j <= mx-1; j++)
                {
                    ra(i,j) = 2*ap::randomreal()-1;
                    ca(i,j).x = 2*ap::randomreal()-1;
                    ca(i,j).y = 2*ap::randomreal()-1;
                    ra(j,i) = ra(i,j);
                    ca(j,i) = ap::conj(ca(i,j));
                }
            }
            for(i = 0; i <= mx-1; i++)
            {
                ca(i,i) = 2*ap::randomreal()-1;
            }
            testrtdproblem(ra, mx, threshold, rtderrors);
            testctdproblem(ca, mx, threshold, ctderrors);
            rmatrixfillsparsea(ra, mx, mx, 0.95);
            cmatrixfillsparsea(ca, mx, mx, 0.95);
            for(i = 0; i <= mx-1; i++)
            {
                for(j = i; j <= mx-1; j++)
                {
                    ra(j,i) = ra(i,j);
                    ca(j,i) = ap::conj(ca(i,j));
                }
            }
            for(i = 0; i <= mx-1; i++)
            {
                ca(i,i) = 2*ap::randomreal()-1;
            }
            testrtdproblem(ra, mx, threshold, rtderrors);
            testctdproblem(ca, mx, threshold, ctderrors);
        }
    }
    
    //
    // report
    //
    waserrors = rqrerrors||rlqerrors||cqrerrors||clqerrors||rbderrors||rhesserrors||rtderrors||ctderrors;
    if( !silent )
    {
        printf("TESTING ORTFAC UNIT\n");
        printf("RQR ERRORS:                              ");
        if( !rqrerrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("RLQ ERRORS:                              ");
        if( !rlqerrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("CQR ERRORS:                              ");
        if( !cqrerrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("CLQ ERRORS:                              ");
        if( !clqerrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("RBD ERRORS:                              ");
        if( !rbderrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("RHESS ERRORS:                            ");
        if( !rhesserrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("RTD ERRORS:                              ");
        if( !rtderrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("CTD ERRORS:                              ");
        if( !ctderrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
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
Diff
*************************************************************************/
static double rmatrixdiff(const ap::real_2d_array& a,
     const ap::real_2d_array& b,
     int m,
     int n)
{
    double result;
    int i;
    int j;

    result = 0;
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            result = ap::maxreal(result, fabs(b(i,j)-a(i,j)));
        }
    }
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
Sparse fill
*************************************************************************/
static void rmatrixfillsparsea(ap::real_2d_array& a,
     int m,
     int n,
     double sparcity)
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


/*************************************************************************
Sparse fill
*************************************************************************/
static void cmatrixfillsparsea(ap::complex_2d_array& a,
     int m,
     int n,
     double sparcity)
{
    int i;
    int j;

    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            if( ap::fp_greater_eq(ap::randomreal(),sparcity) )
            {
                a(i,j).x = 2*ap::randomreal()-1;
                a(i,j).y = 2*ap::randomreal()-1;
            }
            else
            {
                a(i,j) = 0;
            }
        }
    }
}


/*************************************************************************
Matrix multiplication
*************************************************************************/
static void internalmatrixmatrixmultiply(const ap::real_2d_array& a,
     int ai1,
     int ai2,
     int aj1,
     int aj2,
     bool transa,
     const ap::real_2d_array& b,
     int bi1,
     int bi2,
     int bj1,
     int bj2,
     bool transb,
     ap::real_2d_array& c,
     int ci1,
     int ci2,
     int cj1,
     int cj2)
{
    int arows;
    int acols;
    int brows;
    int bcols;
    int crows;
    int ccols;
    int i;
    int j;
    int k;
    int l;
    int r;
    double v;
    ap::real_1d_array work;
    double beta;
    double alpha;

    
    //
    // Pre-setup
    //
    k = ap::maxint(ai2-ai1+1, aj2-aj1+1);
    k = ap::maxint(k, bi2-bi1+1);
    k = ap::maxint(k, bj2-bj1+1);
    work.setbounds(1, k);
    beta = 0;
    alpha = 1;
    
    //
    // Setup
    //
    if( !transa )
    {
        arows = ai2-ai1+1;
        acols = aj2-aj1+1;
    }
    else
    {
        arows = aj2-aj1+1;
        acols = ai2-ai1+1;
    }
    if( !transb )
    {
        brows = bi2-bi1+1;
        bcols = bj2-bj1+1;
    }
    else
    {
        brows = bj2-bj1+1;
        bcols = bi2-bi1+1;
    }
    ap::ap_error::make_assertion(acols==brows, "MatrixMatrixMultiply: incorrect matrix sizes!");
    if( arows<=0||acols<=0||brows<=0||bcols<=0 )
    {
        return;
    }
    crows = arows;
    ccols = bcols;
    
    //
    // Test WORK
    //
    i = ap::maxint(arows, acols);
    i = ap::maxint(brows, i);
    i = ap::maxint(i, bcols);
    work(1) = 0;
    work(i) = 0;
    
    //
    // Prepare C
    //
    if( ap::fp_eq(beta,0) )
    {
        for(i = ci1; i <= ci2; i++)
        {
            for(j = cj1; j <= cj2; j++)
            {
                c(i,j) = 0;
            }
        }
    }
    else
    {
        for(i = ci1; i <= ci2; i++)
        {
            ap::vmul(&c(i, cj1), 1, ap::vlen(cj1,cj2), beta);
        }
    }
    
    //
    // A*B
    //
    if( !transa&&!transb )
    {
        for(l = ai1; l <= ai2; l++)
        {
            for(r = bi1; r <= bi2; r++)
            {
                v = alpha*a(l,aj1+r-bi1);
                k = ci1+l-ai1;
                ap::vadd(&c(k, cj1), 1, &b(r, bj1), 1, ap::vlen(cj1,cj2), v);
            }
        }
        return;
    }
    
    //
    // A*B'
    //
    if( !transa&&transb )
    {
        if( arows*acols<brows*bcols )
        {
            for(r = bi1; r <= bi2; r++)
            {
                for(l = ai1; l <= ai2; l++)
                {
                    v = ap::vdotproduct(&a(l, aj1), 1, &b(r, bj1), 1, ap::vlen(aj1,aj2));
                    c(ci1+l-ai1,cj1+r-bi1) = c(ci1+l-ai1,cj1+r-bi1)+alpha*v;
                }
            }
            return;
        }
        else
        {
            for(l = ai1; l <= ai2; l++)
            {
                for(r = bi1; r <= bi2; r++)
                {
                    v = ap::vdotproduct(&a(l, aj1), 1, &b(r, bj1), 1, ap::vlen(aj1,aj2));
                    c(ci1+l-ai1,cj1+r-bi1) = c(ci1+l-ai1,cj1+r-bi1)+alpha*v;
                }
            }
            return;
        }
    }
    
    //
    // A'*B
    //
    if( transa&&!transb )
    {
        for(l = aj1; l <= aj2; l++)
        {
            for(r = bi1; r <= bi2; r++)
            {
                v = alpha*a(ai1+r-bi1,l);
                k = ci1+l-aj1;
                ap::vadd(&c(k, cj1), 1, &b(r, bj1), 1, ap::vlen(cj1,cj2), v);
            }
        }
        return;
    }
    
    //
    // A'*B'
    //
    if( transa&&transb )
    {
        if( arows*acols<brows*bcols )
        {
            for(r = bi1; r <= bi2; r++)
            {
                for(i = 1; i <= crows; i++)
                {
                    work(i) = 0.0;
                }
                for(l = ai1; l <= ai2; l++)
                {
                    v = alpha*b(r,bj1+l-ai1);
                    k = cj1+r-bi1;
                    ap::vadd(&work(1), 1, &a(l, aj1), 1, ap::vlen(1,crows), v);
                }
                ap::vadd(&c(ci1, k), c.getstride(), &work(1), 1, ap::vlen(ci1,ci2));
            }
            return;
        }
        else
        {
            for(l = aj1; l <= aj2; l++)
            {
                k = ai2-ai1+1;
                ap::vmove(&work(1), 1, &a(ai1, l), a.getstride(), ap::vlen(1,k));
                for(r = bi1; r <= bi2; r++)
                {
                    v = ap::vdotproduct(&work(1), 1, &b(r, bj1), 1, ap::vlen(1,k));
                    c(ci1+l-aj1,cj1+r-bi1) = c(ci1+l-aj1,cj1+r-bi1)+alpha*v;
                }
            }
            return;
        }
    }
}


/*************************************************************************
Problem testing
*************************************************************************/
static void testrqrproblem(const ap::real_2d_array& a,
     int m,
     int n,
     double threshold,
     bool& qrerrors)
{
    int i;
    int j;
    int k;
    ap::real_2d_array b;
    ap::real_1d_array taub;
    ap::real_2d_array q;
    ap::real_2d_array r;
    ap::real_2d_array q2;
    double v;

    
    //
    // Test decompose-and-unpack error
    //
    rmatrixmakeacopy(a, m, n, b);
    rmatrixqr(b, m, n, taub);
    rmatrixqrunpackq(b, m, n, taub, m, q);
    rmatrixqrunpackr(b, m, n, r);
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&q(i, 0), 1, &r(0, j), r.getstride(), ap::vlen(0,m-1));
            qrerrors = qrerrors||ap::fp_greater(fabs(v-a(i,j)),threshold);
        }
    }
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= ap::minint(i, n-1)-1; j++)
        {
            qrerrors = qrerrors||ap::fp_neq(r(i,j),0);
        }
    }
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= m-1; j++)
        {
            v = ap::vdotproduct(&q(i, 0), 1, &q(j, 0), 1, ap::vlen(0,m-1));
            if( i==j )
            {
                v = v-1;
            }
            qrerrors = qrerrors||ap::fp_greater_eq(fabs(v),threshold);
        }
    }
    
    //
    // Test for other errors
    //
    k = 1+ap::randominteger(m);
    rmatrixqrunpackq(b, m, n, taub, k, q2);
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= k-1; j++)
        {
            qrerrors = qrerrors||ap::fp_greater(fabs(q2(i,j)-q(i,j)),10*ap::machineepsilon);
        }
    }
}


/*************************************************************************
Problem testing
*************************************************************************/
static void testcqrproblem(const ap::complex_2d_array& a,
     int m,
     int n,
     double threshold,
     bool& qrerrors)
{
    int i;
    int j;
    int k;
    ap::complex_2d_array b;
    ap::complex_1d_array taub;
    ap::complex_2d_array q;
    ap::complex_2d_array r;
    ap::complex_2d_array q2;
    ap::complex v;

    
    //
    // Test decompose-and-unpack error
    //
    cmatrixmakeacopy(a, m, n, b);
    cmatrixqr(b, m, n, taub);
    cmatrixqrunpackq(b, m, n, taub, m, q);
    cmatrixqrunpackr(b, m, n, r);
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&q(i, 0), 1, "N", &r(0, j), r.getstride(), "N", ap::vlen(0,m-1));
            qrerrors = qrerrors||ap::fp_greater(ap::abscomplex(v-a(i,j)),threshold);
        }
    }
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= ap::minint(i, n-1)-1; j++)
        {
            qrerrors = qrerrors||r(i,j)!=0;
        }
    }
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= m-1; j++)
        {
            v = ap::vdotproduct(&q(i, 0), 1, "N", &q(j, 0), 1, "Conj", ap::vlen(0,m-1));
            if( i==j )
            {
                v = v-1;
            }
            qrerrors = qrerrors||ap::fp_greater_eq(ap::abscomplex(v),threshold);
        }
    }
    
    //
    // Test for other errors
    //
    k = 1+ap::randominteger(m);
    cmatrixqrunpackq(b, m, n, taub, k, q2);
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= k-1; j++)
        {
            qrerrors = qrerrors||ap::fp_greater(ap::abscomplex(q2(i,j)-q(i,j)),10*ap::machineepsilon);
        }
    }
}


/*************************************************************************
Problem testing
*************************************************************************/
static void testrlqproblem(const ap::real_2d_array& a,
     int m,
     int n,
     double threshold,
     bool& lqerrors)
{
    int i;
    int j;
    int k;
    ap::real_2d_array b;
    ap::real_1d_array taub;
    ap::real_2d_array q;
    ap::real_2d_array l;
    ap::real_2d_array q2;
    double v;

    
    //
    // Test decompose-and-unpack error
    //
    rmatrixmakeacopy(a, m, n, b);
    rmatrixlq(b, m, n, taub);
    rmatrixlqunpackq(b, m, n, taub, n, q);
    rmatrixlqunpackl(b, m, n, l);
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&l(i, 0), 1, &q(0, j), q.getstride(), ap::vlen(0,n-1));
            lqerrors = lqerrors||ap::fp_greater_eq(fabs(v-a(i,j)),threshold);
        }
    }
    for(i = 0; i <= m-1; i++)
    {
        for(j = ap::minint(i, n-1)+1; j <= n-1; j++)
        {
            lqerrors = lqerrors||ap::fp_neq(l(i,j),0);
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&q(i, 0), 1, &q(j, 0), 1, ap::vlen(0,n-1));
            if( i==j )
            {
                v = v-1;
            }
            lqerrors = lqerrors||ap::fp_greater_eq(fabs(v),threshold);
        }
    }
    
    //
    // Test for other errors
    //
    k = 1+ap::randominteger(n);
    rmatrixlqunpackq(b, m, n, taub, k, q2);
    for(i = 0; i <= k-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            lqerrors = lqerrors||ap::fp_greater(fabs(q2(i,j)-q(i,j)),10*ap::machineepsilon);
        }
    }
}


/*************************************************************************
Problem testing
*************************************************************************/
static void testclqproblem(const ap::complex_2d_array& a,
     int m,
     int n,
     double threshold,
     bool& lqerrors)
{
    int i;
    int j;
    int k;
    ap::complex_2d_array b;
    ap::complex_1d_array taub;
    ap::complex_2d_array q;
    ap::complex_2d_array l;
    ap::complex_2d_array q2;
    ap::complex v;

    
    //
    // Test decompose-and-unpack error
    //
    cmatrixmakeacopy(a, m, n, b);
    cmatrixlq(b, m, n, taub);
    cmatrixlqunpackq(b, m, n, taub, n, q);
    cmatrixlqunpackl(b, m, n, l);
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&l(i, 0), 1, "N", &q(0, j), q.getstride(), "N", ap::vlen(0,n-1));
            lqerrors = lqerrors||ap::fp_greater_eq(ap::abscomplex(v-a(i,j)),threshold);
        }
    }
    for(i = 0; i <= m-1; i++)
    {
        for(j = ap::minint(i, n-1)+1; j <= n-1; j++)
        {
            lqerrors = lqerrors||l(i,j)!=0;
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&q(i, 0), 1, "N", &q(j, 0), 1, "Conj", ap::vlen(0,n-1));
            if( i==j )
            {
                v = v-1;
            }
            lqerrors = lqerrors||ap::fp_greater_eq(ap::abscomplex(v),threshold);
        }
    }
    
    //
    // Test for other errors
    //
    k = 1+ap::randominteger(n);
    cmatrixlqunpackq(b, m, n, taub, k, q2);
    for(i = 0; i <= k-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            lqerrors = lqerrors||ap::fp_greater(ap::abscomplex(q2(i,j)-q(i,j)),10*ap::machineepsilon);
        }
    }
}


/*************************************************************************
Problem testing
*************************************************************************/
static void testrbdproblem(const ap::real_2d_array& a,
     int m,
     int n,
     double threshold,
     bool& bderrors)
{
    int i;
    int j;
    int k;
    ap::real_2d_array t;
    ap::real_2d_array pt;
    ap::real_2d_array q;
    ap::real_2d_array r;
    ap::real_2d_array bd;
    ap::real_2d_array x;
    ap::real_2d_array r1;
    ap::real_2d_array r2;
    ap::real_1d_array taup;
    ap::real_1d_array tauq;
    ap::real_1d_array d;
    ap::real_1d_array e;
    bool up;
    double v;
    int mtsize;

    
    //
    // Bidiagonal decomposition error
    //
    rmatrixmakeacopy(a, m, n, t);
    rmatrixbd(t, m, n, tauq, taup);
    rmatrixbdunpackq(t, m, n, tauq, m, q);
    rmatrixbdunpackpt(t, m, n, taup, n, pt);
    rmatrixbdunpackdiagonals(t, m, n, up, d, e);
    bd.setlength(m, n);
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            bd(i,j) = 0;
        }
    }
    for(i = 0; i <= ap::minint(m, n)-1; i++)
    {
        bd(i,i) = d(i);
    }
    if( up )
    {
        for(i = 0; i <= ap::minint(m, n)-2; i++)
        {
            bd(i,i+1) = e(i);
        }
    }
    else
    {
        for(i = 0; i <= ap::minint(m, n)-2; i++)
        {
            bd(i+1,i) = e(i);
        }
    }
    r.setlength(m, n);
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&q(i, 0), 1, &bd(0, j), bd.getstride(), ap::vlen(0,m-1));
            r(i,j) = v;
        }
    }
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&r(i, 0), 1, &pt(0, j), pt.getstride(), ap::vlen(0,n-1));
            bderrors = bderrors||ap::fp_greater(fabs(v-a(i,j)),threshold);
        }
    }
    
    //
    // Orthogonality test for Q/PT
    //
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= m-1; j++)
        {
            v = ap::vdotproduct(&q(0, i), q.getstride(), &q(0, j), q.getstride(), ap::vlen(0,m-1));
            if( i==j )
            {
                bderrors = bderrors||ap::fp_greater(fabs(v-1),threshold);
            }
            else
            {
                bderrors = bderrors||ap::fp_greater(fabs(v),threshold);
            }
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&pt(i, 0), 1, &pt(j, 0), 1, ap::vlen(0,n-1));
            if( i==j )
            {
                bderrors = bderrors||ap::fp_greater(fabs(v-1),threshold);
            }
            else
            {
                bderrors = bderrors||ap::fp_greater(fabs(v),threshold);
            }
        }
    }
    
    //
    // Partial unpacking test
    //
    k = 1+ap::randominteger(m);
    rmatrixbdunpackq(t, m, n, tauq, k, r);
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= k-1; j++)
        {
            bderrors = bderrors||ap::fp_greater(fabs(r(i,j)-q(i,j)),10*ap::machineepsilon);
        }
    }
    k = 1+ap::randominteger(n);
    rmatrixbdunpackpt(t, m, n, taup, k, r);
    for(i = 0; i <= k-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            bderrors = bderrors||ap::fp_neq(r(i,j)-pt(i,j),0);
        }
    }
    
    //
    // Multiplication test
    //
    x.setbounds(0, ap::maxint(m, n)-1, 0, ap::maxint(m, n)-1);
    r.setbounds(0, ap::maxint(m, n)-1, 0, ap::maxint(m, n)-1);
    r1.setbounds(0, ap::maxint(m, n)-1, 0, ap::maxint(m, n)-1);
    r2.setbounds(0, ap::maxint(m, n)-1, 0, ap::maxint(m, n)-1);
    for(i = 0; i <= ap::maxint(m, n)-1; i++)
    {
        for(j = 0; j <= ap::maxint(m, n)-1; j++)
        {
            x(i,j) = 2*ap::randomreal()-1;
        }
    }
    mtsize = 1+ap::randominteger(ap::maxint(m, n));
    rmatrixmakeacopy(x, mtsize, m, r);
    internalmatrixmatrixmultiply(r, 0, mtsize-1, 0, m-1, false, q, 0, m-1, 0, m-1, false, r1, 0, mtsize-1, 0, m-1);
    rmatrixmakeacopy(x, mtsize, m, r2);
    rmatrixbdmultiplybyq(t, m, n, tauq, r2, mtsize, m, true, false);
    bderrors = bderrors||ap::fp_greater(rmatrixdiff(r1, r2, mtsize, m),threshold);
    rmatrixmakeacopy(x, mtsize, m, r);
    internalmatrixmatrixmultiply(r, 0, mtsize-1, 0, m-1, false, q, 0, m-1, 0, m-1, true, r1, 0, mtsize-1, 0, m-1);
    rmatrixmakeacopy(x, mtsize, m, r2);
    rmatrixbdmultiplybyq(t, m, n, tauq, r2, mtsize, m, true, true);
    bderrors = bderrors||ap::fp_greater(rmatrixdiff(r1, r2, mtsize, m),threshold);
    rmatrixmakeacopy(x, m, mtsize, r);
    internalmatrixmatrixmultiply(q, 0, m-1, 0, m-1, false, r, 0, m-1, 0, mtsize-1, false, r1, 0, m-1, 0, mtsize-1);
    rmatrixmakeacopy(x, m, mtsize, r2);
    rmatrixbdmultiplybyq(t, m, n, tauq, r2, m, mtsize, false, false);
    bderrors = bderrors||ap::fp_greater(rmatrixdiff(r1, r2, m, mtsize),threshold);
    rmatrixmakeacopy(x, m, mtsize, r);
    internalmatrixmatrixmultiply(q, 0, m-1, 0, m-1, true, r, 0, m-1, 0, mtsize-1, false, r1, 0, m-1, 0, mtsize-1);
    rmatrixmakeacopy(x, m, mtsize, r2);
    rmatrixbdmultiplybyq(t, m, n, tauq, r2, m, mtsize, false, true);
    bderrors = bderrors||ap::fp_greater(rmatrixdiff(r1, r2, m, mtsize),threshold);
    rmatrixmakeacopy(x, mtsize, n, r);
    internalmatrixmatrixmultiply(r, 0, mtsize-1, 0, n-1, false, pt, 0, n-1, 0, n-1, true, r1, 0, mtsize-1, 0, n-1);
    rmatrixmakeacopy(x, mtsize, n, r2);
    rmatrixbdmultiplybyp(t, m, n, taup, r2, mtsize, n, true, false);
    bderrors = bderrors||ap::fp_greater(rmatrixdiff(r1, r2, mtsize, n),threshold);
    rmatrixmakeacopy(x, mtsize, n, r);
    internalmatrixmatrixmultiply(r, 0, mtsize-1, 0, n-1, false, pt, 0, n-1, 0, n-1, false, r1, 0, mtsize-1, 0, n-1);
    rmatrixmakeacopy(x, mtsize, n, r2);
    rmatrixbdmultiplybyp(t, m, n, taup, r2, mtsize, n, true, true);
    bderrors = bderrors||ap::fp_greater(rmatrixdiff(r1, r2, mtsize, n),threshold);
    rmatrixmakeacopy(x, n, mtsize, r);
    internalmatrixmatrixmultiply(pt, 0, n-1, 0, n-1, true, r, 0, n-1, 0, mtsize-1, false, r1, 0, n-1, 0, mtsize-1);
    rmatrixmakeacopy(x, n, mtsize, r2);
    rmatrixbdmultiplybyp(t, m, n, taup, r2, n, mtsize, false, false);
    bderrors = bderrors||ap::fp_greater(rmatrixdiff(r1, r2, n, mtsize),threshold);
    rmatrixmakeacopy(x, n, mtsize, r);
    internalmatrixmatrixmultiply(pt, 0, n-1, 0, n-1, false, r, 0, n-1, 0, mtsize-1, false, r1, 0, n-1, 0, mtsize-1);
    rmatrixmakeacopy(x, n, mtsize, r2);
    rmatrixbdmultiplybyp(t, m, n, taup, r2, n, mtsize, false, true);
    bderrors = bderrors||ap::fp_greater(rmatrixdiff(r1, r2, n, mtsize),threshold);
}


/*************************************************************************
Problem testing
*************************************************************************/
static void testrhessproblem(const ap::real_2d_array& a,
     int n,
     double threshold,
     bool& hesserrors)
{
    ap::real_2d_array b;
    ap::real_2d_array h;
    ap::real_2d_array q;
    ap::real_2d_array t1;
    ap::real_2d_array t2;
    ap::real_1d_array tau;
    int i;
    int j;
    double v;

    rmatrixmakeacopy(a, n, n, b);
    
    //
    // Decomposition
    //
    rmatrixhessenberg(b, n, tau);
    rmatrixhessenbergunpackq(b, n, tau, q);
    rmatrixhessenbergunpackh(b, n, h);
    
    //
    // Matrix properties
    //
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&q(0, i), q.getstride(), &q(0, j), q.getstride(), ap::vlen(0,n-1));
            if( i==j )
            {
                v = v-1;
            }
            hesserrors = hesserrors||ap::fp_greater(fabs(v),threshold);
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= i-2; j++)
        {
            hesserrors = hesserrors||ap::fp_neq(h(i,j),0);
        }
    }
    
    //
    // Decomposition error
    //
    t1.setlength(n, n);
    t2.setlength(n, n);
    internalmatrixmatrixmultiply(q, 0, n-1, 0, n-1, false, h, 0, n-1, 0, n-1, false, t1, 0, n-1, 0, n-1);
    internalmatrixmatrixmultiply(t1, 0, n-1, 0, n-1, false, q, 0, n-1, 0, n-1, true, t2, 0, n-1, 0, n-1);
    hesserrors = hesserrors||ap::fp_greater(rmatrixdiff(t2, a, n, n),threshold);
}


/*************************************************************************
Tridiagonal tester
*************************************************************************/
static void testrtdproblem(const ap::real_2d_array& a,
     int n,
     double threshold,
     bool& tderrors)
{
    int i;
    int j;
    ap::real_2d_array ua;
    ap::real_2d_array la;
    ap::real_2d_array t;
    ap::real_2d_array q;
    ap::real_2d_array t2;
    ap::real_2d_array t3;
    ap::real_1d_array tau;
    ap::real_1d_array d;
    ap::real_1d_array e;
    double v;

    ua.setbounds(0, n-1, 0, n-1);
    la.setbounds(0, n-1, 0, n-1);
    t.setbounds(0, n-1, 0, n-1);
    q.setbounds(0, n-1, 0, n-1);
    t2.setbounds(0, n-1, 0, n-1);
    t3.setbounds(0, n-1, 0, n-1);
    
    //
    // fill
    //
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            ua(i,j) = 0;
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = i; j <= n-1; j++)
        {
            ua(i,j) = a(i,j);
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            la(i,j) = 0;
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= i; j++)
        {
            la(i,j) = a(i,j);
        }
    }
    
    //
    // Test 2tridiagonal: upper
    //
    smatrixtd(ua, n, true, tau, d, e);
    smatrixtdunpackq(ua, n, true, tau, q);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            t(i,j) = 0;
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        t(i,i) = d(i);
    }
    for(i = 0; i <= n-2; i++)
    {
        t(i,i+1) = e(i);
        t(i+1,i) = e(i);
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&q(0, i), q.getstride(), &a(0, j), a.getstride(), ap::vlen(0,n-1));
            t2(i,j) = v;
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&t2(i, 0), 1, &q(0, j), q.getstride(), ap::vlen(0,n-1));
            t3(i,j) = v;
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            tderrors = tderrors||ap::fp_greater(fabs(t3(i,j)-t(i,j)),threshold);
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&q(i, 0), 1, &q(j, 0), 1, ap::vlen(0,n-1));
            if( i==j )
            {
                v = v-1;
            }
            tderrors = tderrors||ap::fp_greater(fabs(v),threshold);
        }
    }
    
    //
    // Test 2tridiagonal: lower
    //
    smatrixtd(la, n, false, tau, d, e);
    smatrixtdunpackq(la, n, false, tau, q);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            t(i,j) = 0;
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        t(i,i) = d(i);
    }
    for(i = 0; i <= n-2; i++)
    {
        t(i,i+1) = e(i);
        t(i+1,i) = e(i);
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&q(0, i), q.getstride(), &a(0, j), a.getstride(), ap::vlen(0,n-1));
            t2(i,j) = v;
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&t2(i, 0), 1, &q(0, j), q.getstride(), ap::vlen(0,n-1));
            t3(i,j) = v;
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            tderrors = tderrors||ap::fp_greater(fabs(t3(i,j)-t(i,j)),threshold);
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&q(i, 0), 1, &q(j, 0), 1, ap::vlen(0,n-1));
            if( i==j )
            {
                v = v-1;
            }
            tderrors = tderrors||ap::fp_greater(fabs(v),threshold);
        }
    }
}


/*************************************************************************
Hermitian problem tester
*************************************************************************/
static void testctdproblem(const ap::complex_2d_array& a,
     int n,
     double threshold,
     bool& tderrors)
{
    int i;
    int j;
    ap::complex_2d_array ua;
    ap::complex_2d_array la;
    ap::complex_2d_array t;
    ap::complex_2d_array q;
    ap::complex_2d_array t2;
    ap::complex_2d_array t3;
    ap::complex_1d_array tau;
    ap::real_1d_array d;
    ap::real_1d_array e;
    ap::complex v;

    ua.setbounds(0, n-1, 0, n-1);
    la.setbounds(0, n-1, 0, n-1);
    t.setbounds(0, n-1, 0, n-1);
    q.setbounds(0, n-1, 0, n-1);
    t2.setbounds(0, n-1, 0, n-1);
    t3.setbounds(0, n-1, 0, n-1);
    
    //
    // fill
    //
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            ua(i,j) = 0;
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = i; j <= n-1; j++)
        {
            ua(i,j) = a(i,j);
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            la(i,j) = 0;
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= i; j++)
        {
            la(i,j) = a(i,j);
        }
    }
    
    //
    // Test 2tridiagonal: upper
    //
    hmatrixtd(ua, n, true, tau, d, e);
    hmatrixtdunpackq(ua, n, true, tau, q);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            t(i,j) = 0;
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        t(i,i) = d(i);
    }
    for(i = 0; i <= n-2; i++)
    {
        t(i,i+1) = e(i);
        t(i+1,i) = e(i);
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&q(0, i), q.getstride(), "Conj", &a(0, j), a.getstride(), "N", ap::vlen(0,n-1));
            t2(i,j) = v;
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&t2(i, 0), 1, "N", &q(0, j), q.getstride(), "N", ap::vlen(0,n-1));
            t3(i,j) = v;
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            tderrors = tderrors||ap::fp_greater(ap::abscomplex(t3(i,j)-t(i,j)),threshold);
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&q(i, 0), 1, "N", &q(j, 0), 1, "Conj", ap::vlen(0,n-1));
            if( i==j )
            {
                v = v-1;
            }
            tderrors = tderrors||ap::fp_greater(ap::abscomplex(v),threshold);
        }
    }
    
    //
    // Test 2tridiagonal: lower
    //
    hmatrixtd(la, n, false, tau, d, e);
    hmatrixtdunpackq(la, n, false, tau, q);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            t(i,j) = 0;
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        t(i,i) = d(i);
    }
    for(i = 0; i <= n-2; i++)
    {
        t(i,i+1) = e(i);
        t(i+1,i) = e(i);
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&q(0, i), q.getstride(), "Conj", &a(0, j), a.getstride(), "N", ap::vlen(0,n-1));
            t2(i,j) = v;
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&t2(i, 0), 1, "N", &q(0, j), q.getstride(), "N", ap::vlen(0,n-1));
            t3(i,j) = v;
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            tderrors = tderrors||ap::fp_greater(ap::abscomplex(t3(i,j)-t(i,j)),threshold);
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&q(i, 0), 1, "N", &q(j, 0), 1, "Conj", ap::vlen(0,n-1));
            if( i==j )
            {
                v = v-1;
            }
            tderrors = tderrors||ap::fp_greater(ap::abscomplex(v),threshold);
        }
    }
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testortfacunit_test_silent()
{
    bool result;

    result = testortfac(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testortfacunit_test()
{
    bool result;

    result = testortfac(false);
    return result;
}




