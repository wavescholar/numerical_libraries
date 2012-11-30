
 
#include <stdio.h>
#include "testevdunit.h"

static void rmatrixfillsparsea(ap::real_2d_array& a,
     int m,
     int n,
     double sparcity);
static void cmatrixfillsparsea(ap::complex_2d_array& a,
     int m,
     int n,
     double sparcity);
static void rmatrixsymmetricsplit(const ap::real_2d_array& a,
     int n,
     ap::real_2d_array& al,
     ap::real_2d_array& au);
static void cmatrixhermitiansplit(const ap::complex_2d_array& a,
     int n,
     ap::complex_2d_array& al,
     ap::complex_2d_array& au);
static void unset2d(ap::real_2d_array& a);
static void cunset2d(ap::complex_2d_array& a);
static void unset1d(ap::real_1d_array& a);
static void cunset1d(ap::complex_1d_array& a);
static double tdtestproduct(const ap::real_1d_array& d,
     const ap::real_1d_array& e,
     int n,
     const ap::real_2d_array& z,
     const ap::real_1d_array& lambda);
static double testproduct(const ap::real_2d_array& a,
     int n,
     const ap::real_2d_array& z,
     const ap::real_1d_array& lambda);
static double testort(const ap::real_2d_array& z, int n);
static double testcproduct(const ap::complex_2d_array& a,
     int n,
     const ap::complex_2d_array& z,
     const ap::real_1d_array& lambda);
static double testcort(const ap::complex_2d_array& z, int n);
static void testsevdproblem(const ap::real_2d_array& a,
     const ap::real_2d_array& al,
     const ap::real_2d_array& au,
     int n,
     double threshold,
     bool& serrors,
     int& failc,
     int& runs);
static void testhevdproblem(const ap::complex_2d_array& a,
     const ap::complex_2d_array& al,
     const ap::complex_2d_array& au,
     int n,
     double threshold,
     bool& herrors,
     int& failc,
     int& runs);
static void testsevdbiproblem(const ap::real_2d_array& afull,
     const ap::real_2d_array& al,
     const ap::real_2d_array& au,
     int n,
     bool distvals,
     double threshold,
     bool& serrors,
     int& failc,
     int& runs);
static void testhevdbiproblem(const ap::complex_2d_array& afull,
     const ap::complex_2d_array& al,
     const ap::complex_2d_array& au,
     int n,
     bool distvals,
     double threshold,
     bool& herrors,
     int& failc,
     int& runs);
static void testtdevdproblem(const ap::real_1d_array& d,
     const ap::real_1d_array& e,
     int n,
     double threshold,
     bool& tderrors,
     int& failc,
     int& runs);
static void testtdevdbiproblem(const ap::real_1d_array& d,
     const ap::real_1d_array& e,
     int n,
     bool distvals,
     double threshold,
     bool& serrors,
     int& failc,
     int& runs);
static void testnsevdproblem(const ap::real_2d_array& a,
     int n,
     double threshold,
     bool& nserrors,
     int& failc,
     int& runs);
static void testevdset(const int& n,
     const double& threshold,
     const double& bithreshold,
     int& failc,
     int& runs,
     bool& nserrors,
     bool& serrors,
     bool& herrors,
     bool& tderrors,
     bool& sbierrors,
     bool& hbierrors,
     bool& tdbierrors);

/*************************************************************************
Testing symmetric EVD subroutine
*************************************************************************/
bool testevd(bool silent)
{
    bool result;
    ap::real_2d_array ra;
    int n;
    int j;
    int failc;
    int runs;
    double failr;
    double failthreshold;
    double threshold;
    double bithreshold;
    bool waserrors;
    bool nserrors;
    bool serrors;
    bool herrors;
    bool tderrors;
    bool sbierrors;
    bool hbierrors;
    bool tdbierrors;
    bool wfailed;

    failthreshold = 0.005;
    threshold = 100000*ap::machineepsilon;
    bithreshold = 1.0E-6;
    nserrors = false;
    serrors = false;
    herrors = false;
    tderrors = false;
    sbierrors = false;
    hbierrors = false;
    tdbierrors = false;
    failc = 0;
    runs = 0;
    
    //
    // Test problems
    //
    for(n = 1; n <= ablasblocksize(ra); n++)
    {
        testevdset(n, threshold, bithreshold, failc, runs, nserrors, serrors, herrors, tderrors, sbierrors, hbierrors, tdbierrors);
    }
    for(j = 2; j <= 3; j++)
    {
        for(n = j*ablasblocksize(ra)-1; n <= j*ablasblocksize(ra)+1; n++)
        {
            testevdset(n, threshold, bithreshold, failc, runs, nserrors, serrors, herrors, tderrors, sbierrors, hbierrors, tdbierrors);
        }
    }
    
    //
    // report
    //
    wfailed = ap::fp_greater(double(failc)/double(runs),failthreshold);
    waserrors = nserrors||serrors||herrors||tderrors||sbierrors||hbierrors||tdbierrors||wfailed;
    if( !silent )
    {
        printf("TESTING EVD UNIT\n");
        printf("NS ERRORS:                               ");
        if( !nserrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("S ERRORS:                                ");
        if( !serrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("H ERRORS:                                ");
        if( !herrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("TD ERRORS:                               ");
        if( !tderrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("SBI ERRORS:                              ");
        if( !sbierrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("HBI ERRORS:                              ");
        if( !hbierrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("TDBI ERRORS:                             ");
        if( !tdbierrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("FAILURE THRESHOLD:                       ");
        if( !wfailed )
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
Copies A to AL (lower half) and AU (upper half), filling unused parts by
random garbage.
*************************************************************************/
static void rmatrixsymmetricsplit(const ap::real_2d_array& a,
     int n,
     ap::real_2d_array& al,
     ap::real_2d_array& au)
{
    int i;
    int j;

    for(i = 0; i <= n-1; i++)
    {
        for(j = i+1; j <= n-1; j++)
        {
            al(i,j) = 2*ap::randomreal()-1;
            al(j,i) = a(i,j);
            au(i,j) = a(i,j);
            au(j,i) = 2*ap::randomreal()-1;
        }
        al(i,i) = a(i,i);
        au(i,i) = a(i,i);
    }
}


/*************************************************************************
Copies A to AL (lower half) and AU (upper half), filling unused parts by
random garbage.
*************************************************************************/
static void cmatrixhermitiansplit(const ap::complex_2d_array& a,
     int n,
     ap::complex_2d_array& al,
     ap::complex_2d_array& au)
{
    int i;
    int j;

    for(i = 0; i <= n-1; i++)
    {
        for(j = i+1; j <= n-1; j++)
        {
            al(i,j) = 2*ap::randomreal()-1;
            al(j,i) = ap::conj(a(i,j));
            au(i,j) = a(i,j);
            au(j,i) = 2*ap::randomreal()-1;
        }
        al(i,i) = a(i,i);
        au(i,i) = a(i,i);
    }
}


/*************************************************************************
Unsets 2D array.
*************************************************************************/
static void unset2d(ap::real_2d_array& a)
{

    a.setbounds(0, 0, 0, 0);
    a(0,0) = 2*ap::randomreal()-1;
}


/*************************************************************************
Unsets 2D array.
*************************************************************************/
static void cunset2d(ap::complex_2d_array& a)
{

    a.setbounds(0, 0, 0, 0);
    a(0,0) = 2*ap::randomreal()-1;
}


/*************************************************************************
Unsets 1D array.
*************************************************************************/
static void unset1d(ap::real_1d_array& a)
{

    a.setbounds(0, 0);
    a(0) = 2*ap::randomreal()-1;
}


/*************************************************************************
Unsets 1D array.
*************************************************************************/
static void cunset1d(ap::complex_1d_array& a)
{

    a.setbounds(0, 0);
    a(0) = 2*ap::randomreal()-1;
}


/*************************************************************************
Tests Z*Lambda*Z' against tridiag(D,E).
Returns relative error.
*************************************************************************/
static double tdtestproduct(const ap::real_1d_array& d,
     const ap::real_1d_array& e,
     int n,
     const ap::real_2d_array& z,
     const ap::real_1d_array& lambda)
{
    double result;
    int i;
    int j;
    int k;
    double v;
    double mx;

    result = 0;
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            
            //
            // Calculate V = A[i,j], A = Z*Lambda*Z'
            //
            v = 0;
            for(k = 0; k <= n-1; k++)
            {
                v = v+z(i,k)*lambda(k)*z(j,k);
            }
            
            //
            // Compare
            //
            if( abs(i-j)==0 )
            {
                result = ap::maxreal(result, fabs(v-d(i)));
            }
            if( abs(i-j)==1 )
            {
                result = ap::maxreal(result, fabs(v-e(ap::minint(i, j))));
            }
            if( abs(i-j)>1 )
            {
                result = ap::maxreal(result, fabs(v));
            }
        }
    }
    mx = 0;
    for(i = 0; i <= n-1; i++)
    {
        mx = ap::maxreal(mx, fabs(d(i)));
    }
    for(i = 0; i <= n-2; i++)
    {
        mx = ap::maxreal(mx, fabs(e(i)));
    }
    if( ap::fp_eq(mx,0) )
    {
        mx = 1;
    }
    result = result/mx;
    return result;
}


/*************************************************************************
Tests Z*Lambda*Z' against A
Returns relative error.
*************************************************************************/
static double testproduct(const ap::real_2d_array& a,
     int n,
     const ap::real_2d_array& z,
     const ap::real_1d_array& lambda)
{
    double result;
    int i;
    int j;
    int k;
    double v;
    double mx;

    result = 0;
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            
            //
            // Calculate V = A[i,j], A = Z*Lambda*Z'
            //
            v = 0;
            for(k = 0; k <= n-1; k++)
            {
                v = v+z(i,k)*lambda(k)*z(j,k);
            }
            
            //
            // Compare
            //
            result = ap::maxreal(result, fabs(v-a(i,j)));
        }
    }
    mx = 0;
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            mx = ap::maxreal(mx, fabs(a(i,j)));
        }
    }
    if( ap::fp_eq(mx,0) )
    {
        mx = 1;
    }
    result = result/mx;
    return result;
}


/*************************************************************************
Tests Z*Z' against diag(1...1)
Returns absolute error.
*************************************************************************/
static double testort(const ap::real_2d_array& z, int n)
{
    double result;
    int i;
    int j;
    double v;

    result = 0;
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&z(0, i), z.getstride(), &z(0, j), z.getstride(), ap::vlen(0,n-1));
            if( i==j )
            {
                v = v-1;
            }
            result = ap::maxreal(result, fabs(v));
        }
    }
    return result;
}


/*************************************************************************
Tests Z*Lambda*Z' against A
Returns relative error.
*************************************************************************/
static double testcproduct(const ap::complex_2d_array& a,
     int n,
     const ap::complex_2d_array& z,
     const ap::real_1d_array& lambda)
{
    double result;
    int i;
    int j;
    int k;
    ap::complex v;
    double mx;

    result = 0;
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            
            //
            // Calculate V = A[i,j], A = Z*Lambda*Z'
            //
            v = 0;
            for(k = 0; k <= n-1; k++)
            {
                v = v+z(i,k)*lambda(k)*ap::conj(z(j,k));
            }
            
            //
            // Compare
            //
            result = ap::maxreal(result, ap::abscomplex(v-a(i,j)));
        }
    }
    mx = 0;
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            mx = ap::maxreal(mx, ap::abscomplex(a(i,j)));
        }
    }
    if( ap::fp_eq(mx,0) )
    {
        mx = 1;
    }
    result = result/mx;
    return result;
}


/*************************************************************************
Tests Z*Z' against diag(1...1)
Returns absolute error.
*************************************************************************/
static double testcort(const ap::complex_2d_array& z, int n)
{
    double result;
    int i;
    int j;
    ap::complex v;

    result = 0;
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&z(0, i), z.getstride(), "N", &z(0, j), z.getstride(), "Conj", ap::vlen(0,n-1));
            if( i==j )
            {
                v = v-1;
            }
            result = ap::maxreal(result, ap::abscomplex(v));
        }
    }
    return result;
}


/*************************************************************************
Tests SEVD problem
*************************************************************************/
static void testsevdproblem(const ap::real_2d_array& a,
     const ap::real_2d_array& al,
     const ap::real_2d_array& au,
     int n,
     double threshold,
     bool& serrors,
     int& failc,
     int& runs)
{
    ap::real_1d_array lambda;
    ap::real_1d_array lambdaref;
    ap::real_2d_array z;
    int i;
    int j;
    double v;

    
    //
    // Test simple EVD: values and full vectors, lower A
    //
    unset1d(lambdaref);
    unset2d(z);
    runs = runs+1;
    if( !smatrixevd(al, n, 1, false, lambdaref, z) )
    {
        failc = failc+1;
        return;
    }
    serrors = serrors||ap::fp_greater(testproduct(a, n, z, lambdaref),threshold);
    serrors = serrors||ap::fp_greater(testort(z, n),threshold);
    for(i = 0; i <= n-2; i++)
    {
        if( ap::fp_less(lambdaref(i+1),lambdaref(i)) )
        {
            serrors = true;
            return;
        }
    }
    
    //
    // Test simple EVD: values and full vectors, upper A
    //
    unset1d(lambda);
    unset2d(z);
    runs = runs+1;
    if( !smatrixevd(au, n, 1, true, lambda, z) )
    {
        failc = failc+1;
        return;
    }
    serrors = serrors||ap::fp_greater(testproduct(a, n, z, lambda),threshold);
    serrors = serrors||ap::fp_greater(testort(z, n),threshold);
    for(i = 0; i <= n-2; i++)
    {
        if( ap::fp_less(lambda(i+1),lambda(i)) )
        {
            serrors = true;
            return;
        }
    }
    
    //
    // Test simple EVD: values only, lower A
    //
    unset1d(lambda);
    unset2d(z);
    runs = runs+1;
    if( !smatrixevd(al, n, 0, false, lambda, z) )
    {
        failc = failc+1;
        return;
    }
    for(i = 0; i <= n-1; i++)
    {
        serrors = serrors||ap::fp_greater(fabs(lambda(i)-lambdaref(i)),threshold);
    }
    
    //
    // Test simple EVD: values only, upper A
    //
    unset1d(lambda);
    unset2d(z);
    runs = runs+1;
    if( !smatrixevd(au, n, 0, true, lambda, z) )
    {
        failc = failc+1;
        return;
    }
    for(i = 0; i <= n-1; i++)
    {
        serrors = serrors||ap::fp_greater(fabs(lambda(i)-lambdaref(i)),threshold);
    }
}


/*************************************************************************
Tests SEVD problem
*************************************************************************/
static void testhevdproblem(const ap::complex_2d_array& a,
     const ap::complex_2d_array& al,
     const ap::complex_2d_array& au,
     int n,
     double threshold,
     bool& herrors,
     int& failc,
     int& runs)
{
    ap::real_1d_array lambda;
    ap::real_1d_array lambdaref;
    ap::complex_2d_array z;
    int i;
    int j;
    ap::complex v;

    
    //
    // Test simple EVD: values and full vectors, lower A
    //
    unset1d(lambdaref);
    cunset2d(z);
    runs = runs+1;
    if( !hmatrixevd(al, n, 1, false, lambdaref, z) )
    {
        failc = failc+1;
        return;
    }
    herrors = herrors||ap::fp_greater(testcproduct(a, n, z, lambdaref),threshold);
    herrors = herrors||ap::fp_greater(testcort(z, n),threshold);
    for(i = 0; i <= n-2; i++)
    {
        if( ap::fp_less(lambdaref(i+1),lambdaref(i)) )
        {
            herrors = true;
            return;
        }
    }
    
    //
    // Test simple EVD: values and full vectors, upper A
    //
    unset1d(lambda);
    cunset2d(z);
    runs = runs+1;
    if( !hmatrixevd(au, n, 1, true, lambda, z) )
    {
        failc = failc+1;
        return;
    }
    herrors = herrors||ap::fp_greater(testcproduct(a, n, z, lambda),threshold);
    herrors = herrors||ap::fp_greater(testcort(z, n),threshold);
    for(i = 0; i <= n-2; i++)
    {
        if( ap::fp_less(lambda(i+1),lambda(i)) )
        {
            herrors = true;
            return;
        }
    }
    
    //
    // Test simple EVD: values only, lower A
    //
    unset1d(lambda);
    cunset2d(z);
    runs = runs+1;
    if( !hmatrixevd(al, n, 0, false, lambda, z) )
    {
        failc = failc+1;
        return;
    }
    for(i = 0; i <= n-1; i++)
    {
        herrors = herrors||ap::fp_greater(fabs(lambda(i)-lambdaref(i)),threshold);
    }
    
    //
    // Test simple EVD: values only, upper A
    //
    unset1d(lambda);
    cunset2d(z);
    runs = runs+1;
    if( !hmatrixevd(au, n, 0, true, lambda, z) )
    {
        failc = failc+1;
        return;
    }
    for(i = 0; i <= n-1; i++)
    {
        herrors = herrors||ap::fp_greater(fabs(lambda(i)-lambdaref(i)),threshold);
    }
}


/*************************************************************************
Tests EVD problem

DistVals    -   is True, when eigenvalues are distinct. Is False, when we
                are solving sparse task with  lots  of  zero  eigenvalues.
                In such cases some tests related to the  eigenvectors  are
                not performed.
*************************************************************************/
static void testsevdbiproblem(const ap::real_2d_array& afull,
     const ap::real_2d_array& al,
     const ap::real_2d_array& au,
     int n,
     bool distvals,
     double threshold,
     bool& serrors,
     int& failc,
     int& runs)
{
    ap::real_1d_array lambda;
    ap::real_1d_array lambdaref;
    ap::real_2d_array z;
    ap::real_2d_array zref;
    ap::real_2d_array a1;
    ap::real_2d_array a2;
    ap::real_2d_array ar;
    bool wsucc;
    int i;
    int j;
    int k;
    int m;
    int i1;
    int i2;
    double v;
    double a;
    double b;

    lambdaref.setbounds(0, n-1);
    zref.setbounds(0, n-1, 0, n-1);
    a1.setbounds(0, n-1, 0, n-1);
    a2.setbounds(0, n-1, 0, n-1);
    
    //
    // Reference EVD
    //
    runs = runs+1;
    if( !smatrixevd(afull, n, 1, true, lambdaref, zref) )
    {
        failc = failc+1;
        return;
    }
    
    //
    // Select random interval boundaries.
    // If there are non-distinct eigenvalues at the boundaries,
    // we move indexes further until values splits. It is done to
    // avoid situations where we can't get definite answer.
    //
    i1 = ap::randominteger(n);
    i2 = i1+ap::randominteger(n-i1);
    while(i1>0)
    {
        if( ap::fp_greater(fabs(lambdaref(i1-1)-lambdaref(i1)),10*threshold) )
        {
            break;
        }
        i1 = i1-1;
    }
    while(i2<n-1)
    {
        if( ap::fp_greater(fabs(lambdaref(i2+1)-lambdaref(i2)),10*threshold) )
        {
            break;
        }
        i2 = i2+1;
    }
    
    //
    // Select A, B
    //
    if( i1>0 )
    {
        a = 0.5*(lambdaref(i1)+lambdaref(i1-1));
    }
    else
    {
        a = lambdaref(0)-1;
    }
    if( i2<n-1 )
    {
        b = 0.5*(lambdaref(i2)+lambdaref(i2+1));
    }
    else
    {
        b = lambdaref(n-1)+1;
    }
    
    //
    // Test interval, no vectors, lower A
    //
    unset1d(lambda);
    unset2d(z);
    runs = runs+1;
    if( !smatrixevdr(al, n, 0, false, a, b, m, lambda, z) )
    {
        failc = failc+1;
        return;
    }
    if( m!=i2-i1+1 )
    {
        failc = failc+1;
        return;
    }
    for(k = 0; k <= m-1; k++)
    {
        serrors = serrors||ap::fp_greater(fabs(lambda(k)-lambdaref(i1+k)),threshold);
    }
    
    //
    // Test interval, no vectors, upper A
    //
    unset1d(lambda);
    unset2d(z);
    runs = runs+1;
    if( !smatrixevdr(au, n, 0, true, a, b, m, lambda, z) )
    {
        failc = failc+1;
        return;
    }
    if( m!=i2-i1+1 )
    {
        failc = failc+1;
        return;
    }
    for(k = 0; k <= m-1; k++)
    {
        serrors = serrors||ap::fp_greater(fabs(lambda(k)-lambdaref(i1+k)),threshold);
    }
    
    //
    // Test indexes, no vectors, lower A
    //
    unset1d(lambda);
    unset2d(z);
    runs = runs+1;
    if( !smatrixevdi(al, n, 0, false, i1, i2, lambda, z) )
    {
        failc = failc+1;
        return;
    }
    m = i2-i1+1;
    for(k = 0; k <= m-1; k++)
    {
        serrors = serrors||ap::fp_greater(fabs(lambda(k)-lambdaref(i1+k)),threshold);
    }
    
    //
    // Test indexes, no vectors, upper A
    //
    unset1d(lambda);
    unset2d(z);
    runs = runs+1;
    if( !smatrixevdi(au, n, 0, true, i1, i2, lambda, z) )
    {
        failc = failc+1;
        return;
    }
    m = i2-i1+1;
    for(k = 0; k <= m-1; k++)
    {
        serrors = serrors||ap::fp_greater(fabs(lambda(k)-lambdaref(i1+k)),threshold);
    }
    
    //
    // Test interval, vectors, lower A
    //
    unset1d(lambda);
    unset2d(z);
    runs = runs+1;
    if( !smatrixevdr(al, n, 1, false, a, b, m, lambda, z) )
    {
        failc = failc+1;
        return;
    }
    if( m!=i2-i1+1 )
    {
        failc = failc+1;
        return;
    }
    for(k = 0; k <= m-1; k++)
    {
        serrors = serrors||ap::fp_greater(fabs(lambda(k)-lambdaref(i1+k)),threshold);
    }
    if( distvals )
    {
        
        //
        // Distinct eigenvalues, test vectors
        //
        for(j = 0; j <= m-1; j++)
        {
            v = ap::vdotproduct(&z(0, j), z.getstride(), &zref(0, i1+j), zref.getstride(), ap::vlen(0,n-1));
            if( ap::fp_less(v,0) )
            {
                ap::vmul(&z(0, j), z.getstride(), ap::vlen(0,n-1), -1);
            }
        }
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                serrors = serrors||ap::fp_greater(fabs(z(i,j)-zref(i,i1+j)),threshold);
            }
        }
    }
    
    //
    // Test interval, vectors, upper A
    //
    unset1d(lambda);
    unset2d(z);
    runs = runs+1;
    if( !smatrixevdr(au, n, 1, true, a, b, m, lambda, z) )
    {
        failc = failc+1;
        return;
    }
    if( m!=i2-i1+1 )
    {
        failc = failc+1;
        return;
    }
    for(k = 0; k <= m-1; k++)
    {
        serrors = serrors||ap::fp_greater(fabs(lambda(k)-lambdaref(i1+k)),threshold);
    }
    if( distvals )
    {
        
        //
        // Distinct eigenvalues, test vectors
        //
        for(j = 0; j <= m-1; j++)
        {
            v = ap::vdotproduct(&z(0, j), z.getstride(), &zref(0, i1+j), zref.getstride(), ap::vlen(0,n-1));
            if( ap::fp_less(v,0) )
            {
                ap::vmul(&z(0, j), z.getstride(), ap::vlen(0,n-1), -1);
            }
        }
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                serrors = serrors||ap::fp_greater(fabs(z(i,j)-zref(i,i1+j)),threshold);
            }
        }
    }
    
    //
    // Test indexes, vectors, lower A
    //
    unset1d(lambda);
    unset2d(z);
    runs = runs+1;
    if( !smatrixevdi(al, n, 1, false, i1, i2, lambda, z) )
    {
        failc = failc+1;
        return;
    }
    m = i2-i1+1;
    for(k = 0; k <= m-1; k++)
    {
        serrors = serrors||ap::fp_greater(fabs(lambda(k)-lambdaref(i1+k)),threshold);
    }
    if( distvals )
    {
        
        //
        // Distinct eigenvalues, test vectors
        //
        for(j = 0; j <= m-1; j++)
        {
            v = ap::vdotproduct(&z(0, j), z.getstride(), &zref(0, i1+j), zref.getstride(), ap::vlen(0,n-1));
            if( ap::fp_less(v,0) )
            {
                ap::vmul(&z(0, j), z.getstride(), ap::vlen(0,n-1), -1);
            }
        }
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                serrors = serrors||ap::fp_greater(fabs(z(i,j)-zref(i,i1+j)),threshold);
            }
        }
    }
    
    //
    // Test indexes, vectors, upper A
    //
    unset1d(lambda);
    unset2d(z);
    runs = runs+1;
    if( !smatrixevdi(au, n, 1, true, i1, i2, lambda, z) )
    {
        failc = failc+1;
        return;
    }
    m = i2-i1+1;
    for(k = 0; k <= m-1; k++)
    {
        serrors = serrors||ap::fp_greater(fabs(lambda(k)-lambdaref(i1+k)),threshold);
    }
    if( distvals )
    {
        
        //
        // Distinct eigenvalues, test vectors
        //
        for(j = 0; j <= m-1; j++)
        {
            v = ap::vdotproduct(&z(0, j), z.getstride(), &zref(0, i1+j), zref.getstride(), ap::vlen(0,n-1));
            if( ap::fp_less(v,0) )
            {
                ap::vmul(&z(0, j), z.getstride(), ap::vlen(0,n-1), -1);
            }
        }
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                serrors = serrors||ap::fp_greater(fabs(z(i,j)-zref(i,i1+j)),threshold);
            }
        }
    }
}


/*************************************************************************
Tests EVD problem

DistVals    -   is True, when eigenvalues are distinct. Is False, when we
                are solving sparse task with  lots  of  zero  eigenvalues.
                In such cases some tests related to the  eigenvectors  are
                not performed.
*************************************************************************/
static void testhevdbiproblem(const ap::complex_2d_array& afull,
     const ap::complex_2d_array& al,
     const ap::complex_2d_array& au,
     int n,
     bool distvals,
     double threshold,
     bool& herrors,
     int& failc,
     int& runs)
{
    ap::real_1d_array lambda;
    ap::real_1d_array lambdaref;
    ap::complex_2d_array z;
    ap::complex_2d_array zref;
    ap::complex_2d_array a1;
    ap::complex_2d_array a2;
    ap::complex_2d_array ar;
    bool wsucc;
    int i;
    int j;
    int k;
    int m;
    int i1;
    int i2;
    ap::complex v;
    double a;
    double b;

    lambdaref.setbounds(0, n-1);
    zref.setbounds(0, n-1, 0, n-1);
    a1.setbounds(0, n-1, 0, n-1);
    a2.setbounds(0, n-1, 0, n-1);
    
    //
    // Reference EVD
    //
    runs = runs+1;
    if( !hmatrixevd(afull, n, 1, true, lambdaref, zref) )
    {
        failc = failc+1;
        return;
    }
    
    //
    // Select random interval boundaries.
    // If there are non-distinct eigenvalues at the boundaries,
    // we move indexes further until values splits. It is done to
    // avoid situations where we can't get definite answer.
    //
    i1 = ap::randominteger(n);
    i2 = i1+ap::randominteger(n-i1);
    while(i1>0)
    {
        if( ap::fp_greater(fabs(lambdaref(i1-1)-lambdaref(i1)),10*threshold) )
        {
            break;
        }
        i1 = i1-1;
    }
    while(i2<n-1)
    {
        if( ap::fp_greater(fabs(lambdaref(i2+1)-lambdaref(i2)),10*threshold) )
        {
            break;
        }
        i2 = i2+1;
    }
    
    //
    // Select A, B
    //
    if( i1>0 )
    {
        a = 0.5*(lambdaref(i1)+lambdaref(i1-1));
    }
    else
    {
        a = lambdaref(0)-1;
    }
    if( i2<n-1 )
    {
        b = 0.5*(lambdaref(i2)+lambdaref(i2+1));
    }
    else
    {
        b = lambdaref(n-1)+1;
    }
    
    //
    // Test interval, no vectors, lower A
    //
    unset1d(lambda);
    cunset2d(z);
    runs = runs+1;
    if( !hmatrixevdr(al, n, 0, false, a, b, m, lambda, z) )
    {
        failc = failc+1;
        return;
    }
    if( m!=i2-i1+1 )
    {
        failc = failc+1;
        return;
    }
    for(k = 0; k <= m-1; k++)
    {
        herrors = herrors||ap::fp_greater(fabs(lambda(k)-lambdaref(i1+k)),threshold);
    }
    
    //
    // Test interval, no vectors, upper A
    //
    unset1d(lambda);
    cunset2d(z);
    runs = runs+1;
    if( !hmatrixevdr(au, n, 0, true, a, b, m, lambda, z) )
    {
        failc = failc+1;
        return;
    }
    if( m!=i2-i1+1 )
    {
        failc = failc+1;
        return;
    }
    for(k = 0; k <= m-1; k++)
    {
        herrors = herrors||ap::fp_greater(fabs(lambda(k)-lambdaref(i1+k)),threshold);
    }
    
    //
    // Test indexes, no vectors, lower A
    //
    unset1d(lambda);
    cunset2d(z);
    runs = runs+1;
    if( !hmatrixevdi(al, n, 0, false, i1, i2, lambda, z) )
    {
        failc = failc+1;
        return;
    }
    m = i2-i1+1;
    for(k = 0; k <= m-1; k++)
    {
        herrors = herrors||ap::fp_greater(fabs(lambda(k)-lambdaref(i1+k)),threshold);
    }
    
    //
    // Test indexes, no vectors, upper A
    //
    unset1d(lambda);
    cunset2d(z);
    runs = runs+1;
    if( !hmatrixevdi(au, n, 0, true, i1, i2, lambda, z) )
    {
        failc = failc+1;
        return;
    }
    m = i2-i1+1;
    for(k = 0; k <= m-1; k++)
    {
        herrors = herrors||ap::fp_greater(fabs(lambda(k)-lambdaref(i1+k)),threshold);
    }
    
    //
    // Test interval, vectors, lower A
    //
    unset1d(lambda);
    cunset2d(z);
    runs = runs+1;
    if( !hmatrixevdr(al, n, 1, false, a, b, m, lambda, z) )
    {
        failc = failc+1;
        return;
    }
    if( m!=i2-i1+1 )
    {
        failc = failc+1;
        return;
    }
    for(k = 0; k <= m-1; k++)
    {
        herrors = herrors||ap::fp_greater(fabs(lambda(k)-lambdaref(i1+k)),threshold);
    }
    if( distvals )
    {
        
        //
        // Distinct eigenvalues, test vectors
        //
        for(j = 0; j <= m-1; j++)
        {
            v = ap::vdotproduct(&z(0, j), z.getstride(), "N", &zref(0, i1+j), zref.getstride(), "Conj", ap::vlen(0,n-1));
            v = ap::conj(v/ap::abscomplex(v));
            ap::vmul(&z(0, j), z.getstride(), ap::vlen(0,n-1), v);
        }
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                herrors = herrors||ap::fp_greater(ap::abscomplex(z(i,j)-zref(i,i1+j)),threshold);
            }
        }
    }
    
    //
    // Test interval, vectors, upper A
    //
    unset1d(lambda);
    cunset2d(z);
    runs = runs+1;
    if( !hmatrixevdr(au, n, 1, true, a, b, m, lambda, z) )
    {
        failc = failc+1;
        return;
    }
    if( m!=i2-i1+1 )
    {
        failc = failc+1;
        return;
    }
    for(k = 0; k <= m-1; k++)
    {
        herrors = herrors||ap::fp_greater(fabs(lambda(k)-lambdaref(i1+k)),threshold);
    }
    if( distvals )
    {
        
        //
        // Distinct eigenvalues, test vectors
        //
        for(j = 0; j <= m-1; j++)
        {
            v = ap::vdotproduct(&z(0, j), z.getstride(), "N", &zref(0, i1+j), zref.getstride(), "Conj", ap::vlen(0,n-1));
            v = ap::conj(v/ap::abscomplex(v));
            ap::vmul(&z(0, j), z.getstride(), ap::vlen(0,n-1), v);
        }
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                herrors = herrors||ap::fp_greater(ap::abscomplex(z(i,j)-zref(i,i1+j)),threshold);
            }
        }
    }
    
    //
    // Test indexes, vectors, lower A
    //
    unset1d(lambda);
    cunset2d(z);
    runs = runs+1;
    if( !hmatrixevdi(al, n, 1, false, i1, i2, lambda, z) )
    {
        failc = failc+1;
        return;
    }
    m = i2-i1+1;
    for(k = 0; k <= m-1; k++)
    {
        herrors = herrors||ap::fp_greater(fabs(lambda(k)-lambdaref(i1+k)),threshold);
    }
    if( distvals )
    {
        
        //
        // Distinct eigenvalues, test vectors
        //
        for(j = 0; j <= m-1; j++)
        {
            v = ap::vdotproduct(&z(0, j), z.getstride(), "N", &zref(0, i1+j), zref.getstride(), "Conj", ap::vlen(0,n-1));
            v = ap::conj(v/ap::abscomplex(v));
            ap::vmul(&z(0, j), z.getstride(), ap::vlen(0,n-1), v);
        }
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                herrors = herrors||ap::fp_greater(ap::abscomplex(z(i,j)-zref(i,i1+j)),threshold);
            }
        }
    }
    
    //
    // Test indexes, vectors, upper A
    //
    unset1d(lambda);
    cunset2d(z);
    runs = runs+1;
    if( !hmatrixevdi(au, n, 1, true, i1, i2, lambda, z) )
    {
        failc = failc+1;
        return;
    }
    m = i2-i1+1;
    for(k = 0; k <= m-1; k++)
    {
        herrors = herrors||ap::fp_greater(fabs(lambda(k)-lambdaref(i1+k)),threshold);
    }
    if( distvals )
    {
        
        //
        // Distinct eigenvalues, test vectors
        //
        for(j = 0; j <= m-1; j++)
        {
            v = ap::vdotproduct(&z(0, j), z.getstride(), "N", &zref(0, i1+j), zref.getstride(), "Conj", ap::vlen(0,n-1));
            v = ap::conj(v/ap::abscomplex(v));
            ap::vmul(&z(0, j), z.getstride(), ap::vlen(0,n-1), v);
        }
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                herrors = herrors||ap::fp_greater(ap::abscomplex(z(i,j)-zref(i,i1+j)),threshold);
            }
        }
    }
}


/*************************************************************************
Tests EVD problem
*************************************************************************/
static void testtdevdproblem(const ap::real_1d_array& d,
     const ap::real_1d_array& e,
     int n,
     double threshold,
     bool& tderrors,
     int& failc,
     int& runs)
{
    ap::real_1d_array lambda;
    ap::real_1d_array ee;
    ap::real_1d_array lambda2;
    ap::real_2d_array z;
    ap::real_2d_array zref;
    ap::real_2d_array a1;
    ap::real_2d_array a2;
    bool wsucc;
    int i;
    int j;
    double v;

    lambda.setbounds(0, n-1);
    lambda2.setbounds(0, n-1);
    zref.setbounds(0, n-1, 0, n-1);
    a1.setbounds(0, n-1, 0, n-1);
    a2.setbounds(0, n-1, 0, n-1);
    if( n>1 )
    {
        ee.setbounds(0, n-2);
    }
    
    //
    // Test simple EVD: values and full vectors
    //
    for(i = 0; i <= n-1; i++)
    {
        lambda(i) = d(i);
    }
    for(i = 0; i <= n-2; i++)
    {
        ee(i) = e(i);
    }
    runs = runs+1;
    wsucc = smatrixtdevd(lambda, ee, n, 2, z);
    if( !wsucc )
    {
        failc = failc+1;
        return;
    }
    tderrors = tderrors||ap::fp_greater(tdtestproduct(d, e, n, z, lambda),threshold);
    tderrors = tderrors||ap::fp_greater(testort(z, n),threshold);
    for(i = 0; i <= n-2; i++)
    {
        if( ap::fp_less(lambda(i+1),lambda(i)) )
        {
            tderrors = true;
            return;
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            zref(i,j) = z(i,j);
        }
    }
    
    //
    // Test values only variant
    //
    for(i = 0; i <= n-1; i++)
    {
        lambda2(i) = d(i);
    }
    for(i = 0; i <= n-2; i++)
    {
        ee(i) = e(i);
    }
    runs = runs+1;
    wsucc = smatrixtdevd(lambda2, ee, n, 0, z);
    if( !wsucc )
    {
        failc = failc+1;
        return;
    }
    for(i = 0; i <= n-1; i++)
    {
        tderrors = tderrors||ap::fp_greater(fabs(lambda2(i)-lambda(i)),threshold);
    }
    
    //
    // Test multiplication variant
    //
    for(i = 0; i <= n-1; i++)
    {
        lambda2(i) = d(i);
    }
    for(i = 0; i <= n-2; i++)
    {
        ee(i) = e(i);
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            a1(i,j) = 2*ap::randomreal()-1;
            a2(i,j) = a1(i,j);
        }
    }
    runs = runs+1;
    wsucc = smatrixtdevd(lambda2, ee, n, 1, a1);
    if( !wsucc )
    {
        failc = failc+1;
        return;
    }
    for(i = 0; i <= n-1; i++)
    {
        tderrors = tderrors||ap::fp_greater(fabs(lambda2(i)-lambda(i)),threshold);
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&a2(i, 0), 1, &zref(0, j), zref.getstride(), ap::vlen(0,n-1));
            
            //
            // next line is a bit complicated because
            // depending on algorithm used we can get either
            // z or -z as eigenvector. so we compare result
            // with both A*ZRef and -A*ZRef
            //
            tderrors = tderrors||ap::fp_greater(fabs(v-a1(i,j)),threshold)&&ap::fp_greater(fabs(v+a1(i,j)),threshold);
        }
    }
    
    //
    // Test first row variant
    //
    for(i = 0; i <= n-1; i++)
    {
        lambda2(i) = d(i);
    }
    for(i = 0; i <= n-2; i++)
    {
        ee(i) = e(i);
    }
    runs = runs+1;
    wsucc = smatrixtdevd(lambda2, ee, n, 3, z);
    if( !wsucc )
    {
        failc = failc+1;
        return;
    }
    for(i = 0; i <= n-1; i++)
    {
        tderrors = tderrors||ap::fp_greater(fabs(lambda2(i)-lambda(i)),threshold);
        
        //
        // next line is a bit complicated because
        // depending on algorithm used we can get either
        // z or -z as eigenvector. so we compare result
        // with both z and -z
        //
        tderrors = tderrors||ap::fp_greater(fabs(z(0,i)-zref(0,i)),threshold)&&ap::fp_greater(fabs(z(0,i)+zref(0,i)),threshold);
    }
}


/*************************************************************************
Tests EVD problem

DistVals    -   is True, when eigenvalues are distinct. Is False, when we
                are solving sparse task with  lots  of  zero  eigenvalues.
                In such cases some tests related to the  eigenvectors  are
                not performed.
*************************************************************************/
static void testtdevdbiproblem(const ap::real_1d_array& d,
     const ap::real_1d_array& e,
     int n,
     bool distvals,
     double threshold,
     bool& serrors,
     int& failc,
     int& runs)
{
    ap::real_1d_array lambda;
    ap::real_1d_array lambdaref;
    ap::real_2d_array z;
    ap::real_2d_array zref;
    ap::real_2d_array a1;
    ap::real_2d_array a2;
    ap::real_2d_array ar;
    bool wsucc;
    int i;
    int j;
    int k;
    int m;
    int i1;
    int i2;
    double v;
    double a;
    double b;

    lambdaref.setbounds(0, n-1);
    zref.setbounds(0, n-1, 0, n-1);
    a1.setbounds(0, n-1, 0, n-1);
    a2.setbounds(0, n-1, 0, n-1);
    
    //
    // Reference EVD
    //
    lambdaref.setlength(n);
    ap::vmove(&lambdaref(0), 1, &d(0), 1, ap::vlen(0,n-1));
    runs = runs+1;
    if( !smatrixtdevd(lambdaref, e, n, 2, zref) )
    {
        failc = failc+1;
        return;
    }
    
    //
    // Select random interval boundaries.
    // If there are non-distinct eigenvalues at the boundaries,
    // we move indexes further until values splits. It is done to
    // avoid situations where we can't get definite answer.
    //
    i1 = ap::randominteger(n);
    i2 = i1+ap::randominteger(n-i1);
    while(i1>0)
    {
        if( ap::fp_greater(fabs(lambdaref(i1-1)-lambdaref(i1)),10*threshold) )
        {
            break;
        }
        i1 = i1-1;
    }
    while(i2<n-1)
    {
        if( ap::fp_greater(fabs(lambdaref(i2+1)-lambdaref(i2)),10*threshold) )
        {
            break;
        }
        i2 = i2+1;
    }
    
    //
    // Test different combinations
    //
    
    //
    // Select A, B
    //
    if( i1>0 )
    {
        a = 0.5*(lambdaref(i1)+lambdaref(i1-1));
    }
    else
    {
        a = lambdaref(0)-1;
    }
    if( i2<n-1 )
    {
        b = 0.5*(lambdaref(i2)+lambdaref(i2+1));
    }
    else
    {
        b = lambdaref(n-1)+1;
    }
    
    //
    // Test interval, no vectors
    //
    lambda.setbounds(0, n-1);
    for(i = 0; i <= n-1; i++)
    {
        lambda(i) = d(i);
    }
    runs = runs+1;
    if( !smatrixtdevdr(lambda, e, n, 0, a, b, m, z) )
    {
        failc = failc+1;
        return;
    }
    if( m!=i2-i1+1 )
    {
        failc = failc+1;
        return;
    }
    for(k = 0; k <= m-1; k++)
    {
        serrors = serrors||ap::fp_greater(fabs(lambda(k)-lambdaref(i1+k)),threshold);
    }
    
    //
    // Test indexes, no vectors
    //
    lambda.setbounds(0, n-1);
    for(i = 0; i <= n-1; i++)
    {
        lambda(i) = d(i);
    }
    runs = runs+1;
    if( !smatrixtdevdi(lambda, e, n, 0, i1, i2, z) )
    {
        failc = failc+1;
        return;
    }
    m = i2-i1+1;
    for(k = 0; k <= m-1; k++)
    {
        serrors = serrors||ap::fp_greater(fabs(lambda(k)-lambdaref(i1+k)),threshold);
    }
    
    //
    // Test interval, transform vectors
    //
    lambda.setbounds(0, n-1);
    for(i = 0; i <= n-1; i++)
    {
        lambda(i) = d(i);
    }
    a1.setbounds(0, n-1, 0, n-1);
    a2.setbounds(0, n-1, 0, n-1);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            a1(i,j) = 2*ap::randomreal()-1;
            a2(i,j) = a1(i,j);
        }
    }
    runs = runs+1;
    if( !smatrixtdevdr(lambda, e, n, 1, a, b, m, a1) )
    {
        failc = failc+1;
        return;
    }
    if( m!=i2-i1+1 )
    {
        failc = failc+1;
        return;
    }
    for(k = 0; k <= m-1; k++)
    {
        serrors = serrors||ap::fp_greater(fabs(lambda(k)-lambdaref(i1+k)),threshold);
    }
    if( distvals )
    {
        ar.setbounds(0, n-1, 0, m-1);
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                v = ap::vdotproduct(&a2(i, 0), 1, &zref(0, i1+j), zref.getstride(), ap::vlen(0,n-1));
                ar(i,j) = v;
            }
        }
        for(j = 0; j <= m-1; j++)
        {
            v = ap::vdotproduct(&a1(0, j), a1.getstride(), &ar(0, j), ar.getstride(), ap::vlen(0,n-1));
            if( ap::fp_less(v,0) )
            {
                ap::vmul(&ar(0, j), ar.getstride(), ap::vlen(0,n-1), -1);
            }
        }
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                serrors = serrors||ap::fp_greater(fabs(a1(i,j)-ar(i,j)),threshold);
            }
        }
    }
    
    //
    // Test indexes, transform vectors
    //
    lambda.setbounds(0, n-1);
    for(i = 0; i <= n-1; i++)
    {
        lambda(i) = d(i);
    }
    a1.setbounds(0, n-1, 0, n-1);
    a2.setbounds(0, n-1, 0, n-1);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            a1(i,j) = 2*ap::randomreal()-1;
            a2(i,j) = a1(i,j);
        }
    }
    runs = runs+1;
    if( !smatrixtdevdi(lambda, e, n, 1, i1, i2, a1) )
    {
        failc = failc+1;
        return;
    }
    m = i2-i1+1;
    for(k = 0; k <= m-1; k++)
    {
        serrors = serrors||ap::fp_greater(fabs(lambda(k)-lambdaref(i1+k)),threshold);
    }
    if( distvals )
    {
        ar.setbounds(0, n-1, 0, m-1);
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                v = ap::vdotproduct(&a2(i, 0), 1, &zref(0, i1+j), zref.getstride(), ap::vlen(0,n-1));
                ar(i,j) = v;
            }
        }
        for(j = 0; j <= m-1; j++)
        {
            v = ap::vdotproduct(&a1(0, j), a1.getstride(), &ar(0, j), ar.getstride(), ap::vlen(0,n-1));
            if( ap::fp_less(v,0) )
            {
                ap::vmul(&ar(0, j), ar.getstride(), ap::vlen(0,n-1), -1);
            }
        }
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                serrors = serrors||ap::fp_greater(fabs(a1(i,j)-ar(i,j)),threshold);
            }
        }
    }
    
    //
    // Test interval, do not transform vectors
    //
    lambda.setbounds(0, n-1);
    for(i = 0; i <= n-1; i++)
    {
        lambda(i) = d(i);
    }
    z.setbounds(0, 0, 0, 0);
    runs = runs+1;
    if( !smatrixtdevdr(lambda, e, n, 2, a, b, m, z) )
    {
        failc = failc+1;
        return;
    }
    if( m!=i2-i1+1 )
    {
        failc = failc+1;
        return;
    }
    for(k = 0; k <= m-1; k++)
    {
        serrors = serrors||ap::fp_greater(fabs(lambda(k)-lambdaref(i1+k)),threshold);
    }
    if( distvals )
    {
        for(j = 0; j <= m-1; j++)
        {
            v = ap::vdotproduct(&z(0, j), z.getstride(), &zref(0, i1+j), zref.getstride(), ap::vlen(0,n-1));
            if( ap::fp_less(v,0) )
            {
                ap::vmul(&z(0, j), z.getstride(), ap::vlen(0,n-1), -1);
            }
        }
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                serrors = serrors||ap::fp_greater(fabs(z(i,j)-zref(i,i1+j)),threshold);
            }
        }
    }
    
    //
    // Test indexes, do not transform vectors
    //
    lambda.setbounds(0, n-1);
    for(i = 0; i <= n-1; i++)
    {
        lambda(i) = d(i);
    }
    z.setbounds(0, 0, 0, 0);
    runs = runs+1;
    if( !smatrixtdevdi(lambda, e, n, 2, i1, i2, z) )
    {
        failc = failc+1;
        return;
    }
    m = i2-i1+1;
    for(k = 0; k <= m-1; k++)
    {
        serrors = serrors||ap::fp_greater(fabs(lambda(k)-lambdaref(i1+k)),threshold);
    }
    if( distvals )
    {
        for(j = 0; j <= m-1; j++)
        {
            v = ap::vdotproduct(&z(0, j), z.getstride(), &zref(0, i1+j), zref.getstride(), ap::vlen(0,n-1));
            if( ap::fp_less(v,0) )
            {
                ap::vmul(&z(0, j), z.getstride(), ap::vlen(0,n-1), -1);
            }
        }
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                serrors = serrors||ap::fp_greater(fabs(z(i,j)-zref(i,i1+j)),threshold);
            }
        }
    }
}


/*************************************************************************
Non-symmetric problem
*************************************************************************/
static void testnsevdproblem(const ap::real_2d_array& a,
     int n,
     double threshold,
     bool& nserrors,
     int& failc,
     int& runs)
{
    double mx;
    int i;
    int j;
    int k;
    int vjob;
    bool needl;
    bool needr;
    ap::real_1d_array wr0;
    ap::real_1d_array wi0;
    ap::real_1d_array wr1;
    ap::real_1d_array wi1;
    ap::real_1d_array wr0s;
    ap::real_1d_array wi0s;
    ap::real_1d_array wr1s;
    ap::real_1d_array wi1s;
    ap::real_2d_array vl;
    ap::real_2d_array vr;
    ap::real_1d_array vec1r;
    ap::real_1d_array vec1i;
    ap::real_1d_array vec2r;
    ap::real_1d_array vec2i;
    ap::real_1d_array vec3r;
    ap::real_1d_array vec3i;
    double curwr;
    double curwi;
    double vt;
    double tmp;

    vec1r.setbounds(0, n-1);
    vec2r.setbounds(0, n-1);
    vec3r.setbounds(0, n-1);
    vec1i.setbounds(0, n-1);
    vec2i.setbounds(0, n-1);
    vec3i.setbounds(0, n-1);
    wr0s.setbounds(0, n-1);
    wr1s.setbounds(0, n-1);
    wi0s.setbounds(0, n-1);
    wi1s.setbounds(0, n-1);
    mx = 0;
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            if( ap::fp_greater(fabs(a(i,j)),mx) )
            {
                mx = fabs(a(i,j));
            }
        }
    }
    if( ap::fp_eq(mx,0) )
    {
        mx = 1;
    }
    
    //
    // Load values-only
    //
    runs = runs+1;
    if( !rmatrixevd(a, n, 0, wr0, wi0, vl, vr) )
    {
        failc = failc+1;
        return;
    }
    
    //
    // Test different jobs
    //
    for(vjob = 1; vjob <= 3; vjob++)
    {
        needr = vjob==1||vjob==3;
        needl = vjob==2||vjob==3;
        runs = runs+1;
        if( !rmatrixevd(a, n, vjob, wr1, wi1, vl, vr) )
        {
            failc = failc+1;
            return;
        }
        
        //
        // Test values:
        // 1. sort by real part
        // 2. test
        //
        ap::vmove(&wr0s(0), 1, &wr0(0), 1, ap::vlen(0,n-1));
        ap::vmove(&wi0s(0), 1, &wi0(0), 1, ap::vlen(0,n-1));
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= n-2-i; j++)
            {
                if( ap::fp_greater(wr0s(j),wr0s(j+1)) )
                {
                    tmp = wr0s(j);
                    wr0s(j) = wr0s(j+1);
                    wr0s(j+1) = tmp;
                    tmp = wi0s(j);
                    wi0s(j) = wi0s(j+1);
                    wi0s(j+1) = tmp;
                }
            }
        }
        ap::vmove(&wr1s(0), 1, &wr1(0), 1, ap::vlen(0,n-1));
        ap::vmove(&wi1s(0), 1, &wi1(0), 1, ap::vlen(0,n-1));
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= n-2-i; j++)
            {
                if( ap::fp_greater(wr1s(j),wr1s(j+1)) )
                {
                    tmp = wr1s(j);
                    wr1s(j) = wr1s(j+1);
                    wr1s(j+1) = tmp;
                    tmp = wi1s(j);
                    wi1s(j) = wi1s(j+1);
                    wi1s(j+1) = tmp;
                }
            }
        }
        for(i = 0; i <= n-1; i++)
        {
            nserrors = nserrors||ap::fp_greater(fabs(wr0s(i)-wr1s(i)),threshold);
            nserrors = nserrors||ap::fp_greater(fabs(wi0s(i)-wi1s(i)),threshold);
        }
        
        //
        // Test right vectors
        //
        if( needr )
        {
            k = 0;
            while(k<=n-1)
            {
                if( ap::fp_eq(wi1(k),0) )
                {
                    ap::vmove(&vec1r(0), 1, &vr(0, k), vr.getstride(), ap::vlen(0,n-1));
                    for(i = 0; i <= n-1; i++)
                    {
                        vec1i(i) = 0;
                    }
                    curwr = wr1(k);
                    curwi = 0;
                }
                if( ap::fp_greater(wi1(k),0) )
                {
                    ap::vmove(&vec1r(0), 1, &vr(0, k), vr.getstride(), ap::vlen(0,n-1));
                    ap::vmove(&vec1i(0), 1, &vr(0, k+1), vr.getstride(), ap::vlen(0,n-1));
                    curwr = wr1(k);
                    curwi = wi1(k);
                }
                if( ap::fp_less(wi1(k),0) )
                {
                    ap::vmove(&vec1r(0), 1, &vr(0, k-1), vr.getstride(), ap::vlen(0,n-1));
                    ap::vmoveneg(&vec1i(0), 1, &vr(0, k), vr.getstride(), ap::vlen(0,n-1));
                    curwr = wr1(k);
                    curwi = wi1(k);
                }
                for(i = 0; i <= n-1; i++)
                {
                    vt = ap::vdotproduct(&a(i, 0), 1, &vec1r(0), 1, ap::vlen(0,n-1));
                    vec2r(i) = vt;
                    vt = ap::vdotproduct(&a(i, 0), 1, &vec1i(0), 1, ap::vlen(0,n-1));
                    vec2i(i) = vt;
                }
                ap::vmove(&vec3r(0), 1, &vec1r(0), 1, ap::vlen(0,n-1), curwr);
                ap::vsub(&vec3r(0), 1, &vec1i(0), 1, ap::vlen(0,n-1), curwi);
                ap::vmove(&vec3i(0), 1, &vec1r(0), 1, ap::vlen(0,n-1), curwi);
                ap::vadd(&vec3i(0), 1, &vec1i(0), 1, ap::vlen(0,n-1), curwr);
                for(i = 0; i <= n-1; i++)
                {
                    nserrors = nserrors||ap::fp_greater(fabs(vec2r(i)-vec3r(i)),threshold);
                    nserrors = nserrors||ap::fp_greater(fabs(vec2i(i)-vec3i(i)),threshold);
                }
                k = k+1;
            }
        }
        
        //
        // Test left vectors
        //
        if( needl )
        {
            k = 0;
            while(k<=n-1)
            {
                if( ap::fp_eq(wi1(k),0) )
                {
                    ap::vmove(&vec1r(0), 1, &vl(0, k), vl.getstride(), ap::vlen(0,n-1));
                    for(i = 0; i <= n-1; i++)
                    {
                        vec1i(i) = 0;
                    }
                    curwr = wr1(k);
                    curwi = 0;
                }
                if( ap::fp_greater(wi1(k),0) )
                {
                    ap::vmove(&vec1r(0), 1, &vl(0, k), vl.getstride(), ap::vlen(0,n-1));
                    ap::vmove(&vec1i(0), 1, &vl(0, k+1), vl.getstride(), ap::vlen(0,n-1));
                    curwr = wr1(k);
                    curwi = wi1(k);
                }
                if( ap::fp_less(wi1(k),0) )
                {
                    ap::vmove(&vec1r(0), 1, &vl(0, k-1), vl.getstride(), ap::vlen(0,n-1));
                    ap::vmoveneg(&vec1i(0), 1, &vl(0, k), vl.getstride(), ap::vlen(0,n-1));
                    curwr = wr1(k);
                    curwi = wi1(k);
                }
                for(j = 0; j <= n-1; j++)
                {
                    vt = ap::vdotproduct(&vec1r(0), 1, &a(0, j), a.getstride(), ap::vlen(0,n-1));
                    vec2r(j) = vt;
                    vt = ap::vdotproduct(&vec1i(0), 1, &a(0, j), a.getstride(), ap::vlen(0,n-1));
                    vec2i(j) = -vt;
                }
                ap::vmove(&vec3r(0), 1, &vec1r(0), 1, ap::vlen(0,n-1), curwr);
                ap::vadd(&vec3r(0), 1, &vec1i(0), 1, ap::vlen(0,n-1), curwi);
                ap::vmove(&vec3i(0), 1, &vec1r(0), 1, ap::vlen(0,n-1), curwi);
                ap::vsub(&vec3i(0), 1, &vec1i(0), 1, ap::vlen(0,n-1), curwr);
                for(i = 0; i <= n-1; i++)
                {
                    nserrors = nserrors||ap::fp_greater(fabs(vec2r(i)-vec3r(i)),threshold);
                    nserrors = nserrors||ap::fp_greater(fabs(vec2i(i)-vec3i(i)),threshold);
                }
                k = k+1;
            }
        }
    }
}


/*************************************************************************
Testing EVD subroutines for one N

NOTES:
* BIThreshold is a threshold for bisection-and-inverse-iteration subroutines.
  special threshold is needed because these subroutines may have much more
  larger error than QR-based algorithms.
*************************************************************************/
static void testevdset(const int& n,
     const double& threshold,
     const double& bithreshold,
     int& failc,
     int& runs,
     bool& nserrors,
     bool& serrors,
     bool& herrors,
     bool& tderrors,
     bool& sbierrors,
     bool& hbierrors,
     bool& tdbierrors)
{
    ap::real_2d_array ra;
    ap::real_2d_array ral;
    ap::real_2d_array rau;
    ap::complex_2d_array ca;
    ap::complex_2d_array cal;
    ap::complex_2d_array cau;
    ap::real_1d_array d;
    ap::real_1d_array e;
    int pass;
    int i;
    int j;
    int mkind;

    
    //
    // Test symmetric problems
    //
    
    //
    // Test symmetric problem: zero, random, sparse matrices.
    //
    ra.setlength(n, n);
    ral.setlength(n, n);
    rau.setlength(n, n);
    ca.setlength(n, n);
    cal.setlength(n, n);
    cau.setlength(n, n);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            ra(i,j) = 0;
            ca(i,j) = 0;
        }
    }
    rmatrixsymmetricsplit(ra, n, ral, rau);
    cmatrixhermitiansplit(ca, n, cal, cau);
    testsevdproblem(ra, ral, rau, n, threshold, serrors, failc, runs);
    testhevdproblem(ca, cal, cau, n, threshold, herrors, failc, runs);
    testsevdbiproblem(ra, ral, rau, n, false, bithreshold, sbierrors, failc, runs);
    testhevdbiproblem(ca, cal, cau, n, false, bithreshold, hbierrors, failc, runs);
    for(i = 0; i <= n-1; i++)
    {
        for(j = i+1; j <= n-1; j++)
        {
            ra(i,j) = 2*ap::randomreal()-1;
            ca(i,j).x = 2*ap::randomreal()-1;
            ca(i,j).y = 2*ap::randomreal()-1;
            ra(j,i) = ra(i,j);
            ca(j,i) = ap::conj(ca(i,j));
        }
        ra(i,i) = 2*ap::randomreal()-1;
        ca(i,i) = 2*ap::randomreal()-1;
    }
    rmatrixsymmetricsplit(ra, n, ral, rau);
    cmatrixhermitiansplit(ca, n, cal, cau);
    testsevdproblem(ra, ral, rau, n, threshold, serrors, failc, runs);
    testhevdproblem(ca, cal, cau, n, threshold, herrors, failc, runs);
    testsevdbiproblem(ra, ral, rau, n, true, bithreshold, sbierrors, failc, runs);
    testhevdbiproblem(ca, cal, cau, n, true, bithreshold, hbierrors, failc, runs);
    rmatrixfillsparsea(ra, n, n, 0.995);
    cmatrixfillsparsea(ca, n, n, 0.995);
    for(i = 0; i <= n-1; i++)
    {
        for(j = i+1; j <= n-1; j++)
        {
            ra(j,i) = ra(i,j);
            ca(j,i) = ap::conj(ca(i,j));
        }
        ca(i,i).y = 0;
    }
    rmatrixsymmetricsplit(ra, n, ral, rau);
    cmatrixhermitiansplit(ca, n, cal, cau);
    testsevdproblem(ra, ral, rau, n, threshold, serrors, failc, runs);
    testhevdproblem(ca, cal, cau, n, threshold, herrors, failc, runs);
    testsevdbiproblem(ra, ral, rau, n, false, bithreshold, sbierrors, failc, runs);
    testhevdbiproblem(ca, cal, cau, n, false, bithreshold, hbierrors, failc, runs);
    
    //
    // testing tridiagonal problems
    //
    for(mkind = 0; mkind <= 4; mkind++)
    {
        d.setlength(n);
        if( n>1 )
        {
            e.setlength(n-1);
        }
        if( mkind==0 )
        {
            
            //
            // Zero matrix
            //
            for(i = 0; i <= n-1; i++)
            {
                d(i) = 0;
            }
            for(i = 0; i <= n-2; i++)
            {
                e(i) = 0;
            }
        }
        if( mkind==1 )
        {
            
            //
            // Diagonal matrix
            //
            for(i = 0; i <= n-1; i++)
            {
                d(i) = 2*ap::randomreal()-1;
            }
            for(i = 0; i <= n-2; i++)
            {
                e(i) = 0;
            }
        }
        if( mkind==2 )
        {
            
            //
            // Off-diagonal matrix
            //
            for(i = 0; i <= n-1; i++)
            {
                d(i) = 0;
            }
            for(i = 0; i <= n-2; i++)
            {
                e(i) = 2*ap::randomreal()-1;
            }
        }
        if( mkind==3 )
        {
            
            //
            // Dense matrix with blocks
            //
            for(i = 0; i <= n-1; i++)
            {
                d(i) = 2*ap::randomreal()-1;
            }
            for(i = 0; i <= n-2; i++)
            {
                e(i) = 2*ap::randomreal()-1;
            }
            j = 1;
            i = 2;
            while(j<=n-2)
            {
                e(j) = 0;
                j = j+i;
                i = i+1;
            }
        }
        if( mkind==4 )
        {
            
            //
            // dense matrix
            //
            for(i = 0; i <= n-1; i++)
            {
                d(i) = 2*ap::randomreal()-1;
            }
            for(i = 0; i <= n-2; i++)
            {
                e(i) = 2*ap::randomreal()-1;
            }
        }
        testtdevdproblem(d, e, n, threshold, tderrors, failc, runs);
        testtdevdbiproblem(d, e, n, mkind==1||mkind==2||mkind==4, bithreshold, tdbierrors, failc, runs);
    }
    
    //
    // Test non-symmetric problems
    //
    
    //
    // Test non-symmetric problems: zero, random, sparse matrices.
    //
    ra.setlength(n, n);
    ca.setlength(n, n);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            ra(i,j) = 0;
            ca(i,j) = 0;
        }
    }
    testnsevdproblem(ra, n, threshold, nserrors, failc, runs);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            ra(i,j) = 2*ap::randomreal()-1;
            ca(i,j).x = 2*ap::randomreal()-1;
            ca(i,j).y = 2*ap::randomreal()-1;
        }
    }
    testnsevdproblem(ra, n, threshold, nserrors, failc, runs);
    rmatrixfillsparsea(ra, n, n, 0.995);
    cmatrixfillsparsea(ca, n, n, 0.995);
    testnsevdproblem(ra, n, threshold, nserrors, failc, runs);
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testevdunit_test_silent()
{
    bool result;

    result = testevd(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testevdunit_test()
{
    bool result;

    result = testevd(false);
    return result;
}




