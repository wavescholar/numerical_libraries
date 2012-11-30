
 
#include <stdio.h>
#include "testrcondunit.h"

static const double threshold50 = 0.25;
static const double threshold90 = 0.10;

static void rmatrixmakeacopy(const ap::real_2d_array& a,
     int m,
     int n,
     ap::real_2d_array& b);
static void rmatrixdrophalf(ap::real_2d_array& a, int n, bool droplower);
static void cmatrixdrophalf(ap::complex_2d_array& a, int n, bool droplower);
static void rmatrixgenzero(ap::real_2d_array& a0, int n);
static bool rmatrixinvmattr(ap::real_2d_array& a,
     int n,
     bool isupper,
     bool isunittriangular);
static bool rmatrixinvmatlu(ap::real_2d_array& a,
     const ap::integer_1d_array& pivots,
     int n);
static bool rmatrixinvmat(ap::real_2d_array& a, int n);
static void rmatrixrefrcond(const ap::real_2d_array& a,
     int n,
     double& rc1,
     double& rcinf);
static void cmatrixmakeacopy(const ap::complex_2d_array& a,
     int m,
     int n,
     ap::complex_2d_array& b);
static void cmatrixgenzero(ap::complex_2d_array& a0, int n);
static bool cmatrixinvmattr(ap::complex_2d_array& a,
     int n,
     bool isupper,
     bool isunittriangular);
static bool cmatrixinvmatlu(ap::complex_2d_array& a,
     const ap::integer_1d_array& pivots,
     int n);
static bool cmatrixinvmat(ap::complex_2d_array& a, int n);
static void cmatrixrefrcond(const ap::complex_2d_array& a,
     int n,
     double& rc1,
     double& rcinf);
static bool testrmatrixtrrcond(int maxn, int passcount);
static bool testcmatrixtrrcond(int maxn, int passcount);
static bool testrmatrixrcond(int maxn, int passcount);
static bool testspdmatrixrcond(int maxn, int passcount);
static bool testcmatrixrcond(int maxn, int passcount);
static bool testhpdmatrixrcond(int maxn, int passcount);

bool testrcond(bool silent)
{
    bool result;
    int maxn;
    int passcount;
    bool waserrors;
    bool rtrerr;
    bool ctrerr;
    bool rerr;
    bool cerr;
    bool spderr;
    bool hpderr;

    maxn = 10;
    passcount = 100;
    
    //
    // report
    //
    rtrerr = !testrmatrixtrrcond(maxn, passcount);
    ctrerr = !testcmatrixtrrcond(maxn, passcount);
    rerr = !testrmatrixrcond(maxn, passcount);
    cerr = !testcmatrixrcond(maxn, passcount);
    spderr = !testspdmatrixrcond(maxn, passcount);
    hpderr = !testhpdmatrixrcond(maxn, passcount);
    waserrors = rtrerr||ctrerr||rerr||cerr||spderr||hpderr;
    if( !silent )
    {
        printf("TESTING RCOND\n");
        printf("REAL TRIANGULAR:                         ");
        if( !rtrerr )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("COMPLEX TRIANGULAR:                      ");
        if( !ctrerr )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("REAL:                                    ");
        if( !rerr )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("SPD:                                     ");
        if( !spderr )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("HPD:                                     ");
        if( !hpderr )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("COMPLEX:                                 ");
        if( !cerr )
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
Generate matrix with given condition number C (2-norm)
*************************************************************************/
static void rmatrixgenzero(ap::real_2d_array& a0, int n)
{
    int i;
    int j;

    a0.setlength(n, n);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            a0(i,j) = 0;
        }
    }
}


/*************************************************************************
triangular inverse
*************************************************************************/
static bool rmatrixinvmattr(ap::real_2d_array& a,
     int n,
     bool isupper,
     bool isunittriangular)
{
    bool result;
    bool nounit;
    int i;
    int j;
    double v;
    double ajj;
    ap::real_1d_array t;

    result = true;
    t.setbounds(0, n-1);
    
    //
    // Test the input parameters.
    //
    nounit = !isunittriangular;
    if( isupper )
    {
        
        //
        // Compute inverse of upper triangular matrix.
        //
        for(j = 0; j <= n-1; j++)
        {
            if( nounit )
            {
                if( ap::fp_eq(a(j,j),0) )
                {
                    result = false;
                    return result;
                }
                a(j,j) = 1/a(j,j);
                ajj = -a(j,j);
            }
            else
            {
                ajj = -1;
            }
            
            //
            // Compute elements 1:j-1 of j-th column.
            //
            if( j>0 )
            {
                ap::vmove(&t(0), 1, &a(0, j), a.getstride(), ap::vlen(0,j-1));
                for(i = 0; i <= j-1; i++)
                {
                    if( i<j-1 )
                    {
                        v = ap::vdotproduct(&a(i, i+1), 1, &t(i+1), 1, ap::vlen(i+1,j-1));
                    }
                    else
                    {
                        v = 0;
                    }
                    if( nounit )
                    {
                        a(i,j) = v+a(i,i)*t(i);
                    }
                    else
                    {
                        a(i,j) = v+t(i);
                    }
                }
                ap::vmul(&a(0, j), a.getstride(), ap::vlen(0,j-1), ajj);
            }
        }
    }
    else
    {
        
        //
        // Compute inverse of lower triangular matrix.
        //
        for(j = n-1; j >= 0; j--)
        {
            if( nounit )
            {
                if( ap::fp_eq(a(j,j),0) )
                {
                    result = false;
                    return result;
                }
                a(j,j) = 1/a(j,j);
                ajj = -a(j,j);
            }
            else
            {
                ajj = -1;
            }
            if( j<n-1 )
            {
                
                //
                // Compute elements j+1:n of j-th column.
                //
                ap::vmove(&t(j+1), 1, &a(j+1, j), a.getstride(), ap::vlen(j+1,n-1));
                for(i = j+1; i <= n-1; i++)
                {
                    if( i>j+1 )
                    {
                        v = ap::vdotproduct(&a(i, j+1), 1, &t(j+1), 1, ap::vlen(j+1,i-1));
                    }
                    else
                    {
                        v = 0;
                    }
                    if( nounit )
                    {
                        a(i,j) = v+a(i,i)*t(i);
                    }
                    else
                    {
                        a(i,j) = v+t(i);
                    }
                }
                ap::vmul(&a(j+1, j), a.getstride(), ap::vlen(j+1,n-1), ajj);
            }
        }
    }
    return result;
}


/*************************************************************************
LU inverse
*************************************************************************/
static bool rmatrixinvmatlu(ap::real_2d_array& a,
     const ap::integer_1d_array& pivots,
     int n)
{
    bool result;
    ap::real_1d_array work;
    int i;
    int iws;
    int j;
    int jb;
    int jj;
    int jp;
    double v;

    result = true;
    
    //
    // Quick return if possible
    //
    if( n==0 )
    {
        return result;
    }
    work.setbounds(0, n-1);
    
    //
    // Form inv(U)
    //
    if( !rmatrixinvmattr(a, n, true, false) )
    {
        result = false;
        return result;
    }
    
    //
    // Solve the equation inv(A)*L = inv(U) for inv(A).
    //
    for(j = n-1; j >= 0; j--)
    {
        
        //
        // Copy current column of L to WORK and replace with zeros.
        //
        for(i = j+1; i <= n-1; i++)
        {
            work(i) = a(i,j);
            a(i,j) = 0;
        }
        
        //
        // Compute current column of inv(A).
        //
        if( j<n-1 )
        {
            for(i = 0; i <= n-1; i++)
            {
                v = ap::vdotproduct(&a(i, j+1), 1, &work(j+1), 1, ap::vlen(j+1,n-1));
                a(i,j) = a(i,j)-v;
            }
        }
    }
    
    //
    // Apply column interchanges.
    //
    for(j = n-2; j >= 0; j--)
    {
        jp = pivots(j);
        if( jp!=j )
        {
            ap::vmove(&work(0), 1, &a(0, j), a.getstride(), ap::vlen(0,n-1));
            ap::vmove(&a(0, j), a.getstride(), &a(0, jp), a.getstride(), ap::vlen(0,n-1));
            ap::vmove(&a(0, jp), a.getstride(), &work(0), 1, ap::vlen(0,n-1));
        }
    }
    return result;
}


/*************************************************************************
Matrix inverse
*************************************************************************/
static bool rmatrixinvmat(ap::real_2d_array& a, int n)
{
    bool result;
    ap::integer_1d_array pivots;

    rmatrixlu(a, n, n, pivots);
    result = rmatrixinvmatlu(a, pivots, n);
    return result;
}


/*************************************************************************
reference RCond
*************************************************************************/
static void rmatrixrefrcond(const ap::real_2d_array& a,
     int n,
     double& rc1,
     double& rcinf)
{
    ap::real_2d_array inva;
    double nrm1a;
    double nrminfa;
    double nrm1inva;
    double nrminfinva;
    double v;
    int k;
    int i;

    
    //
    // inv A
    //
    rmatrixmakeacopy(a, n, n, inva);
    if( !rmatrixinvmat(inva, n) )
    {
        rc1 = 0;
        rcinf = 0;
        return;
    }
    
    //
    // norm A
    //
    nrm1a = 0;
    nrminfa = 0;
    for(k = 0; k <= n-1; k++)
    {
        v = 0;
        for(i = 0; i <= n-1; i++)
        {
            v = v+fabs(a(i,k));
        }
        nrm1a = ap::maxreal(nrm1a, v);
        v = 0;
        for(i = 0; i <= n-1; i++)
        {
            v = v+fabs(a(k,i));
        }
        nrminfa = ap::maxreal(nrminfa, v);
    }
    
    //
    // norm inv A
    //
    nrm1inva = 0;
    nrminfinva = 0;
    for(k = 0; k <= n-1; k++)
    {
        v = 0;
        for(i = 0; i <= n-1; i++)
        {
            v = v+fabs(inva(i,k));
        }
        nrm1inva = ap::maxreal(nrm1inva, v);
        v = 0;
        for(i = 0; i <= n-1; i++)
        {
            v = v+fabs(inva(k,i));
        }
        nrminfinva = ap::maxreal(nrminfinva, v);
    }
    
    //
    // result
    //
    rc1 = nrm1inva*nrm1a;
    rcinf = nrminfinva*nrminfa;
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
Generate matrix with given condition number C (2-norm)
*************************************************************************/
static void cmatrixgenzero(ap::complex_2d_array& a0, int n)
{
    int i;
    int j;

    a0.setlength(n, n);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            a0(i,j) = 0;
        }
    }
}


/*************************************************************************
triangular inverse
*************************************************************************/
static bool cmatrixinvmattr(ap::complex_2d_array& a,
     int n,
     bool isupper,
     bool isunittriangular)
{
    bool result;
    bool nounit;
    int i;
    int j;
    ap::complex v;
    ap::complex ajj;
    ap::complex_1d_array t;

    result = true;
    t.setbounds(0, n-1);
    
    //
    // Test the input parameters.
    //
    nounit = !isunittriangular;
    if( isupper )
    {
        
        //
        // Compute inverse of upper triangular matrix.
        //
        for(j = 0; j <= n-1; j++)
        {
            if( nounit )
            {
                if( a(j,j)==0 )
                {
                    result = false;
                    return result;
                }
                a(j,j) = 1/a(j,j);
                ajj = -a(j,j);
            }
            else
            {
                ajj = -1;
            }
            
            //
            // Compute elements 1:j-1 of j-th column.
            //
            if( j>0 )
            {
                ap::vmove(&t(0), 1, &a(0, j), a.getstride(), "N", ap::vlen(0,j-1));
                for(i = 0; i <= j-1; i++)
                {
                    if( i<j-1 )
                    {
                        v = ap::vdotproduct(&a(i, i+1), 1, "N", &t(i+1), 1, "N", ap::vlen(i+1,j-1));
                    }
                    else
                    {
                        v = 0;
                    }
                    if( nounit )
                    {
                        a(i,j) = v+a(i,i)*t(i);
                    }
                    else
                    {
                        a(i,j) = v+t(i);
                    }
                }
                ap::vmul(&a(0, j), a.getstride(), ap::vlen(0,j-1), ajj);
            }
        }
    }
    else
    {
        
        //
        // Compute inverse of lower triangular matrix.
        //
        for(j = n-1; j >= 0; j--)
        {
            if( nounit )
            {
                if( a(j,j)==0 )
                {
                    result = false;
                    return result;
                }
                a(j,j) = 1/a(j,j);
                ajj = -a(j,j);
            }
            else
            {
                ajj = -1;
            }
            if( j<n-1 )
            {
                
                //
                // Compute elements j+1:n of j-th column.
                //
                ap::vmove(&t(j+1), 1, &a(j+1, j), a.getstride(), "N", ap::vlen(j+1,n-1));
                for(i = j+1; i <= n-1; i++)
                {
                    if( i>j+1 )
                    {
                        v = ap::vdotproduct(&a(i, j+1), 1, "N", &t(j+1), 1, "N", ap::vlen(j+1,i-1));
                    }
                    else
                    {
                        v = 0;
                    }
                    if( nounit )
                    {
                        a(i,j) = v+a(i,i)*t(i);
                    }
                    else
                    {
                        a(i,j) = v+t(i);
                    }
                }
                ap::vmul(&a(j+1, j), a.getstride(), ap::vlen(j+1,n-1), ajj);
            }
        }
    }
    return result;
}


/*************************************************************************
LU inverse
*************************************************************************/
static bool cmatrixinvmatlu(ap::complex_2d_array& a,
     const ap::integer_1d_array& pivots,
     int n)
{
    bool result;
    ap::complex_1d_array work;
    int i;
    int iws;
    int j;
    int jb;
    int jj;
    int jp;
    ap::complex v;

    result = true;
    
    //
    // Quick return if possible
    //
    if( n==0 )
    {
        return result;
    }
    work.setbounds(0, n-1);
    
    //
    // Form inv(U)
    //
    if( !cmatrixinvmattr(a, n, true, false) )
    {
        result = false;
        return result;
    }
    
    //
    // Solve the equation inv(A)*L = inv(U) for inv(A).
    //
    for(j = n-1; j >= 0; j--)
    {
        
        //
        // Copy current column of L to WORK and replace with zeros.
        //
        for(i = j+1; i <= n-1; i++)
        {
            work(i) = a(i,j);
            a(i,j) = 0;
        }
        
        //
        // Compute current column of inv(A).
        //
        if( j<n-1 )
        {
            for(i = 0; i <= n-1; i++)
            {
                v = ap::vdotproduct(&a(i, j+1), 1, "N", &work(j+1), 1, "N", ap::vlen(j+1,n-1));
                a(i,j) = a(i,j)-v;
            }
        }
    }
    
    //
    // Apply column interchanges.
    //
    for(j = n-2; j >= 0; j--)
    {
        jp = pivots(j);
        if( jp!=j )
        {
            ap::vmove(&work(0), 1, &a(0, j), a.getstride(), "N", ap::vlen(0,n-1));
            ap::vmove(&a(0, j), a.getstride(), &a(0, jp), a.getstride(), "N", ap::vlen(0,n-1));
            ap::vmove(&a(0, jp), a.getstride(), &work(0), 1, "N", ap::vlen(0,n-1));
        }
    }
    return result;
}


/*************************************************************************
Matrix inverse
*************************************************************************/
static bool cmatrixinvmat(ap::complex_2d_array& a, int n)
{
    bool result;
    ap::integer_1d_array pivots;

    cmatrixlu(a, n, n, pivots);
    result = cmatrixinvmatlu(a, pivots, n);
    return result;
}


/*************************************************************************
reference RCond
*************************************************************************/
static void cmatrixrefrcond(const ap::complex_2d_array& a,
     int n,
     double& rc1,
     double& rcinf)
{
    ap::complex_2d_array inva;
    double nrm1a;
    double nrminfa;
    double nrm1inva;
    double nrminfinva;
    double v;
    int k;
    int i;

    
    //
    // inv A
    //
    cmatrixmakeacopy(a, n, n, inva);
    if( !cmatrixinvmat(inva, n) )
    {
        rc1 = 0;
        rcinf = 0;
        return;
    }
    
    //
    // norm A
    //
    nrm1a = 0;
    nrminfa = 0;
    for(k = 0; k <= n-1; k++)
    {
        v = 0;
        for(i = 0; i <= n-1; i++)
        {
            v = v+ap::abscomplex(a(i,k));
        }
        nrm1a = ap::maxreal(nrm1a, v);
        v = 0;
        for(i = 0; i <= n-1; i++)
        {
            v = v+ap::abscomplex(a(k,i));
        }
        nrminfa = ap::maxreal(nrminfa, v);
    }
    
    //
    // norm inv A
    //
    nrm1inva = 0;
    nrminfinva = 0;
    for(k = 0; k <= n-1; k++)
    {
        v = 0;
        for(i = 0; i <= n-1; i++)
        {
            v = v+ap::abscomplex(inva(i,k));
        }
        nrm1inva = ap::maxreal(nrm1inva, v);
        v = 0;
        for(i = 0; i <= n-1; i++)
        {
            v = v+ap::abscomplex(inva(k,i));
        }
        nrminfinva = ap::maxreal(nrminfinva, v);
    }
    
    //
    // result
    //
    rc1 = nrm1inva*nrm1a;
    rcinf = nrminfinva*nrminfa;
}


/*************************************************************************
Returns True for successful test, False - for failed test
*************************************************************************/
static bool testrmatrixtrrcond(int maxn, int passcount)
{
    bool result;
    ap::real_2d_array a;
    ap::real_2d_array ea;
    ap::integer_1d_array p;
    int n;
    int i;
    int j;
    int j1;
    int j2;
    int pass;
    bool err50;
    bool err90;
    bool errspec;
    bool errless;
    double erc1;
    double ercinf;
    ap::real_1d_array q50;
    ap::real_1d_array q90;
    double v;
    bool isupper;
    bool isunit;

    err50 = false;
    err90 = false;
    errless = false;
    errspec = false;
    q50.setlength(2);
    q90.setlength(2);
    for(n = 1; n <= maxn; n++)
    {
        
        //
        // special test for zero matrix
        //
        rmatrixgenzero(a, n);
        errspec = errspec||ap::fp_neq(rmatrixtrrcond1(a, n, ap::fp_greater(ap::randomreal(),0.5), false),0);
        errspec = errspec||ap::fp_neq(rmatrixtrrcondinf(a, n, ap::fp_greater(ap::randomreal(),0.5), false),0);
        
        //
        // general test
        //
        a.setlength(n, n);
        for(i = 0; i <= 1; i++)
        {
            q50(i) = 0;
            q90(i) = 0;
        }
        for(pass = 1; pass <= passcount; pass++)
        {
            isupper = ap::fp_greater(ap::randomreal(),0.5);
            isunit = ap::fp_greater(ap::randomreal(),0.5);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    a(i,j) = ap::randomreal()-0.5;
                }
            }
            for(i = 0; i <= n-1; i++)
            {
                a(i,i) = 1+ap::randomreal();
            }
            rmatrixmakeacopy(a, n, n, ea);
            for(i = 0; i <= n-1; i++)
            {
                if( isupper )
                {
                    j1 = 0;
                    j2 = i-1;
                }
                else
                {
                    j1 = i+1;
                    j2 = n-1;
                }
                for(j = j1; j <= j2; j++)
                {
                    ea(i,j) = 0;
                }
                if( isunit )
                {
                    ea(i,i) = 1;
                }
            }
            rmatrixrefrcond(ea, n, erc1, ercinf);
            
            //
            // 1-norm
            //
            v = 1/rmatrixtrrcond1(a, n, isupper, isunit);
            if( ap::fp_greater_eq(v,threshold50*erc1) )
            {
                q50(0) = q50(0)+double(1)/double(passcount);
            }
            if( ap::fp_greater_eq(v,threshold90*erc1) )
            {
                q90(0) = q90(0)+double(1)/double(passcount);
            }
            errless = errless||ap::fp_greater(v,erc1*1.001);
            
            //
            // Inf-norm
            //
            v = 1/rmatrixtrrcondinf(a, n, isupper, isunit);
            if( ap::fp_greater_eq(v,threshold50*ercinf) )
            {
                q50(1) = q50(1)+double(1)/double(passcount);
            }
            if( ap::fp_greater_eq(v,threshold90*ercinf) )
            {
                q90(1) = q90(1)+double(1)/double(passcount);
            }
            errless = errless||ap::fp_greater(v,ercinf*1.001);
        }
        for(i = 0; i <= 1; i++)
        {
            err50 = err50||ap::fp_less(q50(i),0.50);
            err90 = err90||ap::fp_less(q90(i),0.90);
        }
        
        //
        // degenerate matrix test
        //
        if( n>=3 )
        {
            a.setlength(n, n);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    a(i,j) = 0.0;
                }
            }
            a(0,0) = 1;
            a(n-1,n-1) = 1;
            errspec = errspec||ap::fp_neq(rmatrixtrrcond1(a, n, ap::fp_greater(ap::randomreal(),0.5), false),0);
            errspec = errspec||ap::fp_neq(rmatrixtrrcondinf(a, n, ap::fp_greater(ap::randomreal(),0.5), false),0);
        }
        
        //
        // near-degenerate matrix test
        //
        if( n>=2 )
        {
            a.setlength(n, n);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    a(i,j) = 0.0;
                }
            }
            for(i = 0; i <= n-1; i++)
            {
                a(i,i) = 1;
            }
            i = ap::randominteger(n);
            a(i,i) = 0.1*ap::maxrealnumber;
            errspec = errspec||ap::fp_neq(rmatrixtrrcond1(a, n, ap::fp_greater(ap::randomreal(),0.5), false),0);
            errspec = errspec||ap::fp_neq(rmatrixtrrcondinf(a, n, ap::fp_greater(ap::randomreal(),0.5), false),0);
        }
    }
    
    //
    // report
    //
    result = !(err50||err90||errless||errspec);
    return result;
}


/*************************************************************************
Returns True for successful test, False - for failed test
*************************************************************************/
static bool testcmatrixtrrcond(int maxn, int passcount)
{
    bool result;
    ap::complex_2d_array a;
    ap::complex_2d_array ea;
    ap::integer_1d_array p;
    int n;
    int i;
    int j;
    int j1;
    int j2;
    int pass;
    bool err50;
    bool err90;
    bool errspec;
    bool errless;
    double erc1;
    double ercinf;
    ap::real_1d_array q50;
    ap::real_1d_array q90;
    double v;
    bool isupper;
    bool isunit;

    err50 = false;
    err90 = false;
    errless = false;
    errspec = false;
    q50.setlength(2);
    q90.setlength(2);
    for(n = 1; n <= maxn; n++)
    {
        
        //
        // special test for zero matrix
        //
        cmatrixgenzero(a, n);
        errspec = errspec||ap::fp_neq(cmatrixtrrcond1(a, n, ap::fp_greater(ap::randomreal(),0.5), false),0);
        errspec = errspec||ap::fp_neq(cmatrixtrrcondinf(a, n, ap::fp_greater(ap::randomreal(),0.5), false),0);
        
        //
        // general test
        //
        a.setlength(n, n);
        for(i = 0; i <= 1; i++)
        {
            q50(i) = 0;
            q90(i) = 0;
        }
        for(pass = 1; pass <= passcount; pass++)
        {
            isupper = ap::fp_greater(ap::randomreal(),0.5);
            isunit = ap::fp_greater(ap::randomreal(),0.5);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    a(i,j).x = ap::randomreal()-0.5;
                    a(i,j).y = ap::randomreal()-0.5;
                }
            }
            for(i = 0; i <= n-1; i++)
            {
                a(i,i).x = 1+ap::randomreal();
                a(i,i).y = 1+ap::randomreal();
            }
            cmatrixmakeacopy(a, n, n, ea);
            for(i = 0; i <= n-1; i++)
            {
                if( isupper )
                {
                    j1 = 0;
                    j2 = i-1;
                }
                else
                {
                    j1 = i+1;
                    j2 = n-1;
                }
                for(j = j1; j <= j2; j++)
                {
                    ea(i,j) = 0;
                }
                if( isunit )
                {
                    ea(i,i) = 1;
                }
            }
            cmatrixrefrcond(ea, n, erc1, ercinf);
            
            //
            // 1-norm
            //
            v = 1/cmatrixtrrcond1(a, n, isupper, isunit);
            if( ap::fp_greater_eq(v,threshold50*erc1) )
            {
                q50(0) = q50(0)+double(1)/double(passcount);
            }
            if( ap::fp_greater_eq(v,threshold90*erc1) )
            {
                q90(0) = q90(0)+double(1)/double(passcount);
            }
            errless = errless||ap::fp_greater(v,erc1*1.001);
            
            //
            // Inf-norm
            //
            v = 1/cmatrixtrrcondinf(a, n, isupper, isunit);
            if( ap::fp_greater_eq(v,threshold50*ercinf) )
            {
                q50(1) = q50(1)+double(1)/double(passcount);
            }
            if( ap::fp_greater_eq(v,threshold90*ercinf) )
            {
                q90(1) = q90(1)+double(1)/double(passcount);
            }
            errless = errless||ap::fp_greater(v,ercinf*1.001);
        }
        for(i = 0; i <= 1; i++)
        {
            err50 = err50||ap::fp_less(q50(i),0.50);
            err90 = err90||ap::fp_less(q90(i),0.90);
        }
        
        //
        // degenerate matrix test
        //
        if( n>=3 )
        {
            a.setlength(n, n);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    a(i,j) = 0.0;
                }
            }
            a(0,0) = 1;
            a(n-1,n-1) = 1;
            errspec = errspec||ap::fp_neq(cmatrixtrrcond1(a, n, ap::fp_greater(ap::randomreal(),0.5), false),0);
            errspec = errspec||ap::fp_neq(cmatrixtrrcondinf(a, n, ap::fp_greater(ap::randomreal(),0.5), false),0);
        }
        
        //
        // near-degenerate matrix test
        //
        if( n>=2 )
        {
            a.setlength(n, n);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    a(i,j) = 0.0;
                }
            }
            for(i = 0; i <= n-1; i++)
            {
                a(i,i) = 1;
            }
            i = ap::randominteger(n);
            a(i,i) = 0.1*ap::maxrealnumber;
            errspec = errspec||ap::fp_neq(cmatrixtrrcond1(a, n, ap::fp_greater(ap::randomreal(),0.5), false),0);
            errspec = errspec||ap::fp_neq(cmatrixtrrcondinf(a, n, ap::fp_greater(ap::randomreal(),0.5), false),0);
        }
    }
    
    //
    // report
    //
    result = !(err50||err90||errless||errspec);
    return result;
}


/*************************************************************************
Returns True for successful test, False - for failed test
*************************************************************************/
static bool testrmatrixrcond(int maxn, int passcount)
{
    bool result;
    ap::real_2d_array a;
    ap::real_2d_array lua;
    ap::integer_1d_array p;
    int n;
    int i;
    int j;
    int pass;
    bool err50;
    bool err90;
    bool errspec;
    bool errless;
    double erc1;
    double ercinf;
    ap::real_1d_array q50;
    ap::real_1d_array q90;
    double v;

    err50 = false;
    err90 = false;
    errless = false;
    errspec = false;
    q50.setbounds(0, 3);
    q90.setbounds(0, 3);
    for(n = 1; n <= maxn; n++)
    {
        
        //
        // special test for zero matrix
        //
        rmatrixgenzero(a, n);
        rmatrixmakeacopy(a, n, n, lua);
        rmatrixlu(lua, n, n, p);
        errspec = errspec||ap::fp_neq(rmatrixrcond1(a, n),0);
        errspec = errspec||ap::fp_neq(rmatrixrcondinf(a, n),0);
        errspec = errspec||ap::fp_neq(rmatrixlurcond1(lua, n),0);
        errspec = errspec||ap::fp_neq(rmatrixlurcondinf(lua, n),0);
        
        //
        // general test
        //
        a.setbounds(0, n-1, 0, n-1);
        for(i = 0; i <= 3; i++)
        {
            q50(i) = 0;
            q90(i) = 0;
        }
        for(pass = 1; pass <= passcount; pass++)
        {
            rmatrixrndcond(n, exp(ap::randomreal()*log(double(1000))), a);
            rmatrixmakeacopy(a, n, n, lua);
            rmatrixlu(lua, n, n, p);
            rmatrixrefrcond(a, n, erc1, ercinf);
            
            //
            // 1-norm, normal
            //
            v = 1/rmatrixrcond1(a, n);
            if( ap::fp_greater_eq(v,threshold50*erc1) )
            {
                q50(0) = q50(0)+double(1)/double(passcount);
            }
            if( ap::fp_greater_eq(v,threshold90*erc1) )
            {
                q90(0) = q90(0)+double(1)/double(passcount);
            }
            errless = errless||ap::fp_greater(v,erc1*1.001);
            
            //
            // 1-norm, LU
            //
            v = 1/rmatrixlurcond1(lua, n);
            if( ap::fp_greater_eq(v,threshold50*erc1) )
            {
                q50(1) = q50(1)+double(1)/double(passcount);
            }
            if( ap::fp_greater_eq(v,threshold90*erc1) )
            {
                q90(1) = q90(1)+double(1)/double(passcount);
            }
            errless = errless||ap::fp_greater(v,erc1*1.001);
            
            //
            // Inf-norm, normal
            //
            v = 1/rmatrixrcondinf(a, n);
            if( ap::fp_greater_eq(v,threshold50*ercinf) )
            {
                q50(2) = q50(2)+double(1)/double(passcount);
            }
            if( ap::fp_greater_eq(v,threshold90*ercinf) )
            {
                q90(2) = q90(2)+double(1)/double(passcount);
            }
            errless = errless||ap::fp_greater(v,ercinf*1.001);
            
            //
            // Inf-norm, LU
            //
            v = 1/rmatrixlurcondinf(lua, n);
            if( ap::fp_greater_eq(v,threshold50*ercinf) )
            {
                q50(3) = q50(3)+double(1)/double(passcount);
            }
            if( ap::fp_greater_eq(v,threshold90*ercinf) )
            {
                q90(3) = q90(3)+double(1)/double(passcount);
            }
            errless = errless||ap::fp_greater(v,ercinf*1.001);
        }
        for(i = 0; i <= 3; i++)
        {
            err50 = err50||ap::fp_less(q50(i),0.50);
            err90 = err90||ap::fp_less(q90(i),0.90);
        }
        
        //
        // degenerate matrix test
        //
        if( n>=3 )
        {
            a.setlength(n, n);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    a(i,j) = 0.0;
                }
            }
            a(0,0) = 1;
            a(n-1,n-1) = 1;
            errspec = errspec||ap::fp_neq(rmatrixrcond1(a, n),0);
            errspec = errspec||ap::fp_neq(rmatrixrcondinf(a, n),0);
            errspec = errspec||ap::fp_neq(rmatrixlurcond1(a, n),0);
            errspec = errspec||ap::fp_neq(rmatrixlurcondinf(a, n),0);
        }
        
        //
        // near-degenerate matrix test
        //
        if( n>=2 )
        {
            a.setlength(n, n);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    a(i,j) = 0.0;
                }
            }
            for(i = 0; i <= n-1; i++)
            {
                a(i,i) = 1;
            }
            i = ap::randominteger(n);
            a(i,i) = 0.1*ap::maxrealnumber;
            errspec = errspec||ap::fp_neq(rmatrixrcond1(a, n),0);
            errspec = errspec||ap::fp_neq(rmatrixrcondinf(a, n),0);
            errspec = errspec||ap::fp_neq(rmatrixlurcond1(a, n),0);
            errspec = errspec||ap::fp_neq(rmatrixlurcondinf(a, n),0);
        }
    }
    
    //
    // report
    //
    result = !(err50||err90||errless||errspec);
    return result;
}


/*************************************************************************
Returns True for successful test, False - for failed test
*************************************************************************/
static bool testspdmatrixrcond(int maxn, int passcount)
{
    bool result;
    ap::real_2d_array a;
    ap::real_2d_array cha;
    ap::integer_1d_array p;
    int n;
    int i;
    int j;
    int pass;
    bool err50;
    bool err90;
    bool errspec;
    bool errless;
    bool isupper;
    double erc1;
    double ercinf;
    ap::real_1d_array q50;
    ap::real_1d_array q90;
    double v;

    err50 = false;
    err90 = false;
    errless = false;
    errspec = false;
    q50.setlength(2);
    q90.setlength(2);
    for(n = 1; n <= maxn; n++)
    {
        isupper = ap::fp_greater(ap::randomreal(),0.5);
        
        //
        // general test
        //
        a.setlength(n, n);
        for(i = 0; i <= 1; i++)
        {
            q50(i) = 0;
            q90(i) = 0;
        }
        for(pass = 1; pass <= passcount; pass++)
        {
            spdmatrixrndcond(n, exp(ap::randomreal()*log(double(1000))), a);
            rmatrixrefrcond(a, n, erc1, ercinf);
            rmatrixdrophalf(a, n, isupper);
            rmatrixmakeacopy(a, n, n, cha);
            spdmatrixcholesky(cha, n, isupper);
            
            //
            // normal
            //
            v = 1/spdmatrixrcond(a, n, isupper);
            if( ap::fp_greater_eq(v,threshold50*erc1) )
            {
                q50(0) = q50(0)+double(1)/double(passcount);
            }
            if( ap::fp_greater_eq(v,threshold90*erc1) )
            {
                q90(0) = q90(0)+double(1)/double(passcount);
            }
            errless = errless||ap::fp_greater(v,erc1*1.001);
            
            //
            // Cholesky
            //
            v = 1/spdmatrixcholeskyrcond(cha, n, isupper);
            if( ap::fp_greater_eq(v,threshold50*erc1) )
            {
                q50(1) = q50(1)+double(1)/double(passcount);
            }
            if( ap::fp_greater_eq(v,threshold90*erc1) )
            {
                q90(1) = q90(1)+double(1)/double(passcount);
            }
            errless = errless||ap::fp_greater(v,erc1*1.001);
        }
        for(i = 0; i <= 1; i++)
        {
            err50 = err50||ap::fp_less(q50(i),0.50);
            err90 = err90||ap::fp_less(q90(i),0.90);
        }
        
        //
        // degenerate matrix test
        //
        if( n>=3 )
        {
            a.setlength(n, n);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    a(i,j) = 0.0;
                }
            }
            a(0,0) = 1;
            a(n-1,n-1) = 1;
            errspec = errspec||ap::fp_neq(spdmatrixrcond(a, n, isupper),-1);
            errspec = errspec||ap::fp_neq(spdmatrixcholeskyrcond(a, n, isupper),0);
        }
        
        //
        // near-degenerate matrix test
        //
        if( n>=2 )
        {
            a.setlength(n, n);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    a(i,j) = 0.0;
                }
            }
            for(i = 0; i <= n-1; i++)
            {
                a(i,i) = 1;
            }
            i = ap::randominteger(n);
            a(i,i) = 0.1*ap::maxrealnumber;
            errspec = errspec||ap::fp_neq(spdmatrixrcond(a, n, isupper),0);
            errspec = errspec||ap::fp_neq(spdmatrixcholeskyrcond(a, n, isupper),0);
        }
    }
    
    //
    // report
    //
    result = !(err50||err90||errless||errspec);
    return result;
}


/*************************************************************************
Returns True for successful test, False - for failed test
*************************************************************************/
static bool testcmatrixrcond(int maxn, int passcount)
{
    bool result;
    ap::complex_2d_array a;
    ap::complex_2d_array lua;
    ap::integer_1d_array p;
    int n;
    int i;
    int j;
    int pass;
    bool err50;
    bool err90;
    bool errless;
    bool errspec;
    double erc1;
    double ercinf;
    ap::real_1d_array q50;
    ap::real_1d_array q90;
    double v;

    q50.setbounds(0, 3);
    q90.setbounds(0, 3);
    err50 = false;
    err90 = false;
    errless = false;
    errspec = false;
    
    //
    // process
    //
    for(n = 1; n <= maxn; n++)
    {
        
        //
        // special test for zero matrix
        //
        cmatrixgenzero(a, n);
        cmatrixmakeacopy(a, n, n, lua);
        cmatrixlu(lua, n, n, p);
        errspec = errspec||ap::fp_neq(cmatrixrcond1(a, n),0);
        errspec = errspec||ap::fp_neq(cmatrixrcondinf(a, n),0);
        errspec = errspec||ap::fp_neq(cmatrixlurcond1(lua, n),0);
        errspec = errspec||ap::fp_neq(cmatrixlurcondinf(lua, n),0);
        
        //
        // general test
        //
        a.setbounds(0, n-1, 0, n-1);
        for(i = 0; i <= 3; i++)
        {
            q50(i) = 0;
            q90(i) = 0;
        }
        for(pass = 1; pass <= passcount; pass++)
        {
            cmatrixrndcond(n, exp(ap::randomreal()*log(double(1000))), a);
            cmatrixmakeacopy(a, n, n, lua);
            cmatrixlu(lua, n, n, p);
            cmatrixrefrcond(a, n, erc1, ercinf);
            
            //
            // 1-norm, normal
            //
            v = 1/cmatrixrcond1(a, n);
            if( ap::fp_greater_eq(v,threshold50*erc1) )
            {
                q50(0) = q50(0)+double(1)/double(passcount);
            }
            if( ap::fp_greater_eq(v,threshold90*erc1) )
            {
                q90(0) = q90(0)+double(1)/double(passcount);
            }
            errless = errless||ap::fp_greater(v,erc1*1.001);
            
            //
            // 1-norm, LU
            //
            v = 1/cmatrixlurcond1(lua, n);
            if( ap::fp_greater_eq(v,threshold50*erc1) )
            {
                q50(1) = q50(1)+double(1)/double(passcount);
            }
            if( ap::fp_greater_eq(v,threshold90*erc1) )
            {
                q90(1) = q90(1)+double(1)/double(passcount);
            }
            errless = errless||ap::fp_greater(v,erc1*1.001);
            
            //
            // Inf-norm, normal
            //
            v = 1/cmatrixrcondinf(a, n);
            if( ap::fp_greater_eq(v,threshold50*ercinf) )
            {
                q50(2) = q50(2)+double(1)/double(passcount);
            }
            if( ap::fp_greater_eq(v,threshold90*ercinf) )
            {
                q90(2) = q90(2)+double(1)/double(passcount);
            }
            errless = errless||ap::fp_greater(v,ercinf*1.001);
            
            //
            // Inf-norm, LU
            //
            v = 1/cmatrixlurcondinf(lua, n);
            if( ap::fp_greater_eq(v,threshold50*ercinf) )
            {
                q50(3) = q50(3)+double(1)/double(passcount);
            }
            if( ap::fp_greater_eq(v,threshold90*ercinf) )
            {
                q90(3) = q90(3)+double(1)/double(passcount);
            }
            errless = errless||ap::fp_greater(v,ercinf*1.001);
        }
        for(i = 0; i <= 3; i++)
        {
            err50 = err50||ap::fp_less(q50(i),0.50);
            err90 = err90||ap::fp_less(q90(i),0.90);
        }
        
        //
        // degenerate matrix test
        //
        if( n>=3 )
        {
            a.setlength(n, n);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    a(i,j) = 0.0;
                }
            }
            a(0,0) = 1;
            a(n-1,n-1) = 1;
            errspec = errspec||ap::fp_neq(cmatrixrcond1(a, n),0);
            errspec = errspec||ap::fp_neq(cmatrixrcondinf(a, n),0);
            errspec = errspec||ap::fp_neq(cmatrixlurcond1(a, n),0);
            errspec = errspec||ap::fp_neq(cmatrixlurcondinf(a, n),0);
        }
        
        //
        // near-degenerate matrix test
        //
        if( n>=2 )
        {
            a.setlength(n, n);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    a(i,j) = 0.0;
                }
            }
            for(i = 0; i <= n-1; i++)
            {
                a(i,i) = 1;
            }
            i = ap::randominteger(n);
            a(i,i) = 0.1*ap::maxrealnumber;
            errspec = errspec||ap::fp_neq(cmatrixrcond1(a, n),0);
            errspec = errspec||ap::fp_neq(cmatrixrcondinf(a, n),0);
            errspec = errspec||ap::fp_neq(cmatrixlurcond1(a, n),0);
            errspec = errspec||ap::fp_neq(cmatrixlurcondinf(a, n),0);
        }
    }
    
    //
    // report
    //
    result = !(err50||err90||errless||errspec);
    return result;
}


/*************************************************************************
Returns True for successful test, False - for failed test
*************************************************************************/
static bool testhpdmatrixrcond(int maxn, int passcount)
{
    bool result;
    ap::complex_2d_array a;
    ap::complex_2d_array cha;
    ap::integer_1d_array p;
    int n;
    int i;
    int j;
    int pass;
    bool err50;
    bool err90;
    bool errspec;
    bool errless;
    bool isupper;
    double erc1;
    double ercinf;
    ap::real_1d_array q50;
    ap::real_1d_array q90;
    double v;

    err50 = false;
    err90 = false;
    errless = false;
    errspec = false;
    q50.setlength(2);
    q90.setlength(2);
    for(n = 1; n <= maxn; n++)
    {
        isupper = ap::fp_greater(ap::randomreal(),0.5);
        
        //
        // general test
        //
        a.setlength(n, n);
        for(i = 0; i <= 1; i++)
        {
            q50(i) = 0;
            q90(i) = 0;
        }
        for(pass = 1; pass <= passcount; pass++)
        {
            hpdmatrixrndcond(n, exp(ap::randomreal()*log(double(1000))), a);
            cmatrixrefrcond(a, n, erc1, ercinf);
            cmatrixdrophalf(a, n, isupper);
            cmatrixmakeacopy(a, n, n, cha);
            hpdmatrixcholesky(cha, n, isupper);
            
            //
            // normal
            //
            v = 1/hpdmatrixrcond(a, n, isupper);
            if( ap::fp_greater_eq(v,threshold50*erc1) )
            {
                q50(0) = q50(0)+double(1)/double(passcount);
            }
            if( ap::fp_greater_eq(v,threshold90*erc1) )
            {
                q90(0) = q90(0)+double(1)/double(passcount);
            }
            errless = errless||ap::fp_greater(v,erc1*1.001);
            
            //
            // Cholesky
            //
            v = 1/hpdmatrixcholeskyrcond(cha, n, isupper);
            if( ap::fp_greater_eq(v,threshold50*erc1) )
            {
                q50(1) = q50(1)+double(1)/double(passcount);
            }
            if( ap::fp_greater_eq(v,threshold90*erc1) )
            {
                q90(1) = q90(1)+double(1)/double(passcount);
            }
            errless = errless||ap::fp_greater(v,erc1*1.001);
        }
        for(i = 0; i <= 1; i++)
        {
            err50 = err50||ap::fp_less(q50(i),0.50);
            err90 = err90||ap::fp_less(q90(i),0.90);
        }
        
        //
        // degenerate matrix test
        //
        if( n>=3 )
        {
            a.setlength(n, n);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    a(i,j) = 0.0;
                }
            }
            a(0,0) = 1;
            a(n-1,n-1) = 1;
            errspec = errspec||ap::fp_neq(hpdmatrixrcond(a, n, isupper),-1);
            errspec = errspec||ap::fp_neq(hpdmatrixcholeskyrcond(a, n, isupper),0);
        }
        
        //
        // near-degenerate matrix test
        //
        if( n>=2 )
        {
            a.setlength(n, n);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    a(i,j) = 0.0;
                }
            }
            for(i = 0; i <= n-1; i++)
            {
                a(i,i) = 1;
            }
            i = ap::randominteger(n);
            a(i,i) = 0.1*ap::maxrealnumber;
            errspec = errspec||ap::fp_neq(hpdmatrixrcond(a, n, isupper),0);
            errspec = errspec||ap::fp_neq(hpdmatrixcholeskyrcond(a, n, isupper),0);
        }
    }
    
    //
    // report
    //
    result = !(err50||err90||errless||errspec);
    return result;
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testrcondunit_test_silent()
{
    bool result;

    result = testrcond(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testrcondunit_test()
{
    bool result;

    result = testrcond(false);
    return result;
}




