
 
#include <stdio.h>
#include "testmatgenunit.h"

static const int maxsvditerations = 60;

static void unset2d(ap::real_2d_array& a);
static void unset2dc(ap::complex_2d_array& a);
static bool ishpd(ap::complex_2d_array a, int n);
static double svdcond(const ap::real_2d_array& a, int n);
static double extsign(double a, double b);
static double mymax(double a, double b);
static double pythag(double a, double b);

bool testmatgen(bool silent)
{
    bool result;
    ap::real_2d_array a;
    ap::real_2d_array b;
    ap::real_2d_array u;
    ap::real_2d_array v;
    ap::complex_2d_array ca;
    ap::complex_2d_array cb;
    ap::real_2d_array r1;
    ap::real_2d_array r2;
    ap::complex_2d_array c1;
    ap::complex_2d_array c2;
    ap::real_1d_array w;
    int n;
    int maxn;
    int i;
    int j;
    int pass;
    int passcount;
    bool waserrors;
    double cond;
    double threshold;
    double vt;
    ap::complex ct;
    double minw;
    double maxw;
    bool serr;
    bool herr;
    bool spderr;
    bool hpderr;
    bool rerr;
    bool cerr;

    rerr = false;
    cerr = false;
    serr = false;
    herr = false;
    spderr = false;
    hpderr = false;
    waserrors = false;
    maxn = 20;
    passcount = 15;
    threshold = 1000*ap::machineepsilon;
    
    //
    // Testing orthogonal
    //
    for(n = 1; n <= maxn; n++)
    {
        for(pass = 1; pass <= passcount; pass++)
        {
            r1.setbounds(0, n-1, 0, 2*n-1);
            r2.setbounds(0, 2*n-1, 0, n-1);
            c1.setbounds(0, n-1, 0, 2*n-1);
            c2.setbounds(0, 2*n-1, 0, n-1);
            
            //
            // Random orthogonal, real
            //
            unset2d(a);
            unset2d(b);
            rmatrixrndorthogonal(n, a);
            rmatrixrndorthogonal(n, b);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    
                    //
                    // orthogonality test
                    //
                    vt = ap::vdotproduct(&a(i, 0), 1, &a(j, 0), 1, ap::vlen(0,n-1));
                    if( i==j )
                    {
                        rerr = rerr||ap::fp_greater(fabs(vt-1),threshold);
                    }
                    else
                    {
                        rerr = rerr||ap::fp_greater(fabs(vt),threshold);
                    }
                    vt = ap::vdotproduct(&b(i, 0), 1, &b(j, 0), 1, ap::vlen(0,n-1));
                    if( i==j )
                    {
                        rerr = rerr||ap::fp_greater(fabs(vt-1),threshold);
                    }
                    else
                    {
                        rerr = rerr||ap::fp_greater(fabs(vt),threshold);
                    }
                    
                    //
                    // test for difference in A and B
                    //
                    if( n>=2 )
                    {
                        rerr = rerr||ap::fp_eq(a(i,j),b(i,j));
                    }
                }
            }
            
            //
            // Random orthogonal, complex
            //
            unset2dc(ca);
            unset2dc(cb);
            cmatrixrndorthogonal(n, ca);
            cmatrixrndorthogonal(n, cb);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    
                    //
                    // orthogonality test
                    //
                    ct = ap::vdotproduct(&ca(i, 0), 1, "N", &ca(j, 0), 1, "Conj", ap::vlen(0,n-1));
                    if( i==j )
                    {
                        cerr = cerr||ap::fp_greater(ap::abscomplex(ct-1),threshold);
                    }
                    else
                    {
                        cerr = cerr||ap::fp_greater(ap::abscomplex(ct),threshold);
                    }
                    ct = ap::vdotproduct(&cb(i, 0), 1, "N", &cb(j, 0), 1, "Conj", ap::vlen(0,n-1));
                    if( i==j )
                    {
                        cerr = cerr||ap::fp_greater(ap::abscomplex(ct-1),threshold);
                    }
                    else
                    {
                        cerr = cerr||ap::fp_greater(ap::abscomplex(ct),threshold);
                    }
                    
                    //
                    // test for difference in A and B
                    //
                    if( n>=2 )
                    {
                        cerr = cerr||ca(i,j)==cb(i,j);
                    }
                }
            }
            
            //
            // From the right real tests:
            // 1. E*Q is orthogonal
            // 2. Q1<>Q2 (routine result is changing)
            // 3. (E E)'*Q = (Q' Q')' (correct handling of non-square matrices)
            //
            unset2d(a);
            unset2d(b);
            a.setbounds(0, n-1, 0, n-1);
            b.setbounds(0, n-1, 0, n-1);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    a(i,j) = 0;
                    b(i,j) = 0;
                }
                a(i,i) = 1;
                b(i,i) = 1;
            }
            rmatrixrndorthogonalfromtheright(a, n, n);
            rmatrixrndorthogonalfromtheright(b, n, n);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    
                    //
                    // orthogonality test
                    //
                    vt = ap::vdotproduct(&a(i, 0), 1, &a(j, 0), 1, ap::vlen(0,n-1));
                    if( i==j )
                    {
                        rerr = rerr||ap::fp_greater(fabs(vt-1),threshold);
                    }
                    else
                    {
                        rerr = rerr||ap::fp_greater(fabs(vt),threshold);
                    }
                    vt = ap::vdotproduct(&b(i, 0), 1, &b(j, 0), 1, ap::vlen(0,n-1));
                    if( i==j )
                    {
                        rerr = rerr||ap::fp_greater(fabs(vt-1),threshold);
                    }
                    else
                    {
                        rerr = rerr||ap::fp_greater(fabs(vt),threshold);
                    }
                    
                    //
                    // test for difference in A and B
                    //
                    if( n>=2 )
                    {
                        rerr = rerr||ap::fp_eq(a(i,j),b(i,j));
                    }
                }
            }
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    r2(i,j) = 2*ap::randomreal()-1;
                    r2(i+n,j) = r2(i,j);
                }
            }
            rmatrixrndorthogonalfromtheright(r2, 2*n, n);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    rerr = rerr||ap::fp_greater(fabs(r2(i+n,j)-r2(i,j)),threshold);
                }
            }
            
            //
            // From the left real tests:
            // 1. Q*E is orthogonal
            // 2. Q1<>Q2 (routine result is changing)
            // 3. Q*(E E) = (Q Q) (correct handling of non-square matrices)
            //
            unset2d(a);
            unset2d(b);
            a.setbounds(0, n-1, 0, n-1);
            b.setbounds(0, n-1, 0, n-1);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    a(i,j) = 0;
                    b(i,j) = 0;
                }
                a(i,i) = 1;
                b(i,i) = 1;
            }
            rmatrixrndorthogonalfromtheleft(a, n, n);
            rmatrixrndorthogonalfromtheleft(b, n, n);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    
                    //
                    // orthogonality test
                    //
                    vt = ap::vdotproduct(&a(i, 0), 1, &a(j, 0), 1, ap::vlen(0,n-1));
                    if( i==j )
                    {
                        rerr = rerr||ap::fp_greater(fabs(vt-1),threshold);
                    }
                    else
                    {
                        rerr = rerr||ap::fp_greater(fabs(vt),threshold);
                    }
                    vt = ap::vdotproduct(&b(i, 0), 1, &b(j, 0), 1, ap::vlen(0,n-1));
                    if( i==j )
                    {
                        rerr = rerr||ap::fp_greater(fabs(vt-1),threshold);
                    }
                    else
                    {
                        rerr = rerr||ap::fp_greater(fabs(vt),threshold);
                    }
                    
                    //
                    // test for difference in A and B
                    //
                    if( n>=2 )
                    {
                        rerr = rerr||ap::fp_eq(a(i,j),b(i,j));
                    }
                }
            }
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    r1(i,j) = 2*ap::randomreal()-1;
                    r1(i,j+n) = r1(i,j);
                }
            }
            rmatrixrndorthogonalfromtheleft(r1, n, 2*n);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    rerr = rerr||ap::fp_greater(fabs(r1(i,j)-r1(i,j+n)),threshold);
                }
            }
            
            //
            // From the right complex tests:
            // 1. E*Q is orthogonal
            // 2. Q1<>Q2 (routine result is changing)
            // 3. (E E)'*Q = (Q' Q')' (correct handling of non-square matrices)
            //
            unset2dc(ca);
            unset2dc(cb);
            ca.setbounds(0, n-1, 0, n-1);
            cb.setbounds(0, n-1, 0, n-1);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    ca(i,j) = 0;
                    cb(i,j) = 0;
                }
                ca(i,i) = 1;
                cb(i,i) = 1;
            }
            cmatrixrndorthogonalfromtheright(ca, n, n);
            cmatrixrndorthogonalfromtheright(cb, n, n);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    
                    //
                    // orthogonality test
                    //
                    ct = ap::vdotproduct(&ca(i, 0), 1, "N", &ca(j, 0), 1, "Conj", ap::vlen(0,n-1));
                    if( i==j )
                    {
                        cerr = cerr||ap::fp_greater(ap::abscomplex(ct-1),threshold);
                    }
                    else
                    {
                        cerr = cerr||ap::fp_greater(ap::abscomplex(ct),threshold);
                    }
                    ct = ap::vdotproduct(&cb(i, 0), 1, "N", &cb(j, 0), 1, "Conj", ap::vlen(0,n-1));
                    if( i==j )
                    {
                        cerr = cerr||ap::fp_greater(ap::abscomplex(ct-1),threshold);
                    }
                    else
                    {
                        cerr = cerr||ap::fp_greater(ap::abscomplex(ct),threshold);
                    }
                    
                    //
                    // test for difference in A and B
                    //
                    cerr = cerr||ca(i,j)==cb(i,j);
                }
            }
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    c2(i,j) = 2*ap::randomreal()-1;
                    c2(i+n,j) = c2(i,j);
                }
            }
            cmatrixrndorthogonalfromtheright(c2, 2*n, n);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    cerr = cerr||ap::fp_greater(ap::abscomplex(c2(i+n,j)-c2(i,j)),threshold);
                }
            }
            
            //
            // From the left complex tests:
            // 1. Q*E is orthogonal
            // 2. Q1<>Q2 (routine result is changing)
            // 3. Q*(E E) = (Q Q) (correct handling of non-square matrices)
            //
            unset2dc(ca);
            unset2dc(cb);
            ca.setbounds(0, n-1, 0, n-1);
            cb.setbounds(0, n-1, 0, n-1);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    ca(i,j) = 0;
                    cb(i,j) = 0;
                }
                ca(i,i) = 1;
                cb(i,i) = 1;
            }
            cmatrixrndorthogonalfromtheleft(ca, n, n);
            cmatrixrndorthogonalfromtheleft(cb, n, n);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    
                    //
                    // orthogonality test
                    //
                    ct = ap::vdotproduct(&ca(i, 0), 1, "N", &ca(j, 0), 1, "Conj", ap::vlen(0,n-1));
                    if( i==j )
                    {
                        cerr = cerr||ap::fp_greater(ap::abscomplex(ct-1),threshold);
                    }
                    else
                    {
                        cerr = cerr||ap::fp_greater(ap::abscomplex(ct),threshold);
                    }
                    ct = ap::vdotproduct(&cb(i, 0), 1, "N", &cb(j, 0), 1, "Conj", ap::vlen(0,n-1));
                    if( i==j )
                    {
                        cerr = cerr||ap::fp_greater(ap::abscomplex(ct-1),threshold);
                    }
                    else
                    {
                        cerr = cerr||ap::fp_greater(ap::abscomplex(ct),threshold);
                    }
                    
                    //
                    // test for difference in A and B
                    //
                    cerr = cerr||ca(i,j)==cb(i,j);
                }
            }
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    c1(i,j) = 2*ap::randomreal()-1;
                    c1(i,j+n) = c1(i,j);
                }
            }
            cmatrixrndorthogonalfromtheleft(c1, n, 2*n);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    cerr = cerr||ap::fp_greater(ap::abscomplex(c1(i,j)-c1(i,j+n)),threshold);
                }
            }
        }
    }
    
    //
    // Testing GCond
    //
    for(n = 2; n <= maxn; n++)
    {
        for(pass = 1; pass <= passcount; pass++)
        {
            
            //
            // real test
            //
            unset2d(a);
            cond = exp(log(double(1000))*ap::randomreal());
            rmatrixrndcond(n, cond, a);
            b.setbounds(1, n, 1, n);
            for(i = 1; i <= n; i++)
            {
                for(j = 1; j <= n; j++)
                {
                    b(i,j) = a(i-1,j-1);
                }
            }
            if( obsoletesvddecomposition(b, n, n, w, v) )
            {
                maxw = w(1);
                minw = w(1);
                for(i = 2; i <= n; i++)
                {
                    if( ap::fp_greater(w(i),maxw) )
                    {
                        maxw = w(i);
                    }
                    if( ap::fp_less(w(i),minw) )
                    {
                        minw = w(i);
                    }
                }
                vt = maxw/minw/cond;
                if( ap::fp_greater(fabs(log(vt)),log(1+threshold)) )
                {
                    rerr = true;
                }
            }
        }
    }
    
    //
    // Symmetric/SPD
    // N = 2 .. 30
    //
    for(n = 2; n <= maxn; n++)
    {
        
        //
        // SPD matrices
        //
        for(pass = 1; pass <= passcount; pass++)
        {
            
            //
            // Generate A
            //
            unset2d(a);
            cond = exp(log(double(1000))*ap::randomreal());
            spdmatrixrndcond(n, cond, a);
            
            //
            // test condition number
            //
            spderr = spderr||ap::fp_greater(svdcond(a, n)/cond-1,threshold);
            
            //
            // test SPD
            //
            spderr = spderr||!isspd(a, n, true);
            
            //
            // test that A is symmetic
            //
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    spderr = spderr||ap::fp_greater(fabs(a(i,j)-a(j,i)),threshold);
                }
            }
            
            //
            // test for difference between A and B (subsequent matrix)
            //
            unset2d(b);
            spdmatrixrndcond(n, cond, b);
            if( n>=2 )
            {
                for(i = 0; i <= n-1; i++)
                {
                    for(j = 0; j <= n-1; j++)
                    {
                        spderr = spderr||ap::fp_eq(a(i,j),b(i,j));
                    }
                }
            }
        }
        
        //
        // HPD matrices
        //
        for(pass = 1; pass <= passcount; pass++)
        {
            
            //
            // Generate A
            //
            unset2dc(ca);
            cond = exp(log(double(1000))*ap::randomreal());
            hpdmatrixrndcond(n, cond, ca);
            
            //
            // test HPD
            //
            hpderr = hpderr||!ishpd(ca, n);
            
            //
            // test that A is Hermitian
            //
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    hpderr = hpderr||ap::fp_greater(ap::abscomplex(ca(i,j)-ap::conj(ca(j,i))),threshold);
                }
            }
            
            //
            // test for difference between A and B (subsequent matrix)
            //
            unset2dc(cb);
            hpdmatrixrndcond(n, cond, cb);
            if( n>=2 )
            {
                for(i = 0; i <= n-1; i++)
                {
                    for(j = 0; j <= n-1; j++)
                    {
                        hpderr = hpderr||ca(i,j)==cb(i,j);
                    }
                }
            }
        }
        
        //
        // Symmetric matrices
        //
        for(pass = 1; pass <= passcount; pass++)
        {
            
            //
            // test condition number
            //
            unset2d(a);
            cond = exp(log(double(1000))*ap::randomreal());
            smatrixrndcond(n, cond, a);
            serr = serr||ap::fp_greater(svdcond(a, n)/cond-1,threshold);
            
            //
            // test for difference between A and B
            //
            unset2d(b);
            smatrixrndcond(n, cond, b);
            if( n>=2 )
            {
                for(i = 0; i <= n-1; i++)
                {
                    for(j = 0; j <= n-1; j++)
                    {
                        serr = serr||ap::fp_eq(a(i,j),b(i,j));
                    }
                }
            }
        }
        
        //
        // Hermitian matrices
        //
        for(pass = 1; pass <= passcount; pass++)
        {
            
            //
            // Generate A
            //
            unset2dc(ca);
            cond = exp(log(double(1000))*ap::randomreal());
            hmatrixrndcond(n, cond, ca);
            
            //
            // test that A is Hermitian
            //
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    herr = herr||ap::fp_greater(ap::abscomplex(ca(i,j)-ap::conj(ca(j,i))),threshold);
                }
            }
            
            //
            // test for difference between A and B (subsequent matrix)
            //
            unset2dc(cb);
            hmatrixrndcond(n, cond, cb);
            if( n>=2 )
            {
                for(i = 0; i <= n-1; i++)
                {
                    for(j = 0; j <= n-1; j++)
                    {
                        herr = herr||ca(i,j)==cb(i,j);
                    }
                }
            }
        }
    }
    
    //
    // report
    //
    waserrors = rerr||cerr||serr||spderr||herr||hpderr;
    if( !silent )
    {
        printf("TESTING MATRIX GENERATOR\n");
        printf("REAL TEST:                               ");
        if( !rerr )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("COMPLEX TEST:                            ");
        if( !cerr )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("SYMMETRIC TEST:                          ");
        if( !serr )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("HERMITIAN TEST:                          ");
        if( !herr )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("SPD TEST:                                ");
        if( !spderr )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("HPD TEST:                                ");
        if( !hpderr )
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
Test whether matrix is SPD
*************************************************************************/
bool isspd(ap::real_2d_array a, int n, bool isupper)
{
    bool result;
    int i;
    int j;
    double ajj;
    double v;

    
    //
    //     Test the input parameters.
    //
    ap::ap_error::make_assertion(n>=0, "Error in SMatrixCholesky: incorrect function arguments");
    
    //
    //     Quick return if possible
    //
    result = true;
    if( n<=0 )
    {
        return result;
    }
    if( isupper )
    {
        
        //
        // Compute the Cholesky factorization A = U'*U.
        //
        for(j = 0; j <= n-1; j++)
        {
            
            //
            // Compute U(J,J) and test for non-positive-definiteness.
            //
            v = ap::vdotproduct(&a(0, j), a.getstride(), &a(0, j), a.getstride(), ap::vlen(0,j-1));
            ajj = a(j,j)-v;
            if( ap::fp_less_eq(ajj,0) )
            {
                result = false;
                return result;
            }
            ajj = sqrt(ajj);
            a(j,j) = ajj;
            
            //
            // Compute elements J+1:N of row J.
            //
            if( j<n-1 )
            {
                for(i = j+1; i <= n-1; i++)
                {
                    v = ap::vdotproduct(&a(0, i), a.getstride(), &a(0, j), a.getstride(), ap::vlen(0,j-1));
                    a(j,i) = a(j,i)-v;
                }
                v = 1/ajj;
                ap::vmul(&a(j, j+1), 1, ap::vlen(j+1,n-1), v);
            }
        }
    }
    else
    {
        
        //
        // Compute the Cholesky factorization A = L*L'.
        //
        for(j = 0; j <= n-1; j++)
        {
            
            //
            // Compute L(J,J) and test for non-positive-definiteness.
            //
            v = ap::vdotproduct(&a(j, 0), 1, &a(j, 0), 1, ap::vlen(0,j-1));
            ajj = a(j,j)-v;
            if( ap::fp_less_eq(ajj,0) )
            {
                result = false;
                return result;
            }
            ajj = sqrt(ajj);
            a(j,j) = ajj;
            
            //
            // Compute elements J+1:N of column J.
            //
            if( j<n-1 )
            {
                for(i = j+1; i <= n-1; i++)
                {
                    v = ap::vdotproduct(&a(i, 0), 1, &a(j, 0), 1, ap::vlen(0,j-1));
                    a(i,j) = a(i,j)-v;
                }
                v = 1/ajj;
                ap::vmul(&a(j+1, j), a.getstride(), ap::vlen(j+1,n-1), v);
            }
        }
    }
    return result;
}


bool obsoletesvddecomposition(ap::real_2d_array& a,
     int m,
     int n,
     ap::real_1d_array& w,
     ap::real_2d_array& v)
{
    bool result;
    int nm;
    int minmn;
    int l;
    int k;
    int j;
    int jj;
    int its;
    int i;
    double z;
    double y;
    double x;
    double vscale;
    double s;
    double h;
    double g;
    double f;
    double c;
    double anorm;
    ap::real_1d_array rv1;
    bool flag;

    rv1.setbounds(1, n);
    w.setbounds(1, n);
    v.setbounds(1, n, 1, n);
    result = true;
    if( m<n )
    {
        minmn = m;
    }
    else
    {
        minmn = n;
    }
    g = 0.0;
    vscale = 0.0;
    anorm = 0.0;
    for(i = 1; i <= n; i++)
    {
        l = i+1;
        rv1(i) = vscale*g;
        g = 0;
        s = 0;
        vscale = 0;
        if( i<=m )
        {
            for(k = i; k <= m; k++)
            {
                vscale = vscale+fabs(a(k,i));
            }
            if( ap::fp_neq(vscale,0.0) )
            {
                for(k = i; k <= m; k++)
                {
                    a(k,i) = a(k,i)/vscale;
                    s = s+a(k,i)*a(k,i);
                }
                f = a(i,i);
                g = -extsign(sqrt(s), f);
                h = f*g-s;
                a(i,i) = f-g;
                if( i!=n )
                {
                    for(j = l; j <= n; j++)
                    {
                        s = 0.0;
                        for(k = i; k <= m; k++)
                        {
                            s = s+a(k,i)*a(k,j);
                        }
                        f = s/h;
                        for(k = i; k <= m; k++)
                        {
                            a(k,j) = a(k,j)+f*a(k,i);
                        }
                    }
                }
                for(k = i; k <= m; k++)
                {
                    a(k,i) = vscale*a(k,i);
                }
            }
        }
        w(i) = vscale*g;
        g = 0.0;
        s = 0.0;
        vscale = 0.0;
        if( i<=m&&i!=n )
        {
            for(k = l; k <= n; k++)
            {
                vscale = vscale+fabs(a(i,k));
            }
            if( ap::fp_neq(vscale,0.0) )
            {
                for(k = l; k <= n; k++)
                {
                    a(i,k) = a(i,k)/vscale;
                    s = s+a(i,k)*a(i,k);
                }
                f = a(i,l);
                g = -extsign(sqrt(s), f);
                h = f*g-s;
                a(i,l) = f-g;
                for(k = l; k <= n; k++)
                {
                    rv1(k) = a(i,k)/h;
                }
                if( i!=m )
                {
                    for(j = l; j <= m; j++)
                    {
                        s = 0.0;
                        for(k = l; k <= n; k++)
                        {
                            s = s+a(j,k)*a(i,k);
                        }
                        for(k = l; k <= n; k++)
                        {
                            a(j,k) = a(j,k)+s*rv1(k);
                        }
                    }
                }
                for(k = l; k <= n; k++)
                {
                    a(i,k) = vscale*a(i,k);
                }
            }
        }
        anorm = mymax(anorm, fabs(w(i))+fabs(rv1(i)));
    }
    for(i = n; i >= 1; i--)
    {
        if( i<n )
        {
            if( ap::fp_neq(g,0.0) )
            {
                for(j = l; j <= n; j++)
                {
                    v(j,i) = a(i,j)/a(i,l)/g;
                }
                for(j = l; j <= n; j++)
                {
                    s = 0.0;
                    for(k = l; k <= n; k++)
                    {
                        s = s+a(i,k)*v(k,j);
                    }
                    for(k = l; k <= n; k++)
                    {
                        v(k,j) = v(k,j)+s*v(k,i);
                    }
                }
            }
            for(j = l; j <= n; j++)
            {
                v(i,j) = 0.0;
                v(j,i) = 0.0;
            }
        }
        v(i,i) = 1.0;
        g = rv1(i);
        l = i;
    }
    for(i = minmn; i >= 1; i--)
    {
        l = i+1;
        g = w(i);
        if( i<n )
        {
            for(j = l; j <= n; j++)
            {
                a(i,j) = 0.0;
            }
        }
        if( ap::fp_neq(g,0.0) )
        {
            g = 1.0/g;
            if( i!=n )
            {
                for(j = l; j <= n; j++)
                {
                    s = 0.0;
                    for(k = l; k <= m; k++)
                    {
                        s = s+a(k,i)*a(k,j);
                    }
                    f = s/a(i,i)*g;
                    for(k = i; k <= m; k++)
                    {
                        a(k,j) = a(k,j)+f*a(k,i);
                    }
                }
            }
            for(j = i; j <= m; j++)
            {
                a(j,i) = a(j,i)*g;
            }
        }
        else
        {
            for(j = i; j <= m; j++)
            {
                a(j,i) = 0.0;
            }
        }
        a(i,i) = a(i,i)+1.0;
    }
    for(k = n; k >= 1; k--)
    {
        for(its = 1; its <= maxsvditerations; its++)
        {
            flag = true;
            for(l = k; l >= 1; l--)
            {
                nm = l-1;
                if( ap::fp_eq(fabs(rv1(l))+anorm,anorm) )
                {
                    flag = false;
                    break;
                }
                if( ap::fp_eq(fabs(w(nm))+anorm,anorm) )
                {
                    break;
                }
            }
            if( flag )
            {
                c = 0.0;
                s = 1.0;
                for(i = l; i <= k; i++)
                {
                    f = s*rv1(i);
                    if( ap::fp_neq(fabs(f)+anorm,anorm) )
                    {
                        g = w(i);
                        h = pythag(f, g);
                        w(i) = h;
                        h = 1.0/h;
                        c = g*h;
                        s = -f*h;
                        for(j = 1; j <= m; j++)
                        {
                            y = a(j,nm);
                            z = a(j,i);
                            a(j,nm) = y*c+z*s;
                            a(j,i) = -y*s+z*c;
                        }
                    }
                }
            }
            z = w(k);
            if( l==k )
            {
                if( ap::fp_less(z,0.0) )
                {
                    w(k) = -z;
                    for(j = 1; j <= n; j++)
                    {
                        v(j,k) = -v(j,k);
                    }
                }
                break;
            }
            if( its==maxsvditerations )
            {
                result = false;
                return result;
            }
            x = w(l);
            nm = k-1;
            y = w(nm);
            g = rv1(nm);
            h = rv1(k);
            f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g = pythag(f, double(1));
            f = ((x-z)*(x+z)+h*(y/(f+extsign(g, f))-h))/x;
            c = 1.0;
            s = 1.0;
            for(j = l; j <= nm; j++)
            {
                i = j+1;
                g = rv1(i);
                y = w(i);
                h = s*g;
                g = c*g;
                z = pythag(f, h);
                rv1(j) = z;
                c = f/z;
                s = h/z;
                f = x*c+g*s;
                g = -x*s+g*c;
                h = y*s;
                y = y*c;
                for(jj = 1; jj <= n; jj++)
                {
                    x = v(jj,j);
                    z = v(jj,i);
                    v(jj,j) = x*c+z*s;
                    v(jj,i) = -x*s+z*c;
                }
                z = pythag(f, h);
                w(j) = z;
                if( ap::fp_neq(z,0.0) )
                {
                    z = 1.0/z;
                    c = f*z;
                    s = h*z;
                }
                f = c*g+s*y;
                x = -s*g+c*y;
                for(jj = 1; jj <= m; jj++)
                {
                    y = a(jj,j);
                    z = a(jj,i);
                    a(jj,j) = y*c+z*s;
                    a(jj,i) = -y*s+z*c;
                }
            }
            rv1(l) = 0.0;
            rv1(k) = f;
            w(k) = x;
        }
    }
    return result;
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
static void unset2dc(ap::complex_2d_array& a)
{

    a.setbounds(0, 0, 0, 0);
    a(0,0) = 2*ap::randomreal()-1;
}


/*************************************************************************
Tests whether A is HPD
*************************************************************************/
static bool ishpd(ap::complex_2d_array a, int n)
{
    bool result;
    int j;
    double ajj;
    ap::complex v;
    double r;
    ap::complex_1d_array t;
    ap::complex_1d_array t2;
    ap::complex_1d_array t3;
    int i;
    ap::complex_2d_array a1;

    t.setbounds(0, n-1);
    t2.setbounds(0, n-1);
    t3.setbounds(0, n-1);
    result = true;
    
    //
    // Compute the Cholesky factorization A = U'*U.
    //
    for(j = 0; j <= n-1; j++)
    {
        
        //
        // Compute U(J,J) and test for non-positive-definiteness.
        //
        v = ap::vdotproduct(&a(0, j), a.getstride(), "Conj", &a(0, j), a.getstride(), "N", ap::vlen(0,j-1));
        ajj = (a(j,j)-v).x;
        if( ap::fp_less_eq(ajj,0) )
        {
            a(j,j) = ajj;
            result = false;
            return result;
        }
        ajj = sqrt(ajj);
        a(j,j) = ajj;
        
        //
        // Compute elements J+1:N-1 of row J.
        //
        if( j<n-1 )
        {
            ap::vmove(&t2(0), 1, &a(0, j), a.getstride(), "Conj", ap::vlen(0,j-1));
            ap::vmove(&t3(j+1), 1, &a(j, j+1), 1, "N", ap::vlen(j+1,n-1));
            for(i = j+1; i <= n-1; i++)
            {
                v = ap::vdotproduct(&a(0, i), a.getstride(), "N", &t2(0), 1, "N", ap::vlen(0,j-1));
                t3(i) = t3(i)-v;
            }
            ap::vmove(&a(j, j+1), 1, &t3(j+1), 1, "N", ap::vlen(j+1,n-1));
            r = 1/ajj;
            ap::vmul(&a(j, j+1), 1, ap::vlen(j+1,n-1), r);
        }
    }
    return result;
}


/*************************************************************************
SVD condition number
*************************************************************************/
static double svdcond(const ap::real_2d_array& a, int n)
{
    double result;
    ap::real_2d_array a1;
    ap::real_2d_array v;
    ap::real_1d_array w;
    int i;
    int j;
    double minw;
    double maxw;

    a1.setbounds(1, n, 1, n);
    for(i = 1; i <= n; i++)
    {
        for(j = 1; j <= n; j++)
        {
            a1(i,j) = a(i-1,j-1);
        }
    }
    if( !obsoletesvddecomposition(a1, n, n, w, v) )
    {
        result = 0;
        return result;
    }
    minw = w(1);
    maxw = w(1);
    for(i = 2; i <= n; i++)
    {
        if( ap::fp_less(w(i),minw) )
        {
            minw = w(i);
        }
        if( ap::fp_greater(w(i),maxw) )
        {
            maxw = w(i);
        }
    }
    result = maxw/minw;
    return result;
}


static double extsign(double a, double b)
{
    double result;

    if( ap::fp_greater_eq(b,0) )
    {
        result = fabs(a);
    }
    else
    {
        result = -fabs(a);
    }
    return result;
}


static double mymax(double a, double b)
{
    double result;

    if( ap::fp_greater(a,b) )
    {
        result = a;
    }
    else
    {
        result = b;
    }
    return result;
}


static double pythag(double a, double b)
{
    double result;

    if( ap::fp_less(fabs(a),fabs(b)) )
    {
        result = fabs(b)*sqrt(1+ap::sqr(a/b));
    }
    else
    {
        result = fabs(a)*sqrt(1+ap::sqr(b/a));
    }
    return result;
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testmatgenunit_test_silent()
{
    bool result;

    result = testmatgen(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testmatgenunit_test()
{
    bool result;

    result = testmatgen(false);
    return result;
}




