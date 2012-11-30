
 
#include <stdio.h>
#include "testtrfacunit.h"

static void testcluproblem(const ap::complex_2d_array& a,
     int m,
     int n,
     double threshold,
     bool& err,
     bool& properr);
static void testrluproblem(const ap::real_2d_array& a,
     int m,
     int n,
     double threshold,
     bool& err,
     bool& properr);

bool testtrfac(bool silent)
{
    bool result;
    ap::real_2d_array ra;
    ap::real_2d_array ral;
    ap::real_2d_array rau;
    ap::complex_2d_array ca;
    ap::complex_2d_array cal;
    ap::complex_2d_array cau;
    int m;
    int n;
    int mx;
    int maxmn;
    int i;
    int j;
    int minij;
    int pass;
    ap::complex vc;
    double vr;
    bool waserrors;
    bool spderr;
    bool hpderr;
    bool rerr;
    bool cerr;
    bool properr;
    double threshold;

    rerr = false;
    spderr = false;
    cerr = false;
    hpderr = false;
    properr = false;
    waserrors = false;
    maxmn = 4*ablasblocksize(ra)+1;
    threshold = 1000*ap::machineepsilon*maxmn;
    
    //
    // test LU
    //
    for(mx = 1; mx <= maxmn; mx++)
    {
        
        //
        // Initialize N/M, both are <=MX,
        // at least one of them is exactly equal to MX
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
        
        //
        // First, test on zero matrix
        //
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
        testcluproblem(ca, m, n, threshold, cerr, properr);
        testrluproblem(ra, m, n, threshold, rerr, properr);
        
        //
        // Second, random matrix with moderate condition number
        //
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
        for(i = 0; i <= ap::minint(m, n)-1; i++)
        {
            ra(i,i) = 1+10*ap::randomreal();
            ca(i,i) = 1+10*ap::randomreal();
        }
        cmatrixrndorthogonalfromtheleft(ca, m, n);
        cmatrixrndorthogonalfromtheright(ca, m, n);
        rmatrixrndorthogonalfromtheleft(ra, m, n);
        rmatrixrndorthogonalfromtheright(ra, m, n);
        testcluproblem(ca, m, n, threshold, cerr, properr);
        testrluproblem(ra, m, n, threshold, rerr, properr);
    }
    
    //
    // Test Cholesky
    //
    for(n = 1; n <= maxmn; n++)
    {
        
        //
        // Load CA (HPD matrix with low condition number),
        //      CAL and CAU - its lower and upper triangles
        //
        hpdmatrixrndcond(n, 1+50*ap::randomreal(), ca);
        cal.setlength(n, n);
        cau.setlength(n, n);
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                cal(i,j) = i;
                cau(i,j) = j;
            }
        }
        for(i = 0; i <= n-1; i++)
        {
            ap::vmove(&cal(i, 0), 1, &ca(i, 0), 1, "N", ap::vlen(0,i));
            ap::vmove(&cau(i, i), 1, &ca(i, i), 1, "N", ap::vlen(i,n-1));
        }
        
        //
        // Test HPDMatrixCholesky:
        // 1. it must leave upper (lower) part unchanged
        // 2. max(A-L*L^H) must be small
        //
        if( hpdmatrixcholesky(cal, n, false) )
        {
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    if( j>i )
                    {
                        hpderr = hpderr||cal(i,j)!=i;
                    }
                    else
                    {
                        vc = ap::vdotproduct(&cal(i, 0), 1, "N", &cal(j, 0), 1, "Conj", ap::vlen(0,j));
                        hpderr = hpderr||ap::fp_greater(ap::abscomplex(ca(i,j)-vc),threshold);
                    }
                }
            }
        }
        else
        {
            hpderr = true;
        }
        if( hpdmatrixcholesky(cau, n, true) )
        {
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    if( j<i )
                    {
                        hpderr = hpderr||cau(i,j)!=j;
                    }
                    else
                    {
                        vc = ap::vdotproduct(&cau(0, i), cau.getstride(), "Conj", &cau(0, j), cau.getstride(), "N", ap::vlen(0,i));
                        hpderr = hpderr||ap::fp_greater(ap::abscomplex(ca(i,j)-vc),threshold);
                    }
                }
            }
        }
        else
        {
            hpderr = true;
        }
        
        //
        // Load RA (SPD matrix with low condition number),
        //      RAL and RAU - its lower and upper triangles
        //
        spdmatrixrndcond(n, 1+50*ap::randomreal(), ra);
        ral.setlength(n, n);
        rau.setlength(n, n);
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                ral(i,j) = i;
                rau(i,j) = j;
            }
        }
        for(i = 0; i <= n-1; i++)
        {
            ap::vmove(&ral(i, 0), 1, &ra(i, 0), 1, ap::vlen(0,i));
            ap::vmove(&rau(i, i), 1, &ra(i, i), 1, ap::vlen(i,n-1));
        }
        
        //
        // Test SPDMatrixCholesky:
        // 1. it must leave upper (lower) part unchanged
        // 2. max(A-L*L^H) must be small
        //
        if( spdmatrixcholesky(ral, n, false) )
        {
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    if( j>i )
                    {
                        spderr = spderr||ap::fp_neq(ral(i,j),i);
                    }
                    else
                    {
                        vr = ap::vdotproduct(&ral(i, 0), 1, &ral(j, 0), 1, ap::vlen(0,j));
                        spderr = spderr||ap::fp_greater(fabs(ra(i,j)-vr),threshold);
                    }
                }
            }
        }
        else
        {
            spderr = true;
        }
        if( spdmatrixcholesky(rau, n, true) )
        {
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    if( j<i )
                    {
                        spderr = spderr||ap::fp_neq(rau(i,j),j);
                    }
                    else
                    {
                        vr = ap::vdotproduct(&rau(0, i), rau.getstride(), &rau(0, j), rau.getstride(), ap::vlen(0,i));
                        spderr = spderr||ap::fp_greater(fabs(ra(i,j)-vr),threshold);
                    }
                }
            }
        }
        else
        {
            spderr = true;
        }
    }
    
    //
    // report
    //
    waserrors = rerr||spderr||cerr||hpderr||properr;
    if( !silent )
    {
        printf("TESTING TRIANGULAR FACTORIZATIONS\n");
        printf("* REAL:                                  ");
        if( rerr )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* SPD:                                   ");
        if( spderr )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* COMPLEX:                               ");
        if( cerr )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* HPD:                                   ");
        if( hpderr )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* OTHER PROPERTIES:                      ");
        if( properr )
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
        printf("\n\n");
    }
    result = !waserrors;
    return result;
}


static void testcluproblem(const ap::complex_2d_array& a,
     int m,
     int n,
     double threshold,
     bool& err,
     bool& properr)
{
    ap::complex_2d_array ca;
    ap::complex_2d_array cl;
    ap::complex_2d_array cu;
    ap::complex_2d_array ca2;
    ap::complex_1d_array ct;
    int i;
    int j;
    int minmn;
    ap::complex v;
    ap::integer_1d_array p;

    minmn = ap::minint(m, n);
    
    //
    // PLU test
    //
    ca.setlength(m, n);
    for(i = 0; i <= m-1; i++)
    {
        ap::vmove(&ca(i, 0), 1, &a(i, 0), 1, "N", ap::vlen(0,n-1));
    }
    cmatrixplu(ca, m, n, p);
    for(i = 0; i <= minmn-1; i++)
    {
        if( p(i)<i||p(i)>=m )
        {
            properr = false;
            return;
        }
    }
    cl.setlength(m, minmn);
    for(j = 0; j <= minmn-1; j++)
    {
        for(i = 0; i <= j-1; i++)
        {
            cl(i,j) = 0.0;
        }
        cl(j,j) = 1.0;
        for(i = j+1; i <= m-1; i++)
        {
            cl(i,j) = ca(i,j);
        }
    }
    cu.setlength(minmn, n);
    for(i = 0; i <= minmn-1; i++)
    {
        for(j = 0; j <= i-1; j++)
        {
            cu(i,j) = 0.0;
        }
        for(j = i; j <= n-1; j++)
        {
            cu(i,j) = ca(i,j);
        }
    }
    ca2.setlength(m, n);
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&cl(i, 0), 1, "N", &cu(0, j), cu.getstride(), "N", ap::vlen(0,minmn-1));
            ca2(i,j) = v;
        }
    }
    ct.setlength(n);
    for(i = minmn-1; i >= 0; i--)
    {
        if( i!=p(i) )
        {
            ap::vmove(&ct(0), 1, &ca2(i, 0), 1, "N", ap::vlen(0,n-1));
            ap::vmove(&ca2(i, 0), 1, &ca2(p(i), 0), 1, "N", ap::vlen(0,n-1));
            ap::vmove(&ca2(p(i), 0), 1, &ct(0), 1, "N", ap::vlen(0,n-1));
        }
    }
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            err = err||ap::fp_greater(ap::abscomplex(a(i,j)-ca2(i,j)),threshold);
        }
    }
    
    //
    // LUP test
    //
    ca.setlength(m, n);
    for(i = 0; i <= m-1; i++)
    {
        ap::vmove(&ca(i, 0), 1, &a(i, 0), 1, "N", ap::vlen(0,n-1));
    }
    cmatrixlup(ca, m, n, p);
    for(i = 0; i <= minmn-1; i++)
    {
        if( p(i)<i||p(i)>=n )
        {
            properr = false;
            return;
        }
    }
    cl.setlength(m, minmn);
    for(j = 0; j <= minmn-1; j++)
    {
        for(i = 0; i <= j-1; i++)
        {
            cl(i,j) = 0.0;
        }
        for(i = j; i <= m-1; i++)
        {
            cl(i,j) = ca(i,j);
        }
    }
    cu.setlength(minmn, n);
    for(i = 0; i <= minmn-1; i++)
    {
        for(j = 0; j <= i-1; j++)
        {
            cu(i,j) = 0.0;
        }
        cu(i,i) = 1.0;
        for(j = i+1; j <= n-1; j++)
        {
            cu(i,j) = ca(i,j);
        }
    }
    ca2.setlength(m, n);
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&cl(i, 0), 1, "N", &cu(0, j), cu.getstride(), "N", ap::vlen(0,minmn-1));
            ca2(i,j) = v;
        }
    }
    ct.setlength(m);
    for(i = minmn-1; i >= 0; i--)
    {
        if( i!=p(i) )
        {
            ap::vmove(&ct(0), 1, &ca2(0, i), ca2.getstride(), "N", ap::vlen(0,m-1));
            ap::vmove(&ca2(0, i), ca2.getstride(), &ca2(0, p(i)), ca2.getstride(), "N", ap::vlen(0,m-1));
            ap::vmove(&ca2(0, p(i)), ca2.getstride(), &ct(0), 1, "N", ap::vlen(0,m-1));
        }
    }
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            err = err||ap::fp_greater(ap::abscomplex(a(i,j)-ca2(i,j)),threshold);
        }
    }
}


static void testrluproblem(const ap::real_2d_array& a,
     int m,
     int n,
     double threshold,
     bool& err,
     bool& properr)
{
    ap::real_2d_array ca;
    ap::real_2d_array cl;
    ap::real_2d_array cu;
    ap::real_2d_array ca2;
    ap::real_1d_array ct;
    int i;
    int j;
    int minmn;
    double v;
    ap::integer_1d_array p;

    minmn = ap::minint(m, n);
    
    //
    // PLU test
    //
    ca.setlength(m, n);
    for(i = 0; i <= m-1; i++)
    {
        ap::vmove(&ca(i, 0), 1, &a(i, 0), 1, ap::vlen(0,n-1));
    }
    rmatrixplu(ca, m, n, p);
    for(i = 0; i <= minmn-1; i++)
    {
        if( p(i)<i||p(i)>=m )
        {
            properr = false;
            return;
        }
    }
    cl.setlength(m, minmn);
    for(j = 0; j <= minmn-1; j++)
    {
        for(i = 0; i <= j-1; i++)
        {
            cl(i,j) = 0.0;
        }
        cl(j,j) = 1.0;
        for(i = j+1; i <= m-1; i++)
        {
            cl(i,j) = ca(i,j);
        }
    }
    cu.setlength(minmn, n);
    for(i = 0; i <= minmn-1; i++)
    {
        for(j = 0; j <= i-1; j++)
        {
            cu(i,j) = 0.0;
        }
        for(j = i; j <= n-1; j++)
        {
            cu(i,j) = ca(i,j);
        }
    }
    ca2.setlength(m, n);
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&cl(i, 0), 1, &cu(0, j), cu.getstride(), ap::vlen(0,minmn-1));
            ca2(i,j) = v;
        }
    }
    ct.setlength(n);
    for(i = minmn-1; i >= 0; i--)
    {
        if( i!=p(i) )
        {
            ap::vmove(&ct(0), 1, &ca2(i, 0), 1, ap::vlen(0,n-1));
            ap::vmove(&ca2(i, 0), 1, &ca2(p(i), 0), 1, ap::vlen(0,n-1));
            ap::vmove(&ca2(p(i), 0), 1, &ct(0), 1, ap::vlen(0,n-1));
        }
    }
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            err = err||ap::fp_greater(fabs(a(i,j)-ca2(i,j)),threshold);
        }
    }
    
    //
    // LUP test
    //
    ca.setlength(m, n);
    for(i = 0; i <= m-1; i++)
    {
        ap::vmove(&ca(i, 0), 1, &a(i, 0), 1, ap::vlen(0,n-1));
    }
    rmatrixlup(ca, m, n, p);
    for(i = 0; i <= minmn-1; i++)
    {
        if( p(i)<i||p(i)>=n )
        {
            properr = false;
            return;
        }
    }
    cl.setlength(m, minmn);
    for(j = 0; j <= minmn-1; j++)
    {
        for(i = 0; i <= j-1; i++)
        {
            cl(i,j) = 0.0;
        }
        for(i = j; i <= m-1; i++)
        {
            cl(i,j) = ca(i,j);
        }
    }
    cu.setlength(minmn, n);
    for(i = 0; i <= minmn-1; i++)
    {
        for(j = 0; j <= i-1; j++)
        {
            cu(i,j) = 0.0;
        }
        cu(i,i) = 1.0;
        for(j = i+1; j <= n-1; j++)
        {
            cu(i,j) = ca(i,j);
        }
    }
    ca2.setlength(m, n);
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            v = ap::vdotproduct(&cl(i, 0), 1, &cu(0, j), cu.getstride(), ap::vlen(0,minmn-1));
            ca2(i,j) = v;
        }
    }
    ct.setlength(m);
    for(i = minmn-1; i >= 0; i--)
    {
        if( i!=p(i) )
        {
            ap::vmove(&ct(0), 1, &ca2(0, i), ca2.getstride(), ap::vlen(0,m-1));
            ap::vmove(&ca2(0, i), ca2.getstride(), &ca2(0, p(i)), ca2.getstride(), ap::vlen(0,m-1));
            ap::vmove(&ca2(0, p(i)), ca2.getstride(), &ct(0), 1, ap::vlen(0,m-1));
        }
    }
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            err = err||ap::fp_greater(fabs(a(i,j)-ca2(i,j)),threshold);
        }
    }
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testtrfacunit_test_silent()
{
    bool result;

    result = testtrfac(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testtrfacunit_test()
{
    bool result;

    result = testtrfac(false);
    return result;
}




