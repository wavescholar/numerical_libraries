
 
#include <stdio.h>
#include "testmatinvunit.h"

static void rmatrixmakeacopy(const ap::real_2d_array& a,
     int m,
     int n,
     ap::real_2d_array& b);
static void cmatrixmakeacopy(const ap::complex_2d_array& a,
     int m,
     int n,
     ap::complex_2d_array& b);
static bool rmatrixcheckinverse(const ap::real_2d_array& a,
     const ap::real_2d_array& inva,
     int n,
     double threshold,
     int info,
     const matinvreport& rep);
static bool spdmatrixcheckinverse(ap::real_2d_array a,
     ap::real_2d_array inva,
     bool isupper,
     int n,
     double threshold,
     int info,
     const matinvreport& rep);
static bool hpdmatrixcheckinverse(ap::complex_2d_array a,
     ap::complex_2d_array inva,
     bool isupper,
     int n,
     double threshold,
     int info,
     const matinvreport& rep);
static bool rmatrixcheckinversesingular(const ap::real_2d_array& inva,
     int n,
     double threshold,
     int info,
     const matinvreport& rep);
static bool cmatrixcheckinverse(const ap::complex_2d_array& a,
     const ap::complex_2d_array& inva,
     int n,
     double threshold,
     int info,
     const matinvreport& rep);
static bool cmatrixcheckinversesingular(const ap::complex_2d_array& inva,
     int n,
     double threshold,
     int info,
     const matinvreport& rep);
static void rmatrixdrophalf(ap::real_2d_array& a, int n, bool droplower);
static void cmatrixdrophalf(ap::complex_2d_array& a, int n, bool droplower);
static void testrtrinv(int maxn,
     int passcount,
     double threshold,
     bool& rtrerrors);
static void testctrinv(int maxn,
     int passcount,
     double threshold,
     bool& ctrerrors);
static void testrinv(int maxn, int passcount, double threshold, bool& rerrors);
static void testcinv(int maxn, int passcount, double threshold, bool& cerrors);
static void testspdinv(int maxn,
     int passcount,
     double threshold,
     bool& spderrors);
static void testhpdinv(int maxn,
     int passcount,
     double threshold,
     bool& hpderrors);
static void unset2d(ap::real_2d_array& x);
static void unset1d(ap::real_1d_array& x);
static void cunset2d(ap::complex_2d_array& x);
static void cunset1d(ap::complex_1d_array& x);
static void unsetrep(matinvreport& r);

/*************************************************************************
Test
*************************************************************************/
bool testmatinv(bool silent)
{
    bool result;
    int maxrn;
    int maxcn;
    int passcount;
    double threshold;
    double rcondtol;
    bool rtrerrors;
    bool ctrerrors;
    bool rerrors;
    bool cerrors;
    bool spderrors;
    bool hpderrors;
    bool waserrors;
    ap::real_2d_array emptyra;
    ap::real_2d_array emptyca;

    maxrn = 3*ablasblocksize(emptyra)+1;
    maxcn = 3*ablasblocksize(emptyca)+1;
    passcount = 1;
    threshold = 10000*ap::machineepsilon;
    rcondtol = 0.01;
    rtrerrors = false;
    ctrerrors = false;
    rerrors = false;
    cerrors = false;
    spderrors = false;
    hpderrors = false;
    testrtrinv(maxrn, passcount, threshold, rtrerrors);
    testctrinv(maxcn, passcount, threshold, ctrerrors);
    testrinv(maxrn, passcount, threshold, rerrors);
    testspdinv(maxrn, passcount, threshold, spderrors);
    testcinv(maxcn, passcount, threshold, cerrors);
    testhpdinv(maxcn, passcount, threshold, hpderrors);
    waserrors = rtrerrors||ctrerrors||rerrors||cerrors||spderrors||hpderrors;
    if( !silent )
    {
        printf("TESTING MATINV\n");
        printf("* REAL TRIANGULAR:                        ");
        if( rtrerrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* COMPLEX TRIANGULAR:                     ");
        if( ctrerrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
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
Checks whether inverse is correct
Returns True on success.
*************************************************************************/
static bool rmatrixcheckinverse(const ap::real_2d_array& a,
     const ap::real_2d_array& inva,
     int n,
     double threshold,
     int info,
     const matinvreport& rep)
{
    bool result;
    int i;
    int j;
    double v;

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
            for(j = 0; j <= n-1; j++)
            {
                v = ap::vdotproduct(&a(i, 0), 1, &inva(0, j), inva.getstride(), ap::vlen(0,n-1));
                if( i==j )
                {
                    v = v-1;
                }
                result = result&&ap::fp_less_eq(fabs(v),threshold);
            }
        }
    }
    return result;
}


/*************************************************************************
Checks whether inverse is correct
Returns True on success.
*************************************************************************/
static bool spdmatrixcheckinverse(ap::real_2d_array a,
     ap::real_2d_array inva,
     bool isupper,
     int n,
     double threshold,
     int info,
     const matinvreport& rep)
{
    bool result;
    int i;
    int j;
    double v;

    for(i = 0; i <= n-2; i++)
    {
        if( isupper )
        {
            ap::vmove(&a(i+1, i), a.getstride(), &a(i, i+1), 1, ap::vlen(i+1,n-1));
            ap::vmove(&inva(i+1, i), inva.getstride(), &inva(i, i+1), 1, ap::vlen(i+1,n-1));
        }
        else
        {
            ap::vmove(&a(i, i+1), 1, &a(i+1, i), a.getstride(), ap::vlen(i+1,n-1));
            ap::vmove(&inva(i, i+1), 1, &inva(i+1, i), inva.getstride(), ap::vlen(i+1,n-1));
        }
    }
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
            for(j = 0; j <= n-1; j++)
            {
                v = ap::vdotproduct(&a(i, 0), 1, &inva(0, j), inva.getstride(), ap::vlen(0,n-1));
                if( i==j )
                {
                    v = v-1;
                }
                result = result&&ap::fp_less_eq(fabs(v),threshold);
            }
        }
    }
    return result;
}


/*************************************************************************
Checks whether inverse is correct
Returns True on success.
*************************************************************************/
static bool hpdmatrixcheckinverse(ap::complex_2d_array a,
     ap::complex_2d_array inva,
     bool isupper,
     int n,
     double threshold,
     int info,
     const matinvreport& rep)
{
    bool result;
    int i;
    int j;
    ap::complex v;

    for(i = 0; i <= n-2; i++)
    {
        if( isupper )
        {
            ap::vmove(&a(i+1, i), a.getstride(), &a(i, i+1), 1, "Conj", ap::vlen(i+1,n-1));
            ap::vmove(&inva(i+1, i), inva.getstride(), &inva(i, i+1), 1, "Conj", ap::vlen(i+1,n-1));
        }
        else
        {
            ap::vmove(&a(i, i+1), 1, &a(i+1, i), a.getstride(), "Conj", ap::vlen(i+1,n-1));
            ap::vmove(&inva(i, i+1), 1, &inva(i+1, i), inva.getstride(), "Conj", ap::vlen(i+1,n-1));
        }
    }
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
            for(j = 0; j <= n-1; j++)
            {
                v = ap::vdotproduct(&a(i, 0), 1, "N", &inva(0, j), inva.getstride(), "N", ap::vlen(0,n-1));
                if( i==j )
                {
                    v = v-1;
                }
                result = result&&ap::fp_less_eq(ap::abscomplex(v),threshold);
            }
        }
    }
    return result;
}


/*************************************************************************
Checks whether inversion result indicate singular matrix
Returns True on success.
*************************************************************************/
static bool rmatrixcheckinversesingular(const ap::real_2d_array& inva,
     int n,
     double threshold,
     int info,
     const matinvreport& rep)
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
                for(j = 0; j <= n-1; j++)
                {
                    result = result&&ap::fp_eq(inva(i,j),0);
                }
            }
        }
    }
    return result;
}


/*************************************************************************
Checks whether inverse is correct
Returns True on success.
*************************************************************************/
static bool cmatrixcheckinverse(const ap::complex_2d_array& a,
     const ap::complex_2d_array& inva,
     int n,
     double threshold,
     int info,
     const matinvreport& rep)
{
    bool result;
    int i;
    int j;
    ap::complex v;

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
            for(j = 0; j <= n-1; j++)
            {
                v = ap::vdotproduct(&a(i, 0), 1, "N", &inva(0, j), inva.getstride(), "N", ap::vlen(0,n-1));
                if( i==j )
                {
                    v = v-1;
                }
                result = result&&ap::fp_less_eq(ap::abscomplex(v),threshold);
            }
        }
    }
    return result;
}


/*************************************************************************
Checks whether inversion result indicate singular matrix
Returns True on success.
*************************************************************************/
static bool cmatrixcheckinversesingular(const ap::complex_2d_array& inva,
     int n,
     double threshold,
     int info,
     const matinvreport& rep)
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
                for(j = 0; j <= n-1; j++)
                {
                    result = result&&inva(i,j)==0;
                }
            }
        }
    }
    return result;
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
Real TR inverse
*************************************************************************/
static void testrtrinv(int maxn,
     int passcount,
     double threshold,
     bool& rtrerrors)
{
    ap::real_2d_array a;
    ap::real_2d_array b;
    int n;
    int pass;
    int i;
    int j;
    int task;
    bool isupper;
    bool isunit;
    double v;
    bool waserrors;
    int info;
    matinvreport rep;

    waserrors = false;
    
    //
    // Test
    //
    for(n = 1; n <= maxn; n++)
    {
        a.setlength(n, n);
        b.setlength(n, n);
        for(task = 0; task <= 3; task++)
        {
            for(pass = 1; pass <= passcount; pass++)
            {
                
                //
                // Determine task
                //
                isupper = task%2==0;
                isunit = task/2%2==0;
                
                //
                // Generate matrix
                //
                for(i = 0; i <= n-1; i++)
                {
                    for(j = 0; j <= n-1; j++)
                    {
                        if( i==j )
                        {
                            a(i,i) = 1+ap::randomreal();
                        }
                        else
                        {
                            a(i,j) = 0.2*ap::randomreal()-0.1;
                        }
                        b(i,j) = a(i,j);
                    }
                }
                
                //
                // Inverse
                //
                rmatrixtrinverse(b, n, isupper, isunit, info, rep);
                if( info<=0 )
                {
                    rtrerrors = true;
                    return;
                }
                
                //
                // Structural test
                //
                if( isunit )
                {
                    for(i = 0; i <= n-1; i++)
                    {
                        rtrerrors = rtrerrors||ap::fp_neq(a(i,i),b(i,i));
                    }
                }
                if( isupper )
                {
                    for(i = 0; i <= n-1; i++)
                    {
                        for(j = 0; j <= i-1; j++)
                        {
                            rtrerrors = rtrerrors||ap::fp_neq(a(i,j),b(i,j));
                        }
                    }
                }
                else
                {
                    for(i = 0; i <= n-1; i++)
                    {
                        for(j = i+1; j <= n-1; j++)
                        {
                            rtrerrors = rtrerrors||ap::fp_neq(a(i,j),b(i,j));
                        }
                    }
                }
                
                //
                // Inverse test
                //
                for(i = 0; i <= n-1; i++)
                {
                    for(j = 0; j <= n-1; j++)
                    {
                        if( j<i&&isupper||j>i&&!isupper )
                        {
                            a(i,j) = 0;
                            b(i,j) = 0;
                        }
                    }
                }
                if( isunit )
                {
                    for(i = 0; i <= n-1; i++)
                    {
                        a(i,i) = 1;
                        b(i,i) = 1;
                    }
                }
                for(i = 0; i <= n-1; i++)
                {
                    for(j = 0; j <= n-1; j++)
                    {
                        v = ap::vdotproduct(&a(i, 0), 1, &b(0, j), b.getstride(), ap::vlen(0,n-1));
                        if( j!=i )
                        {
                            rtrerrors = rtrerrors||ap::fp_greater(fabs(v),threshold);
                        }
                        else
                        {
                            rtrerrors = rtrerrors||ap::fp_greater(fabs(v-1),threshold);
                        }
                    }
                }
            }
        }
    }
}


/*************************************************************************
Complex TR inverse
*************************************************************************/
static void testctrinv(int maxn,
     int passcount,
     double threshold,
     bool& ctrerrors)
{
    ap::complex_2d_array a;
    ap::complex_2d_array b;
    int n;
    int pass;
    int i;
    int j;
    int task;
    bool isupper;
    bool isunit;
    ap::complex v;
    bool waserrors;
    int info;
    matinvreport rep;

    waserrors = false;
    
    //
    // Test
    //
    for(n = 1; n <= maxn; n++)
    {
        a.setlength(n, n);
        b.setlength(n, n);
        for(task = 0; task <= 3; task++)
        {
            for(pass = 1; pass <= passcount; pass++)
            {
                
                //
                // Determine task
                //
                isupper = task%2==0;
                isunit = task/2%2==0;
                
                //
                // Generate matrix
                //
                for(i = 0; i <= n-1; i++)
                {
                    for(j = 0; j <= n-1; j++)
                    {
                        if( i==j )
                        {
                            a(i,i).x = 1+ap::randomreal();
                            a(i,i).y = 1+ap::randomreal();
                        }
                        else
                        {
                            a(i,j).x = 0.2*ap::randomreal()-0.1;
                            a(i,j).y = 0.2*ap::randomreal()-0.1;
                        }
                        b(i,j) = a(i,j);
                    }
                }
                
                //
                // Inverse
                //
                cmatrixtrinverse(b, n, isupper, isunit, info, rep);
                if( info<=0 )
                {
                    ctrerrors = true;
                    return;
                }
                
                //
                // Structural test
                //
                if( isunit )
                {
                    for(i = 0; i <= n-1; i++)
                    {
                        ctrerrors = ctrerrors||a(i,i)!=b(i,i);
                    }
                }
                if( isupper )
                {
                    for(i = 0; i <= n-1; i++)
                    {
                        for(j = 0; j <= i-1; j++)
                        {
                            ctrerrors = ctrerrors||a(i,j)!=b(i,j);
                        }
                    }
                }
                else
                {
                    for(i = 0; i <= n-1; i++)
                    {
                        for(j = i+1; j <= n-1; j++)
                        {
                            ctrerrors = ctrerrors||a(i,j)!=b(i,j);
                        }
                    }
                }
                
                //
                // Inverse test
                //
                for(i = 0; i <= n-1; i++)
                {
                    for(j = 0; j <= n-1; j++)
                    {
                        if( j<i&&isupper||j>i&&!isupper )
                        {
                            a(i,j) = 0;
                            b(i,j) = 0;
                        }
                    }
                }
                if( isunit )
                {
                    for(i = 0; i <= n-1; i++)
                    {
                        a(i,i) = 1;
                        b(i,i) = 1;
                    }
                }
                for(i = 0; i <= n-1; i++)
                {
                    for(j = 0; j <= n-1; j++)
                    {
                        v = ap::vdotproduct(&a(i, 0), 1, "N", &b(0, j), b.getstride(), "N", ap::vlen(0,n-1));
                        if( j!=i )
                        {
                            ctrerrors = ctrerrors||ap::fp_greater(ap::abscomplex(v),threshold);
                        }
                        else
                        {
                            ctrerrors = ctrerrors||ap::fp_greater(ap::abscomplex(v-1),threshold);
                        }
                    }
                }
            }
        }
    }
}


/*************************************************************************
Real test
*************************************************************************/
static void testrinv(int maxn,
     int passcount,
     double threshold,
     bool& rerrors)
{
    ap::real_2d_array a;
    ap::real_2d_array lua;
    ap::real_2d_array inva;
    ap::real_2d_array invlua;
    ap::integer_1d_array p;
    int i;
    int j;
    int k;
    int n;
    int pass;
    int taskkind;
    double v;
    int info;
    matinvreport rep;

    
    //
    // General square matrices:
    // * test general solvers
    // * test least squares solver
    //
    for(pass = 1; pass <= passcount; pass++)
    {
        for(n = 1; n <= maxn; n++)
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
            rmatrixmakeacopy(a, n, n, inva);
            rmatrixmakeacopy(lua, n, n, invlua);
            info = 0;
            unsetrep(rep);
            rmatrixinverse(inva, n, info, rep);
            rerrors = rerrors||!rmatrixcheckinverse(a, inva, n, threshold, info, rep);
            info = 0;
            unsetrep(rep);
            rmatrixluinverse(invlua, p, n, info, rep);
            rerrors = rerrors||!rmatrixcheckinverse(a, invlua, n, threshold, info, rep);
            
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
            // 2. test different methods
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
                rmatrixmakeacopy(a, n, n, lua);
                rmatrixlu(lua, n, n, p);
                info = 0;
                unsetrep(rep);
                rmatrixinverse(a, n, info, rep);
                rerrors = rerrors||!rmatrixcheckinversesingular(a, n, threshold, info, rep);
                info = 0;
                unsetrep(rep);
                rmatrixluinverse(lua, p, n, info, rep);
                rerrors = rerrors||!rmatrixcheckinversesingular(lua, n, threshold, info, rep);
            }
        }
    }
}


/*************************************************************************
Complex test
*************************************************************************/
static void testcinv(int maxn,
     int passcount,
     double threshold,
     bool& cerrors)
{
    ap::complex_2d_array a;
    ap::complex_2d_array lua;
    ap::complex_2d_array inva;
    ap::complex_2d_array invlua;
    ap::integer_1d_array p;
    int i;
    int j;
    int k;
    int n;
    int pass;
    int taskkind;
    double v;
    int info;
    matinvreport rep;

    
    //
    // General square matrices:
    // * test general solvers
    // * test least squares solver
    //
    for(pass = 1; pass <= passcount; pass++)
    {
        for(n = 1; n <= maxn; n++)
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
            cmatrixmakeacopy(a, n, n, inva);
            cmatrixmakeacopy(lua, n, n, invlua);
            info = 0;
            unsetrep(rep);
            cmatrixinverse(inva, n, info, rep);
            cerrors = cerrors||!cmatrixcheckinverse(a, inva, n, threshold, info, rep);
            info = 0;
            unsetrep(rep);
            cmatrixluinverse(invlua, p, n, info, rep);
            cerrors = cerrors||!cmatrixcheckinverse(a, invlua, n, threshold, info, rep);
            
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
            // 2. test different methods
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
                cmatrixmakeacopy(a, n, n, lua);
                cmatrixlu(lua, n, n, p);
                info = 0;
                unsetrep(rep);
                cmatrixinverse(a, n, info, rep);
                cerrors = cerrors||!cmatrixcheckinversesingular(a, n, threshold, info, rep);
                info = 0;
                unsetrep(rep);
                cmatrixluinverse(lua, p, n, info, rep);
                cerrors = cerrors||!cmatrixcheckinversesingular(lua, n, threshold, info, rep);
            }
        }
    }
}


/*************************************************************************
SPD test
*************************************************************************/
static void testspdinv(int maxn,
     int passcount,
     double threshold,
     bool& spderrors)
{
    ap::real_2d_array a;
    ap::real_2d_array cha;
    ap::real_2d_array inva;
    ap::real_2d_array invcha;
    bool isupper;
    int i;
    int j;
    int k;
    int n;
    int pass;
    int taskkind;
    double v;
    int info;
    matinvreport rep;

    
    //
    // General square matrices:
    // * test general solvers
    // * test least squares solver
    //
    for(pass = 1; pass <= passcount; pass++)
    {
        for(n = 1; n <= maxn; n++)
        {
            isupper = ap::fp_greater(ap::randomreal(),0.5);
            
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
            spdmatrixrndcond(n, double(1000), a);
            rmatrixdrophalf(a, n, isupper);
            rmatrixmakeacopy(a, n, n, cha);
            if( !spdmatrixcholesky(cha, n, isupper) )
            {
                continue;
            }
            rmatrixmakeacopy(a, n, n, inva);
            rmatrixmakeacopy(cha, n, n, invcha);
            info = 0;
            unsetrep(rep);
            spdmatrixinverse(inva, n, isupper, info, rep);
            spderrors = spderrors||!spdmatrixcheckinverse(a, inva, isupper, n, threshold, info, rep);
            info = 0;
            unsetrep(rep);
            spdmatrixcholeskyinverse(invcha, n, isupper, info, rep);
            spderrors = spderrors||!spdmatrixcheckinverse(a, invcha, isupper, n, threshold, info, rep);
            
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
            // 2. test different methods
            //
            for(taskkind = 0; taskkind <= 2; taskkind++)
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
                info = 0;
                unsetrep(rep);
                spdmatrixcholeskyinverse(a, n, isupper, info, rep);
                if( info!=-3&&info!=1 )
                {
                    spderrors = true;
                }
                else
                {
                    spderrors = spderrors||ap::fp_less(rep.r1,0)||ap::fp_greater(rep.r1,1000*ap::machineepsilon);
                    spderrors = spderrors||ap::fp_less(rep.rinf,0)||ap::fp_greater(rep.rinf,1000*ap::machineepsilon);
                }
            }
        }
    }
}


/*************************************************************************
HPD test
*************************************************************************/
static void testhpdinv(int maxn,
     int passcount,
     double threshold,
     bool& hpderrors)
{
    ap::complex_2d_array a;
    ap::complex_2d_array cha;
    ap::complex_2d_array inva;
    ap::complex_2d_array invcha;
    bool isupper;
    int i;
    int j;
    int k;
    int n;
    int pass;
    int taskkind;
    ap::complex v;
    int info;
    matinvreport rep;

    
    //
    // General square matrices:
    // * test general solvers
    // * test least squares solver
    //
    for(pass = 1; pass <= passcount; pass++)
    {
        for(n = 1; n <= maxn; n++)
        {
            isupper = ap::fp_greater(ap::randomreal(),0.5);
            
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
            hpdmatrixrndcond(n, double(1000), a);
            cmatrixdrophalf(a, n, isupper);
            cmatrixmakeacopy(a, n, n, cha);
            if( !hpdmatrixcholesky(cha, n, isupper) )
            {
                continue;
            }
            cmatrixmakeacopy(a, n, n, inva);
            cmatrixmakeacopy(cha, n, n, invcha);
            info = 0;
            unsetrep(rep);
            hpdmatrixinverse(inva, n, isupper, info, rep);
            hpderrors = hpderrors||!hpdmatrixcheckinverse(a, inva, isupper, n, threshold, info, rep);
            info = 0;
            unsetrep(rep);
            hpdmatrixcholeskyinverse(invcha, n, isupper, info, rep);
            hpderrors = hpderrors||!hpdmatrixcheckinverse(a, invcha, isupper, n, threshold, info, rep);
            
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
            // 2. test different methods
            //
            for(taskkind = 0; taskkind <= 2; taskkind++)
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
                        for(j = 0; j <= n-1; j++)
                        {
                            a(i,j).x = 2*ap::randomreal()-1;
                            a(i,j).y = 2*ap::randomreal()-1;
                        }
                    }
                    k = ap::randominteger(n);
                    ap::vmul(&a(k, 0), 1, ap::vlen(0,n-1), 0);
                    ap::vmul(&a(0, k), a.getstride(), ap::vlen(0,n-1), 0);
                }
                info = 0;
                unsetrep(rep);
                hpdmatrixcholeskyinverse(a, n, isupper, info, rep);
                if( info!=-3&&info!=1 )
                {
                    hpderrors = true;
                }
                else
                {
                    hpderrors = hpderrors||ap::fp_less(rep.r1,0)||ap::fp_greater(rep.r1,1000*ap::machineepsilon);
                    hpderrors = hpderrors||ap::fp_less(rep.rinf,0)||ap::fp_greater(rep.rinf,1000*ap::machineepsilon);
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
static void unsetrep(matinvreport& r)
{

    r.r1 = -1;
    r.rinf = -1;
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testmatinvunit_test_silent()
{
    bool result;

    result = testmatinv(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testmatinvunit_test()
{
    bool result;

    result = testmatinv(false);
    return result;
}




