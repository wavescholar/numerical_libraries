
 
#include <stdio.h>
#include "testinverseupdateunit.h"

static void makeacopy(const ap::real_2d_array& a,
     int m,
     int n,
     ap::real_2d_array& b);
static void matlu(ap::real_2d_array& a,
     int m,
     int n,
     ap::integer_1d_array& pivots);
static void generaterandomorthogonalmatrix(ap::real_2d_array& a0, int n);
static void generaterandommatrixcond(ap::real_2d_array& a0, int n, double c);
static bool invmattr(ap::real_2d_array& a,
     int n,
     bool isupper,
     bool isunittriangular);
static bool invmatlu(ap::real_2d_array& a,
     const ap::integer_1d_array& pivots,
     int n);
static bool invmat(ap::real_2d_array& a, int n);
static double matrixdiff(const ap::real_2d_array& a,
     const ap::real_2d_array& b,
     int m,
     int n);
static bool updandinv(ap::real_2d_array& a,
     const ap::real_1d_array& u,
     const ap::real_1d_array& v,
     int n);

bool testinverseupdate(bool silent)
{
    bool result;
    ap::real_2d_array a;
    ap::real_2d_array inva;
    ap::real_2d_array b1;
    ap::real_2d_array b2;
    ap::real_1d_array u;
    ap::real_1d_array v;
    int n;
    int maxn;
    int i;
    int j;
    int updrow;
    int updcol;
    double val;
    int pass;
    int passcount;
    bool waserrors;
    double threshold;
    double c;

    waserrors = false;
    maxn = 10;
    passcount = 100;
    threshold = 1.0E-6;
    
    //
    // process
    //
    for(n = 1; n <= maxn; n++)
    {
        a.setbounds(0, n-1, 0, n-1);
        b1.setbounds(0, n-1, 0, n-1);
        b2.setbounds(0, n-1, 0, n-1);
        u.setbounds(0, n-1);
        v.setbounds(0, n-1);
        for(pass = 1; pass <= passcount; pass++)
        {
            c = exp(ap::randomreal()*log(double(10)));
            generaterandommatrixcond(a, n, c);
            makeacopy(a, n, n, inva);
            if( !invmat(inva, n) )
            {
                waserrors = true;
                break;
            }
            
            //
            // Test simple update
            //
            updrow = ap::randominteger(n);
            updcol = ap::randominteger(n);
            val = 0.1*(2*ap::randomreal()-1);
            for(i = 0; i <= n-1; i++)
            {
                if( i==updrow )
                {
                    u(i) = val;
                }
                else
                {
                    u(i) = 0;
                }
                if( i==updcol )
                {
                    v(i) = 1;
                }
                else
                {
                    v(i) = 0;
                }
            }
            makeacopy(a, n, n, b1);
            if( !updandinv(b1, u, v, n) )
            {
                waserrors = true;
                break;
            }
            makeacopy(inva, n, n, b2);
            rmatrixinvupdatesimple(b2, n, updrow, updcol, val);
            waserrors = waserrors||ap::fp_greater(matrixdiff(b1, b2, n, n),threshold);
            
            //
            // Test row update
            //
            updrow = ap::randominteger(n);
            for(i = 0; i <= n-1; i++)
            {
                if( i==updrow )
                {
                    u(i) = 1;
                }
                else
                {
                    u(i) = 0;
                }
                v(i) = 0.1*(2*ap::randomreal()-1);
            }
            makeacopy(a, n, n, b1);
            if( !updandinv(b1, u, v, n) )
            {
                waserrors = true;
                break;
            }
            makeacopy(inva, n, n, b2);
            rmatrixinvupdaterow(b2, n, updrow, v);
            waserrors = waserrors||ap::fp_greater(matrixdiff(b1, b2, n, n),threshold);
            
            //
            // Test column update
            //
            updcol = ap::randominteger(n);
            for(i = 0; i <= n-1; i++)
            {
                if( i==updcol )
                {
                    v(i) = 1;
                }
                else
                {
                    v(i) = 0;
                }
                u(i) = 0.1*(2*ap::randomreal()-1);
            }
            makeacopy(a, n, n, b1);
            if( !updandinv(b1, u, v, n) )
            {
                waserrors = true;
                break;
            }
            makeacopy(inva, n, n, b2);
            rmatrixinvupdatecolumn(b2, n, updcol, u);
            waserrors = waserrors||ap::fp_greater(matrixdiff(b1, b2, n, n),threshold);
            
            //
            // Test full update
            //
            for(i = 0; i <= n-1; i++)
            {
                v(i) = 0.1*(2*ap::randomreal()-1);
                u(i) = 0.1*(2*ap::randomreal()-1);
            }
            makeacopy(a, n, n, b1);
            if( !updandinv(b1, u, v, n) )
            {
                waserrors = true;
                break;
            }
            makeacopy(inva, n, n, b2);
            rmatrixinvupdateuv(b2, n, u, v);
            waserrors = waserrors||ap::fp_greater(matrixdiff(b1, b2, n, n),threshold);
        }
    }
    
    //
    // report
    //
    if( !silent )
    {
        printf("TESTING INVERSE UPDATE (REAL)\n");
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
static void makeacopy(const ap::real_2d_array& a,
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
LU decomposition
*************************************************************************/
static void matlu(ap::real_2d_array& a,
     int m,
     int n,
     ap::integer_1d_array& pivots)
{
    int i;
    int j;
    int jp;
    ap::real_1d_array t1;
    double s;

    pivots.setbounds(0, ap::minint(m-1, n-1));
    t1.setbounds(0, ap::maxint(m-1, n-1));
    ap::ap_error::make_assertion(m>=0&&n>=0, "Error in LUDecomposition: incorrect function arguments");
    
    //
    // Quick return if possible
    //
    if( m==0||n==0 )
    {
        return;
    }
    for(j = 0; j <= ap::minint(m-1, n-1); j++)
    {
        
        //
        // Find pivot and test for singularity.
        //
        jp = j;
        for(i = j+1; i <= m-1; i++)
        {
            if( ap::fp_greater(fabs(a(i,j)),fabs(a(jp,j))) )
            {
                jp = i;
            }
        }
        pivots(j) = jp;
        if( ap::fp_neq(a(jp,j),0) )
        {
            
            //
            //Apply the interchange to rows
            //
            if( jp!=j )
            {
                ap::vmove(&t1(0), 1, &a(j, 0), 1, ap::vlen(0,n-1));
                ap::vmove(&a(j, 0), 1, &a(jp, 0), 1, ap::vlen(0,n-1));
                ap::vmove(&a(jp, 0), 1, &t1(0), 1, ap::vlen(0,n-1));
            }
            
            //
            //Compute elements J+1:M of J-th column.
            //
            if( j<m )
            {
                jp = j+1;
                s = 1/a(j,j);
                ap::vmul(&a(jp, j), a.getstride(), ap::vlen(jp,m-1), s);
            }
        }
        if( j<ap::minint(m, n)-1 )
        {
            
            //
            //Update trailing submatrix.
            //
            jp = j+1;
            for(i = j+1; i <= m-1; i++)
            {
                s = a(i,j);
                ap::vsub(&a(i, jp), 1, &a(j, jp), 1, ap::vlen(jp,n-1), s);
            }
        }
    }
}


/*************************************************************************
Generate matrix with given condition number C (2-norm)
*************************************************************************/
static void generaterandomorthogonalmatrix(ap::real_2d_array& a0, int n)
{
    double t;
    double lambda;
    int s;
    int i;
    int j;
    double u1;
    double u2;
    ap::real_1d_array w;
    ap::real_1d_array v;
    ap::real_2d_array a;
    double sm;

    if( n<=0 )
    {
        return;
    }
    w.setbounds(1, n);
    v.setbounds(1, n);
    a.setbounds(1, n, 1, n);
    a0.setbounds(0, n-1, 0, n-1);
    
    //
    // Prepare A
    //
    for(i = 1; i <= n; i++)
    {
        for(j = 1; j <= n; j++)
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
    
    //
    // Calculate A using Stewart algorithm
    //
    for(s = 2; s <= n; s++)
    {
        
        //
        // Prepare v and Lambda = v'*v
        //
        do
        {
            i = 1;
            while(i<=s)
            {
                u1 = 2*ap::randomreal()-1;
                u2 = 2*ap::randomreal()-1;
                sm = u1*u1+u2*u2;
                if( ap::fp_eq(sm,0)||ap::fp_greater(sm,1) )
                {
                    continue;
                }
                sm = sqrt(-2*log(sm)/sm);
                v(i) = u1*sm;
                if( i+1<=s )
                {
                    v(i+1) = u2*sm;
                }
                i = i+2;
            }
            lambda = ap::vdotproduct(&v(1), 1, &v(1), 1, ap::vlen(1,s));
        }
        while(ap::fp_eq(lambda,0));
        lambda = 2/lambda;
        
        //
        // A * (I - 2 vv'/v'v ) =
        //   = A - (2/v'v) * A * v * v' =
        //   = A - (2/v'v) * w * v'
        //  where w = Av
        //
        for(i = 1; i <= s; i++)
        {
            t = ap::vdotproduct(&a(i, 1), 1, &v(1), 1, ap::vlen(1,s));
            w(i) = t;
        }
        for(i = 1; i <= s; i++)
        {
            t = w(i)*lambda;
            ap::vsub(&a(i, 1), 1, &v(1), 1, ap::vlen(1,s), t);
        }
    }
    
    //
    //
    //
    for(i = 1; i <= n; i++)
    {
        for(j = 1; j <= n; j++)
        {
            a0(i-1,j-1) = a(i,j);
        }
    }
}


static void generaterandommatrixcond(ap::real_2d_array& a0, int n, double c)
{
    double l1;
    double l2;
    ap::real_2d_array q1;
    ap::real_2d_array q2;
    ap::real_1d_array cc;
    int i;
    int j;
    int k;

    generaterandomorthogonalmatrix(q1, n);
    generaterandomorthogonalmatrix(q2, n);
    cc.setbounds(0, n-1);
    l1 = 0;
    l2 = log(1/c);
    cc(0) = exp(l1);
    for(i = 1; i <= n-2; i++)
    {
        cc(i) = exp(ap::randomreal()*(l2-l1)+l1);
    }
    cc(n-1) = exp(l2);
    a0.setbounds(0, n-1, 0, n-1);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            a0(i,j) = 0;
            for(k = 0; k <= n-1; k++)
            {
                a0(i,j) = a0(i,j)+q1(i,k)*cc(k)*q2(j,k);
            }
        }
    }
}


/*************************************************************************
triangular inverse
*************************************************************************/
static bool invmattr(ap::real_2d_array& a,
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
static bool invmatlu(ap::real_2d_array& a,
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
    if( !invmattr(a, n, true, false) )
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
static bool invmat(ap::real_2d_array& a, int n)
{
    bool result;
    ap::integer_1d_array pivots;

    matlu(a, n, n, pivots);
    result = invmatlu(a, pivots, n);
    return result;
}


/*************************************************************************
Diff
*************************************************************************/
static double matrixdiff(const ap::real_2d_array& a,
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
Update and inverse
*************************************************************************/
static bool updandinv(ap::real_2d_array& a,
     const ap::real_1d_array& u,
     const ap::real_1d_array& v,
     int n)
{
    bool result;
    ap::integer_1d_array pivots;
    int i;
    double r;

    for(i = 0; i <= n-1; i++)
    {
        r = u(i);
        ap::vadd(&a(i, 0), 1, &v(0), 1, ap::vlen(0,n-1), r);
    }
    matlu(a, n, n, pivots);
    result = invmatlu(a, pivots, n);
    return result;
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testinverseupdateunit_test_silent()
{
    bool result;

    result = testinverseupdate(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testinverseupdateunit_test()
{
    bool result;

    result = testinverseupdate(false);
    return result;
}




