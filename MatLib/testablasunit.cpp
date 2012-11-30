
 
#include <stdio.h>
#include "testablasunit.h"

static void naivematrixmatrixmultiply(const ap::real_2d_array& a,
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
     double alpha,
     ap::real_2d_array& c,
     int ci1,
     int ci2,
     int cj1,
     int cj2,
     double beta);
static bool testtrsm(int minn, int maxn);
static bool testsyrk(int minn, int maxn);
static bool testgemm(int minn, int maxn);
static bool testtrans(int minn, int maxn);
static bool testrank1(int minn, int maxn);
static bool testmv(int minn, int maxn);
static bool testcopy(int minn, int maxn);

bool testablas(bool silent)
{
    bool result;
    double threshold;
    bool trsmerrors;
    bool syrkerrors;
    bool gemmerrors;
    bool transerrors;
    bool rank1errors;
    bool mverrors;
    bool copyerrors;
    bool waserrors;
    ap::real_2d_array ra;

    trsmerrors = false;
    syrkerrors = false;
    gemmerrors = false;
    transerrors = false;
    rank1errors = false;
    mverrors = false;
    copyerrors = false;
    waserrors = false;
    threshold = 10000*ap::machineepsilon;
    trsmerrors = trsmerrors||testtrsm(1, 3*ablasblocksize(ra)+1);
    syrkerrors = syrkerrors||testsyrk(1, 3*ablasblocksize(ra)+1);
    gemmerrors = gemmerrors||testgemm(1, 3*ablasblocksize(ra)+1);
    transerrors = transerrors||testtrans(1, 3*ablasblocksize(ra)+1);
    rank1errors = rank1errors||testrank1(1, 3*ablasblocksize(ra)+1);
    mverrors = mverrors||testmv(1, 3*ablasblocksize(ra)+1);
    copyerrors = copyerrors||testcopy(1, 3*ablasblocksize(ra)+1);
    
    //
    // report
    //
    waserrors = trsmerrors||syrkerrors||gemmerrors||transerrors||rank1errors||mverrors||copyerrors;
    if( !silent )
    {
        printf("TESTING ABLAS\n");
        printf("* TRSM:                                  ");
        if( trsmerrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* SYRK:                                  ");
        if( syrkerrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* GEMM:                                  ");
        if( gemmerrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* TRANS:                                 ");
        if( transerrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* RANK1:                                 ");
        if( rank1errors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* MV:                                    ");
        if( mverrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* COPY:                                  ");
        if( copyerrors )
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


/*************************************************************************
Reference implementation

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************/
void refcmatrixrighttrsm(int m,
     int n,
     const ap::complex_2d_array& a,
     int i1,
     int j1,
     bool isupper,
     bool isunit,
     int optype,
     ap::complex_2d_array& x,
     int i2,
     int j2)
{
    ap::complex_2d_array a1;
    ap::complex_2d_array a2;
    ap::complex_1d_array tx;
    int i;
    int j;
    ap::complex vc;
    bool rupper;

    if( n*m==0 )
    {
        return;
    }
    a1.setlength(n, n);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            a1(i,j) = 0;
        }
    }
    if( isupper )
    {
        for(i = 0; i <= n-1; i++)
        {
            for(j = i; j <= n-1; j++)
            {
                a1(i,j) = a(i1+i,j1+j);
            }
        }
    }
    else
    {
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= i; j++)
            {
                a1(i,j) = a(i1+i,j1+j);
            }
        }
    }
    rupper = isupper;
    if( isunit )
    {
        for(i = 0; i <= n-1; i++)
        {
            a1(i,i) = 1;
        }
    }
    a2.setlength(n, n);
    if( optype==0 )
    {
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                a2(i,j) = a1(i,j);
            }
        }
    }
    if( optype==1 )
    {
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                a2(i,j) = a1(j,i);
            }
        }
        rupper = !rupper;
    }
    if( optype==2 )
    {
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                a2(i,j) = ap::conj(a1(j,i));
            }
        }
        rupper = !rupper;
    }
    internalcmatrixtrinverse(a2, n, rupper, false);
    tx.setlength(n);
    for(i = 0; i <= m-1; i++)
    {
        ap::vmove(&tx(0), 1, &x(i2+i, j2), 1, "N", ap::vlen(0,n-1));
        for(j = 0; j <= n-1; j++)
        {
            vc = ap::vdotproduct(&tx(0), 1, "N", &a2(0, j), a2.getstride(), "N", ap::vlen(0,n-1));
            x(i2+i,j2+j) = vc;
        }
    }
}


/*************************************************************************
Reference implementation

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************/
void refcmatrixlefttrsm(int m,
     int n,
     const ap::complex_2d_array& a,
     int i1,
     int j1,
     bool isupper,
     bool isunit,
     int optype,
     ap::complex_2d_array& x,
     int i2,
     int j2)
{
    ap::complex_2d_array a1;
    ap::complex_2d_array a2;
    ap::complex_1d_array tx;
    int i;
    int j;
    ap::complex vc;
    bool rupper;

    if( n*m==0 )
    {
        return;
    }
    a1.setlength(m, m);
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= m-1; j++)
        {
            a1(i,j) = 0;
        }
    }
    if( isupper )
    {
        for(i = 0; i <= m-1; i++)
        {
            for(j = i; j <= m-1; j++)
            {
                a1(i,j) = a(i1+i,j1+j);
            }
        }
    }
    else
    {
        for(i = 0; i <= m-1; i++)
        {
            for(j = 0; j <= i; j++)
            {
                a1(i,j) = a(i1+i,j1+j);
            }
        }
    }
    rupper = isupper;
    if( isunit )
    {
        for(i = 0; i <= m-1; i++)
        {
            a1(i,i) = 1;
        }
    }
    a2.setlength(m, m);
    if( optype==0 )
    {
        for(i = 0; i <= m-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                a2(i,j) = a1(i,j);
            }
        }
    }
    if( optype==1 )
    {
        for(i = 0; i <= m-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                a2(i,j) = a1(j,i);
            }
        }
        rupper = !rupper;
    }
    if( optype==2 )
    {
        for(i = 0; i <= m-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                a2(i,j) = ap::conj(a1(j,i));
            }
        }
        rupper = !rupper;
    }
    internalcmatrixtrinverse(a2, m, rupper, false);
    tx.setlength(m);
    for(j = 0; j <= n-1; j++)
    {
        ap::vmove(&tx(0), 1, &x(i2, j2+j), x.getstride(), "N", ap::vlen(0,m-1));
        for(i = 0; i <= m-1; i++)
        {
            vc = ap::vdotproduct(&a2(i, 0), 1, "N", &tx(0), 1, "N", ap::vlen(0,m-1));
            x(i2+i,j2+j) = vc;
        }
    }
}


/*************************************************************************
Reference implementation

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************/
void refrmatrixrighttrsm(int m,
     int n,
     const ap::real_2d_array& a,
     int i1,
     int j1,
     bool isupper,
     bool isunit,
     int optype,
     ap::real_2d_array& x,
     int i2,
     int j2)
{
    ap::real_2d_array a1;
    ap::real_2d_array a2;
    ap::real_1d_array tx;
    int i;
    int j;
    double vr;
    bool rupper;

    if( n*m==0 )
    {
        return;
    }
    a1.setlength(n, n);
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            a1(i,j) = 0;
        }
    }
    if( isupper )
    {
        for(i = 0; i <= n-1; i++)
        {
            for(j = i; j <= n-1; j++)
            {
                a1(i,j) = a(i1+i,j1+j);
            }
        }
    }
    else
    {
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= i; j++)
            {
                a1(i,j) = a(i1+i,j1+j);
            }
        }
    }
    rupper = isupper;
    if( isunit )
    {
        for(i = 0; i <= n-1; i++)
        {
            a1(i,i) = 1;
        }
    }
    a2.setlength(n, n);
    if( optype==0 )
    {
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                a2(i,j) = a1(i,j);
            }
        }
    }
    if( optype==1 )
    {
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                a2(i,j) = a1(j,i);
            }
        }
        rupper = !rupper;
    }
    internalrmatrixtrinverse(a2, n, rupper, false);
    tx.setlength(n);
    for(i = 0; i <= m-1; i++)
    {
        ap::vmove(&tx(0), 1, &x(i2+i, j2), 1, ap::vlen(0,n-1));
        for(j = 0; j <= n-1; j++)
        {
            vr = ap::vdotproduct(&tx(0), 1, &a2(0, j), a2.getstride(), ap::vlen(0,n-1));
            x(i2+i,j2+j) = vr;
        }
    }
}


/*************************************************************************
Reference implementation

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************/
void refrmatrixlefttrsm(int m,
     int n,
     const ap::real_2d_array& a,
     int i1,
     int j1,
     bool isupper,
     bool isunit,
     int optype,
     ap::real_2d_array& x,
     int i2,
     int j2)
{
    ap::real_2d_array a1;
    ap::real_2d_array a2;
    ap::real_1d_array tx;
    int i;
    int j;
    double vr;
    bool rupper;

    if( n*m==0 )
    {
        return;
    }
    a1.setlength(m, m);
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= m-1; j++)
        {
            a1(i,j) = 0;
        }
    }
    if( isupper )
    {
        for(i = 0; i <= m-1; i++)
        {
            for(j = i; j <= m-1; j++)
            {
                a1(i,j) = a(i1+i,j1+j);
            }
        }
    }
    else
    {
        for(i = 0; i <= m-1; i++)
        {
            for(j = 0; j <= i; j++)
            {
                a1(i,j) = a(i1+i,j1+j);
            }
        }
    }
    rupper = isupper;
    if( isunit )
    {
        for(i = 0; i <= m-1; i++)
        {
            a1(i,i) = 1;
        }
    }
    a2.setlength(m, m);
    if( optype==0 )
    {
        for(i = 0; i <= m-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                a2(i,j) = a1(i,j);
            }
        }
    }
    if( optype==1 )
    {
        for(i = 0; i <= m-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                a2(i,j) = a1(j,i);
            }
        }
        rupper = !rupper;
    }
    internalrmatrixtrinverse(a2, m, rupper, false);
    tx.setlength(m);
    for(j = 0; j <= n-1; j++)
    {
        ap::vmove(&tx(0), 1, &x(i2, j2+j), x.getstride(), ap::vlen(0,m-1));
        for(i = 0; i <= m-1; i++)
        {
            vr = ap::vdotproduct(&a2(i, 0), 1, &tx(0), 1, ap::vlen(0,m-1));
            x(i2+i,j2+j) = vr;
        }
    }
}


/*************************************************************************
Internal subroutine.
Triangular matrix inversion

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************/
bool internalcmatrixtrinverse(ap::complex_2d_array& a,
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
                    if( i+1<j )
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
            if( j+1<n )
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
Internal subroutine.
Triangular matrix inversion

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************/
bool internalrmatrixtrinverse(ap::real_2d_array& a,
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
Reference SYRK subroutine.

  -- ALGLIB routine --
     16.12.2009
     Bochkanov Sergey
*************************************************************************/
void refcmatrixsyrk(int n,
     int k,
     double alpha,
     const ap::complex_2d_array& a,
     int ia,
     int ja,
     int optypea,
     double beta,
     ap::complex_2d_array& c,
     int ic,
     int jc,
     bool isupper)
{
    ap::complex_2d_array ae;
    int i;
    int j;
    ap::complex vc;

    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            if( isupper&&j>=i||!isupper&&j<=i )
            {
                if( ap::fp_eq(beta,0) )
                {
                    c(i+ic,j+jc) = 0;
                }
                else
                {
                    c(i+ic,j+jc) = c(i+ic,j+jc)*beta;
                }
            }
        }
    }
    if( ap::fp_eq(alpha,0) )
    {
        return;
    }
    if( n*k>0 )
    {
        ae.setlength(n, k);
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= k-1; j++)
        {
            if( optypea==0 )
            {
                ae(i,j) = a(ia+i,ja+j);
            }
            if( optypea==2 )
            {
                ae(i,j) = ap::conj(a(ia+j,ja+i));
            }
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            vc = 0;
            if( k>0 )
            {
                vc = ap::vdotproduct(&ae(i, 0), 1, "N", &ae(j, 0), 1, "Conj", ap::vlen(0,k-1));
            }
            vc = alpha*vc;
            if( isupper&&j>=i )
            {
                c(ic+i,jc+j) = vc+c(ic+i,jc+j);
            }
            if( !isupper&&j<=i )
            {
                c(ic+i,jc+j) = vc+c(ic+i,jc+j);
            }
        }
    }
}


/*************************************************************************
Reference SYRK subroutine.

  -- ALGLIB routine --
     16.12.2009
     Bochkanov Sergey
*************************************************************************/
void refrmatrixsyrk(int n,
     int k,
     double alpha,
     const ap::real_2d_array& a,
     int ia,
     int ja,
     int optypea,
     double beta,
     ap::real_2d_array& c,
     int ic,
     int jc,
     bool isupper)
{
    ap::real_2d_array ae;
    int i;
    int j;
    double vr;

    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            if( isupper&&j>=i||!isupper&&j<=i )
            {
                if( ap::fp_eq(beta,0) )
                {
                    c(i+ic,j+jc) = 0;
                }
                else
                {
                    c(i+ic,j+jc) = c(i+ic,j+jc)*beta;
                }
            }
        }
    }
    if( ap::fp_eq(alpha,0) )
    {
        return;
    }
    if( n*k>0 )
    {
        ae.setlength(n, k);
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= k-1; j++)
        {
            if( optypea==0 )
            {
                ae(i,j) = a(ia+i,ja+j);
            }
            if( optypea==1 )
            {
                ae(i,j) = a(ia+j,ja+i);
            }
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            vr = 0;
            if( k>0 )
            {
                vr = ap::vdotproduct(&ae(i, 0), 1, &ae(j, 0), 1, ap::vlen(0,k-1));
            }
            vr = alpha*vr;
            if( isupper&&j>=i )
            {
                c(ic+i,jc+j) = vr+c(ic+i,jc+j);
            }
            if( !isupper&&j<=i )
            {
                c(ic+i,jc+j) = vr+c(ic+i,jc+j);
            }
        }
    }
}


/*************************************************************************
Reference GEMM,
ALGLIB subroutine
*************************************************************************/
void refcmatrixgemm(int m,
     int n,
     int k,
     ap::complex alpha,
     const ap::complex_2d_array& a,
     int ia,
     int ja,
     int optypea,
     const ap::complex_2d_array& b,
     int ib,
     int jb,
     int optypeb,
     ap::complex beta,
     ap::complex_2d_array& c,
     int ic,
     int jc)
{
    ap::complex_2d_array ae;
    ap::complex_2d_array be;
    int i;
    int j;
    ap::complex vc;

    ae.setlength(m, k);
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= k-1; j++)
        {
            if( optypea==0 )
            {
                ae(i,j) = a(ia+i,ja+j);
            }
            if( optypea==1 )
            {
                ae(i,j) = a(ia+j,ja+i);
            }
            if( optypea==2 )
            {
                ae(i,j) = ap::conj(a(ia+j,ja+i));
            }
        }
    }
    be.setlength(k, n);
    for(i = 0; i <= k-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            if( optypeb==0 )
            {
                be(i,j) = b(ib+i,jb+j);
            }
            if( optypeb==1 )
            {
                be(i,j) = b(ib+j,jb+i);
            }
            if( optypeb==2 )
            {
                be(i,j) = ap::conj(b(ib+j,jb+i));
            }
        }
    }
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            vc = ap::vdotproduct(&ae(i, 0), 1, "N", &be(0, j), be.getstride(), "N", ap::vlen(0,k-1));
            vc = alpha*vc;
            if( beta!=0 )
            {
                vc = vc+beta*c(ic+i,jc+j);
            }
            c(ic+i,jc+j) = vc;
        }
    }
}


/*************************************************************************
Reference GEMM,
ALGLIB subroutine
*************************************************************************/
void refrmatrixgemm(int m,
     int n,
     int k,
     double alpha,
     const ap::real_2d_array& a,
     int ia,
     int ja,
     int optypea,
     const ap::real_2d_array& b,
     int ib,
     int jb,
     int optypeb,
     double beta,
     ap::real_2d_array& c,
     int ic,
     int jc)
{
    ap::real_2d_array ae;
    ap::real_2d_array be;
    int i;
    int j;
    double vc;

    ae.setlength(m, k);
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= k-1; j++)
        {
            if( optypea==0 )
            {
                ae(i,j) = a(ia+i,ja+j);
            }
            if( optypea==1 )
            {
                ae(i,j) = a(ia+j,ja+i);
            }
        }
    }
    be.setlength(k, n);
    for(i = 0; i <= k-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            if( optypeb==0 )
            {
                be(i,j) = b(ib+i,jb+j);
            }
            if( optypeb==1 )
            {
                be(i,j) = b(ib+j,jb+i);
            }
        }
    }
    for(i = 0; i <= m-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            vc = ap::vdotproduct(&ae(i, 0), 1, &be(0, j), be.getstride(), ap::vlen(0,k-1));
            vc = alpha*vc;
            if( ap::fp_neq(beta,0) )
            {
                vc = vc+beta*c(ic+i,jc+j);
            }
            c(ic+i,jc+j) = vc;
        }
    }
}


static void naivematrixmatrixmultiply(const ap::real_2d_array& a,
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
     double alpha,
     ap::real_2d_array& c,
     int ci1,
     int ci2,
     int cj1,
     int cj2,
     double beta)
{
    int arows;
    int acols;
    int brows;
    int bcols;
    int i;
    int j;
    int k;
    int l;
    int r;
    double v;
    ap::real_1d_array x1;
    ap::real_1d_array x2;

    
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
    ap::ap_error::make_assertion(acols==brows, "NaiveMatrixMatrixMultiply: incorrect matrix sizes!");
    if( arows<=0||acols<=0||brows<=0||bcols<=0 )
    {
        return;
    }
    l = arows;
    r = bcols;
    k = acols;
    x1.setbounds(1, k);
    x2.setbounds(1, k);
    for(i = 1; i <= l; i++)
    {
        for(j = 1; j <= r; j++)
        {
            if( !transa )
            {
                if( !transb )
                {
                    v = ap::vdotproduct(&b(bi1, bj1+j-1), b.getstride(), &a(ai1+i-1, aj1), 1, ap::vlen(bi1,bi2));
                }
                else
                {
                    v = ap::vdotproduct(&b(bi1+j-1, bj1), 1, &a(ai1+i-1, aj1), 1, ap::vlen(bj1,bj2));
                }
            }
            else
            {
                if( !transb )
                {
                    v = ap::vdotproduct(&b(bi1, bj1+j-1), b.getstride(), &a(ai1, aj1+i-1), a.getstride(), ap::vlen(bi1,bi2));
                }
                else
                {
                    v = ap::vdotproduct(&b(bi1+j-1, bj1), 1, &a(ai1, aj1+i-1), a.getstride(), ap::vlen(bj1,bj2));
                }
            }
            if( ap::fp_eq(beta,0) )
            {
                c(ci1+i-1,cj1+j-1) = alpha*v;
            }
            else
            {
                c(ci1+i-1,cj1+j-1) = beta*c(ci1+i-1,cj1+j-1)+alpha*v;
            }
        }
    }
}


/*************************************************************************
?Matrix????TRSM tests

Returns False for passed test, True - for failed
*************************************************************************/
static bool testtrsm(int minn, int maxn)
{
    bool result;
    int n;
    int m;
    int mx;
    int i;
    int j;
    int optype;
    int uppertype;
    int unittype;
    int xoffsi;
    int xoffsj;
    int aoffsitype;
    int aoffsjtype;
    int aoffsi;
    int aoffsj;
    ap::real_2d_array refra;
    ap::real_2d_array refrxl;
    ap::real_2d_array refrxr;
    ap::complex_2d_array refca;
    ap::complex_2d_array refcxl;
    ap::complex_2d_array refcxr;
    ap::real_2d_array ra;
    ap::complex_2d_array ca;
    ap::real_2d_array rxr1;
    ap::real_2d_array rxl1;
    ap::complex_2d_array cxr1;
    ap::complex_2d_array cxl1;
    ap::real_2d_array rxr2;
    ap::real_2d_array rxl2;
    ap::complex_2d_array cxr2;
    ap::complex_2d_array cxl2;
    double threshold;

    threshold = ap::sqr(double(maxn))*100*ap::machineepsilon;
    result = false;
    for(mx = minn; mx <= maxn; mx++)
    {
        
        //
        // Select random M/N in [1,MX] such that max(M,N)=MX
        //
        m = 1+ap::randominteger(mx);
        n = 1+ap::randominteger(mx);
        if( ap::fp_greater(ap::randomreal(),0.5) )
        {
            m = mx;
        }
        else
        {
            n = mx;
        }
        
        //
        // Initialize RefRA/RefCA by random matrices whose upper
        // and lower triangle submatrices are non-degenerate
        // well-conditioned matrices.
        //
        // Matrix size is 2Mx2M (four copies of same MxM matrix
        // to test different offsets)
        //
        refra.setlength(2*m, 2*m);
        for(i = 0; i <= m-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                refra(i,j) = 0.2*ap::randomreal()-0.1;
            }
        }
        for(i = 0; i <= m-1; i++)
        {
            refra(i,i) = (2*ap::randominteger(1)-1)*(2*m+ap::randomreal());
        }
        for(i = 0; i <= m-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                refra(i+m,j) = refra(i,j);
                refra(i,j+m) = refra(i,j);
                refra(i+m,j+m) = refra(i,j);
            }
        }
        refca.setlength(2*m, 2*m);
        for(i = 0; i <= m-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                refca(i,j).x = 0.2*ap::randomreal()-0.1;
                refca(i,j).y = 0.2*ap::randomreal()-0.1;
            }
        }
        for(i = 0; i <= m-1; i++)
        {
            refca(i,i).x = (2*ap::randominteger(2)-1)*(2*m+ap::randomreal());
            refca(i,i).y = (2*ap::randominteger(2)-1)*(2*m+ap::randomreal());
        }
        for(i = 0; i <= m-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                refca(i+m,j) = refca(i,j);
                refca(i,j+m) = refca(i,j);
                refca(i+m,j+m) = refca(i,j);
            }
        }
        
        //
        // Generate random XL/XR.
        //
        // XR is NxM matrix (matrix for 'Right' subroutines)
        // XL is MxN matrix (matrix for 'Left' subroutines)
        //
        refrxr.setlength(n, m);
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                refrxr(i,j) = 2*ap::randomreal()-1;
            }
        }
        refrxl.setlength(m, n);
        for(i = 0; i <= m-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                refrxl(i,j) = 2*ap::randomreal()-1;
            }
        }
        refcxr.setlength(n, m);
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                refcxr(i,j).x = 2*ap::randomreal()-1;
                refcxr(i,j).y = 2*ap::randomreal()-1;
            }
        }
        refcxl.setlength(m, n);
        for(i = 0; i <= m-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                refcxl(i,j).x = 2*ap::randomreal()-1;
                refcxl(i,j).y = 2*ap::randomreal()-1;
            }
        }
        
        //
        // test different types of operations, offsets, and so on...
        //
        // to avoid unnecessary slowdown we don't test ALL possible
        // combinations of operation types. We just generate one random
        // set of parameters and test it.
        //
        ra.setlength(2*m, 2*m);
        rxr1.setlength(n, m);
        rxr2.setlength(n, m);
        rxl1.setlength(m, n);
        rxl2.setlength(m, n);
        ca.setlength(2*m, 2*m);
        cxr1.setlength(n, m);
        cxr2.setlength(n, m);
        cxl1.setlength(m, n);
        cxl2.setlength(m, n);
        optype = ap::randominteger(3);
        uppertype = ap::randominteger(2);
        unittype = ap::randominteger(2);
        xoffsi = ap::randominteger(2);
        xoffsj = ap::randominteger(2);
        aoffsitype = ap::randominteger(2);
        aoffsjtype = ap::randominteger(2);
        aoffsi = m*aoffsitype;
        aoffsj = m*aoffsjtype;
        
        //
        // copy A, XR, XL (fill unused parts with random garbage)
        //
        for(i = 0; i <= 2*m-1; i++)
        {
            for(j = 0; j <= 2*m-1; j++)
            {
                if( i>=aoffsi&&i<aoffsi+m&&j>=aoffsj&&j<aoffsj+m )
                {
                    ca(i,j) = refca(i,j);
                    ra(i,j) = refra(i,j);
                }
                else
                {
                    ca(i,j) = ap::randomreal();
                    ra(i,j) = ap::randomreal();
                }
            }
        }
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                if( i>=xoffsi&&j>=xoffsj )
                {
                    cxr1(i,j) = refcxr(i,j);
                    cxr2(i,j) = refcxr(i,j);
                    rxr1(i,j) = refrxr(i,j);
                    rxr2(i,j) = refrxr(i,j);
                }
                else
                {
                    cxr1(i,j) = ap::randomreal();
                    cxr2(i,j) = cxr1(i,j);
                    rxr1(i,j) = ap::randomreal();
                    rxr2(i,j) = rxr1(i,j);
                }
            }
        }
        for(i = 0; i <= m-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                if( i>=xoffsi&&j>=xoffsj )
                {
                    cxl1(i,j) = refcxl(i,j);
                    cxl2(i,j) = refcxl(i,j);
                    rxl1(i,j) = refrxl(i,j);
                    rxl2(i,j) = refrxl(i,j);
                }
                else
                {
                    cxl1(i,j) = ap::randomreal();
                    cxl2(i,j) = cxl1(i,j);
                    rxl1(i,j) = ap::randomreal();
                    rxl2(i,j) = rxl1(i,j);
                }
            }
        }
        
        //
        // Test CXR
        //
        cmatrixrighttrsm(n-xoffsi, m-xoffsj, ca, aoffsi, aoffsj, uppertype==0, unittype==0, optype, cxr1, xoffsi, xoffsj);
        refcmatrixrighttrsm(n-xoffsi, m-xoffsj, ca, aoffsi, aoffsj, uppertype==0, unittype==0, optype, cxr2, xoffsi, xoffsj);
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                result = result||ap::fp_greater(ap::abscomplex(cxr1(i,j)-cxr2(i,j)),threshold);
            }
        }
        
        //
        // Test CXL
        //
        cmatrixlefttrsm(m-xoffsi, n-xoffsj, ca, aoffsi, aoffsj, uppertype==0, unittype==0, optype, cxl1, xoffsi, xoffsj);
        refcmatrixlefttrsm(m-xoffsi, n-xoffsj, ca, aoffsi, aoffsj, uppertype==0, unittype==0, optype, cxl2, xoffsi, xoffsj);
        for(i = 0; i <= m-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                result = result||ap::fp_greater(ap::abscomplex(cxl1(i,j)-cxl2(i,j)),threshold);
            }
        }
        if( optype<2 )
        {
            
            //
            // Test RXR
            //
            rmatrixrighttrsm(n-xoffsi, m-xoffsj, ra, aoffsi, aoffsj, uppertype==0, unittype==0, optype, rxr1, xoffsi, xoffsj);
            refrmatrixrighttrsm(n-xoffsi, m-xoffsj, ra, aoffsi, aoffsj, uppertype==0, unittype==0, optype, rxr2, xoffsi, xoffsj);
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= m-1; j++)
                {
                    result = result||ap::fp_greater(fabs(rxr1(i,j)-rxr2(i,j)),threshold);
                }
            }
            
            //
            // Test RXL
            //
            rmatrixlefttrsm(m-xoffsi, n-xoffsj, ra, aoffsi, aoffsj, uppertype==0, unittype==0, optype, rxl1, xoffsi, xoffsj);
            refrmatrixlefttrsm(m-xoffsi, n-xoffsj, ra, aoffsi, aoffsj, uppertype==0, unittype==0, optype, rxl2, xoffsi, xoffsj);
            for(i = 0; i <= m-1; i++)
            {
                for(j = 0; j <= n-1; j++)
                {
                    result = result||ap::fp_greater(fabs(rxl1(i,j)-rxl2(i,j)),threshold);
                }
            }
        }
    }
    return result;
}


/*************************************************************************
SYRK tests

Returns False for passed test, True - for failed
*************************************************************************/
static bool testsyrk(int minn, int maxn)
{
    bool result;
    int n;
    int k;
    int mx;
    int i;
    int j;
    int uppertype;
    int xoffsi;
    int xoffsj;
    int aoffsitype;
    int aoffsjtype;
    int aoffsi;
    int aoffsj;
    int alphatype;
    int betatype;
    ap::real_2d_array refra;
    ap::real_2d_array refrc;
    ap::complex_2d_array refca;
    ap::complex_2d_array refcc;
    double alpha;
    double beta;
    ap::real_2d_array ra1;
    ap::real_2d_array ra2;
    ap::complex_2d_array ca1;
    ap::complex_2d_array ca2;
    ap::real_2d_array rc;
    ap::real_2d_array rct;
    ap::complex_2d_array cc;
    ap::complex_2d_array cct;
    double threshold;

    threshold = maxn*100*ap::machineepsilon;
    result = false;
    for(mx = minn; mx <= maxn; mx++)
    {
        
        //
        // Select random M/N in [1,MX] such that max(M,N)=MX
        //
        k = 1+ap::randominteger(mx);
        n = 1+ap::randominteger(mx);
        if( ap::fp_greater(ap::randomreal(),0.5) )
        {
            k = mx;
        }
        else
        {
            n = mx;
        }
        
        //
        // Initialize RefRA/RefCA by random Hermitian matrices,
        // RefRC/RefCC by random matrices
        //
        // RA/CA size is 2Nx2N (four copies of same NxN matrix
        // to test different offsets)
        //
        refra.setlength(2*n, 2*n);
        refca.setlength(2*n, 2*n);
        for(i = 0; i <= n-1; i++)
        {
            refra(i,i) = 2*ap::randomreal()-1;
            refca(i,i) = 2*ap::randomreal()-1;
            for(j = i+1; j <= n-1; j++)
            {
                refra(i,j) = 2*ap::randomreal()-1;
                refca(i,j).x = 2*ap::randomreal()-1;
                refca(i,j).y = 2*ap::randomreal()-1;
                refra(j,i) = refra(i,j);
                refca(j,i) = ap::conj(refca(i,j));
            }
        }
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                refra(i+n,j) = refra(i,j);
                refra(i,j+n) = refra(i,j);
                refra(i+n,j+n) = refra(i,j);
                refca(i+n,j) = refca(i,j);
                refca(i,j+n) = refca(i,j);
                refca(i+n,j+n) = refca(i,j);
            }
        }
        refrc.setlength(n, k);
        refcc.setlength(n, k);
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= k-1; j++)
            {
                refrc(i,j) = 2*ap::randomreal()-1;
                refcc(i,j).x = 2*ap::randomreal()-1;
                refcc(i,j).y = 2*ap::randomreal()-1;
            }
        }
        
        //
        // test different types of operations, offsets, and so on...
        //
        // to avoid unnecessary slowdown we don't test ALL possible
        // combinations of operation types. We just generate one random
        // set of parameters and test it.
        //
        ra1.setlength(2*n, 2*n);
        ra2.setlength(2*n, 2*n);
        ca1.setlength(2*n, 2*n);
        ca2.setlength(2*n, 2*n);
        rc.setlength(n, k);
        rct.setlength(k, n);
        cc.setlength(n, k);
        cct.setlength(k, n);
        uppertype = ap::randominteger(2);
        xoffsi = ap::randominteger(2);
        xoffsj = ap::randominteger(2);
        aoffsitype = ap::randominteger(2);
        aoffsjtype = ap::randominteger(2);
        alphatype = ap::randominteger(2);
        betatype = ap::randominteger(2);
        aoffsi = n*aoffsitype;
        aoffsj = n*aoffsjtype;
        alpha = alphatype*(2*ap::randomreal()-1);
        beta = betatype*(2*ap::randomreal()-1);
        
        //
        // copy A, C (fill unused parts with random garbage)
        //
        for(i = 0; i <= 2*n-1; i++)
        {
            for(j = 0; j <= 2*n-1; j++)
            {
                if( i>=aoffsi&&i<aoffsi+n&&j>=aoffsj&&j<aoffsj+n )
                {
                    ca1(i,j) = refca(i,j);
                    ca2(i,j) = refca(i,j);
                    ra1(i,j) = refra(i,j);
                    ra2(i,j) = refra(i,j);
                }
                else
                {
                    ca1(i,j) = ap::randomreal();
                    ca2(i,j) = ca1(i,j);
                    ra1(i,j) = ap::randomreal();
                    ra2(i,j) = ra1(i,j);
                }
            }
        }
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= k-1; j++)
            {
                if( i>=xoffsi&&j>=xoffsj )
                {
                    rc(i,j) = refrc(i,j);
                    rct(j,i) = refrc(i,j);
                    cc(i,j) = refcc(i,j);
                    cct(j,i) = refcc(i,j);
                }
                else
                {
                    rc(i,j) = ap::randomreal();
                    rct(j,i) = rc(i,j);
                    cc(i,j) = ap::randomreal();
                    cct(j,i) = cct(j,i);
                }
            }
        }
        
        //
        // Test complex
        // Only one of transform types is selected and tested
        //
        if( ap::fp_greater(ap::randomreal(),0.5) )
        {
            cmatrixsyrk(n-xoffsi, k-xoffsj, alpha, cc, xoffsi, xoffsj, 0, beta, ca1, aoffsi, aoffsj, uppertype==0);
            refcmatrixsyrk(n-xoffsi, k-xoffsj, alpha, cc, xoffsi, xoffsj, 0, beta, ca2, aoffsi, aoffsj, uppertype==0);
        }
        else
        {
            cmatrixsyrk(n-xoffsi, k-xoffsj, alpha, cct, xoffsj, xoffsi, 2, beta, ca1, aoffsi, aoffsj, uppertype==0);
            refcmatrixsyrk(n-xoffsi, k-xoffsj, alpha, cct, xoffsj, xoffsi, 2, beta, ca2, aoffsi, aoffsj, uppertype==0);
        }
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                result = result||ap::fp_greater(ap::abscomplex(ca1(i,j)-ca2(i,j)),threshold);
            }
        }
        
        //
        // Test real
        // Only one of transform types is selected and tested
        //
        if( ap::fp_greater(ap::randomreal(),0.5) )
        {
            rmatrixsyrk(n-xoffsi, k-xoffsj, alpha, rc, xoffsi, xoffsj, 0, beta, ra1, aoffsi, aoffsj, uppertype==0);
            refrmatrixsyrk(n-xoffsi, k-xoffsj, alpha, rc, xoffsi, xoffsj, 0, beta, ra2, aoffsi, aoffsj, uppertype==0);
        }
        else
        {
            rmatrixsyrk(n-xoffsi, k-xoffsj, alpha, rct, xoffsj, xoffsi, 1, beta, ra1, aoffsi, aoffsj, uppertype==0);
            refrmatrixsyrk(n-xoffsi, k-xoffsj, alpha, rct, xoffsj, xoffsi, 1, beta, ra2, aoffsi, aoffsj, uppertype==0);
        }
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                result = result||ap::fp_greater(fabs(ra1(i,j)-ra2(i,j)),threshold);
            }
        }
    }
    return result;
}


/*************************************************************************
GEMM tests

Returns False for passed test, True - for failed
*************************************************************************/
static bool testgemm(int minn, int maxn)
{
    bool result;
    int m;
    int n;
    int k;
    int mx;
    int i;
    int j;
    int aoffsi;
    int aoffsj;
    int aoptype;
    int aoptyper;
    int boffsi;
    int boffsj;
    int boptype;
    int boptyper;
    int coffsi;
    int coffsj;
    ap::real_2d_array refra;
    ap::real_2d_array refrb;
    ap::real_2d_array refrc;
    ap::complex_2d_array refca;
    ap::complex_2d_array refcb;
    ap::complex_2d_array refcc;
    double alphar;
    double betar;
    ap::complex alphac;
    ap::complex betac;
    ap::real_2d_array rc1;
    ap::real_2d_array rc2;
    ap::complex_2d_array cc1;
    ap::complex_2d_array cc2;
    double threshold;

    threshold = maxn*100*ap::machineepsilon;
    result = false;
    for(mx = minn; mx <= maxn; mx++)
    {
        
        //
        // Select random M/N/K in [1,MX] such that max(M,N,K)=MX
        //
        m = 1+ap::randominteger(mx);
        n = 1+ap::randominteger(mx);
        k = 1+ap::randominteger(mx);
        i = ap::randominteger(3);
        if( i==0 )
        {
            m = mx;
        }
        if( i==1 )
        {
            n = mx;
        }
        if( i==2 )
        {
            k = mx;
        }
        
        //
        // Initialize A/B/C by random matrices with size (MaxN+1)*(MaxN+1)
        //
        refra.setlength(maxn+1, maxn+1);
        refrb.setlength(maxn+1, maxn+1);
        refrc.setlength(maxn+1, maxn+1);
        refca.setlength(maxn+1, maxn+1);
        refcb.setlength(maxn+1, maxn+1);
        refcc.setlength(maxn+1, maxn+1);
        for(i = 0; i <= maxn; i++)
        {
            for(j = 0; j <= maxn; j++)
            {
                refra(i,j) = 2*ap::randomreal()-1;
                refrb(i,j) = 2*ap::randomreal()-1;
                refrc(i,j) = 2*ap::randomreal()-1;
                refca(i,j).x = 2*ap::randomreal()-1;
                refca(i,j).y = 2*ap::randomreal()-1;
                refcb(i,j).x = 2*ap::randomreal()-1;
                refcb(i,j).y = 2*ap::randomreal()-1;
                refcc(i,j).x = 2*ap::randomreal()-1;
                refcc(i,j).y = 2*ap::randomreal()-1;
            }
        }
        
        //
        // test different types of operations, offsets, and so on...
        //
        // to avoid unnecessary slowdown we don't test ALL possible
        // combinations of operation types. We just generate one random
        // set of parameters and test it.
        //
        rc1.setlength(maxn+1, maxn+1);
        rc2.setlength(maxn+1, maxn+1);
        cc1.setlength(maxn+1, maxn+1);
        cc2.setlength(maxn+1, maxn+1);
        aoffsi = ap::randominteger(2);
        aoffsj = ap::randominteger(2);
        aoptype = ap::randominteger(3);
        aoptyper = ap::randominteger(2);
        boffsi = ap::randominteger(2);
        boffsj = ap::randominteger(2);
        boptype = ap::randominteger(3);
        boptyper = ap::randominteger(2);
        coffsi = ap::randominteger(2);
        coffsj = ap::randominteger(2);
        alphar = ap::randominteger(2)*(2*ap::randomreal()-1);
        betar = ap::randominteger(2)*(2*ap::randomreal()-1);
        if( ap::fp_greater(ap::randomreal(),0.5) )
        {
            alphac.x = 2*ap::randomreal()-1;
            alphac.y = 2*ap::randomreal()-1;
        }
        else
        {
            alphac = 0;
        }
        if( ap::fp_greater(ap::randomreal(),0.5) )
        {
            betac.x = 2*ap::randomreal()-1;
            betac.y = 2*ap::randomreal()-1;
        }
        else
        {
            betac = 0;
        }
        
        //
        // copy C
        //
        for(i = 0; i <= maxn; i++)
        {
            for(j = 0; j <= maxn; j++)
            {
                rc1(i,j) = refrc(i,j);
                rc2(i,j) = refrc(i,j);
                cc1(i,j) = refcc(i,j);
                cc2(i,j) = refcc(i,j);
            }
        }
        
        //
        // Test complex
        //
        cmatrixgemm(m, n, k, alphac, refca, aoffsi, aoffsj, aoptype, refcb, boffsi, boffsj, boptype, betac, cc1, coffsi, coffsj);
        refcmatrixgemm(m, n, k, alphac, refca, aoffsi, aoffsj, aoptype, refcb, boffsi, boffsj, boptype, betac, cc2, coffsi, coffsj);
        for(i = 0; i <= maxn; i++)
        {
            for(j = 0; j <= maxn; j++)
            {
                result = result||ap::fp_greater(ap::abscomplex(cc1(i,j)-cc2(i,j)),threshold);
            }
        }
        
        //
        // Test real
        //
        rmatrixgemm(m, n, k, alphar, refra, aoffsi, aoffsj, aoptyper, refrb, boffsi, boffsj, boptyper, betar, rc1, coffsi, coffsj);
        refrmatrixgemm(m, n, k, alphar, refra, aoffsi, aoffsj, aoptyper, refrb, boffsi, boffsj, boptyper, betar, rc2, coffsi, coffsj);
        for(i = 0; i <= maxn; i++)
        {
            for(j = 0; j <= maxn; j++)
            {
                result = result||ap::fp_greater(fabs(rc1(i,j)-rc2(i,j)),threshold);
            }
        }
    }
    return result;
}


/*************************************************************************
transpose tests

Returns False for passed test, True - for failed
*************************************************************************/
static bool testtrans(int minn, int maxn)
{
    bool result;
    int m;
    int n;
    int mx;
    int i;
    int j;
    int aoffsi;
    int aoffsj;
    int boffsi;
    int boffsj;
    double v1;
    double v2;
    double threshold;
    ap::real_2d_array refra;
    ap::real_2d_array refrb;
    ap::complex_2d_array refca;
    ap::complex_2d_array refcb;

    result = false;
    threshold = 1000*ap::machineepsilon;
    for(mx = minn; mx <= maxn; mx++)
    {
        
        //
        // Select random M/N in [1,MX] such that max(M,N)=MX
        // Generate random V1 and V2 which are used to fill
        // RefRB/RefCB with control values.
        //
        m = 1+ap::randominteger(mx);
        n = 1+ap::randominteger(mx);
        if( ap::randominteger(2)==0 )
        {
            m = mx;
        }
        else
        {
            n = mx;
        }
        v1 = ap::randomreal();
        v2 = ap::randomreal();
        
        //
        // Initialize A by random matrix with size (MaxN+1)*(MaxN+1)
        // Fill B with control values
        //
        refra.setlength(maxn+1, maxn+1);
        refrb.setlength(maxn+1, maxn+1);
        refca.setlength(maxn+1, maxn+1);
        refcb.setlength(maxn+1, maxn+1);
        for(i = 0; i <= maxn; i++)
        {
            for(j = 0; j <= maxn; j++)
            {
                refra(i,j) = 2*ap::randomreal()-1;
                refca(i,j).x = 2*ap::randomreal()-1;
                refca(i,j).y = 2*ap::randomreal()-1;
                refrb(i,j) = i*v1+j*v2;
                refcb(i,j) = i*v1+j*v2;
            }
        }
        
        //
        // test different offsets (zero or one)
        //
        // to avoid unnecessary slowdown we don't test ALL possible
        // combinations of operation types. We just generate one random
        // set of parameters and test it.
        //
        aoffsi = ap::randominteger(2);
        aoffsj = ap::randominteger(2);
        boffsi = ap::randominteger(2);
        boffsj = ap::randominteger(2);
        rmatrixtranspose(m, n, refra, aoffsi, aoffsj, refrb, boffsi, boffsj);
        for(i = 0; i <= maxn; i++)
        {
            for(j = 0; j <= maxn; j++)
            {
                if( i<boffsi||i>=boffsi+n||j<boffsj||j>=boffsj+m )
                {
                    result = result||ap::fp_greater(fabs(refrb(i,j)-(v1*i+v2*j)),threshold);
                }
                else
                {
                    result = result||ap::fp_greater(fabs(refrb(i,j)-refra(aoffsi+j-boffsj,aoffsj+i-boffsi)),threshold);
                }
            }
        }
        cmatrixtranspose(m, n, refca, aoffsi, aoffsj, refcb, boffsi, boffsj);
        for(i = 0; i <= maxn; i++)
        {
            for(j = 0; j <= maxn; j++)
            {
                if( i<boffsi||i>=boffsi+n||j<boffsj||j>=boffsj+m )
                {
                    result = result||ap::fp_greater(ap::abscomplex(refcb(i,j)-(v1*i+v2*j)),threshold);
                }
                else
                {
                    result = result||ap::fp_greater(ap::abscomplex(refcb(i,j)-refca(aoffsi+j-boffsj,aoffsj+i-boffsi)),threshold);
                }
            }
        }
    }
    return result;
}


/*************************************************************************
rank-1tests

Returns False for passed test, True - for failed
*************************************************************************/
static bool testrank1(int minn, int maxn)
{
    bool result;
    int m;
    int n;
    int mx;
    int i;
    int j;
    int aoffsi;
    int aoffsj;
    int uoffs;
    int voffs;
    double threshold;
    ap::real_2d_array refra;
    ap::real_2d_array refrb;
    ap::complex_2d_array refca;
    ap::complex_2d_array refcb;
    ap::real_1d_array ru;
    ap::real_1d_array rv;
    ap::complex_1d_array cu;
    ap::complex_1d_array cv;

    result = false;
    threshold = 1000*ap::machineepsilon;
    for(mx = minn; mx <= maxn; mx++)
    {
        
        //
        // Select random M/N in [1,MX] such that max(M,N)=MX
        //
        m = 1+ap::randominteger(mx);
        n = 1+ap::randominteger(mx);
        if( ap::randominteger(2)==0 )
        {
            m = mx;
        }
        else
        {
            n = mx;
        }
        
        //
        // Initialize A by random matrix with size (MaxN+1)*(MaxN+1)
        // Fill B with control values
        //
        refra.setlength(maxn+maxn, maxn+maxn);
        refrb.setlength(maxn+maxn, maxn+maxn);
        refca.setlength(maxn+maxn, maxn+maxn);
        refcb.setlength(maxn+maxn, maxn+maxn);
        for(i = 0; i <= 2*maxn-1; i++)
        {
            for(j = 0; j <= 2*maxn-1; j++)
            {
                refra(i,j) = 2*ap::randomreal()-1;
                refca(i,j).x = 2*ap::randomreal()-1;
                refca(i,j).y = 2*ap::randomreal()-1;
                refrb(i,j) = refra(i,j);
                refcb(i,j) = refca(i,j);
            }
        }
        ru.setlength(2*m);
        cu.setlength(2*m);
        for(i = 0; i <= 2*m-1; i++)
        {
            ru(i) = 2*ap::randomreal()-1;
            cu(i).x = 2*ap::randomreal()-1;
            cu(i).y = 2*ap::randomreal()-1;
        }
        rv.setlength(2*n);
        cv.setlength(2*n);
        for(i = 0; i <= 2*n-1; i++)
        {
            rv(i) = 2*ap::randomreal()-1;
            cv(i).x = 2*ap::randomreal()-1;
            cv(i).y = 2*ap::randomreal()-1;
        }
        
        //
        // test different offsets (zero or one)
        //
        // to avoid unnecessary slowdown we don't test ALL possible
        // combinations of operation types. We just generate one random
        // set of parameters and test it.
        //
        aoffsi = ap::randominteger(maxn);
        aoffsj = ap::randominteger(maxn);
        uoffs = ap::randominteger(m);
        voffs = ap::randominteger(n);
        cmatrixrank1(m, n, refca, aoffsi, aoffsj, cu, uoffs, cv, voffs);
        for(i = 0; i <= 2*maxn-1; i++)
        {
            for(j = 0; j <= 2*maxn-1; j++)
            {
                if( i<aoffsi||i>=aoffsi+m||j<aoffsj||j>=aoffsj+n )
                {
                    result = result||ap::fp_greater(ap::abscomplex(refca(i,j)-refcb(i,j)),threshold);
                }
                else
                {
                    result = result||ap::fp_greater(ap::abscomplex(refca(i,j)-(refcb(i,j)+cu(i-aoffsi+uoffs)*cv(j-aoffsj+voffs))),threshold);
                }
            }
        }
        rmatrixrank1(m, n, refra, aoffsi, aoffsj, ru, uoffs, rv, voffs);
        for(i = 0; i <= 2*maxn-1; i++)
        {
            for(j = 0; j <= 2*maxn-1; j++)
            {
                if( i<aoffsi||i>=aoffsi+m||j<aoffsj||j>=aoffsj+n )
                {
                    result = result||ap::fp_greater(fabs(refra(i,j)-refrb(i,j)),threshold);
                }
                else
                {
                    result = result||ap::fp_greater(fabs(refra(i,j)-(refrb(i,j)+ru(i-aoffsi+uoffs)*rv(j-aoffsj+voffs))),threshold);
                }
            }
        }
    }
    return result;
}


/*************************************************************************
MV tests

Returns False for passed test, True - for failed
*************************************************************************/
static bool testmv(int minn, int maxn)
{
    bool result;
    int m;
    int n;
    int mx;
    int i;
    int j;
    int aoffsi;
    int aoffsj;
    int xoffs;
    int yoffs;
    int opca;
    int opra;
    double threshold;
    double rv1;
    double rv2;
    ap::complex cv1;
    ap::complex cv2;
    ap::real_2d_array refra;
    ap::complex_2d_array refca;
    ap::real_1d_array rx;
    ap::real_1d_array ry;
    ap::complex_1d_array cx;
    ap::complex_1d_array cy;

    result = false;
    threshold = 1000*ap::machineepsilon;
    for(mx = minn; mx <= maxn; mx++)
    {
        
        //
        // Select random M/N in [1,MX] such that max(M,N)=MX
        //
        m = 1+ap::randominteger(mx);
        n = 1+ap::randominteger(mx);
        if( ap::randominteger(2)==0 )
        {
            m = mx;
        }
        else
        {
            n = mx;
        }
        
        //
        // Initialize A by random matrix with size (MaxN+MaxN)*(MaxN+MaxN)
        // Initialize X by random vector with size (MaxN+MaxN)
        // Fill Y by control values
        //
        refra.setlength(maxn+maxn, maxn+maxn);
        refca.setlength(maxn+maxn, maxn+maxn);
        for(i = 0; i <= 2*maxn-1; i++)
        {
            for(j = 0; j <= 2*maxn-1; j++)
            {
                refra(i,j) = 2*ap::randomreal()-1;
                refca(i,j).x = 2*ap::randomreal()-1;
                refca(i,j).y = 2*ap::randomreal()-1;
            }
        }
        rx.setlength(2*maxn);
        cx.setlength(2*maxn);
        ry.setlength(2*maxn);
        cy.setlength(2*maxn);
        for(i = 0; i <= 2*maxn-1; i++)
        {
            rx(i) = 2*ap::randomreal()-1;
            cx(i).x = 2*ap::randomreal()-1;
            cx(i).y = 2*ap::randomreal()-1;
            ry(i) = i;
            cy(i) = i;
        }
        
        //
        // test different offsets (zero or one)
        //
        // to avoid unnecessary slowdown we don't test ALL possible
        // combinations of operation types. We just generate one random
        // set of parameters and test it.
        //
        aoffsi = ap::randominteger(maxn);
        aoffsj = ap::randominteger(maxn);
        xoffs = ap::randominteger(maxn);
        yoffs = ap::randominteger(maxn);
        opca = ap::randominteger(3);
        opra = ap::randominteger(2);
        cmatrixmv(m, n, refca, aoffsi, aoffsj, opca, cx, xoffs, cy, yoffs);
        for(i = 0; i <= 2*maxn-1; i++)
        {
            if( i<yoffs||i>=yoffs+m )
            {
                result = result||cy(i)!=i;
            }
            else
            {
                cv1 = cy(i);
                if( opca==0 )
                {
                    cv2 = ap::vdotproduct(&refca(aoffsi+i-yoffs, aoffsj), 1, "N", &cx(xoffs), 1, "N", ap::vlen(aoffsj,aoffsj+n-1));
                }
                if( opca==1 )
                {
                    cv2 = ap::vdotproduct(&refca(aoffsi, aoffsj+i-yoffs), refca.getstride(), "N", &cx(xoffs), 1, "N", ap::vlen(aoffsi,aoffsi+n-1));
                }
                if( opca==2 )
                {
                    cv2 = ap::vdotproduct(&refca(aoffsi, aoffsj+i-yoffs), refca.getstride(), "Conj", &cx(xoffs), 1, "N", ap::vlen(aoffsi,aoffsi+n-1));
                }
                result = result||ap::fp_greater(ap::abscomplex(cv1-cv2),threshold);
            }
        }
        rmatrixmv(m, n, refra, aoffsi, aoffsj, opra, rx, xoffs, ry, yoffs);
        for(i = 0; i <= 2*maxn-1; i++)
        {
            if( i<yoffs||i>=yoffs+m )
            {
                result = result||ap::fp_neq(ry(i),i);
            }
            else
            {
                rv1 = ry(i);
                if( opra==0 )
                {
                    rv2 = ap::vdotproduct(&refra(aoffsi+i-yoffs, aoffsj), 1, &rx(xoffs), 1, ap::vlen(aoffsj,aoffsj+n-1));
                }
                if( opra==1 )
                {
                    rv2 = ap::vdotproduct(&refra(aoffsi, aoffsj+i-yoffs), refra.getstride(), &rx(xoffs), 1, ap::vlen(aoffsi,aoffsi+n-1));
                }
                result = result||ap::fp_greater(fabs(rv1-rv2),threshold);
            }
        }
    }
    return result;
}


/*************************************************************************
COPY tests

Returns False for passed test, True - for failed
*************************************************************************/
static bool testcopy(int minn, int maxn)
{
    bool result;
    int m;
    int n;
    int mx;
    int i;
    int j;
    int aoffsi;
    int aoffsj;
    int boffsi;
    int boffsj;
    double threshold;
    double rv1;
    double rv2;
    ap::complex cv1;
    ap::complex cv2;
    ap::real_2d_array ra;
    ap::real_2d_array rb;
    ap::complex_2d_array ca;
    ap::complex_2d_array cb;

    result = false;
    threshold = 1000*ap::machineepsilon;
    for(mx = minn; mx <= maxn; mx++)
    {
        
        //
        // Select random M/N in [1,MX] such that max(M,N)=MX
        //
        m = 1+ap::randominteger(mx);
        n = 1+ap::randominteger(mx);
        if( ap::randominteger(2)==0 )
        {
            m = mx;
        }
        else
        {
            n = mx;
        }
        
        //
        // Initialize A by random matrix with size (MaxN+MaxN)*(MaxN+MaxN)
        // Initialize X by random vector with size (MaxN+MaxN)
        // Fill Y by control values
        //
        ra.setlength(maxn+maxn, maxn+maxn);
        ca.setlength(maxn+maxn, maxn+maxn);
        rb.setlength(maxn+maxn, maxn+maxn);
        cb.setlength(maxn+maxn, maxn+maxn);
        for(i = 0; i <= 2*maxn-1; i++)
        {
            for(j = 0; j <= 2*maxn-1; j++)
            {
                ra(i,j) = 2*ap::randomreal()-1;
                ca(i,j).x = 2*ap::randomreal()-1;
                ca(i,j).y = 2*ap::randomreal()-1;
                rb(i,j) = 1+2*i+3*j;
                cb(i,j) = 1+2*i+3*j;
            }
        }
        
        //
        // test different offsets (zero or one)
        //
        // to avoid unnecessary slowdown we don't test ALL possible
        // combinations of operation types. We just generate one random
        // set of parameters and test it.
        //
        aoffsi = ap::randominteger(maxn);
        aoffsj = ap::randominteger(maxn);
        boffsi = ap::randominteger(maxn);
        boffsj = ap::randominteger(maxn);
        cmatrixcopy(m, n, ca, aoffsi, aoffsj, cb, boffsi, boffsj);
        for(i = 0; i <= 2*maxn-1; i++)
        {
            for(j = 0; j <= 2*maxn-1; j++)
            {
                if( i<boffsi||i>=boffsi+m||j<boffsj||j>=boffsj+n )
                {
                    result = result||cb(i,j)!=1+2*i+3*j;
                }
                else
                {
                    result = result||ap::fp_greater(ap::abscomplex(ca(aoffsi+i-boffsi,aoffsj+j-boffsj)-cb(i,j)),threshold);
                }
            }
        }
        rmatrixcopy(m, n, ra, aoffsi, aoffsj, rb, boffsi, boffsj);
        for(i = 0; i <= 2*maxn-1; i++)
        {
            for(j = 0; j <= 2*maxn-1; j++)
            {
                if( i<boffsi||i>=boffsi+m||j<boffsj||j>=boffsj+n )
                {
                    result = result||ap::fp_neq(rb(i,j),1+2*i+3*j);
                }
                else
                {
                    result = result||ap::fp_greater(fabs(ra(aoffsi+i-boffsi,aoffsj+j-boffsj)-rb(i,j)),threshold);
                }
            }
        }
    }
    return result;
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testablasunit_test_silent()
{
    bool result;

    result = testablas(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testablasunit_test()
{
    bool result;

    result = testablas(false);
    return result;
}




