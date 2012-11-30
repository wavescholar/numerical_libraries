
 
#include <stdio.h>
#include "testblasunit.h"

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

bool testblas(bool silent)
{
    bool result;
    int pass;
    int passcount;
    int n;
    int i;
    int i1;
    int i2;
    int j;
    int j1;
    int j2;
    int l;
    int k;
    int r;
    int i3;
    int j3;
    int col1;
    int col2;
    int row1;
    int row2;
    ap::real_1d_array x1;
    ap::real_1d_array x2;
    ap::real_2d_array a;
    ap::real_2d_array b;
    ap::real_2d_array c1;
    ap::real_2d_array c2;
    double err;
    double e1;
    double e2;
    double e3;
    double v;
    double scl1;
    double scl2;
    double scl3;
    bool was1;
    bool was2;
    bool trans1;
    bool trans2;
    double threshold;
    bool n2errors;
    bool hsnerrors;
    bool amaxerrors;
    bool mverrors;
    bool iterrors;
    bool cterrors;
    bool mmerrors;
    bool waserrors;

    n2errors = false;
    amaxerrors = false;
    hsnerrors = false;
    mverrors = false;
    iterrors = false;
    cterrors = false;
    mmerrors = false;
    waserrors = false;
    threshold = 10000*ap::machineepsilon;
    
    //
    // Test Norm2
    //
    passcount = 1000;
    e1 = 0;
    e2 = 0;
    e3 = 0;
    scl2 = 0.5*ap::maxrealnumber;
    scl3 = 2*ap::minrealnumber;
    for(pass = 1; pass <= passcount; pass++)
    {
        n = 1+ap::randominteger(1000);
        i1 = ap::randominteger(10);
        i2 = n+i1-1;
        x1.setbounds(i1, i2);
        x2.setbounds(i1, i2);
        for(i = i1; i <= i2; i++)
        {
            x1(i) = 2*ap::randomreal()-1;
        }
        v = 0;
        for(i = i1; i <= i2; i++)
        {
            v = v+ap::sqr(x1(i));
        }
        v = sqrt(v);
        e1 = ap::maxreal(e1, fabs(v-vectornorm2(x1, i1, i2)));
        for(i = i1; i <= i2; i++)
        {
            x2(i) = scl2*x1(i);
        }
        e2 = ap::maxreal(e2, fabs(v*scl2-vectornorm2(x2, i1, i2)));
        for(i = i1; i <= i2; i++)
        {
            x2(i) = scl3*x1(i);
        }
        e3 = ap::maxreal(e3, fabs(v*scl3-vectornorm2(x2, i1, i2)));
    }
    e2 = e2/scl2;
    e3 = e3/scl3;
    n2errors = ap::fp_greater_eq(e1,threshold)||ap::fp_greater_eq(e2,threshold)||ap::fp_greater_eq(e3,threshold);
    
    //
    // Testing VectorAbsMax, Column/Row AbsMax
    //
    x1.setbounds(1, 5);
    x1(1) = 2.0;
    x1(2) = 0.2;
    x1(3) = -1.3;
    x1(4) = 0.7;
    x1(5) = -3.0;
    amaxerrors = vectoridxabsmax(x1, 1, 5)!=5||vectoridxabsmax(x1, 1, 4)!=1||vectoridxabsmax(x1, 2, 4)!=3;
    n = 30;
    x1.setbounds(1, n);
    a.setbounds(1, n, 1, n);
    for(i = 1; i <= n; i++)
    {
        for(j = 1; j <= n; j++)
        {
            a(i,j) = 2*ap::randomreal()-1;
        }
    }
    was1 = false;
    was2 = false;
    for(pass = 1; pass <= 1000; pass++)
    {
        j = 1+ap::randominteger(n);
        i1 = 1+ap::randominteger(n);
        i2 = i1+ap::randominteger(n+1-i1);
        ap::vmove(&x1(i1), 1, &a(i1, j), a.getstride(), ap::vlen(i1,i2));
        if( vectoridxabsmax(x1, i1, i2)!=columnidxabsmax(a, i1, i2, j) )
        {
            was1 = true;
        }
        i = 1+ap::randominteger(n);
        j1 = 1+ap::randominteger(n);
        j2 = j1+ap::randominteger(n+1-j1);
        ap::vmove(&x1(j1), 1, &a(i, j1), 1, ap::vlen(j1,j2));
        if( vectoridxabsmax(x1, j1, j2)!=rowidxabsmax(a, j1, j2, i) )
        {
            was2 = true;
        }
    }
    amaxerrors = amaxerrors||was1||was2;
    
    //
    // Testing upper Hessenberg 1-norm
    //
    a.setbounds(1, 3, 1, 3);
    x1.setbounds(1, 3);
    a(1,1) = 2;
    a(1,2) = 3;
    a(1,3) = 1;
    a(2,1) = 4;
    a(2,2) = -5;
    a(2,3) = 8;
    a(3,1) = 99;
    a(3,2) = 3;
    a(3,3) = 1;
    hsnerrors = ap::fp_greater(fabs(upperhessenberg1norm(a, 1, 3, 1, 3, x1)-11),threshold);
    
    //
    // Testing MatrixVectorMultiply
    //
    a.setbounds(2, 3, 3, 5);
    x1.setbounds(1, 3);
    x2.setbounds(1, 2);
    a(2,3) = 2;
    a(2,4) = -1;
    a(2,5) = -1;
    a(3,3) = 1;
    a(3,4) = -2;
    a(3,5) = 2;
    x1(1) = 1;
    x1(2) = 2;
    x1(3) = 1;
    x2(1) = -1;
    x2(2) = -1;
    matrixvectormultiply(a, 2, 3, 3, 5, false, x1, 1, 3, 1.0, x2, 1, 2, 1.0);
    matrixvectormultiply(a, 2, 3, 3, 5, true, x2, 1, 2, 1.0, x1, 1, 3, 1.0);
    e1 = fabs(x1(1)+5)+fabs(x1(2)-8)+fabs(x1(3)+1)+fabs(x2(1)+2)+fabs(x2(2)+2);
    x1(1) = 1;
    x1(2) = 2;
    x1(3) = 1;
    x2(1) = -1;
    x2(2) = -1;
    matrixvectormultiply(a, 2, 3, 3, 5, false, x1, 1, 3, 1.0, x2, 1, 2, 0.0);
    matrixvectormultiply(a, 2, 3, 3, 5, true, x2, 1, 2, 1.0, x1, 1, 3, 0.0);
    e2 = fabs(x1(1)+3)+fabs(x1(2)-3)+fabs(x1(3)+1)+fabs(x2(1)+1)+fabs(x2(2)+1);
    mverrors = ap::fp_greater_eq(e1+e2,threshold);
    
    //
    // testing inplace transpose
    //
    n = 10;
    a.setbounds(1, n, 1, n);
    b.setbounds(1, n, 1, n);
    x1.setbounds(1, n-1);
    for(i = 1; i <= n; i++)
    {
        for(j = 1; j <= n; j++)
        {
            a(i,j) = ap::randomreal();
        }
    }
    passcount = 10000;
    was1 = false;
    for(pass = 1; pass <= passcount; pass++)
    {
        i1 = 1+ap::randominteger(n);
        i2 = i1+ap::randominteger(n-i1+1);
        j1 = 1+ap::randominteger(n-(i2-i1));
        j2 = j1+(i2-i1);
        copymatrix(a, i1, i2, j1, j2, b, i1, i2, j1, j2);
        inplacetranspose(b, i1, i2, j1, j2, x1);
        for(i = i1; i <= i2; i++)
        {
            for(j = j1; j <= j2; j++)
            {
                if( ap::fp_neq(a(i,j),b(i1+(j-j1),j1+(i-i1))) )
                {
                    was1 = true;
                }
            }
        }
    }
    iterrors = was1;
    
    //
    // testing copy and transpose
    //
    n = 10;
    a.setbounds(1, n, 1, n);
    b.setbounds(1, n, 1, n);
    for(i = 1; i <= n; i++)
    {
        for(j = 1; j <= n; j++)
        {
            a(i,j) = ap::randomreal();
        }
    }
    passcount = 10000;
    was1 = false;
    for(pass = 1; pass <= passcount; pass++)
    {
        i1 = 1+ap::randominteger(n);
        i2 = i1+ap::randominteger(n-i1+1);
        j1 = 1+ap::randominteger(n);
        j2 = j1+ap::randominteger(n-j1+1);
        copyandtranspose(a, i1, i2, j1, j2, b, j1, j2, i1, i2);
        for(i = i1; i <= i2; i++)
        {
            for(j = j1; j <= j2; j++)
            {
                if( ap::fp_neq(a(i,j),b(j,i)) )
                {
                    was1 = true;
                }
            }
        }
    }
    cterrors = was1;
    
    //
    // Testing MatrixMatrixMultiply
    //
    n = 10;
    a.setbounds(1, 2*n, 1, 2*n);
    b.setbounds(1, 2*n, 1, 2*n);
    c1.setbounds(1, 2*n, 1, 2*n);
    c2.setbounds(1, 2*n, 1, 2*n);
    x1.setbounds(1, n);
    x2.setbounds(1, n);
    for(i = 1; i <= 2*n; i++)
    {
        for(j = 1; j <= 2*n; j++)
        {
            a(i,j) = ap::randomreal();
            b(i,j) = ap::randomreal();
        }
    }
    passcount = 1000;
    was1 = false;
    for(pass = 1; pass <= passcount; pass++)
    {
        for(i = 1; i <= 2*n; i++)
        {
            for(j = 1; j <= 2*n; j++)
            {
                c1(i,j) = 2.1*i+3.1*j;
                c2(i,j) = c1(i,j);
            }
        }
        l = 1+ap::randominteger(n);
        k = 1+ap::randominteger(n);
        r = 1+ap::randominteger(n);
        i1 = 1+ap::randominteger(n);
        j1 = 1+ap::randominteger(n);
        i2 = 1+ap::randominteger(n);
        j2 = 1+ap::randominteger(n);
        i3 = 1+ap::randominteger(n);
        j3 = 1+ap::randominteger(n);
        trans1 = ap::fp_greater(ap::randomreal(),0.5);
        trans2 = ap::fp_greater(ap::randomreal(),0.5);
        if( trans1 )
        {
            col1 = l;
            row1 = k;
        }
        else
        {
            col1 = k;
            row1 = l;
        }
        if( trans2 )
        {
            col2 = k;
            row2 = r;
        }
        else
        {
            col2 = r;
            row2 = k;
        }
        scl1 = ap::randomreal();
        scl2 = ap::randomreal();
        matrixmatrixmultiply(a, i1, i1+row1-1, j1, j1+col1-1, trans1, b, i2, i2+row2-1, j2, j2+col2-1, trans2, scl1, c1, i3, i3+l-1, j3, j3+r-1, scl2, x1);
        naivematrixmatrixmultiply(a, i1, i1+row1-1, j1, j1+col1-1, trans1, b, i2, i2+row2-1, j2, j2+col2-1, trans2, scl1, c2, i3, i3+l-1, j3, j3+r-1, scl2);
        err = 0;
        for(i = 1; i <= l; i++)
        {
            for(j = 1; j <= r; j++)
            {
                err = ap::maxreal(err, fabs(c1(i3+i-1,j3+j-1)-c2(i3+i-1,j3+j-1)));
            }
        }
        if( ap::fp_greater(err,threshold) )
        {
            was1 = true;
            break;
        }
    }
    mmerrors = was1;
    
    //
    // report
    //
    waserrors = n2errors||amaxerrors||hsnerrors||mverrors||iterrors||cterrors||mmerrors;
    if( !silent )
    {
        printf("TESTING BLAS\n");
        printf("VectorNorm2:                             ");
        if( n2errors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("AbsMax (vector/row/column):              ");
        if( amaxerrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("UpperHessenberg1Norm:                    ");
        if( hsnerrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("MatrixVectorMultiply:                    ");
        if( mverrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("InplaceTranspose:                        ");
        if( iterrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("CopyAndTranspose:                        ");
        if( cterrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("MatrixMatrixMultiply:                    ");
        if( mmerrors )
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
Silent unit test
*************************************************************************/
bool testblasunit_test_silent()
{
    bool result;

    result = testblas(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testblasunit_test()
{
    bool result;

    result = testblas(false);
    return result;
}




