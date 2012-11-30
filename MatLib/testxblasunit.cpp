
 
#include <stdio.h>
#include "testxblasunit.h"

static const double xchunk = 1048576;
static const int xchunkcount = 4;

bool testxblas(bool silent)
{
    bool result;
    bool approxerrors;
    bool exactnesserrors;
    bool waserrors;
    double approxthreshold;
    int maxn;
    int passcount;
    int n;
    int i;
    int pass;
    double rv1;
    double rv2;
    double rv2err;
    ap::complex cv1;
    ap::complex cv2;
    double cv2err;
    double cv2errx;
    double cv2erry;
    ap::real_1d_array rx;
    ap::real_1d_array ry;
    ap::complex_1d_array cx;
    ap::complex_1d_array cy;
    ap::real_1d_array temp;
    double b;
    double s;

    approxerrors = false;
    exactnesserrors = false;
    waserrors = false;
    approxthreshold = 1000*ap::machineepsilon;
    maxn = 1000;
    passcount = 10;
    
    //
    // tests:
    // 1. ability to calculate dot product
    // 2. higher precision
    //
    for(n = 1; n <= maxn; n++)
    {
        for(pass = 1; pass <= passcount; pass++)
        {
            
            //
            //  ability to approximately calculate real dot product
            //
            rx.setlength(n);
            ry.setlength(n);
            temp.setlength(n);
            for(i = 0; i <= n-1; i++)
            {
                if( ap::fp_greater(ap::randomreal(),0.2) )
                {
                    rx(i) = 2*ap::randomreal()-1;
                }
                else
                {
                    rx(i) = 0;
                }
                if( ap::fp_greater(ap::randomreal(),0.2) )
                {
                    ry(i) = 2*ap::randomreal()-1;
                }
                else
                {
                    ry(i) = 0;
                }
            }
            rv1 = ap::vdotproduct(&rx(0), 1, &ry(0), 1, ap::vlen(0,n-1));
            xdot(rx, ry, n, temp, rv2, rv2err);
            approxerrors = approxerrors||ap::fp_greater(fabs(rv1-rv2),approxthreshold);
            
            //
            //  ability to approximately calculate complex dot product
            //
            cx.setlength(n);
            cy.setlength(n);
            temp.setlength(2*n);
            for(i = 0; i <= n-1; i++)
            {
                if( ap::fp_greater(ap::randomreal(),0.2) )
                {
                    cx(i).x = 2*ap::randomreal()-1;
                    cx(i).y = 2*ap::randomreal()-1;
                }
                else
                {
                    cx(i) = 0;
                }
                if( ap::fp_greater(ap::randomreal(),0.2) )
                {
                    cy(i).x = 2*ap::randomreal()-1;
                    cy(i).y = 2*ap::randomreal()-1;
                }
                else
                {
                    cy(i) = 0;
                }
            }
            cv1 = ap::vdotproduct(&cx(0), 1, "N", &cy(0), 1, "N", ap::vlen(0,n-1));
            xcdot(cx, cy, n, temp, cv2, cv2err);
            approxerrors = approxerrors||ap::fp_greater(ap::abscomplex(cv1-cv2),approxthreshold);
        }
    }
    
    //
    // test of precision: real
    //
    n = 50000;
    rx.setlength(n);
    ry.setlength(n);
    temp.setlength(n);
    for(pass = 0; pass <= passcount-1; pass++)
    {
        ap::ap_error::make_assertion(n%2==0, "");
        
        //
        // First test: X + X + ... + X - X - X - ... - X = 1*X
        //
        s = exp(double(ap::maxint(pass, 50)));
        if( pass==passcount-1&&pass>1 )
        {
            s = ap::maxrealnumber;
        }
        ry(0) = (2*ap::randomreal()-1)*s*sqrt(2*ap::randomreal());
        for(i = 1; i <= n-1; i++)
        {
            ry(i) = ry(0);
        }
        for(i = 0; i <= n/2-1; i++)
        {
            rx(i) = 1;
        }
        for(i = n/2; i <= n-2; i++)
        {
            rx(i) = -1;
        }
        rx(n-1) = 0;
        xdot(rx, ry, n, temp, rv2, rv2err);
        exactnesserrors = exactnesserrors||ap::fp_less(rv2err,0);
        exactnesserrors = exactnesserrors||ap::fp_greater(rv2err,4*ap::machineepsilon*fabs(ry(0)));
        exactnesserrors = exactnesserrors||ap::fp_greater(fabs(rv2-ry(0)),rv2err);
        
        //
        // First test: X + X + ... + X = N*X
        //
        s = exp(double(ap::maxint(pass, 50)));
        if( pass==passcount-1&&pass>1 )
        {
            s = ap::maxrealnumber;
        }
        ry(0) = (2*ap::randomreal()-1)*s*sqrt(2*ap::randomreal());
        for(i = 1; i <= n-1; i++)
        {
            ry(i) = ry(0);
        }
        for(i = 0; i <= n-1; i++)
        {
            rx(i) = 1;
        }
        xdot(rx, ry, n, temp, rv2, rv2err);
        exactnesserrors = exactnesserrors||ap::fp_less(rv2err,0);
        exactnesserrors = exactnesserrors||ap::fp_greater(rv2err,4*ap::machineepsilon*fabs(ry(0))*n);
        exactnesserrors = exactnesserrors||ap::fp_greater(fabs(rv2-n*ry(0)),rv2err);
    }
    
    //
    // test of precision: complex
    //
    n = 50000;
    cx.setlength(n);
    cy.setlength(n);
    temp.setlength(2*n);
    for(pass = 0; pass <= passcount-1; pass++)
    {
        ap::ap_error::make_assertion(n%2==0, "");
        
        //
        // First test: X + X + ... + X - X - X - ... - X = 1*X
        //
        s = exp(double(ap::maxint(pass, 50)));
        if( pass==passcount-1&&pass>1 )
        {
            s = ap::maxrealnumber;
        }
        cy(0).x = (2*ap::randomreal()-1)*s*sqrt(2*ap::randomreal());
        cy(0).y = (2*ap::randomreal()-1)*s*sqrt(2*ap::randomreal());
        for(i = 1; i <= n-1; i++)
        {
            cy(i) = cy(0);
        }
        for(i = 0; i <= n/2-1; i++)
        {
            cx(i) = 1;
        }
        for(i = n/2; i <= n-2; i++)
        {
            cx(i) = -1;
        }
        cx(n-1) = 0;
        xcdot(cx, cy, n, temp, cv2, cv2err);
        exactnesserrors = exactnesserrors||ap::fp_less(cv2err,0);
        exactnesserrors = exactnesserrors||ap::fp_greater(cv2err,4*ap::machineepsilon*ap::abscomplex(cy(0)));
        exactnesserrors = exactnesserrors||ap::fp_greater(ap::abscomplex(cv2-cy(0)),cv2err);
        
        //
        // First test: X + X + ... + X = N*X
        //
        s = exp(double(ap::maxint(pass, 50)));
        if( pass==passcount-1&&pass>1 )
        {
            s = ap::maxrealnumber;
        }
        cy(0) = (2*ap::randomreal()-1)*s*sqrt(2*ap::randomreal());
        for(i = 1; i <= n-1; i++)
        {
            cy(i) = cy(0);
        }
        for(i = 0; i <= n-1; i++)
        {
            cx(i) = 1;
        }
        xcdot(cx, cy, n, temp, cv2, cv2err);
        exactnesserrors = exactnesserrors||ap::fp_less(cv2err,0);
        exactnesserrors = exactnesserrors||ap::fp_greater(cv2err,4*ap::machineepsilon*ap::abscomplex(cy(0))*n);
        exactnesserrors = exactnesserrors||ap::fp_greater(ap::abscomplex(cv2-n*cy(0)),cv2err);
    }
    
    //
    // report
    //
    waserrors = approxerrors||exactnesserrors;
    if( !silent )
    {
        printf("TESTING XBLAS\n");
        printf("APPROX.TESTS:                            ");
        if( approxerrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("EXACT TESTS:                             ");
        if( exactnesserrors )
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
    
    //
    // end
    //
    result = !waserrors;
    return result;
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testxblasunit_test_silent()
{
    bool result;

    result = testxblas(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testxblasunit_test()
{
    bool result;

    result = testxblas(false);
    return result;
}




