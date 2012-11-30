
 
#include <stdio.h>
#include "testkmeansunit.h"

static void simpletest1(int nvars,
     int nc,
     int passcount,
     bool& converrors,
     bool& othererrors,
     bool& simpleerrors);
static double rnormal();
static double rsphere(ap::real_2d_array& xy, int n, int i);

bool testkmeans(bool silent)
{
    bool result;
    int nf;
    int maxnf;
    int nc;
    int maxnc;
    int passcount;
    int pass;
    bool waserrors;
    bool converrors;
    bool simpleerrors;
    bool complexerrors;
    bool othererrors;

    
    //
    // Primary settings
    //
    maxnf = 5;
    maxnc = 5;
    passcount = 10;
    waserrors = false;
    converrors = false;
    othererrors = false;
    simpleerrors = false;
    complexerrors = false;
    
    //
    //
    //
    for(nf = 1; nf <= maxnf; nf++)
    {
        for(nc = 1; nc <= maxnc; nc++)
        {
            simpletest1(nf, nc, passcount, converrors, othererrors, simpleerrors);
        }
    }
    
    //
    // Final report
    //
    waserrors = converrors||othererrors||simpleerrors||complexerrors;
    if( !silent )
    {
        printf("K-MEANS TEST\n");
        printf("TOTAL RESULTS:                           ");
        if( !waserrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("* CONVERGENCE:                           ");
        if( !converrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("* SIMPLE TASKS:                          ");
        if( !simpleerrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("* COMPLEX TASKS:                         ");
        if( !complexerrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("* OTHER PROPERTIES:                      ");
        if( !othererrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        if( waserrors )
        {
            printf("TEST SUMMARY: FAILED\n");
        }
        else
        {
            printf("TEST SUMMARY: PASSED\n");
        }
        printf("\n\n");
    }
    result = !waserrors;
    return result;
}


/*************************************************************************
Simple test 1: ellipsoid in NF-dimensional space.
compare k-means centers with random centers
*************************************************************************/
static void simpletest1(int nvars,
     int nc,
     int passcount,
     bool& converrors,
     bool& othererrors,
     bool& simpleerrors)
{
    int npoints;
    int majoraxis;
    ap::real_2d_array xy;
    ap::real_1d_array tmp;
    double v;
    int i;
    int j;
    int info;
    ap::real_2d_array c;
    ap::integer_1d_array xyc;
    int pass;
    int restarts;
    double ekmeans;
    double erandom;
    double dclosest;
    int cclosest;

    npoints = nc*100;
    restarts = 5;
    passcount = 10;
    tmp.setbounds(0, nvars-1);
    for(pass = 1; pass <= passcount; pass++)
    {
        
        //
        // Fill
        //
        xy.setbounds(0, npoints-1, 0, nvars-1);
        majoraxis = ap::randominteger(nvars);
        for(i = 0; i <= npoints-1; i++)
        {
            rsphere(xy, nvars, i);
            xy(i,majoraxis) = nc*xy(i,majoraxis);
        }
        
        //
        // Test
        //
        kmeansgenerate(xy, npoints, nvars, nc, restarts, info, c, xyc);
        if( info<0 )
        {
            converrors = true;
            return;
        }
        
        //
        // Test that XYC is correct mapping to cluster centers
        //
        for(i = 0; i <= npoints-1; i++)
        {
            cclosest = -1;
            dclosest = ap::maxrealnumber;
            for(j = 0; j <= nc-1; j++)
            {
                ap::vmove(&tmp(0), 1, &xy(i, 0), 1, ap::vlen(0,nvars-1));
                ap::vsub(&tmp(0), 1, &c(0, j), c.getstride(), ap::vlen(0,nvars-1));
                v = ap::vdotproduct(&tmp(0), 1, &tmp(0), 1, ap::vlen(0,nvars-1));
                if( ap::fp_less(v,dclosest) )
                {
                    cclosest = j;
                    dclosest = v;
                }
            }
            if( cclosest!=xyc(i) )
            {
                othererrors = true;
                return;
            }
        }
        
        //
        // Use first NC rows of XY as random centers
        // (XY is totally random, so it is as good as any other choice).
        //
        // Compare potential functions.
        //
        ekmeans = 0;
        for(i = 0; i <= npoints-1; i++)
        {
            ap::vmove(&tmp(0), 1, &xy(i, 0), 1, ap::vlen(0,nvars-1));
            ap::vsub(&tmp(0), 1, &c(0, xyc(i)), c.getstride(), ap::vlen(0,nvars-1));
            v = ap::vdotproduct(&tmp(0), 1, &tmp(0), 1, ap::vlen(0,nvars-1));
            ekmeans = ekmeans+v;
        }
        erandom = 0;
        for(i = 0; i <= npoints-1; i++)
        {
            dclosest = ap::maxrealnumber;
            for(j = 0; j <= nc-1; j++)
            {
                ap::vmove(&tmp(0), 1, &xy(i, 0), 1, ap::vlen(0,nvars-1));
                ap::vsub(&tmp(0), 1, &xy(j, 0), 1, ap::vlen(0,nvars-1));
                v = ap::vdotproduct(&tmp(0), 1, &tmp(0), 1, ap::vlen(0,nvars-1));
                if( ap::fp_less(v,dclosest) )
                {
                    dclosest = v;
                }
            }
            erandom = erandom+v;
        }
        if( ap::fp_less(erandom,ekmeans) )
        {
            simpleerrors = true;
            return;
        }
    }
}


/*************************************************************************
Random normal number
*************************************************************************/
static double rnormal()
{
    double result;
    double u;
    double v;
    double s;
    double x1;
    double x2;

    while(true)
    {
        u = 2*ap::randomreal()-1;
        v = 2*ap::randomreal()-1;
        s = ap::sqr(u)+ap::sqr(v);
        if( ap::fp_greater(s,0)&&ap::fp_less(s,1) )
        {
            s = sqrt(-2*log(s)/s);
            x1 = u*s;
            x2 = v*s;
            break;
        }
    }
    result = x1;
    return result;
}


/*************************************************************************
Random point from sphere
*************************************************************************/
static double rsphere(ap::real_2d_array& xy, int n, int i)
{
    double result;
    int j;
    double v;

    for(j = 0; j <= n-1; j++)
    {
        xy(i,j) = rnormal();
    }
    v = ap::vdotproduct(&xy(i, 0), 1, &xy(i, 0), 1, ap::vlen(0,n-1));
    v = ap::randomreal()/sqrt(v);
    ap::vmul(&xy(i, 0), 1, ap::vlen(0,n-1), v);
    return result;
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testkmeansunit_test_silent()
{
    bool result;

    result = testkmeans(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testkmeansunit_test()
{
    bool result;

    result = testkmeans(false);
    return result;
}




