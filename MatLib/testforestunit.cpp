
 
#include <stdio.h>
#include "testforestunit.h"

static void testprocessing(bool& err);
static void basictest1(int nvars, int nclasses, int passcount, bool& err);
static void basictest2(bool& err);
static void basictest3(bool& err);
static void basictest4(bool& err);
static void basictest5(bool& err);
static double rnormal();
static double rsphere(ap::real_2d_array& xy, int n, int i);
static void unsetdf(decisionforest& df);

bool testforest(bool silent)
{
    bool result;
    int ncmax;
    int nvmax;
    int passcount;
    int nvars;
    int nclasses;
    bool waserrors;
    bool basicerrors;
    bool procerrors;
    int i;
    int j;

    
    //
    // Primary settings
    //
    nvmax = 4;
    ncmax = 3;
    passcount = 10;
    basicerrors = false;
    procerrors = false;
    waserrors = false;
    
    //
    // Tests
    //
    testprocessing(procerrors);
    for(nvars = 1; nvars <= nvmax; nvars++)
    {
        for(nclasses = 1; nclasses <= ncmax; nclasses++)
        {
            basictest1(nvars, nclasses, passcount, basicerrors);
        }
    }
    basictest2(basicerrors);
    basictest3(basicerrors);
    basictest4(basicerrors);
    basictest5(basicerrors);
    
    //
    // Final report
    //
    waserrors = basicerrors||procerrors;
    if( !silent )
    {
        printf("RANDOM FOREST TEST\n");
        printf("TOTAL RESULTS:                           ");
        if( !waserrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("* PROCESSING FUNCTIONS:                  ");
        if( !procerrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("* BASIC TESTS:                           ");
        if( !basicerrors )
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
Processing functions test
*************************************************************************/
static void testprocessing(bool& err)
{
    int nvars;
    int nclasses;
    int nsample;
    int ntrees;
    int nfeatures;
    int flags;
    decisionforest df1;
    decisionforest df2;
    int npoints;
    ap::real_2d_array xy;
    int pass;
    int passcount;
    int i;
    int j;
    bool allsame;
    int rlen;
    int info;
    dfreport rep;
    ap::real_1d_array x1;
    ap::real_1d_array x2;
    ap::real_1d_array y1;
    ap::real_1d_array y2;
    ap::real_1d_array ra;
    ap::real_1d_array ra2;
    double v;

    passcount = 100;
    
    //
    // Main cycle
    //
    for(pass = 1; pass <= passcount; pass++)
    {
        
        //
        // initialize parameters
        //
        nvars = 1+ap::randominteger(5);
        nclasses = 1+ap::randominteger(3);
        ntrees = 1+ap::randominteger(4);
        nfeatures = 1+ap::randominteger(nvars);
        flags = 0;
        if( ap::fp_greater(ap::randomreal(),0.5) )
        {
            flags = flags+2;
        }
        
        //
        // Initialize arrays and data
        //
        npoints = 10+ap::randominteger(50);
        nsample = ap::maxint(10, ap::randominteger(npoints));
        x1.setbounds(0, nvars-1);
        x2.setbounds(0, nvars-1);
        y1.setbounds(0, nclasses-1);
        y2.setbounds(0, nclasses-1);
        xy.setbounds(0, npoints-1, 0, nvars);
        for(i = 0; i <= npoints-1; i++)
        {
            for(j = 0; j <= nvars-1; j++)
            {
                if( j%2==0 )
                {
                    xy(i,j) = 2*ap::randomreal()-1;
                }
                else
                {
                    xy(i,j) = ap::randominteger(2);
                }
            }
            if( nclasses==1 )
            {
                xy(i,nvars) = 2*ap::randomreal()-1;
            }
            else
            {
                xy(i,nvars) = ap::randominteger(nclasses);
            }
        }
        
        //
        // create forest
        //
        dfbuildinternal(xy, npoints, nvars, nclasses, ntrees, nsample, nfeatures, flags, info, df1, rep);
        if( info<=0 )
        {
            err = true;
            return;
        }
        
        //
        // Same inputs leads to same outputs
        //
        for(i = 0; i <= nvars-1; i++)
        {
            x1(i) = 2*ap::randomreal()-1;
            x2(i) = x1(i);
        }
        for(i = 0; i <= nclasses-1; i++)
        {
            y1(i) = 2*ap::randomreal()-1;
            y2(i) = 2*ap::randomreal()-1;
        }
        dfprocess(df1, x1, y1);
        dfprocess(df1, x2, y2);
        allsame = true;
        for(i = 0; i <= nclasses-1; i++)
        {
            allsame = allsame&&ap::fp_eq(y1(i),y2(i));
        }
        err = err||!allsame;
        
        //
        // Same inputs on original forest leads to same outputs
        // on copy created using DFCopy
        //
        unsetdf(df2);
        dfcopy(df1, df2);
        for(i = 0; i <= nvars-1; i++)
        {
            x1(i) = 2*ap::randomreal()-1;
            x2(i) = x1(i);
        }
        for(i = 0; i <= nclasses-1; i++)
        {
            y1(i) = 2*ap::randomreal()-1;
            y2(i) = 2*ap::randomreal()-1;
        }
        dfprocess(df1, x1, y1);
        dfprocess(df2, x2, y2);
        allsame = true;
        for(i = 0; i <= nclasses-1; i++)
        {
            allsame = allsame&&ap::fp_eq(y1(i),y2(i));
        }
        err = err||!allsame;
        
        //
        // Same inputs on original forest leads to same outputs
        // on copy created using DFSerialize
        //
        unsetdf(df2);
        ra.setbounds(0, 0);
        ra(0) = 0;
        rlen = 0;
        dfserialize(df1, ra, rlen);
        ra2.setbounds(0, rlen-1);
        for(i = 0; i <= rlen-1; i++)
        {
            ra2(i) = ra(i);
        }
        dfunserialize(ra2, df2);
        for(i = 0; i <= nvars-1; i++)
        {
            x1(i) = 2*ap::randomreal()-1;
            x2(i) = x1(i);
        }
        for(i = 0; i <= nclasses-1; i++)
        {
            y1(i) = 2*ap::randomreal()-1;
            y2(i) = 2*ap::randomreal()-1;
        }
        dfprocess(df1, x1, y1);
        dfprocess(df2, x2, y2);
        allsame = true;
        for(i = 0; i <= nclasses-1; i++)
        {
            allsame = allsame&&ap::fp_eq(y1(i),y2(i));
        }
        err = err||!allsame;
        
        //
        // Normalization properties
        //
        if( nclasses>1 )
        {
            for(i = 0; i <= nvars-1; i++)
            {
                x1(i) = 2*ap::randomreal()-1;
            }
            dfprocess(df1, x1, y1);
            v = 0;
            for(i = 0; i <= nclasses-1; i++)
            {
                v = v+y1(i);
                err = err||ap::fp_less(y1(i),0);
            }
            err = err||ap::fp_greater(fabs(v-1),1000*ap::machineepsilon);
        }
    }
}


/*************************************************************************
Basic test:  one-tree forest built using full sample must remember all the
training cases
*************************************************************************/
static void basictest1(int nvars, int nclasses, int passcount, bool& err)
{
    int pass;
    ap::real_2d_array xy;
    int npoints;
    int i;
    int j;
    int k;
    double s;
    int info;
    decisionforest df;
    ap::real_1d_array x;
    ap::real_1d_array y;
    dfreport rep;
    bool hassame;

    if( nclasses==1 )
    {
        
        //
        // only classification tasks
        //
        return;
    }
    for(pass = 1; pass <= passcount; pass++)
    {
        
        //
        // select number of points
        //
        if( pass<=3&&passcount>3 )
        {
            npoints = pass;
        }
        else
        {
            npoints = 100+ap::randominteger(100);
        }
        
        //
        // Prepare task
        //
        xy.setbounds(0, npoints-1, 0, nvars);
        x.setbounds(0, nvars-1);
        y.setbounds(0, nclasses-1);
        for(i = 0; i <= npoints-1; i++)
        {
            for(j = 0; j <= nvars-1; j++)
            {
                xy(i,j) = 2*ap::randomreal()-1;
            }
            xy(i,nvars) = ap::randominteger(nclasses);
        }
        
        //
        // Test
        //
        dfbuildinternal(xy, npoints, nvars, nclasses, 1, npoints, 1, 1, info, df, rep);
        if( info<=0 )
        {
            err = true;
            return;
        }
        for(i = 0; i <= npoints-1; i++)
        {
            ap::vmove(&x(0), 1, &xy(i, 0), 1, ap::vlen(0,nvars-1));
            dfprocess(df, x, y);
            s = 0;
            for(j = 0; j <= nclasses-1; j++)
            {
                if( ap::fp_less(y(j),0) )
                {
                    err = true;
                    return;
                }
                s = s+y(j);
            }
            if( ap::fp_greater(fabs(s-1),1000*ap::machineepsilon) )
            {
                err = true;
                return;
            }
            if( ap::fp_greater(fabs(y(ap::round(xy(i,nvars)))-1),1000*ap::machineepsilon) )
            {
                
                //
                // not an error if there exists such K,J that XY[K,J]=XY[I,J]
                // (may be we just can't distinguish two tied values).
                //
                // definitely error otherwise.
                //
                hassame = false;
                for(k = 0; k <= npoints-1; k++)
                {
                    if( k!=i )
                    {
                        for(j = 0; j <= nvars-1; j++)
                        {
                            if( ap::fp_eq(xy(k,j),xy(i,j)) )
                            {
                                hassame = true;
                            }
                        }
                    }
                }
                if( !hassame )
                {
                    err = true;
                    return;
                }
            }
        }
    }
}


/*************************************************************************
Basic test:  tests generalization ability on a simple noisy classification
task:
* 0<x<1 - P(class=0)=1
* 1<x<2 - P(class=0)=2-x
* 2<x<3 - P(class=0)=0
*************************************************************************/
static void basictest2(bool& err)
{
    int pass;
    int passcount;
    ap::real_2d_array xy;
    int npoints;
    int ntrees;
    int i;
    int j;
    double s;
    int info;
    decisionforest df;
    ap::real_1d_array x;
    ap::real_1d_array y;
    dfreport rep;
    bool hassame;

    passcount = 1;
    for(pass = 1; pass <= passcount; pass++)
    {
        
        //
        // select npoints and ntrees
        //
        npoints = 3000;
        ntrees = 50;
        
        //
        // Prepare task
        //
        xy.setbounds(0, npoints-1, 0, 1);
        x.setbounds(0, 0);
        y.setbounds(0, 1);
        for(i = 0; i <= npoints-1; i++)
        {
            xy(i,0) = 3*ap::randomreal();
            if( ap::fp_less_eq(xy(i,0),1) )
            {
                xy(i,1) = 0;
            }
            else
            {
                if( ap::fp_less_eq(xy(i,0),2) )
                {
                    if( ap::fp_less(ap::randomreal(),xy(i,0)-1) )
                    {
                        xy(i,1) = 1;
                    }
                    else
                    {
                        xy(i,1) = 0;
                    }
                }
                else
                {
                    xy(i,1) = 1;
                }
            }
        }
        
        //
        // Test
        //
        dfbuildinternal(xy, npoints, 1, 2, ntrees, ap::round(0.05*npoints), 1, 0, info, df, rep);
        if( info<=0 )
        {
            err = true;
            return;
        }
        x(0) = 0.0;
        while(ap::fp_less_eq(x(0),3.0))
        {
            dfprocess(df, x, y);
            
            //
            // Test for basic properties
            //
            s = 0;
            for(j = 0; j <= 1; j++)
            {
                if( ap::fp_less(y(j),0) )
                {
                    err = true;
                    return;
                }
                s = s+y(j);
            }
            if( ap::fp_greater(fabs(s-1),1000*ap::machineepsilon) )
            {
                err = true;
                return;
            }
            
            //
            // test for good correlation with results
            //
            if( ap::fp_less(x(0),1) )
            {
                err = err||ap::fp_less(y(0),0.8);
            }
            if( ap::fp_greater_eq(x(0),1)&&ap::fp_less_eq(x(0),2) )
            {
                err = err||ap::fp_greater(fabs(y(1)-(x(0)-1)),0.5);
            }
            if( ap::fp_greater(x(0),2) )
            {
                err = err||ap::fp_less(y(1),0.8);
            }
            x(0) = x(0)+0.01;
        }
    }
}


/*************************************************************************
Basic test:  tests  generalization ability on a simple classification task
(no noise):
* |x|<1, |y|<1
* x^2+y^2<=0.25 - P(class=0)=1
* x^2+y^2>0.25  - P(class=0)=0
*************************************************************************/
static void basictest3(bool& err)
{
    int pass;
    int passcount;
    ap::real_2d_array xy;
    int npoints;
    int ntrees;
    int i;
    int j;
    int k;
    double s;
    int info;
    decisionforest df;
    ap::real_1d_array x;
    ap::real_1d_array y;
    dfreport rep;
    int testgridsize;
    double r;

    passcount = 1;
    testgridsize = 50;
    for(pass = 1; pass <= passcount; pass++)
    {
        
        //
        // select npoints and ntrees
        //
        npoints = 2000;
        ntrees = 100;
        
        //
        // Prepare task
        //
        xy.setbounds(0, npoints-1, 0, 2);
        x.setbounds(0, 1);
        y.setbounds(0, 1);
        for(i = 0; i <= npoints-1; i++)
        {
            xy(i,0) = 2*ap::randomreal()-1;
            xy(i,1) = 2*ap::randomreal()-1;
            if( ap::fp_less_eq(ap::sqr(xy(i,0))+ap::sqr(xy(i,1)),0.25) )
            {
                xy(i,2) = 0;
            }
            else
            {
                xy(i,2) = 1;
            }
        }
        
        //
        // Test
        //
        dfbuildinternal(xy, npoints, 2, 2, ntrees, ap::round(0.1*npoints), 1, 0, info, df, rep);
        if( info<=0 )
        {
            err = true;
            return;
        }
        for(i = -testgridsize/2; i <= testgridsize/2; i++)
        {
            for(j = -testgridsize/2; j <= testgridsize/2; j++)
            {
                x(0) = double(i)/double(testgridsize/2);
                x(1) = double(j)/double(testgridsize/2);
                dfprocess(df, x, y);
                
                //
                // Test for basic properties
                //
                s = 0;
                for(k = 0; k <= 1; k++)
                {
                    if( ap::fp_less(y(k),0) )
                    {
                        err = true;
                        return;
                    }
                    s = s+y(k);
                }
                if( ap::fp_greater(fabs(s-1),1000*ap::machineepsilon) )
                {
                    err = true;
                    return;
                }
                
                //
                // test for good correlation with results
                //
                r = sqrt(ap::sqr(x(0))+ap::sqr(x(1)));
                if( ap::fp_less(r,0.5*0.5) )
                {
                    err = err||ap::fp_less(y(0),0.6);
                }
                if( ap::fp_greater(r,0.5*1.5) )
                {
                    err = err||ap::fp_less(y(1),0.6);
                }
            }
        }
    }
}


/*************************************************************************
Basic test: simple regression task without noise:
* |x|<1, |y|<1
* F(x,y) = x^2+y
*************************************************************************/
static void basictest4(bool& err)
{
    int pass;
    int passcount;
    ap::real_2d_array xy;
    int npoints;
    int ntrees;
    int ns;
    int strongc;
    int i;
    int j;
    int k;
    double s;
    int info;
    decisionforest df;
    decisionforest df2;
    ap::real_1d_array x;
    ap::real_1d_array y;
    dfreport rep;
    dfreport rep2;
    int testgridsize;
    double maxerr;
    double maxerr2;
    double avgerr;
    double avgerr2;
    int cnt;
    double ey;

    passcount = 1;
    testgridsize = 50;
    for(pass = 1; pass <= passcount; pass++)
    {
        
        //
        // select npoints and ntrees
        //
        npoints = 5000;
        ntrees = 100;
        ns = ap::round(0.1*npoints);
        strongc = 1;
        
        //
        // Prepare task
        //
        xy.setbounds(0, npoints-1, 0, 2);
        x.setbounds(0, 1);
        y.setbounds(0, 0);
        for(i = 0; i <= npoints-1; i++)
        {
            xy(i,0) = 2*ap::randomreal()-1;
            xy(i,1) = 2*ap::randomreal()-1;
            xy(i,2) = ap::sqr(xy(i,0))+xy(i,1);
        }
        
        //
        // Test
        //
        dfbuildinternal(xy, npoints, 2, 1, ntrees, ns, 1, 0, info, df, rep);
        if( info<=0 )
        {
            err = true;
            return;
        }
        dfbuildinternal(xy, npoints, 2, 1, ntrees, ns, 1, strongc, info, df2, rep2);
        if( info<=0 )
        {
            err = true;
            return;
        }
        maxerr = 0;
        maxerr2 = 0;
        avgerr = 0;
        avgerr2 = 0;
        cnt = 0;
        for(i = ap::round(-0.7*testgridsize/2); i <= ap::round(0.7*testgridsize/2); i++)
        {
            for(j = ap::round(-0.7*testgridsize/2); j <= ap::round(0.7*testgridsize/2); j++)
            {
                x(0) = double(i)/double(testgridsize/2);
                x(1) = double(j)/double(testgridsize/2);
                ey = ap::sqr(x(0))+x(1);
                dfprocess(df, x, y);
                maxerr = ap::maxreal(maxerr, fabs(y(0)-ey));
                avgerr = avgerr+fabs(y(0)-ey);
                dfprocess(df2, x, y);
                maxerr2 = ap::maxreal(maxerr2, fabs(y(0)-ey));
                avgerr2 = avgerr2+fabs(y(0)-ey);
                cnt = cnt+1;
            }
        }
        avgerr = avgerr/cnt;
        avgerr2 = avgerr2/cnt;
        err = err||ap::fp_greater(maxerr,0.2);
        err = err||ap::fp_greater(maxerr2,0.2);
        err = err||ap::fp_greater(avgerr,0.1);
        err = err||ap::fp_greater(avgerr2,0.1);
    }
}


/*************************************************************************
Basic test: extended variable selection leads to better results.

Next task CAN be solved without EVS but it is very unlikely. With EVS
it can be easily and exactly solved.

Task matrix:
    1 0 0 0 ... 0   0
    0 1 0 0 ... 0   1
    0 0 1 0 ... 0   2
    0 0 0 1 ... 0   3
    0 0 0 0 ... 1   N-1
*************************************************************************/
static void basictest5(bool& err)
{
    ap::real_2d_array xy;
    int nvars;
    int npoints;
    int nfeatures;
    int nsample;
    int ntrees;
    int evs;
    int i;
    int j;
    bool eflag;
    int info;
    decisionforest df;
    ap::real_1d_array x;
    ap::real_1d_array y;
    dfreport rep;

    
    //
    // select npoints and ntrees
    //
    npoints = 50;
    nvars = npoints;
    ntrees = 1;
    nsample = npoints;
    evs = 2;
    nfeatures = 1;
    
    //
    // Prepare task
    //
    xy.setbounds(0, npoints-1, 0, nvars);
    x.setbounds(0, nvars-1);
    y.setbounds(0, 0);
    for(i = 0; i <= npoints-1; i++)
    {
        for(j = 0; j <= nvars-1; j++)
        {
            xy(i,j) = 0;
        }
        xy(i,i) = 1;
        xy(i,nvars) = i;
    }
    
    //
    // Without EVS
    //
    dfbuildinternal(xy, npoints, nvars, 1, ntrees, nsample, nfeatures, 0, info, df, rep);
    if( info<=0 )
    {
        err = true;
        return;
    }
    eflag = false;
    for(i = 0; i <= npoints-1; i++)
    {
        ap::vmove(&x(0), 1, &xy(i, 0), 1, ap::vlen(0,nvars-1));
        dfprocess(df, x, y);
        if( ap::fp_greater(fabs(y(0)-xy(i,nvars)),1000*ap::machineepsilon) )
        {
            eflag = true;
        }
    }
    if( !eflag )
    {
        err = true;
        return;
    }
    
    //
    // With EVS
    //
    dfbuildinternal(xy, npoints, nvars, 1, ntrees, nsample, nfeatures, evs, info, df, rep);
    if( info<=0 )
    {
        err = true;
        return;
    }
    eflag = false;
    for(i = 0; i <= npoints-1; i++)
    {
        ap::vmove(&x(0), 1, &xy(i, 0), 1, ap::vlen(0,nvars-1));
        dfprocess(df, x, y);
        if( ap::fp_greater(fabs(y(0)-xy(i,nvars)),1000*ap::machineepsilon) )
        {
            eflag = true;
        }
    }
    if( eflag )
    {
        err = true;
        return;
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
Unsets DF
*************************************************************************/
static void unsetdf(decisionforest& df)
{
    ap::real_2d_array xy;
    int info;
    dfreport rep;

    xy.setbounds(0, 0, 0, 1);
    xy(0,0) = 0;
    xy(0,1) = 0;
    dfbuildinternal(xy, 1, 1, 1, 1, 1, 1, 0, info, df, rep);
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testforestunit_test_silent()
{
    bool result;

    result = testforest(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testforestunit_test()
{
    bool result;

    result = testforest(false);
    return result;
}




