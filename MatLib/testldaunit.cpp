
 
#include <stdio.h>
#include "testldaunit.h"

static void gensimpleset(int nfeatures,
     int nclasses,
     int nsamples,
     int axis,
     ap::real_2d_array& xy);
static void gendeg1set(int nfeatures,
     int nclasses,
     int nsamples,
     int axis,
     ap::real_2d_array& xy);
static double generatenormal(double mean, double sigma);
static bool testwn(const ap::real_2d_array& xy,
     const ap::real_2d_array& wn,
     int ns,
     int nf,
     int nc,
     int ndeg);
static double calcj(int nf,
     const ap::real_2d_array& st,
     const ap::real_2d_array& sw,
     const ap::real_1d_array& w,
     double& p,
     double& q);
static void fishers(const ap::real_2d_array& xy,
     int npoints,
     int nfeatures,
     int nclasses,
     ap::real_2d_array& st,
     ap::real_2d_array& sw);

bool testlda(bool silent)
{
    bool result;
    int maxnf;
    int maxns;
    int maxnc;
    int passcount;
    bool ldanerrors;
    bool lda1errors;
    bool waserrors;
    int nf;
    int nc;
    int ns;
    int i;
    int info;
    int pass;
    int axis;
    ap::real_2d_array xy;
    ap::real_2d_array wn;
    ap::real_1d_array w1;

    
    //
    // Primary settings
    //
    maxnf = 10;
    maxns = 1000;
    maxnc = 5;
    passcount = 1;
    waserrors = false;
    ldanerrors = false;
    lda1errors = false;
    
    //
    // General tests
    //
    for(nf = 1; nf <= maxnf; nf++)
    {
        for(nc = 2; nc <= maxnc; nc++)
        {
            for(pass = 1; pass <= passcount; pass++)
            {
                
                //
                // Simple test for LDA-N/LDA-1
                //
                axis = ap::randominteger(nf);
                ns = maxns/2+ap::randominteger(maxns/2);
                gensimpleset(nf, nc, ns, axis, xy);
                fisherldan(xy, ns, nf, nc, info, wn);
                if( info!=1 )
                {
                    ldanerrors = true;
                    continue;
                }
                ldanerrors = ldanerrors||!testwn(xy, wn, ns, nf, nc, 0);
                ldanerrors = ldanerrors||ap::fp_less_eq(fabs(wn(axis,0)),0.75);
                fisherlda(xy, ns, nf, nc, info, w1);
                for(i = 0; i <= nf-1; i++)
                {
                    lda1errors = lda1errors||ap::fp_neq(w1(i),wn(i,0));
                }
                
                //
                // Degenerate test for LDA-N
                //
                if( nf>=3 )
                {
                    ns = maxns/2+ap::randominteger(maxns/2);
                    
                    //
                    // there are two duplicate features,
                    // axis is oriented along non-duplicate feature
                    //
                    axis = ap::randominteger(nf-2);
                    gendeg1set(nf, nc, ns, axis, xy);
                    fisherldan(xy, ns, nf, nc, info, wn);
                    if( info!=2 )
                    {
                        ldanerrors = true;
                        continue;
                    }
                    ldanerrors = ldanerrors||ap::fp_less_eq(wn(axis,0),0.75);
                    fisherlda(xy, ns, nf, nc, info, w1);
                    for(i = 0; i <= nf-1; i++)
                    {
                        lda1errors = lda1errors||ap::fp_neq(w1(i),wn(i,0));
                    }
                }
            }
        }
    }
    
    //
    // Final report
    //
    waserrors = ldanerrors||lda1errors;
    if( !silent )
    {
        printf("LDA TEST\n");
        printf("FISHER LDA-N:                            ");
        if( !ldanerrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("FISHER LDA-1:                            ");
        if( !lda1errors )
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
Generates 'simple' set - a sequence of unit 'balls' at (0,0), (1,0), (2,0)
and so on.
*************************************************************************/
static void gensimpleset(int nfeatures,
     int nclasses,
     int nsamples,
     int axis,
     ap::real_2d_array& xy)
{
    int i;
    int j;
    int c;
    double v;

    ap::ap_error::make_assertion(axis>=0&&axis<nfeatures, "GenSimpleSet: wrong Axis!");
    xy.setbounds(0, nsamples-1, 0, nfeatures);
    for(i = 0; i <= nsamples-1; i++)
    {
        for(j = 0; j <= nfeatures-1; j++)
        {
            xy(i,j) = generatenormal(0.0, 1.0);
        }
        c = i%nclasses;
        xy(i,axis) = xy(i,axis)+c;
        xy(i,nfeatures) = c;
    }
}


/*************************************************************************
Generates 'degenerate' set #1.
NFeatures>=3.
*************************************************************************/
static void gendeg1set(int nfeatures,
     int nclasses,
     int nsamples,
     int axis,
     ap::real_2d_array& xy)
{
    int i;
    int j;
    int c;
    double v;

    ap::ap_error::make_assertion(axis>=0&&axis<nfeatures, "GenDeg1Set: wrong Axis!");
    ap::ap_error::make_assertion(nfeatures>=3, "GenDeg1Set: wrong NFeatures!");
    xy.setbounds(0, nsamples-1, 0, nfeatures);
    if( axis>=nfeatures-2 )
    {
        axis = nfeatures-3;
    }
    for(i = 0; i <= nsamples-1; i++)
    {
        for(j = 0; j <= nfeatures-2; j++)
        {
            xy(i,j) = generatenormal(0.0, 1.0);
        }
        xy(i,nfeatures-1) = xy(i,nfeatures-2);
        c = i%nclasses;
        xy(i,axis) = xy(i,axis)+c;
        xy(i,nfeatures) = c;
    }
}


/*************************************************************************
Normal random number
*************************************************************************/
static double generatenormal(double mean, double sigma)
{
    double result;
    double u;
    double v;
    double s;
    double sum;

    result = mean;
    while(true)
    {
        u = (2*ap::randominteger(2)-1)*ap::randomreal();
        v = (2*ap::randominteger(2)-1)*ap::randomreal();
        sum = u*u+v*v;
        if( ap::fp_less(sum,1)&&ap::fp_greater(sum,0) )
        {
            sum = sqrt(-2*log(sum)/sum);
            result = sigma*u*sum+mean;
            return result;
        }
    }
    return result;
}


/*************************************************************************
Tests WN for correctness
*************************************************************************/
static bool testwn(const ap::real_2d_array& xy,
     const ap::real_2d_array& wn,
     int ns,
     int nf,
     int nc,
     int ndeg)
{
    bool result;
    ap::real_2d_array st;
    ap::real_2d_array sw;
    ap::real_2d_array a;
    ap::real_2d_array z;
    ap::real_1d_array tx;
    ap::real_1d_array jp;
    ap::real_1d_array jq;
    ap::real_1d_array work;
    int i;
    int j;
    double v;
    double wprev;
    double tol;
    double p;
    double q;

    tol = 10000;
    result = true;
    fishers(xy, ns, nf, nc, st, sw);
    
    //
    // Test for decreasing of J
    //
    tx.setbounds(0, nf-1);
    jp.setbounds(0, nf-1);
    jq.setbounds(0, nf-1);
    for(j = 0; j <= nf-1; j++)
    {
        ap::vmove(&tx(0), 1, &wn(0, j), wn.getstride(), ap::vlen(0,nf-1));
        v = calcj(nf, st, sw, tx, p, q);
        jp(j) = p;
        jq(j) = q;
    }
    for(i = 1; i <= nf-1-ndeg; i++)
    {
        result = result&&ap::fp_greater_eq(jp(i-1)/jq(i-1),(1-tol*ap::machineepsilon)*jp(i)/jq(i));
    }
    for(i = nf-1-ndeg+1; i <= nf-1; i++)
    {
        result = result&&ap::fp_less_eq(jp(i),tol*ap::machineepsilon*jp(0));
    }
    
    //
    // Test for J optimality
    //
    ap::vmove(&tx(0), 1, &wn(0, 0), wn.getstride(), ap::vlen(0,nf-1));
    v = calcj(nf, st, sw, tx, p, q);
    for(i = 0; i <= nf-1; i++)
    {
        wprev = tx(i);
        tx(i) = wprev+0.01;
        result = result&&ap::fp_greater_eq(v,(1-tol*ap::machineepsilon)*calcj(nf, st, sw, tx, p, q));
        tx(i) = wprev-0.01;
        result = result&&ap::fp_greater_eq(v,(1-tol*ap::machineepsilon)*calcj(nf, st, sw, tx, p, q));
        tx(i) = wprev;
    }
    
    //
    // Test for linear independence of W
    //
    work.setbounds(1, nf);
    a.setbounds(0, nf-1, 0, nf-1);
    matrixmatrixmultiply(wn, 0, nf-1, 0, nf-1, false, wn, 0, nf-1, 0, nf-1, true, 1.0, a, 0, nf-1, 0, nf-1, 0.0, work);
    if( smatrixevd(a, nf, 1, true, tx, z) )
    {
        result = result&&ap::fp_greater(tx(0),tx(nf-1)*1000*ap::machineepsilon);
    }
    
    //
    // Test for other properties
    //
    for(j = 0; j <= nf-1; j++)
    {
        v = ap::vdotproduct(&wn(0, j), wn.getstride(), &wn(0, j), wn.getstride(), ap::vlen(0,nf-1));
        v = sqrt(v);
        result = result&&ap::fp_less_eq(fabs(v-1),1000*ap::machineepsilon);
        v = 0;
        for(i = 0; i <= nf-1; i++)
        {
            v = v+wn(i,j);
        }
        result = result&&ap::fp_greater_eq(v,0);
    }
    return result;
}


/*************************************************************************
Calculates J
*************************************************************************/
static double calcj(int nf,
     const ap::real_2d_array& st,
     const ap::real_2d_array& sw,
     const ap::real_1d_array& w,
     double& p,
     double& q)
{
    double result;
    ap::real_1d_array tx;
    int i;
    int j;
    double v;

    tx.setbounds(0, nf-1);
    for(i = 0; i <= nf-1; i++)
    {
        v = ap::vdotproduct(&st(i, 0), 1, &w(0), 1, ap::vlen(0,nf-1));
        tx(i) = v;
    }
    v = ap::vdotproduct(&w(0), 1, &tx(0), 1, ap::vlen(0,nf-1));
    p = v;
    for(i = 0; i <= nf-1; i++)
    {
        v = ap::vdotproduct(&sw(i, 0), 1, &w(0), 1, ap::vlen(0,nf-1));
        tx(i) = v;
    }
    v = ap::vdotproduct(&w(0), 1, &tx(0), 1, ap::vlen(0,nf-1));
    q = v;
    result = p/q;
    return result;
}


/*************************************************************************
Calculates ST/SW
*************************************************************************/
static void fishers(const ap::real_2d_array& xy,
     int npoints,
     int nfeatures,
     int nclasses,
     ap::real_2d_array& st,
     ap::real_2d_array& sw)
{
    int i;
    int j;
    int k;
    double v;
    ap::integer_1d_array c;
    ap::real_1d_array mu;
    ap::real_2d_array muc;
    ap::integer_1d_array nc;
    ap::real_1d_array tf;
    ap::real_1d_array work;

    
    //
    // Prepare temporaries
    //
    tf.setbounds(0, nfeatures-1);
    work.setbounds(1, nfeatures);
    
    //
    // Convert class labels from reals to integers (just for convenience)
    //
    c.setbounds(0, npoints-1);
    for(i = 0; i <= npoints-1; i++)
    {
        c(i) = ap::round(xy(i,nfeatures));
    }
    
    //
    // Calculate class sizes and means
    //
    mu.setbounds(0, nfeatures-1);
    muc.setbounds(0, nclasses-1, 0, nfeatures-1);
    nc.setbounds(0, nclasses-1);
    for(j = 0; j <= nfeatures-1; j++)
    {
        mu(j) = 0;
    }
    for(i = 0; i <= nclasses-1; i++)
    {
        nc(i) = 0;
        for(j = 0; j <= nfeatures-1; j++)
        {
            muc(i,j) = 0;
        }
    }
    for(i = 0; i <= npoints-1; i++)
    {
        ap::vadd(&mu(0), 1, &xy(i, 0), 1, ap::vlen(0,nfeatures-1));
        ap::vadd(&muc(c(i), 0), 1, &xy(i, 0), 1, ap::vlen(0,nfeatures-1));
        nc(c(i)) = nc(c(i))+1;
    }
    for(i = 0; i <= nclasses-1; i++)
    {
        v = double(1)/double(nc(i));
        ap::vmul(&muc(i, 0), 1, ap::vlen(0,nfeatures-1), v);
    }
    v = double(1)/double(npoints);
    ap::vmul(&mu(0), 1, ap::vlen(0,nfeatures-1), v);
    
    //
    // Create ST matrix
    //
    st.setbounds(0, nfeatures-1, 0, nfeatures-1);
    for(i = 0; i <= nfeatures-1; i++)
    {
        for(j = 0; j <= nfeatures-1; j++)
        {
            st(i,j) = 0;
        }
    }
    for(k = 0; k <= npoints-1; k++)
    {
        ap::vmove(&tf(0), 1, &xy(k, 0), 1, ap::vlen(0,nfeatures-1));
        ap::vsub(&tf(0), 1, &mu(0), 1, ap::vlen(0,nfeatures-1));
        for(i = 0; i <= nfeatures-1; i++)
        {
            v = tf(i);
            ap::vadd(&st(i, 0), 1, &tf(0), 1, ap::vlen(0,nfeatures-1), v);
        }
    }
    
    //
    // Create SW matrix
    //
    sw.setbounds(0, nfeatures-1, 0, nfeatures-1);
    for(i = 0; i <= nfeatures-1; i++)
    {
        for(j = 0; j <= nfeatures-1; j++)
        {
            sw(i,j) = 0;
        }
    }
    for(k = 0; k <= npoints-1; k++)
    {
        ap::vmove(&tf(0), 1, &xy(k, 0), 1, ap::vlen(0,nfeatures-1));
        ap::vsub(&tf(0), 1, &muc(c(k), 0), 1, ap::vlen(0,nfeatures-1));
        for(i = 0; i <= nfeatures-1; i++)
        {
            v = tf(i);
            ap::vadd(&sw(i, 0), 1, &tf(0), 1, ap::vlen(0,nfeatures-1), v);
        }
    }
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testldaunit_test_silent()
{
    bool result;

    result = testlda(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testldaunit_test()
{
    bool result;

    result = testlda(false);
    return result;
}




