
 
#include <stdio.h>
#include "testhqrndunit.h"

static void unsetstate(hqrndstate& state);

void calculatemv(const ap::real_1d_array& x,
     int n,
     double& mean,
     double& means,
     double& stddev,
     double& stddevs)
{
    int i;
    double v;
    double v1;
    double v2;
    double variance;

    mean = 0;
    means = 1;
    stddev = 0;
    stddevs = 1;
    variance = 0;
    if( n<=1 )
    {
        return;
    }
    
    //
    // Mean
    //
    for(i = 0; i <= n-1; i++)
    {
        mean = mean+x(i);
    }
    mean = mean/n;
    
    //
    // Variance (using corrected two-pass algorithm)
    //
    if( n!=1 )
    {
        v1 = 0;
        for(i = 0; i <= n-1; i++)
        {
            v1 = v1+ap::sqr(x(i)-mean);
        }
        v2 = 0;
        for(i = 0; i <= n-1; i++)
        {
            v2 = v2+(x(i)-mean);
        }
        v2 = ap::sqr(v2)/n;
        variance = (v1-v2)/(n-1);
        if( ap::fp_less(variance,0) )
        {
            variance = 0;
        }
        stddev = sqrt(variance);
    }
    
    //
    // Errors
    //
    means = stddev/sqrt(double(n));
    stddevs = stddev*sqrt(double(2))/sqrt(double(n-1));
}


bool testhqrnd(bool silent)
{
    bool result;
    bool waserrors;
    int samplesize;
    double sigmathreshold;
    int passcount;
    int n;
    int i;
    int pass;
    int s1;
    int s2;
    int i1;
    int i2;
    double r1;
    double r2;
    ap::real_1d_array x;
    double mean;
    double means;
    double stddev;
    double stddevs;
    double lambda;
    bool seederrors;
    bool urerrors;
    double ursigmaerr;
    bool uierrors;
    double uisigmaerr;
    bool normerrors;
    double normsigmaerr;
    bool experrors;
    double expsigmaerr;
    hqrndstate state;

    waserrors = false;
    sigmathreshold = 7;
    samplesize = 100000;
    passcount = 50;
    x.setbounds(0, samplesize-1);
    
    //
    // Test seed errors
    //
    seederrors = false;
    for(pass = 1; pass <= passcount; pass++)
    {
        s1 = 1+ap::randominteger(32000);
        s2 = 1+ap::randominteger(32000);
        unsetstate(state);
        hqrndseed(s1, s2, state);
        i1 = hqrnduniformi(100, state);
        unsetstate(state);
        hqrndseed(s1, s2, state);
        i2 = hqrnduniformi(100, state);
        seederrors = seederrors||i1!=i2;
        unsetstate(state);
        hqrndseed(s1, s2, state);
        r1 = hqrnduniformr(state);
        unsetstate(state);
        hqrndseed(s1, s2, state);
        r2 = hqrnduniformr(state);
        seederrors = seederrors||ap::fp_neq(r1,r2);
    }
    
    //
    // Test HQRNDRandomize() and real uniform generator
    //
    unsetstate(state);
    hqrndrandomize(state);
    urerrors = false;
    ursigmaerr = 0;
    for(i = 0; i <= samplesize-1; i++)
    {
        x(i) = hqrnduniformr(state);
    }
    for(i = 0; i <= samplesize-1; i++)
    {
        urerrors = urerrors||ap::fp_less_eq(x(i),0)||ap::fp_greater_eq(x(i),1);
    }
    calculatemv(x, samplesize, mean, means, stddev, stddevs);
    if( ap::fp_neq(means,0) )
    {
        ursigmaerr = ap::maxreal(ursigmaerr, fabs((mean-0.5)/means));
    }
    else
    {
        urerrors = true;
    }
    if( ap::fp_neq(stddevs,0) )
    {
        ursigmaerr = ap::maxreal(ursigmaerr, fabs((stddev-sqrt(double(1)/double(12)))/stddevs));
    }
    else
    {
        urerrors = true;
    }
    urerrors = urerrors||ap::fp_greater(ursigmaerr,sigmathreshold);
    
    //
    // Test HQRNDRandomize() and integer uniform
    //
    unsetstate(state);
    hqrndrandomize(state);
    uierrors = false;
    uisigmaerr = 0;
    for(n = 2; n <= 10; n++)
    {
        for(i = 0; i <= samplesize-1; i++)
        {
            x(i) = hqrnduniformi(n, state);
        }
        for(i = 0; i <= samplesize-1; i++)
        {
            uierrors = uierrors||ap::fp_less(x(i),0)||ap::fp_greater_eq(x(i),n);
        }
        calculatemv(x, samplesize, mean, means, stddev, stddevs);
        if( ap::fp_neq(means,0) )
        {
            uisigmaerr = ap::maxreal(uisigmaerr, fabs((mean-0.5*(n-1))/means));
        }
        else
        {
            uierrors = true;
        }
        if( ap::fp_neq(stddevs,0) )
        {
            uisigmaerr = ap::maxreal(uisigmaerr, fabs((stddev-sqrt((ap::sqr(double(n))-1)/12))/stddevs));
        }
        else
        {
            uierrors = true;
        }
    }
    uierrors = uierrors||ap::fp_greater(uisigmaerr,sigmathreshold);
    
    //
    // Special 'close-to-limit' test on uniformity of integers
    // (straightforward implementation like 'RND mod N' will return
    //  non-uniform numbers for N=2/3*LIMIT)
    //
    unsetstate(state);
    hqrndrandomize(state);
    uierrors = false;
    uisigmaerr = 0;
    n = ap::round(2.0/3.0*2147483563.0);
    for(i = 0; i <= samplesize-1; i++)
    {
        x(i) = hqrnduniformi(n, state);
    }
    for(i = 0; i <= samplesize-1; i++)
    {
        uierrors = uierrors||ap::fp_less(x(i),0)||ap::fp_greater_eq(x(i),n);
    }
    calculatemv(x, samplesize, mean, means, stddev, stddevs);
    if( ap::fp_neq(means,0) )
    {
        uisigmaerr = ap::maxreal(uisigmaerr, fabs((mean-0.5*(n-1))/means));
    }
    else
    {
        uierrors = true;
    }
    if( ap::fp_neq(stddevs,0) )
    {
        uisigmaerr = ap::maxreal(uisigmaerr, fabs((stddev-sqrt((ap::sqr(double(n))-1)/12))/stddevs));
    }
    else
    {
        uierrors = true;
    }
    uierrors = uierrors||ap::fp_greater(uisigmaerr,sigmathreshold);
    
    //
    // Test normal
    //
    unsetstate(state);
    hqrndrandomize(state);
    normerrors = false;
    normsigmaerr = 0;
    i = 0;
    while(i<samplesize)
    {
        hqrndnormal2(state, r1, r2);
        x(i) = r1;
        if( i+1<samplesize )
        {
            x(i+1) = r2;
        }
        i = i+2;
    }
    calculatemv(x, samplesize, mean, means, stddev, stddevs);
    if( ap::fp_neq(means,0) )
    {
        normsigmaerr = ap::maxreal(normsigmaerr, fabs((mean-0)/means));
    }
    else
    {
        normerrors = true;
    }
    if( ap::fp_neq(stddevs,0) )
    {
        normsigmaerr = ap::maxreal(normsigmaerr, fabs((stddev-1)/stddevs));
    }
    else
    {
        normerrors = true;
    }
    normerrors = normerrors||ap::fp_greater(normsigmaerr,sigmathreshold);
    
    //
    // Test exponential
    //
    unsetstate(state);
    hqrndrandomize(state);
    experrors = false;
    expsigmaerr = 0;
    lambda = 2+5*ap::randomreal();
    for(i = 0; i <= samplesize-1; i++)
    {
        x(i) = hqrndexponential(lambda, state);
    }
    for(i = 0; i <= samplesize-1; i++)
    {
        uierrors = uierrors||ap::fp_less(x(i),0);
    }
    calculatemv(x, samplesize, mean, means, stddev, stddevs);
    if( ap::fp_neq(means,0) )
    {
        expsigmaerr = ap::maxreal(expsigmaerr, fabs((mean-1.0/lambda)/means));
    }
    else
    {
        experrors = true;
    }
    if( ap::fp_neq(stddevs,0) )
    {
        expsigmaerr = ap::maxreal(expsigmaerr, fabs((stddev-1.0/lambda)/stddevs));
    }
    else
    {
        experrors = true;
    }
    experrors = experrors||ap::fp_greater(expsigmaerr,sigmathreshold);
    
    //
    // Final report
    //
    waserrors = seederrors||urerrors||uierrors||normerrors||experrors;
    if( !silent )
    {
        printf("RNG TEST\n");
        printf("SEED TEST:                               ");
        if( !seederrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("UNIFORM CONTINUOUS:                      ");
        if( !urerrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("UNIFORM INTEGER:                         ");
        if( !uierrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("NORMAL:                                  ");
        if( !normerrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("EXPONENTIAL:                             ");
        if( !experrors )
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
Unsets HQRNDState structure
*************************************************************************/
static void unsetstate(hqrndstate& state)
{

    state.s1 = 0;
    state.s2 = 0;
    state.v = 0;
    state.magicv = 0;
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testhqrndunit_test_silent()
{
    bool result;

    result = testhqrnd(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testhqrndunit_test()
{
    bool result;

    result = testhqrnd(false);
    return result;
}




