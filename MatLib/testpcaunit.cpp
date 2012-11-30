
 
#include <stdio.h>
#include "testpcaunit.h"

static void calculatemv(const ap::real_1d_array& x,
     int n,
     double& mean,
     double& means,
     double& stddev,
     double& stddevs);

bool testpca(bool silent)
{
    bool result;
    int passcount;
    int maxn;
    int maxm;
    double threshold;
    int m;
    int n;
    int i;
    int j;
    int k;
    int info;
    ap::real_1d_array means;
    ap::real_1d_array s;
    ap::real_1d_array t2;
    ap::real_1d_array t3;
    ap::real_2d_array v;
    ap::real_2d_array x;
    double t;
    double h;
    double tmean;
    double tmeans;
    double tstddev;
    double tstddevs;
    double tmean2;
    double tmeans2;
    double tstddev2;
    double tstddevs2;
    bool pcaconverrors;
    bool pcaorterrors;
    bool pcavarerrors;
    bool pcaopterrors;
    bool waserrors;

    
    //
    // Primary settings
    //
    maxm = 10;
    maxn = 100;
    passcount = 1;
    threshold = 1000*ap::machineepsilon;
    waserrors = false;
    pcaconverrors = false;
    pcaorterrors = false;
    pcavarerrors = false;
    pcaopterrors = false;
    
    //
    // Test 1: N random points in M-dimensional space
    //
    for(m = 1; m <= maxm; m++)
    {
        for(n = 1; n <= maxn; n++)
        {
            
            //
            // Generate task
            //
            x.setbounds(0, n-1, 0, m-1);
            means.setbounds(0, m-1);
            for(j = 0; j <= m-1; j++)
            {
                means(j) = 1.5*ap::randomreal()-0.75;
            }
            for(i = 0; i <= n-1; i++)
            {
                for(j = 0; j <= m-1; j++)
                {
                    x(i,j) = means(j)+(2*ap::randomreal()-1);
                }
            }
            
            //
            // Solve
            //
            pcabuildbasis(x, n, m, info, s, v);
            if( info!=1 )
            {
                pcaconverrors = true;
                continue;
            }
            
            //
            // Orthogonality test
            //
            for(i = 0; i <= m-1; i++)
            {
                for(j = 0; j <= m-1; j++)
                {
                    t = ap::vdotproduct(&v(0, i), v.getstride(), &v(0, j), v.getstride(), ap::vlen(0,m-1));
                    if( i==j )
                    {
                        t = t-1;
                    }
                    pcaorterrors = pcaorterrors||ap::fp_greater(fabs(t),threshold);
                }
            }
            
            //
            // Variance test
            //
            t2.setbounds(0, n-1);
            for(k = 0; k <= m-1; k++)
            {
                for(i = 0; i <= n-1; i++)
                {
                    t = ap::vdotproduct(&x(i, 0), 1, &v(0, k), v.getstride(), ap::vlen(0,m-1));
                    t2(i) = t;
                }
                calculatemv(t2, n, tmean, tmeans, tstddev, tstddevs);
                if( n!=1 )
                {
                    t = ap::sqr(tstddev)*n/(n-1);
                }
                else
                {
                    t = 0;
                }
                pcavarerrors = pcavarerrors||ap::fp_greater(fabs(t-s(k)),threshold);
            }
            for(k = 0; k <= m-2; k++)
            {
                pcavarerrors = pcavarerrors||ap::fp_less(s(k),s(k+1));
            }
            
            //
            // Optimality: different perturbations in V[..,0] can't
            // increase variance of projection - can only decrease.
            //
            t2.setbounds(0, n-1);
            t3.setbounds(0, n-1);
            for(i = 0; i <= n-1; i++)
            {
                t = ap::vdotproduct(&x(i, 0), 1, &v(0, 0), v.getstride(), ap::vlen(0,m-1));
                t2(i) = t;
            }
            calculatemv(t2, n, tmean, tmeans, tstddev, tstddevs);
            for(k = 0; k <= 2*m-1; k++)
            {
                h = 0.001;
                if( k%2!=0 )
                {
                    h = -h;
                }
                ap::vmove(&t3(0), 1, &t2(0), 1, ap::vlen(0,n-1));
                ap::vadd(&t3(0), 1, &x(0, k/2), x.getstride(), ap::vlen(0,n-1), h);
                t = 0;
                for(j = 0; j <= m-1; j++)
                {
                    if( j!=k/2 )
                    {
                        t = t+ap::sqr(v(j,0));
                    }
                    else
                    {
                        t = t+ap::sqr(v(j,0)+h);
                    }
                }
                t = 1/sqrt(t);
                ap::vmul(&t3(0), 1, ap::vlen(0,n-1), t);
                calculatemv(t3, n, tmean2, tmeans2, tstddev2, tstddevs2);
                pcaopterrors = pcaopterrors||ap::fp_greater(tstddev2,tstddev+threshold);
            }
        }
    }
    
    //
    // Special test for N=0
    //
    for(m = 1; m <= maxm; m++)
    {
        
        //
        // Solve
        //
        pcabuildbasis(x, 0, m, info, s, v);
        if( info!=1 )
        {
            pcaconverrors = true;
            continue;
        }
        
        //
        // Orthogonality test
        //
        for(i = 0; i <= m-1; i++)
        {
            for(j = 0; j <= m-1; j++)
            {
                t = ap::vdotproduct(&v(0, i), v.getstride(), &v(0, j), v.getstride(), ap::vlen(0,m-1));
                if( i==j )
                {
                    t = t-1;
                }
                pcaorterrors = pcaorterrors||ap::fp_greater(fabs(t),threshold);
            }
        }
    }
    
    //
    // Final report
    //
    waserrors = pcaconverrors||pcaorterrors||pcavarerrors||pcaopterrors;
    if( !silent )
    {
        printf("PCA TEST\n");
        printf("TOTAL RESULTS:                           ");
        if( !waserrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("* CONVERGENCE                            ");
        if( !pcaconverrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("* ORTOGONALITY                           ");
        if( !pcaorterrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("* VARIANCE REPORT                        ");
        if( !pcavarerrors )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("* OPTIMALITY                             ");
        if( !pcaopterrors )
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
Moments estimates and their errors
*************************************************************************/
static void calculatemv(const ap::real_1d_array& x,
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
        variance = (v1-v2)/n;
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


/*************************************************************************
Silent unit test
*************************************************************************/
bool testpcaunit_test_silent()
{
    bool result;

    result = testpca(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testpcaunit_test()
{
    bool result;

    result = testpca(false);
    return result;
}




