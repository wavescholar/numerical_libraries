
 
#include <stdio.h>
#include "testspdgevdunit.h"

/*************************************************************************
Testing bidiagonal SVD decomposition subroutine
*************************************************************************/
bool testspdgevd(bool silent)
{
    bool result;
    int pass;
    int n;
    int passcount;
    int maxn;
    int atask;
    int btask;
    ap::real_1d_array d;
    ap::real_1d_array t1;
    ap::real_2d_array a;
    ap::real_2d_array b;
    ap::real_2d_array afull;
    ap::real_2d_array bfull;
    ap::real_2d_array l;
    ap::real_2d_array z;
    bool isuppera;
    bool isupperb;
    int i;
    int j;
    int minij;
    double v;
    double v1;
    double v2;
    bool cw;
    double err;
    double valerr;
    double threshold;
    bool waserrors;
    bool wfailed;
    bool wnsorted;

    threshold = 10000*ap::machineepsilon;
    valerr = 0;
    wfailed = false;
    wnsorted = false;
    maxn = 20;
    passcount = 5;
    
    //
    // Main cycle
    //
    for(n = 1; n <= maxn; n++)
    {
        for(pass = 1; pass <= passcount; pass++)
        {
            for(atask = 0; atask <= 1; atask++)
            {
                for(btask = 0; btask <= 1; btask++)
                {
                    isuppera = atask==0;
                    isupperb = btask==0;
                    
                    //
                    // Initialize A, B, AFull, BFull
                    //
                    t1.setbounds(0, n-1);
                    a.setbounds(0, n-1, 0, n-1);
                    b.setbounds(0, n-1, 0, n-1);
                    afull.setbounds(0, n-1, 0, n-1);
                    bfull.setbounds(0, n-1, 0, n-1);
                    l.setbounds(0, n-1, 0, n-1);
                    for(i = 0; i <= n-1; i++)
                    {
                        for(j = 0; j <= n-1; j++)
                        {
                            a(i,j) = 2*ap::randomreal()-1;
                            a(j,i) = a(i,j);
                            afull(i,j) = a(i,j);
                            afull(j,i) = a(i,j);
                        }
                    }
                    for(i = 0; i <= n-1; i++)
                    {
                        for(j = i+1; j <= n-1; j++)
                        {
                            l(i,j) = ap::randomreal();
                            l(j,i) = l(i,j);
                        }
                        l(i,i) = 1.5+ap::randomreal();
                    }
                    for(i = 0; i <= n-1; i++)
                    {
                        for(j = 0; j <= n-1; j++)
                        {
                            minij = ap::minint(i, j);
                            v = ap::vdotproduct(&l(i, 0), 1, &l(0, j), l.getstride(), ap::vlen(0,minij));
                            b(i,j) = v;
                            b(j,i) = v;
                            bfull(i,j) = v;
                            bfull(j,i) = v;
                        }
                    }
                    for(i = 0; i <= n-1; i++)
                    {
                        for(j = 0; j <= n-1; j++)
                        {
                            if( isuppera )
                            {
                                if( j<i )
                                {
                                    a(i,j) = 2*ap::randomreal()-1;
                                }
                            }
                            else
                            {
                                if( i<j )
                                {
                                    a(i,j) = 2*ap::randomreal()-1;
                                }
                            }
                            if( isupperb )
                            {
                                if( j<i )
                                {
                                    b(i,j) = 2*ap::randomreal()-1;
                                }
                            }
                            else
                            {
                                if( i<j )
                                {
                                    b(i,j) = 2*ap::randomreal()-1;
                                }
                            }
                        }
                    }
                    
                    //
                    // Problem 1
                    //
                    if( !smatrixgevd(a, n, isuppera, b, isupperb, 1, 1, d, z) )
                    {
                        wfailed = true;
                        continue;
                    }
                    err = 0;
                    for(j = 0; j <= n-1; j++)
                    {
                        for(i = 0; i <= n-1; i++)
                        {
                            v1 = ap::vdotproduct(&afull(i, 0), 1, &z(0, j), z.getstride(), ap::vlen(0,n-1));
                            v2 = ap::vdotproduct(&bfull(i, 0), 1, &z(0, j), z.getstride(), ap::vlen(0,n-1));
                            err = ap::maxreal(err, fabs(v1-d(j)*v2));
                        }
                    }
                    valerr = ap::maxreal(err, valerr);
                    
                    //
                    // Problem 2
                    //
                    if( !smatrixgevd(a, n, isuppera, b, isupperb, 1, 2, d, z) )
                    {
                        wfailed = true;
                        continue;
                    }
                    err = 0;
                    for(j = 0; j <= n-1; j++)
                    {
                        for(i = 0; i <= n-1; i++)
                        {
                            v1 = ap::vdotproduct(&bfull(i, 0), 1, &z(0, j), z.getstride(), ap::vlen(0,n-1));
                            t1(i) = v1;
                        }
                        for(i = 0; i <= n-1; i++)
                        {
                            v2 = ap::vdotproduct(&afull(i, 0), 1, &t1(0), 1, ap::vlen(0,n-1));
                            err = ap::maxreal(err, fabs(v2-d(j)*z(i,j)));
                        }
                    }
                    valerr = ap::maxreal(err, valerr);
                    
                    //
                    // Test problem 3
                    //
                    if( !smatrixgevd(a, n, isuppera, b, isupperb, 1, 3, d, z) )
                    {
                        wfailed = true;
                        continue;
                    }
                    err = 0;
                    for(j = 0; j <= n-1; j++)
                    {
                        for(i = 0; i <= n-1; i++)
                        {
                            v1 = ap::vdotproduct(&afull(i, 0), 1, &z(0, j), z.getstride(), ap::vlen(0,n-1));
                            t1(i) = v1;
                        }
                        for(i = 0; i <= n-1; i++)
                        {
                            v2 = ap::vdotproduct(&bfull(i, 0), 1, &t1(0), 1, ap::vlen(0,n-1));
                            err = ap::maxreal(err, fabs(v2-d(j)*z(i,j)));
                        }
                    }
                    valerr = ap::maxreal(err, valerr);
                }
            }
        }
    }
    
    //
    // report
    //
    waserrors = ap::fp_greater(valerr,threshold)||wfailed||wnsorted;
    if( !silent )
    {
        printf("TESTING SYMMETRIC GEVD\n");
        printf("Av-lambdav error (generalized):          %5.3le\n",
            double(valerr));
        printf("Eigen values order:                      ");
        if( !wnsorted )
        {
            printf("OK\n");
        }
        else
        {
            printf("FAILED\n");
        }
        printf("Always converged:                        ");
        if( !wfailed )
        {
            printf("YES\n");
        }
        else
        {
            printf("NO\n");
        }
        printf("Threshold:                               %5.3le\n",
            double(threshold));
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
Silent unit test
*************************************************************************/
bool testspdgevdunit_test_silent()
{
    bool result;

    result = testspdgevd(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testspdgevdunit_test()
{
    bool result;

    result = testspdgevd(false);
    return result;
}




