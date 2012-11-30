
 
#include <stdio.h>
#include "testinvldltunit.h"

static void generatematrix(ap::real_2d_array& a, int n, int task);
static void restorematrix(ap::real_2d_array& a, int n, bool upperin);

bool testinvldlt(bool silent)
{
    bool result;
    ap::real_2d_array a;
    ap::real_2d_array a2;
    ap::real_2d_array a3;
    ap::integer_1d_array p;
    int n;
    int pass;
    int mtask;
    int i;
    int j;
    int k;
    int minij;
    bool upperin;
    bool cr;
    double v;
    double err;
    bool waserrors;
    int passcount;
    int maxn;
    int htask;
    double threshold;

    err = 0;
    passcount = 10;
    maxn = 20;
    threshold = 10000000*ap::machineepsilon;
    waserrors = false;
    
    //
    // Test
    //
    for(n = 1; n <= maxn; n++)
    {
        a.setbounds(0, n-1, 0, n-1);
        a2.setbounds(0, n-1, 0, n-1);
        a3.setbounds(0, n-1, 0, n-1);
        for(mtask = 2; mtask <= 2; mtask++)
        {
            for(htask = 0; htask <= 1; htask++)
            {
                for(pass = 1; pass <= passcount; pass++)
                {
                    upperin = htask==0;
                    
                    //
                    // Prepare task:
                    // * A contains symmetric matrix
                    // * A2, A3 contains its upper (or lower) half
                    //
                    generatematrix(a, n, mtask);
                    for(i = 0; i <= n-1; i++)
                    {
                        for(j = 0; j <= n-1; j++)
                        {
                            a2(i,j) = a(i,j);
                            a3(i,j) = a(i,j);
                        }
                    }
                    for(i = 0; i <= n-1; i++)
                    {
                        for(j = 0; j <= n-1; j++)
                        {
                            if( upperin )
                            {
                                if( j<i )
                                {
                                    a2(i,j) = 0;
                                    a3(i,j) = 0;
                                }
                            }
                            else
                            {
                                if( i<j )
                                {
                                    a2(i,j) = 0;
                                    a3(i,j) = 0;
                                }
                            }
                        }
                    }
                    
                    //
                    // Test 1: inv(A2)
                    //
                    smatrixinverse(a2, n, upperin);
                    restorematrix(a2, n, upperin);
                    for(i = 0; i <= n-1; i++)
                    {
                        for(j = 0; j <= n-1; j++)
                        {
                            v = ap::vdotproduct(&a(i, 0), 1, &a2(0, j), a2.getstride(), ap::vlen(0,n-1));
                            if( i==j )
                            {
                                v = v-1;
                            }
                            err = ap::maxreal(err, fabs(v));
                        }
                    }
                    
                    //
                    // Test 2: inv(LDLt(A3))
                    //
                    smatrixldlt(a3, n, upperin, p);
                    smatrixldltinverse(a3, p, n, upperin);
                    restorematrix(a3, n, upperin);
                    for(i = 0; i <= n-1; i++)
                    {
                        for(j = 0; j <= n-1; j++)
                        {
                            v = ap::vdotproduct(&a(i, 0), 1, &a3(0, j), a3.getstride(), ap::vlen(0,n-1));
                            if( i==j )
                            {
                                v = v-1;
                            }
                            err = ap::maxreal(err, fabs(v));
                        }
                    }
                }
            }
        }
    }
    
    //
    // report
    //
    waserrors = ap::fp_greater(err,threshold);
    if( !silent )
    {
        printf("TESTING LDLT INVERSE\n");
        printf("ERROR:                                   %5.3le\n",
            double(err));
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


static void generatematrix(ap::real_2d_array& a, int n, int task)
{
    int i;
    int j;

    if( task==0 )
    {
        
        //
        // Zero matrix
        //
        for(i = 0; i <= n-1; i++)
        {
            for(j = 0; j <= n-1; j++)
            {
                a(i,j) = 0;
            }
        }
    }
    if( task==1 )
    {
        
        //
        // Sparse matrix
        //
        for(i = 0; i <= n-1; i++)
        {
            for(j = i+1; j <= n-1; j++)
            {
                if( ap::fp_greater(ap::randomreal(),0.95) )
                {
                    a(i,j) = 2*ap::randomreal()-1;
                }
                else
                {
                    a(i,j) = 0;
                }
                a(j,i) = a(i,j);
            }
            if( ap::fp_greater(ap::randomreal(),0.95) )
            {
                a(i,i) = (2*ap::randominteger(2)-1)*(0.8+ap::randomreal());
            }
            else
            {
                a(i,i) = 0;
            }
        }
    }
    if( task==2 )
    {
        
        //
        // Dense matrix
        //
        for(i = 0; i <= n-1; i++)
        {
            for(j = i+1; j <= n-1; j++)
            {
                a(i,j) = 2*ap::randomreal()-1;
                a(j,i) = a(i,j);
            }
            a(i,i) = (2*ap::randominteger(2)-1)*(0.7+ap::randomreal());
        }
    }
}


static void restorematrix(ap::real_2d_array& a, int n, bool upperin)
{
    int i;
    int j;

    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-1; j++)
        {
            if( upperin )
            {
                if( j<i )
                {
                    a(i,j) = a(j,i);
                }
            }
            else
            {
                if( i<j )
                {
                    a(i,j) = a(j,i);
                }
            }
        }
    }
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testinvldltunit_test_silent()
{
    bool result;

    result = testinvldlt(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testinvldltunit_test()
{
    bool result;

    result = testinvldlt(false);
    return result;
}




