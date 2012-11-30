
 
#include <stdio.h>
#include "testldltunit.h"

static void generatematrix(ap::real_2d_array& a, int n, int task);

bool testldlt(bool silent)
{
    bool result;
    ap::real_2d_array a;
    ap::real_2d_array a2;
    ap::real_2d_array l;
    ap::real_2d_array d;
    ap::real_2d_array u;
    ap::real_2d_array t;
    ap::real_2d_array t2;
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
    threshold = 1000*ap::machineepsilon;
    waserrors = false;
    
    //
    // Test
    //
    for(n = 1; n <= maxn; n++)
    {
        a.setbounds(0, n-1, 0, n-1);
        a2.setbounds(0, n-1, 0, n-1);
        l.setbounds(0, n-1, 0, n-1);
        u.setbounds(0, n-1, 0, n-1);
        d.setbounds(0, n-1, 0, n-1);
        t.setbounds(0, n-1, 0, n-1);
        t2.setbounds(0, n-1, 0, n-1);
        for(mtask = 0; mtask <= 2; mtask++)
        {
            for(htask = 0; htask <= 1; htask++)
            {
                for(pass = 1; pass <= passcount; pass++)
                {
                    upperin = htask==0;
                    
                    //
                    // Prepare task:
                    // * A contains symmetric matrix
                    // * A2 contains its upper (or lower) half
                    //
                    generatematrix(a, n, mtask);
                    for(i = 0; i <= n-1; i++)
                    {
                        for(j = 0; j <= n-1; j++)
                        {
                            a2(i,j) = a(i,j);
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
                                }
                            }
                            else
                            {
                                if( i<j )
                                {
                                    a2(i,j) = 0;
                                }
                            }
                        }
                    }
                    
                    //
                    // LDLt
                    //
                    smatrixldlt(a2, n, upperin, p);
                    
                    //
                    // Test (upper or lower)
                    //
                    if( upperin )
                    {
                        
                        //
                        // Unpack D
                        //
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = 0; j <= n-1; j++)
                            {
                                d(i,j) = 0;
                            }
                        }
                        k = 0;
                        while(k<n)
                        {
                            if( p(k)>=0 )
                            {
                                d(k,k) = a2(k,k);
                                k = k+1;
                            }
                            else
                            {
                                d(k,k) = a2(k,k);
                                d(k,k+1) = a2(k,k+1);
                                d(k+1,k) = a2(k,k+1);
                                d(k+1,k+1) = a2(k+1,k+1);
                                k = k+2;
                            }
                        }
                        
                        //
                        // Unpack U
                        //
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = 0; j <= n-1; j++)
                            {
                                u(i,j) = 0;
                            }
                            u(i,i) = 1;
                        }
                        k = 0;
                        while(k<n)
                        {
                            
                            //
                            // unpack Uk
                            //
                            for(i = 0; i <= n-1; i++)
                            {
                                for(j = 0; j <= n-1; j++)
                                {
                                    t(i,j) = 0;
                                }
                                t(i,i) = 1;
                            }
                            if( p(k)>=0 )
                            {
                                for(i = 0; i <= k-1; i++)
                                {
                                    t(i,k) = a2(i,k);
                                }
                            }
                            else
                            {
                                for(i = 0; i <= k-1; i++)
                                {
                                    t(i,k) = a2(i,k);
                                    t(i,k+1) = a2(i,k+1);
                                }
                            }
                            
                            //
                            // multiply U
                            //
                            for(i = 0; i <= n-1; i++)
                            {
                                for(j = 0; j <= n-1; j++)
                                {
                                    v = ap::vdotproduct(&t(i, 0), 1, &u(0, j), u.getstride(), ap::vlen(0,n-1));
                                    t2(i,j) = v;
                                }
                            }
                            for(i = 0; i <= n-1; i++)
                            {
                                for(j = 0; j <= n-1; j++)
                                {
                                    u(i,j) = t2(i,j);
                                }
                            }
                            
                            //
                            // permutations
                            //
                            if( p(k)>=0 )
                            {
                                for(j = 0; j <= n-1; j++)
                                {
                                    v = u(k,j);
                                    u(k,j) = u(p(k),j);
                                    u(p(k),j) = v;
                                }
                                k = k+1;
                            }
                            else
                            {
                                for(j = 0; j <= n-1; j++)
                                {
                                    v = u(k,j);
                                    u(k,j) = u(n+p(k),j);
                                    u(n+p(k),j) = v;
                                }
                                k = k+2;
                            }
                        }
                        
                        //
                        // Calculate U*D*U', store result in T2
                        //
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = 0; j <= n-1; j++)
                            {
                                v = ap::vdotproduct(&u(i, 0), 1, &d(0, j), d.getstride(), ap::vlen(0,n-1));
                                t(i,j) = v;
                            }
                        }
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = 0; j <= n-1; j++)
                            {
                                v = ap::vdotproduct(&t(i, 0), 1, &u(j, 0), 1, ap::vlen(0,n-1));
                                t2(i,j) = v;
                            }
                        }
                    }
                    else
                    {
                        
                        //
                        // Unpack D
                        //
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = 0; j <= n-1; j++)
                            {
                                d(i,j) = 0;
                            }
                        }
                        k = 0;
                        while(k<n)
                        {
                            if( p(k)>=0 )
                            {
                                d(k,k) = a2(k,k);
                                k = k+1;
                            }
                            else
                            {
                                d(k,k) = a2(k,k);
                                d(k,k+1) = a2(k+1,k);
                                d(k+1,k) = a2(k+1,k);
                                d(k+1,k+1) = a2(k+1,k+1);
                                k = k+2;
                            }
                        }
                        
                        //
                        // Unpack L
                        //
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = 0; j <= n-1; j++)
                            {
                                l(i,j) = 0;
                            }
                            l(i,i) = 1;
                        }
                        k = 0;
                        while(k<n)
                        {
                            
                            //
                            // permutations
                            //
                            if( p(k)>=0 )
                            {
                                for(i = 0; i <= n-1; i++)
                                {
                                    v = l(i,k);
                                    l(i,k) = l(i,p(k));
                                    l(i,p(k)) = v;
                                }
                            }
                            else
                            {
                                for(i = 0; i <= n-1; i++)
                                {
                                    v = l(i,k+1);
                                    l(i,k+1) = l(i,n+p(k+1));
                                    l(i,n+p(k+1)) = v;
                                }
                            }
                            
                            //
                            // unpack Lk
                            //
                            for(i = 0; i <= n-1; i++)
                            {
                                for(j = 0; j <= n-1; j++)
                                {
                                    t(i,j) = 0;
                                }
                                t(i,i) = 1;
                            }
                            if( p(k)>=0 )
                            {
                                for(i = k+1; i <= n-1; i++)
                                {
                                    t(i,k) = a2(i,k);
                                }
                            }
                            else
                            {
                                for(i = k+2; i <= n-1; i++)
                                {
                                    t(i,k) = a2(i,k);
                                    t(i,k+1) = a2(i,k+1);
                                }
                            }
                            
                            //
                            // multiply L
                            //
                            for(i = 0; i <= n-1; i++)
                            {
                                for(j = 0; j <= n-1; j++)
                                {
                                    v = ap::vdotproduct(&l(i, 0), 1, &t(0, j), t.getstride(), ap::vlen(0,n-1));
                                    t2(i,j) = v;
                                }
                            }
                            for(i = 0; i <= n-1; i++)
                            {
                                for(j = 0; j <= n-1; j++)
                                {
                                    l(i,j) = t2(i,j);
                                }
                            }
                            
                            //
                            // Next K
                            //
                            if( p(k)>=0 )
                            {
                                k = k+1;
                            }
                            else
                            {
                                k = k+2;
                            }
                        }
                        
                        //
                        // Calculate L*D*L', store result in T2
                        //
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = 0; j <= n-1; j++)
                            {
                                v = ap::vdotproduct(&l(i, 0), 1, &d(0, j), d.getstride(), ap::vlen(0,n-1));
                                t(i,j) = v;
                            }
                        }
                        for(i = 0; i <= n-1; i++)
                        {
                            for(j = 0; j <= n-1; j++)
                            {
                                v = ap::vdotproduct(&t(i, 0), 1, &l(j, 0), 1, ap::vlen(0,n-1));
                                t2(i,j) = v;
                            }
                        }
                    }
                    
                    //
                    // Test
                    //
                    for(i = 0; i <= n-1; i++)
                    {
                        for(j = 0; j <= n-1; j++)
                        {
                            err = ap::maxreal(err, fabs(a(i,j)-t2(i,j)));
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
        printf("TESTING LDLT DECOMPOSITION\n");
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
            a(i,i) = (2*ap::randominteger(2)-1)*(0.8+ap::randomreal());
        }
    }
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testldltunit_test_silent()
{
    bool result;

    result = testldlt(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testldltunit_test()
{
    bool result;

    result = testldlt(false);
    return result;
}




