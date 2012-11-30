
 
#include <stdio.h>
#include "testtsortunit.h"

static void unset2d(ap::complex_2d_array& a);
static void unset1d(ap::real_1d_array& a);
static void unset1di(ap::integer_1d_array& a);
static void testsortresults(const ap::real_1d_array& asorted,
     const ap::integer_1d_array& p1,
     const ap::integer_1d_array& p2,
     const ap::real_1d_array& aoriginal,
     int n,
     bool& waserrors);

/*************************************************************************
Testing tag sort
*************************************************************************/
bool testsort(bool silent)
{
    bool result;
    bool waserrors;
    int n;
    int i;
    int pass;
    int passcount;
    int maxn;
    ap::real_1d_array a;
    ap::real_1d_array a0;
    ap::real_1d_array a2;
    ap::integer_1d_array p1;
    ap::integer_1d_array p2;

    waserrors = false;
    maxn = 100;
    passcount = 10;
    
    //
    // Test
    //
    for(n = 1; n <= maxn; n++)
    {
        for(pass = 1; pass <= passcount; pass++)
        {
            
            //
            // (probably) distinct sort
            //
            unset1di(p1);
            unset1di(p2);
            a.setbounds(0, n-1);
            a0.setbounds(0, n-1);
            for(i = 0; i <= n-1; i++)
            {
                a(i) = 2*ap::randomreal()-1;
                a0(i) = a(i);
            }
            tagsort(a0, n, p1, p2);
            testsortresults(a0, p1, p2, a, n, waserrors);
            
            //
            // non-distinct sort
            //
            unset1di(p1);
            unset1di(p2);
            a.setbounds(0, n-1);
            a0.setbounds(0, n-1);
            for(i = 0; i <= n-1; i++)
            {
                a(i) = i/2;
                a0(i) = a(i);
            }
            tagsort(a0, n, p1, p2);
            testsortresults(a0, p1, p2, a, n, waserrors);
        }
    }
    
    //
    // report
    //
    if( !silent )
    {
        printf("TESTING TAGSORT\n");
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
Unsets 2D array.
*************************************************************************/
static void unset2d(ap::complex_2d_array& a)
{

    a.setbounds(0, 0, 0, 0);
    a(0,0) = 2*ap::randomreal()-1;
}


/*************************************************************************
Unsets 1D array.
*************************************************************************/
static void unset1d(ap::real_1d_array& a)
{

    a.setbounds(0, 0);
    a(0) = 2*ap::randomreal()-1;
}


/*************************************************************************
Unsets 1D array.
*************************************************************************/
static void unset1di(ap::integer_1d_array& a)
{

    a.setbounds(0, 0);
    a(0) = ap::randominteger(3)-1;
}


static void testsortresults(const ap::real_1d_array& asorted,
     const ap::integer_1d_array& p1,
     const ap::integer_1d_array& p2,
     const ap::real_1d_array& aoriginal,
     int n,
     bool& waserrors)
{
    int i;
    ap::real_1d_array a2;
    double t;
    ap::integer_1d_array f;

    a2.setbounds(0, n-1);
    f.setbounds(0, n-1);
    
    //
    // is set ordered?
    //
    for(i = 0; i <= n-2; i++)
    {
        waserrors = waserrors||ap::fp_greater(asorted(i),asorted(i+1));
    }
    
    //
    // P1 correctness
    //
    for(i = 0; i <= n-1; i++)
    {
        waserrors = waserrors||ap::fp_neq(asorted(i),aoriginal(p1(i)));
    }
    for(i = 0; i <= n-1; i++)
    {
        f(i) = 0;
    }
    for(i = 0; i <= n-1; i++)
    {
        f(p1(i)) = f(p1(i))+1;
    }
    for(i = 0; i <= n-1; i++)
    {
        waserrors = waserrors||f(i)!=1;
    }
    
    //
    // P2 correctness
    //
    for(i = 0; i <= n-1; i++)
    {
        a2(i) = aoriginal(i);
    }
    for(i = 0; i <= n-1; i++)
    {
        if( p2(i)!=i )
        {
            t = a2(i);
            a2(i) = a2(p2(i));
            a2(p2(i)) = t;
        }
    }
    for(i = 0; i <= n-1; i++)
    {
        waserrors = waserrors||ap::fp_neq(asorted(i),a2(i));
    }
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testtsortunit_test_silent()
{
    bool result;

    result = testsort(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testtsortunit_test()
{
    bool result;

    result = testsort(false);
    return result;
}




