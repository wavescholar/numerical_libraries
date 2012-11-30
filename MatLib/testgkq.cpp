
 
#include <stdio.h>
#include "testgkq.h"

static double mapkind(int k);

/*************************************************************************
Test
*************************************************************************/
bool testgkqunit(bool silent)
{
    bool result;
    int pkind;
    double errtol;
    double eps;
    double nonstricterrtol;
    int n;
    int i;
    int k;
    int info;
    double err;
    int akind;
    int bkind;
    double alphac;
    double betac;
    ap::real_1d_array x1;
    ap::real_1d_array wg1;
    ap::real_1d_array wk1;
    ap::real_1d_array x2;
    ap::real_1d_array wg2;
    ap::real_1d_array wk2;
    int info1;
    int info2;
    bool successatleastonce;
    bool intblerrors;
    bool vstblerrors;
    bool generrors;
    bool waserrors;

    intblerrors = false;
    vstblerrors = false;
    generrors = false;
    waserrors = false;
    errtol = 10000*ap::machineepsilon;
    nonstricterrtol = 1000*errtol;
    
    //
    // test recurrence-based Legendre nodes against the precalculated table
    //
    for(pkind = 0; pkind <= 5; pkind++)
    {
        n = 0;
        if( pkind==0 )
        {
            n = 15;
        }
        if( pkind==1 )
        {
            n = 21;
        }
        if( pkind==2 )
        {
            n = 31;
        }
        if( pkind==3 )
        {
            n = 41;
        }
        if( pkind==4 )
        {
            n = 51;
        }
        if( pkind==5 )
        {
            n = 61;
        }
        gkqlegendrecalc(n, info, x1, wk1, wg1);
        gkqlegendretbl(n, x2, wk2, wg2, eps);
        if( info<=0 )
        {
            generrors = true;
            break;
        }
        for(i = 0; i <= n-1; i++)
        {
            vstblerrors = vstblerrors||ap::fp_greater(fabs(x1(i)-x2(i)),errtol);
            vstblerrors = vstblerrors||ap::fp_greater(fabs(wk1(i)-wk2(i)),errtol);
            vstblerrors = vstblerrors||ap::fp_greater(fabs(wg1(i)-wg2(i)),errtol);
        }
    }
    
    //
    // Test recurrence-baced Gauss-Kronrod nodes against Gauss-only nodes
    // calculated with subroutines from GQ unit.
    //
    for(k = 1; k <= 30; k++)
    {
        n = 2*k+1;
        
        //
        // Gauss-Legendre
        //
        err = 0;
        gkqgenerategausslegendre(n, info1, x1, wk1, wg1);
        gqgenerategausslegendre(k, info2, x2, wg2);
        if( info1>0&&info2>0 )
        {
            for(i = 0; i <= k-1; i++)
            {
                err = ap::maxreal(err, fabs(x1(2*i+1)-x2(i)));
                err = ap::maxreal(err, fabs(wg1(2*i+1)-wg2(i)));
            }
        }
        else
        {
            generrors = true;
        }
        generrors = generrors||ap::fp_greater(err,errtol);
    }
    for(k = 1; k <= 15; k++)
    {
        n = 2*k+1;
        
        //
        // Gauss-Jacobi
        //
        successatleastonce = false;
        err = 0;
        for(akind = 0; akind <= 9; akind++)
        {
            for(bkind = 0; bkind <= 9; bkind++)
            {
                alphac = mapkind(akind);
                betac = mapkind(bkind);
                gkqgenerategaussjacobi(n, alphac, betac, info1, x1, wk1, wg1);
                gqgenerategaussjacobi(k, alphac, betac, info2, x2, wg2);
                if( info1>0&&info2>0 )
                {
                    successatleastonce = true;
                    for(i = 0; i <= k-1; i++)
                    {
                        err = ap::maxreal(err, fabs(x1(2*i+1)-x2(i)));
                        err = ap::maxreal(err, fabs(wg1(2*i+1)-wg2(i)));
                    }
                }
                else
                {
                    generrors = generrors||info1!=-5;
                }
            }
        }
        generrors = generrors||ap::fp_greater(err,errtol)||!successatleastonce;
    }
    
    //
    // end
    //
    waserrors = intblerrors||vstblerrors||generrors;
    if( !silent )
    {
        printf("TESTING GAUSS-KRONROD QUADRATURES\n");
        printf("FINAL RESULT:                             ");
        if( waserrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* PRE-CALCULATED TABLE:                   ");
        if( intblerrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* CALCULATED AGAINST THE TABLE:           ");
        if( vstblerrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* GENERAL PROPERTIES:                     ");
        if( generrors )
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
    result = !waserrors;
    return result;
}


/*************************************************************************
Maps:
    0   =>  -0.9
    1   =>  -0.5
    2   =>  -0.1
    3   =>   0.0
    4   =>  +0.1
    5   =>  +0.5
    6   =>  +0.9
    7   =>  +1.0
    8   =>  +1.5
    9   =>  +2.0
*************************************************************************/
static double mapkind(int k)
{
    double result;

    result = 0;
    if( k==0 )
    {
        result = -0.9;
    }
    if( k==1 )
    {
        result = -0.5;
    }
    if( k==2 )
    {
        result = -0.1;
    }
    if( k==3 )
    {
        result = 0.0;
    }
    if( k==4 )
    {
        result = +0.1;
    }
    if( k==5 )
    {
        result = +0.5;
    }
    if( k==6 )
    {
        result = +0.9;
    }
    if( k==7 )
    {
        result = +1.0;
    }
    if( k==8 )
    {
        result = +1.5;
    }
    if( k==9 )
    {
        result = +2.0;
    }
    return result;
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testgkq_test_silent()
{
    bool result;

    result = testgkqunit(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testgkq_test()
{
    bool result;

    result = testgkqunit(false);
    return result;
}




