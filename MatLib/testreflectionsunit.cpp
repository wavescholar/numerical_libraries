
 
#include <stdio.h>
#include "testreflectionsunit.h"

bool testreflections(bool silent)
{
    bool result;
    int i;
    int j;
    int n;
    int m;
    int maxmn;
    ap::real_1d_array x;
    ap::real_1d_array v;
    ap::real_1d_array work;
    ap::real_2d_array h;
    ap::real_2d_array a;
    ap::real_2d_array b;
    ap::real_2d_array c;
    double tmp;
    double beta;
    double tau;
    double err;
    double mer;
    double mel;
    double meg;
    int pass;
    int passcount;
    double threshold;
    int tasktype;
    double xscale;

    passcount = 10;
    threshold = 100*ap::machineepsilon;
    mer = 0;
    mel = 0;
    meg = 0;
    for(pass = 1; pass <= passcount; pass++)
    {
        for(n = 1; n <= 10; n++)
        {
            for(m = 1; m <= 10; m++)
            {
                
                //
                // Task
                //
                n = 1+ap::randominteger(10);
                m = 1+ap::randominteger(10);
                maxmn = ap::maxint(m, n);
                
                //
                // Initialize
                //
                x.setbounds(1, maxmn);
                v.setbounds(1, maxmn);
                work.setbounds(1, maxmn);
                h.setbounds(1, maxmn, 1, maxmn);
                a.setbounds(1, maxmn, 1, maxmn);
                b.setbounds(1, maxmn, 1, maxmn);
                c.setbounds(1, maxmn, 1, maxmn);
                
                //
                // GenerateReflection, three tasks are possible:
                // * random X
                // * zero X
                // * non-zero X[1], all other are zeros
                // * random X, near underflow scale
                // * random X, near overflow scale
                //
                for(tasktype = 0; tasktype <= 4; tasktype++)
                {
                    xscale = 1;
                    if( tasktype==0 )
                    {
                        for(i = 1; i <= n; i++)
                        {
                            x(i) = 2*ap::randomreal()-1;
                        }
                    }
                    if( tasktype==1 )
                    {
                        for(i = 1; i <= n; i++)
                        {
                            x(i) = 0;
                        }
                    }
                    if( tasktype==2 )
                    {
                        x(1) = 2*ap::randomreal()-1;
                        for(i = 2; i <= n; i++)
                        {
                            x(i) = 0;
                        }
                    }
                    if( tasktype==3 )
                    {
                        for(i = 1; i <= n; i++)
                        {
                            x(i) = (ap::randominteger(21)-10)*ap::minrealnumber;
                        }
                        xscale = 10*ap::minrealnumber;
                    }
                    if( tasktype==4 )
                    {
                        for(i = 1; i <= n; i++)
                        {
                            x(i) = (2*ap::randomreal()-1)*ap::maxrealnumber;
                        }
                        xscale = ap::maxrealnumber;
                    }
                    ap::vmove(&v(1), 1, &x(1), 1, ap::vlen(1,n));
                    generatereflection(v, n, tau);
                    beta = v(1);
                    v(1) = 1;
                    for(i = 1; i <= n; i++)
                    {
                        for(j = 1; j <= n; j++)
                        {
                            if( i==j )
                            {
                                h(i,j) = 1-tau*v(i)*v(j);
                            }
                            else
                            {
                                h(i,j) = -tau*v(i)*v(j);
                            }
                        }
                    }
                    err = 0;
                    for(i = 1; i <= n; i++)
                    {
                        tmp = ap::vdotproduct(&h(i, 1), 1, &x(1), 1, ap::vlen(1,n));
                        if( i==1 )
                        {
                            err = ap::maxreal(err, fabs(tmp-beta));
                        }
                        else
                        {
                            err = ap::maxreal(err, fabs(tmp));
                        }
                    }
                    meg = ap::maxreal(meg, err/xscale);
                }
                
                //
                // ApplyReflectionFromTheLeft
                //
                for(i = 1; i <= m; i++)
                {
                    x(i) = 2*ap::randomreal()-1;
                    v(i) = x(i);
                }
                for(i = 1; i <= m; i++)
                {
                    for(j = 1; j <= n; j++)
                    {
                        a(i,j) = 2*ap::randomreal()-1;
                        b(i,j) = a(i,j);
                    }
                }
                generatereflection(v, m, tau);
                beta = v(1);
                v(1) = 1;
                applyreflectionfromtheleft(b, tau, v, 1, m, 1, n, work);
                for(i = 1; i <= m; i++)
                {
                    for(j = 1; j <= m; j++)
                    {
                        if( i==j )
                        {
                            h(i,j) = 1-tau*v(i)*v(j);
                        }
                        else
                        {
                            h(i,j) = -tau*v(i)*v(j);
                        }
                    }
                }
                for(i = 1; i <= m; i++)
                {
                    for(j = 1; j <= n; j++)
                    {
                        tmp = ap::vdotproduct(&h(i, 1), 1, &a(1, j), a.getstride(), ap::vlen(1,m));
                        c(i,j) = tmp;
                    }
                }
                err = 0;
                for(i = 1; i <= m; i++)
                {
                    for(j = 1; j <= n; j++)
                    {
                        err = ap::maxreal(err, fabs(b(i,j)-c(i,j)));
                    }
                }
                mel = ap::maxreal(mel, err);
                
                //
                // ApplyReflectionFromTheRight
                //
                for(i = 1; i <= n; i++)
                {
                    x(i) = 2*ap::randomreal()-1;
                    v(i) = x(i);
                }
                for(i = 1; i <= m; i++)
                {
                    for(j = 1; j <= n; j++)
                    {
                        a(i,j) = 2*ap::randomreal()-1;
                        b(i,j) = a(i,j);
                    }
                }
                generatereflection(v, n, tau);
                beta = v(1);
                v(1) = 1;
                applyreflectionfromtheright(b, tau, v, 1, m, 1, n, work);
                for(i = 1; i <= n; i++)
                {
                    for(j = 1; j <= n; j++)
                    {
                        if( i==j )
                        {
                            h(i,j) = 1-tau*v(i)*v(j);
                        }
                        else
                        {
                            h(i,j) = -tau*v(i)*v(j);
                        }
                    }
                }
                for(i = 1; i <= m; i++)
                {
                    for(j = 1; j <= n; j++)
                    {
                        tmp = ap::vdotproduct(&a(i, 1), 1, &h(1, j), h.getstride(), ap::vlen(1,n));
                        c(i,j) = tmp;
                    }
                }
                err = 0;
                for(i = 1; i <= m; i++)
                {
                    for(j = 1; j <= n; j++)
                    {
                        err = ap::maxreal(err, fabs(b(i,j)-c(i,j)));
                    }
                }
                mer = ap::maxreal(mer, err);
            }
        }
    }
    
    //
    // Overflow crash test
    //
    x.setbounds(1, 10);
    v.setbounds(1, 10);
    for(i = 1; i <= 10; i++)
    {
        v(i) = ap::maxrealnumber*0.01*(2*ap::randomreal()-1);
    }
    generatereflection(v, 10, tau);
    result = ap::fp_less_eq(meg,threshold)&&ap::fp_less_eq(mel,threshold)&&ap::fp_less_eq(mer,threshold);
    if( !silent )
    {
        printf("TESTING REFLECTIONS\n");
        printf("Pass count is %0ld\n",
            long(passcount));
        printf("Generate     absolute error is       %5.3le\n",
            double(meg));
        printf("Apply(Left)  absolute error is       %5.3le\n",
            double(mel));
        printf("Apply(Right) absolute error is       %5.3le\n",
            double(mer));
        printf("Overflow crash test passed\n");
        if( result )
        {
            printf("TEST PASSED\n");
        }
        else
        {
            printf("TEST FAILED\n");
        }
    }
    return result;
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testreflectionsunit_test_silent()
{
    bool result;

    result = testreflections(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testreflectionsunit_test()
{
    bool result;

    result = testreflections(false);
    return result;
}




