
 
#include <stdio.h>
#include "testcreflunit.h"

bool testcrefl(bool silent)
{
    bool result;
    int i;
    int j;
    int k;
    int n;
    int m;
    int maxmn;
    ap::complex_1d_array x;
    ap::complex_1d_array v;
    ap::complex_1d_array work;
    ap::complex_2d_array h;
    ap::complex_2d_array a;
    ap::complex_2d_array b;
    ap::complex_2d_array c;
    ap::complex tmp;
    ap::complex beta;
    ap::complex tau;
    double err;
    double mer;
    double mel;
    double meg;
    int pass;
    int passcount;
    bool waserrors;
    double threshold;

    threshold = 1000*ap::machineepsilon;
    passcount = 1000;
    mer = 0;
    mel = 0;
    meg = 0;
    for(pass = 1; pass <= passcount; pass++)
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
        // GenerateReflection
        //
        for(i = 1; i <= n; i++)
        {
            x(i).x = 2*ap::randomreal()-1;
            x(i).y = 2*ap::randomreal()-1;
            v(i) = x(i);
        }
        complexgeneratereflection(v, n, tau);
        beta = v(1);
        v(1) = 1;
        for(i = 1; i <= n; i++)
        {
            for(j = 1; j <= n; j++)
            {
                if( i==j )
                {
                    h(i,j) = 1-tau*v(i)*ap::conj(v(j));
                }
                else
                {
                    h(i,j) = -tau*v(i)*ap::conj(v(j));
                }
            }
        }
        err = 0;
        for(i = 1; i <= n; i++)
        {
            tmp = ap::vdotproduct(&h(1, i), h.getstride(), "Conj", &x(1), 1, "N", ap::vlen(1,n));
            if( i==1 )
            {
                err = ap::maxreal(err, ap::abscomplex(tmp-beta));
            }
            else
            {
                err = ap::maxreal(err, ap::abscomplex(tmp));
            }
        }
        err = ap::maxreal(err, fabs(beta.y));
        meg = ap::maxreal(meg, err);
        
        //
        // ApplyReflectionFromTheLeft
        //
        for(i = 1; i <= m; i++)
        {
            x(i).x = 2*ap::randomreal()-1;
            x(i).y = 2*ap::randomreal()-1;
            v(i) = x(i);
        }
        for(i = 1; i <= m; i++)
        {
            for(j = 1; j <= n; j++)
            {
                a(i,j).x = 2*ap::randomreal()-1;
                a(i,j).y = 2*ap::randomreal()-1;
                b(i,j) = a(i,j);
            }
        }
        complexgeneratereflection(v, m, tau);
        beta = v(1);
        v(1) = 1;
        complexapplyreflectionfromtheleft(b, tau, v, 1, m, 1, n, work);
        for(i = 1; i <= m; i++)
        {
            for(j = 1; j <= m; j++)
            {
                if( i==j )
                {
                    h(i,j) = 1-tau*v(i)*ap::conj(v(j));
                }
                else
                {
                    h(i,j) = -tau*v(i)*ap::conj(v(j));
                }
            }
        }
        for(i = 1; i <= m; i++)
        {
            for(j = 1; j <= n; j++)
            {
                tmp = ap::vdotproduct(&h(i, 1), 1, "N", &a(1, j), a.getstride(), "N", ap::vlen(1,m));
                c(i,j) = tmp;
            }
        }
        err = 0;
        for(i = 1; i <= m; i++)
        {
            for(j = 1; j <= n; j++)
            {
                err = ap::maxreal(err, ap::abscomplex(b(i,j)-c(i,j)));
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
        complexgeneratereflection(v, n, tau);
        beta = v(1);
        v(1) = 1;
        complexapplyreflectionfromtheright(b, tau, v, 1, m, 1, n, work);
        for(i = 1; i <= n; i++)
        {
            for(j = 1; j <= n; j++)
            {
                if( i==j )
                {
                    h(i,j) = 1-tau*v(i)*ap::conj(v(j));
                }
                else
                {
                    h(i,j) = -tau*v(i)*ap::conj(v(j));
                }
            }
        }
        for(i = 1; i <= m; i++)
        {
            for(j = 1; j <= n; j++)
            {
                tmp = ap::vdotproduct(&a(i, 1), 1, "N", &h(1, j), h.getstride(), "N", ap::vlen(1,n));
                c(i,j) = tmp;
            }
        }
        err = 0;
        for(i = 1; i <= m; i++)
        {
            for(j = 1; j <= n; j++)
            {
                err = ap::maxreal(err, ap::abscomplex(b(i,j)-c(i,j)));
            }
        }
        mer = ap::maxreal(mer, err);
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
    complexgeneratereflection(v, 10, tau);
    
    //
    // report
    //
    waserrors = ap::fp_greater(meg,threshold)||ap::fp_greater(mel,threshold)||ap::fp_greater(mer,threshold);
    if( !silent )
    {
        printf("TESTING COMPLEX REFLECTIONS\n");
        printf("Generate error:                          %5.3le\n",
            double(meg));
        printf("Apply(L) error:                          %5.3le\n",
            double(mel));
        printf("Apply(R) error:                          %5.3le\n",
            double(mer));
        printf("Threshold:                               %5.3le\n",
            double(threshold));
        printf("Overflow crash test:                     PASSED\n");
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
bool testcreflunit_test_silent()
{
    bool result;

    result = testcrefl(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testcreflunit_test()
{
    bool result;

    result = testcrefl(false);
    return result;
}




