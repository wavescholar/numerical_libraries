
 
#include <stdio.h>
#include "testgq.h"

static double mapkind(int k);
static void buildgausslegendrequadrature(int n,
     ap::real_1d_array& x,
     ap::real_1d_array& w);
static void buildgaussjacobiquadrature(int n,
     double alpha,
     double beta,
     ap::real_1d_array& x,
     ap::real_1d_array& w);
static void buildgausslaguerrequadrature(int n,
     double alpha,
     ap::real_1d_array& x,
     ap::real_1d_array& w);

/*************************************************************************
Test
*************************************************************************/
bool testgqunit(bool silent)
{
    bool result;
    ap::real_1d_array alpha;
    ap::real_1d_array beta;
    ap::real_1d_array x;
    ap::real_1d_array w;
    ap::real_1d_array x2;
    ap::real_1d_array w2;
    double err;
    int n;
    int i;
    int info;
    int akind;
    int bkind;
    double alphac;
    double betac;
    double errtol;
    double nonstricterrtol;
    double stricterrtol;
    bool recerrors;
    bool specerrors;
    bool waserrors;

    recerrors = false;
    specerrors = false;
    waserrors = false;
    errtol = 1.0E-12;
    nonstricterrtol = 1.0E-6;
    stricterrtol = 1000*ap::machineepsilon;
    
    //
    // Three tests for rec-based Gauss quadratures with known weights/nodes:
    // 1. Gauss-Legendre with N=2
    // 2. Gauss-Legendre with N=5
    // 3. Gauss-Chebyshev with N=1, 2, 4, 8, ..., 512
    //
    err = 0;
    alpha.setlength(2);
    beta.setlength(2);
    alpha(0) = 0;
    alpha(1) = 0;
    beta(1) = double(1)/double(4*1*1-1);
    gqgeneraterec(alpha, beta, 2.0, 2, info, x, w);
    if( info>0 )
    {
        err = ap::maxreal(err, fabs(x(0)+sqrt(double(3))/3));
        err = ap::maxreal(err, fabs(x(1)-sqrt(double(3))/3));
        err = ap::maxreal(err, fabs(w(0)-1));
        err = ap::maxreal(err, fabs(w(1)-1));
        for(i = 0; i <= 0; i++)
        {
            recerrors = recerrors||ap::fp_greater_eq(x(i),x(i+1));
        }
    }
    else
    {
        recerrors = true;
    }
    alpha.setlength(5);
    beta.setlength(5);
    alpha(0) = 0;
    for(i = 1; i <= 4; i++)
    {
        alpha(i) = 0;
        beta(i) = ap::sqr(double(i))/(4*ap::sqr(double(i))-1);
    }
    gqgeneraterec(alpha, beta, 2.0, 5, info, x, w);
    if( info>0 )
    {
        err = ap::maxreal(err, fabs(x(0)+sqrt(245+14*sqrt(double(70)))/21));
        err = ap::maxreal(err, fabs(x(0)+x(4)));
        err = ap::maxreal(err, fabs(x(1)+sqrt(245-14*sqrt(double(70)))/21));
        err = ap::maxreal(err, fabs(x(1)+x(3)));
        err = ap::maxreal(err, fabs(x(2)));
        err = ap::maxreal(err, fabs(w(0)-(322-13*sqrt(double(70)))/900));
        err = ap::maxreal(err, fabs(w(0)-w(4)));
        err = ap::maxreal(err, fabs(w(1)-(322+13*sqrt(double(70)))/900));
        err = ap::maxreal(err, fabs(w(1)-w(3)));
        err = ap::maxreal(err, fabs(w(2)-double(128)/double(225)));
        for(i = 0; i <= 3; i++)
        {
            recerrors = recerrors||ap::fp_greater_eq(x(i),x(i+1));
        }
    }
    else
    {
        recerrors = true;
    }
    n = 1;
    while(n<=512)
    {
        alpha.setlength(n);
        beta.setlength(n);
        for(i = 0; i <= n-1; i++)
        {
            alpha(i) = 0;
            if( i==0 )
            {
                beta(i) = 0;
            }
            if( i==1 )
            {
                beta(i) = double(1)/double(2);
            }
            if( i>1 )
            {
                beta(i) = double(1)/double(4);
            }
        }
        gqgeneraterec(alpha, beta, ap::pi(), n, info, x, w);
        if( info>0 )
        {
            for(i = 0; i <= n-1; i++)
            {
                err = ap::maxreal(err, fabs(x(i)-cos(ap::pi()*(n-i-0.5)/n)));
                err = ap::maxreal(err, fabs(w(i)-ap::pi()/n));
            }
            for(i = 0; i <= n-2; i++)
            {
                recerrors = recerrors||ap::fp_greater_eq(x(i),x(i+1));
            }
        }
        else
        {
            recerrors = true;
        }
        n = n*2;
    }
    recerrors = recerrors||ap::fp_greater(err,errtol);
    
    //
    // Three tests for rec-based Gauss-Lobatto quadratures with known weights/nodes:
    // 1. Gauss-Lobatto with N=3
    // 2. Gauss-Lobatto with N=4
    // 3. Gauss-Lobatto with N=6
    //
    err = 0;
    alpha.setlength(2);
    beta.setlength(2);
    alpha(0) = 0;
    alpha(1) = 0;
    beta(0) = 0;
    beta(1) = double(1*1)/double(4*1*1-1);
    gqgenerategausslobattorec(alpha, beta, 2.0, double(-1), double(+1), 3, info, x, w);
    if( info>0 )
    {
        err = ap::maxreal(err, fabs(x(0)+1));
        err = ap::maxreal(err, fabs(x(1)));
        err = ap::maxreal(err, fabs(x(2)-1));
        err = ap::maxreal(err, fabs(w(0)-double(1)/double(3)));
        err = ap::maxreal(err, fabs(w(1)-double(4)/double(3)));
        err = ap::maxreal(err, fabs(w(2)-double(1)/double(3)));
        for(i = 0; i <= 1; i++)
        {
            recerrors = recerrors||ap::fp_greater_eq(x(i),x(i+1));
        }
    }
    else
    {
        recerrors = true;
    }
    alpha.setlength(3);
    beta.setlength(3);
    alpha(0) = 0;
    alpha(1) = 0;
    alpha(2) = 0;
    beta(0) = 0;
    beta(1) = double(1*1)/double(4*1*1-1);
    beta(2) = double(2*2)/double(4*2*2-1);
    gqgenerategausslobattorec(alpha, beta, 2.0, double(-1), double(+1), 4, info, x, w);
    if( info>0 )
    {
        err = ap::maxreal(err, fabs(x(0)+1));
        err = ap::maxreal(err, fabs(x(1)+sqrt(double(5))/5));
        err = ap::maxreal(err, fabs(x(2)-sqrt(double(5))/5));
        err = ap::maxreal(err, fabs(x(3)-1));
        err = ap::maxreal(err, fabs(w(0)-double(1)/double(6)));
        err = ap::maxreal(err, fabs(w(1)-double(5)/double(6)));
        err = ap::maxreal(err, fabs(w(2)-double(5)/double(6)));
        err = ap::maxreal(err, fabs(w(3)-double(1)/double(6)));
        for(i = 0; i <= 2; i++)
        {
            recerrors = recerrors||ap::fp_greater_eq(x(i),x(i+1));
        }
    }
    else
    {
        recerrors = true;
    }
    alpha.setlength(5);
    beta.setlength(5);
    alpha(0) = 0;
    alpha(1) = 0;
    alpha(2) = 0;
    alpha(3) = 0;
    alpha(4) = 0;
    beta(0) = 0;
    beta(1) = double(1*1)/double(4*1*1-1);
    beta(2) = double(2*2)/double(4*2*2-1);
    beta(3) = double(3*3)/double(4*3*3-1);
    beta(4) = double(4*4)/double(4*4*4-1);
    gqgenerategausslobattorec(alpha, beta, 2.0, double(-1), double(+1), 6, info, x, w);
    if( info>0 )
    {
        err = ap::maxreal(err, fabs(x(0)+1));
        err = ap::maxreal(err, fabs(x(1)+sqrt((7+2*sqrt(double(7)))/21)));
        err = ap::maxreal(err, fabs(x(2)+sqrt((7-2*sqrt(double(7)))/21)));
        err = ap::maxreal(err, fabs(x(3)-sqrt((7-2*sqrt(double(7)))/21)));
        err = ap::maxreal(err, fabs(x(4)-sqrt((7+2*sqrt(double(7)))/21)));
        err = ap::maxreal(err, fabs(x(5)-1));
        err = ap::maxreal(err, fabs(w(0)-double(1)/double(15)));
        err = ap::maxreal(err, fabs(w(1)-(14-sqrt(double(7)))/30));
        err = ap::maxreal(err, fabs(w(2)-(14+sqrt(double(7)))/30));
        err = ap::maxreal(err, fabs(w(3)-(14+sqrt(double(7)))/30));
        err = ap::maxreal(err, fabs(w(4)-(14-sqrt(double(7)))/30));
        err = ap::maxreal(err, fabs(w(5)-double(1)/double(15)));
        for(i = 0; i <= 4; i++)
        {
            recerrors = recerrors||ap::fp_greater_eq(x(i),x(i+1));
        }
    }
    else
    {
        recerrors = true;
    }
    recerrors = recerrors||ap::fp_greater(err,errtol);
    
    //
    // Three tests for rec-based Gauss-Radau quadratures with known weights/nodes:
    // 1. Gauss-Radau with N=2
    // 2. Gauss-Radau with N=3
    // 3. Gauss-Radau with N=3 (another case)
    //
    err = 0;
    alpha.setlength(1);
    beta.setlength(2);
    alpha(0) = 0;
    beta(0) = 0;
    beta(1) = double(1*1)/double(4*1*1-1);
    gqgenerategaussradaurec(alpha, beta, 2.0, double(-1), 2, info, x, w);
    if( info>0 )
    {
        err = ap::maxreal(err, fabs(x(0)+1));
        err = ap::maxreal(err, fabs(x(1)-double(1)/double(3)));
        err = ap::maxreal(err, fabs(w(0)-0.5));
        err = ap::maxreal(err, fabs(w(1)-1.5));
        for(i = 0; i <= 0; i++)
        {
            recerrors = recerrors||ap::fp_greater_eq(x(i),x(i+1));
        }
    }
    else
    {
        recerrors = true;
    }
    alpha.setlength(2);
    beta.setlength(3);
    alpha(0) = 0;
    alpha(1) = 0;
    for(i = 0; i <= 2; i++)
    {
        beta(i) = ap::sqr(double(i))/(4*ap::sqr(double(i))-1);
    }
    gqgenerategaussradaurec(alpha, beta, 2.0, double(-1), 3, info, x, w);
    if( info>0 )
    {
        err = ap::maxreal(err, fabs(x(0)+1));
        err = ap::maxreal(err, fabs(x(1)-(1-sqrt(double(6)))/5));
        err = ap::maxreal(err, fabs(x(2)-(1+sqrt(double(6)))/5));
        err = ap::maxreal(err, fabs(w(0)-double(2)/double(9)));
        err = ap::maxreal(err, fabs(w(1)-(16+sqrt(double(6)))/18));
        err = ap::maxreal(err, fabs(w(2)-(16-sqrt(double(6)))/18));
        for(i = 0; i <= 1; i++)
        {
            recerrors = recerrors||ap::fp_greater_eq(x(i),x(i+1));
        }
    }
    else
    {
        recerrors = true;
    }
    alpha.setlength(2);
    beta.setlength(3);
    alpha(0) = 0;
    alpha(1) = 0;
    for(i = 0; i <= 2; i++)
    {
        beta(i) = ap::sqr(double(i))/(4*ap::sqr(double(i))-1);
    }
    gqgenerategaussradaurec(alpha, beta, 2.0, double(+1), 3, info, x, w);
    if( info>0 )
    {
        err = ap::maxreal(err, fabs(x(2)-1));
        err = ap::maxreal(err, fabs(x(1)+(1-sqrt(double(6)))/5));
        err = ap::maxreal(err, fabs(x(0)+(1+sqrt(double(6)))/5));
        err = ap::maxreal(err, fabs(w(2)-double(2)/double(9)));
        err = ap::maxreal(err, fabs(w(1)-(16+sqrt(double(6)))/18));
        err = ap::maxreal(err, fabs(w(0)-(16-sqrt(double(6)))/18));
        for(i = 0; i <= 1; i++)
        {
            recerrors = recerrors||ap::fp_greater_eq(x(i),x(i+1));
        }
    }
    else
    {
        recerrors = true;
    }
    recerrors = recerrors||ap::fp_greater(err,errtol);
    
    //
    // test recurrence-based special cases (Legendre, Jacobi, Hermite, ...)
    // against another implementation (polynomial root-finder)
    //
    for(n = 1; n <= 20; n++)
    {
        
        //
        // test gauss-legendre
        //
        err = 0;
        gqgenerategausslegendre(n, info, x, w);
        if( info>0 )
        {
            buildgausslegendrequadrature(n, x2, w2);
            for(i = 0; i <= n-1; i++)
            {
                err = ap::maxreal(err, fabs(x(i)-x2(i)));
                err = ap::maxreal(err, fabs(w(i)-w2(i)));
            }
        }
        else
        {
            specerrors = true;
        }
        specerrors = specerrors||ap::fp_greater(err,errtol);
        
        //
        // Test Gauss-Jacobi.
        // Since task is much more difficult we will use less strict
        // threshold.
        //
        err = 0;
        for(akind = 0; akind <= 9; akind++)
        {
            for(bkind = 0; bkind <= 9; bkind++)
            {
                alphac = mapkind(akind);
                betac = mapkind(bkind);
                gqgenerategaussjacobi(n, alphac, betac, info, x, w);
                if( info>0 )
                {
                    buildgaussjacobiquadrature(n, alphac, betac, x2, w2);
                    for(i = 0; i <= n-1; i++)
                    {
                        err = ap::maxreal(err, fabs(x(i)-x2(i)));
                        err = ap::maxreal(err, fabs(w(i)-w2(i)));
                    }
                }
                else
                {
                    specerrors = true;
                }
            }
        }
        specerrors = specerrors||ap::fp_greater(err,nonstricterrtol);
        
        //
        // special test for Gauss-Jacobi (Chebyshev weight
        // function with analytically known nodes/weights)
        //
        err = 0;
        gqgenerategaussjacobi(n, -0.5, -0.5, info, x, w);
        if( info>0 )
        {
            for(i = 0; i <= n-1; i++)
            {
                err = ap::maxreal(err, fabs(x(i)+cos(ap::pi()*(i+0.5)/n)));
                err = ap::maxreal(err, fabs(w(i)-ap::pi()/n));
            }
        }
        else
        {
            specerrors = true;
        }
        specerrors = specerrors||ap::fp_greater(err,stricterrtol);
        
        //
        // Test Gauss-Laguerre
        //
        err = 0;
        for(akind = 0; akind <= 9; akind++)
        {
            alphac = mapkind(akind);
            gqgenerategausslaguerre(n, alphac, info, x, w);
            if( info>0 )
            {
                buildgausslaguerrequadrature(n, alphac, x2, w2);
                for(i = 0; i <= n-1; i++)
                {
                    err = ap::maxreal(err, fabs(x(i)-x2(i)));
                    err = ap::maxreal(err, fabs(w(i)-w2(i)));
                }
            }
            else
            {
                specerrors = true;
            }
        }
        specerrors = specerrors||ap::fp_greater(err,nonstricterrtol);
        
        //
        // Test Gauss-Hermite
        //
        err = 0;
        gqgenerategausshermite(n, info, x, w);
        if( info>0 )
        {
            buildgausshermitequadrature(n, x2, w2);
            for(i = 0; i <= n-1; i++)
            {
                err = ap::maxreal(err, fabs(x(i)-x2(i)));
                err = ap::maxreal(err, fabs(w(i)-w2(i)));
            }
        }
        else
        {
            specerrors = true;
        }
        specerrors = specerrors||ap::fp_greater(err,nonstricterrtol);
    }
    
    //
    // end
    //
    waserrors = recerrors||specerrors;
    if( !silent )
    {
        printf("TESTING GAUSS QUADRATURES\n");
        printf("FINAL RESULT:                             ");
        if( waserrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* SPECIAL CASES (LEGENDRE/JACOBI/..)      ");
        if( specerrors )
        {
            printf("FAILED\n");
        }
        else
        {
            printf("OK\n");
        }
        printf("* RECURRENCE-BASED:                       ");
        if( recerrors )
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
Gauss-Hermite, another variant
*************************************************************************/
void buildgausshermitequadrature(int n,
     ap::real_1d_array& x,
     ap::real_1d_array& w)
{
    int i;
    int j;
    double r;
    double r1;
    double p1;
    double p2;
    double p3;
    double dp3;
    double pipm4;
    double tmp;

    x.setbounds(0, n-1);
    w.setbounds(0, n-1);
    pipm4 = pow(ap::pi(), -0.25);
    for(i = 0; i <= (n+1)/2-1; i++)
    {
        if( i==0 )
        {
            r = sqrt(double(2*n+1))-1.85575*pow(double(2*n+1), -double(1)/double(6));
        }
        else
        {
            if( i==1 )
            {
                r = r-1.14*pow(double(n), 0.426)/r;
            }
            else
            {
                if( i==2 )
                {
                    r = 1.86*r-0.86*x(0);
                }
                else
                {
                    if( i==3 )
                    {
                        r = 1.91*r-0.91*x(1);
                    }
                    else
                    {
                        r = 2*r-x(i-2);
                    }
                }
            }
        }
        do
        {
            p2 = 0;
            p3 = pipm4;
            for(j = 0; j <= n-1; j++)
            {
                p1 = p2;
                p2 = p3;
                p3 = p2*r*sqrt(double(2)/double(j+1))-p1*sqrt(double(j)/double(j+1));
            }
            dp3 = sqrt(double(2*j))*p2;
            r1 = r;
            r = r-p3/dp3;
        }
        while(ap::fp_greater_eq(fabs(r-r1),ap::machineepsilon*(1+fabs(r))*100));
        x(i) = r;
        w(i) = 2/(dp3*dp3);
        x(n-1-i) = -x(i);
        w(n-1-i) = w(i);
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-2-i; j++)
        {
            if( ap::fp_greater_eq(x(j),x(j+1)) )
            {
                tmp = x(j);
                x(j) = x(j+1);
                x(j+1) = tmp;
                tmp = w(j);
                w(j) = w(j+1);
                w(j+1) = tmp;
            }
        }
    }
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
Gauss-Legendre, another variant
*************************************************************************/
static void buildgausslegendrequadrature(int n,
     ap::real_1d_array& x,
     ap::real_1d_array& w)
{
    int i;
    int j;
    double r;
    double r1;
    double p1;
    double p2;
    double p3;
    double dp3;
    double tmp;

    x.setbounds(0, n-1);
    w.setbounds(0, n-1);
    for(i = 0; i <= (n+1)/2-1; i++)
    {
        r = cos(ap::pi()*(4*i+3)/(4*n+2));
        do
        {
            p2 = 0;
            p3 = 1;
            for(j = 0; j <= n-1; j++)
            {
                p1 = p2;
                p2 = p3;
                p3 = ((2*j+1)*r*p2-j*p1)/(j+1);
            }
            dp3 = n*(r*p3-p2)/(r*r-1);
            r1 = r;
            r = r-p3/dp3;
        }
        while(ap::fp_greater_eq(fabs(r-r1),ap::machineepsilon*(1+fabs(r))*100));
        x(i) = r;
        x(n-1-i) = -r;
        w(i) = 2/((1-r*r)*dp3*dp3);
        w(n-1-i) = 2/((1-r*r)*dp3*dp3);
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-2-i; j++)
        {
            if( ap::fp_greater_eq(x(j),x(j+1)) )
            {
                tmp = x(j);
                x(j) = x(j+1);
                x(j+1) = tmp;
                tmp = w(j);
                w(j) = w(j+1);
                w(j+1) = tmp;
            }
        }
    }
}


/*************************************************************************
Gauss-Jacobi, another variant
*************************************************************************/
static void buildgaussjacobiquadrature(int n,
     double alpha,
     double beta,
     ap::real_1d_array& x,
     ap::real_1d_array& w)
{
    int i;
    int j;
    double r;
    double r1;
    double t1;
    double t2;
    double t3;
    double p1;
    double p2;
    double p3;
    double pp;
    double an;
    double bn;
    double a;
    double b;
    double c;
    double tmpsgn;
    double tmp;
    double alfbet;
    double temp;
    int its;

    x.setbounds(0, n-1);
    w.setbounds(0, n-1);
    for(i = 0; i <= n-1; i++)
    {
        if( i==0 )
        {
            an = alpha/n;
            bn = beta/n;
            t1 = (1+alpha)*(2.78/(4+n*n)+0.768*an/n);
            t2 = 1+1.48*an+0.96*bn+0.452*an*an+0.83*an*bn;
            r = (t2-t1)/t2;
        }
        else
        {
            if( i==1 )
            {
                t1 = (4.1+alpha)/((1+alpha)*(1+0.156*alpha));
                t2 = 1+0.06*(n-8)*(1+0.12*alpha)/n;
                t3 = 1+0.012*beta*(1+0.25*fabs(alpha))/n;
                r = r-t1*t2*t3*(1-r);
            }
            else
            {
                if( i==2 )
                {
                    t1 = (1.67+0.28*alpha)/(1+0.37*alpha);
                    t2 = 1+0.22*(n-8)/n;
                    t3 = 1+8*beta/((6.28+beta)*n*n);
                    r = r-t1*t2*t3*(x(0)-r);
                }
                else
                {
                    if( i<n-2 )
                    {
                        r = 3*x(i-1)-3*x(i-2)+x(i-3);
                    }
                    else
                    {
                        if( i==n-2 )
                        {
                            t1 = (1+0.235*beta)/(0.766+0.119*beta);
                            t2 = 1/(1+0.639*(n-4)/(1+0.71*(n-4)));
                            t3 = 1/(1+20*alpha/((7.5+alpha)*n*n));
                            r = r+t1*t2*t3*(r-x(i-2));
                        }
                        else
                        {
                            if( i==n-1 )
                            {
                                t1 = (1+0.37*beta)/(1.67+0.28*beta);
                                t2 = 1/(1+0.22*(n-8)/n);
                                t3 = 1/(1+8*alpha/((6.28+alpha)*n*n));
                                r = r+t1*t2*t3*(r-x(i-2));
                            }
                        }
                    }
                }
            }
        }
        alfbet = alpha+beta;
        do
        {
            temp = 2+alfbet;
            p1 = (alpha-beta+temp*r)*0.5;
            p2 = 1;
            for(j = 2; j <= n; j++)
            {
                p3 = p2;
                p2 = p1;
                temp = 2*j+alfbet;
                a = 2*j*(j+alfbet)*(temp-2);
                b = (temp-1)*(alpha*alpha-beta*beta+temp*(temp-2)*r);
                c = 2*(j-1+alpha)*(j-1+beta)*temp;
                p1 = (b*p2-c*p3)/a;
            }
            pp = (n*(alpha-beta-temp*r)*p1+2*(n+alpha)*(n+beta)*p2)/(temp*(1-r*r));
            r1 = r;
            r = r1-p1/pp;
        }
        while(ap::fp_greater_eq(fabs(r-r1),ap::machineepsilon*(1+fabs(r))*100));
        x(i) = r;
        w(i) = exp(lngamma(alpha+n, tmpsgn)+lngamma(beta+n, tmpsgn)-lngamma(double(n+1), tmpsgn)-lngamma(n+alfbet+1, tmpsgn))*temp*pow(double(2), alfbet)/(pp*p2);
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-2-i; j++)
        {
            if( ap::fp_greater_eq(x(j),x(j+1)) )
            {
                tmp = x(j);
                x(j) = x(j+1);
                x(j+1) = tmp;
                tmp = w(j);
                w(j) = w(j+1);
                w(j+1) = tmp;
            }
        }
    }
}


/*************************************************************************
Gauss-Laguerre, another variant
*************************************************************************/
static void buildgausslaguerrequadrature(int n,
     double alpha,
     ap::real_1d_array& x,
     ap::real_1d_array& w)
{
    int i;
    int j;
    double r;
    double r1;
    double p1;
    double p2;
    double p3;
    double dp3;
    double tsg;
    double tmp;

    x.setbounds(0, n-1);
    w.setbounds(0, n-1);
    for(i = 0; i <= n-1; i++)
    {
        if( i==0 )
        {
            r = (1+alpha)*(3+0.92*alpha)/(1+2.4*n+1.8*alpha);
        }
        else
        {
            if( i==1 )
            {
                r = r+(15+6.25*alpha)/(1+0.9*alpha+2.5*n);
            }
            else
            {
                r = r+((1+2.55*(i-1))/(1.9*(i-1))+1.26*(i-1)*alpha/(1+3.5*(i-1)))/(1+0.3*alpha)*(r-x(i-2));
            }
        }
        do
        {
            p2 = 0;
            p3 = 1;
            for(j = 0; j <= n-1; j++)
            {
                p1 = p2;
                p2 = p3;
                p3 = ((-r+2*j+alpha+1)*p2-(j+alpha)*p1)/(j+1);
            }
            dp3 = (n*p3-(n+alpha)*p2)/r;
            r1 = r;
            r = r-p3/dp3;
        }
        while(ap::fp_greater_eq(fabs(r-r1),ap::machineepsilon*(1+fabs(r))*100));
        x(i) = r;
        w(i) = -exp(lngamma(alpha+n, tsg)-lngamma(double(n), tsg))/(dp3*n*p2);
    }
    for(i = 0; i <= n-1; i++)
    {
        for(j = 0; j <= n-2-i; j++)
        {
            if( ap::fp_greater_eq(x(j),x(j+1)) )
            {
                tmp = x(j);
                x(j) = x(j+1);
                x(j+1) = tmp;
                tmp = w(j);
                w(j) = w(j+1);
                w(j+1) = tmp;
            }
        }
    }
}


/*************************************************************************
Silent unit test
*************************************************************************/
bool testgq_test_silent()
{
    bool result;

    result = testgqunit(true);
    return result;
}


/*************************************************************************
Unit test
*************************************************************************/
bool testgq_test()
{
    bool result;

    result = testgqunit(false);
    return result;
}




