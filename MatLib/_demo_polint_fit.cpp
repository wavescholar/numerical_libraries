#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "polint.h"

int main(int argc, char **argv)
{
    int m;
    int n;
    ap::real_1d_array x;
    ap::real_1d_array y;
    ap::real_1d_array w;
    ap::real_1d_array xc;
    ap::real_1d_array yc;
    ap::integer_1d_array dc;
    polynomialfitreport rep;
    int info;
    barycentricinterpolant p;
    int i;
    int j;
    double a;
    double b;
    double v;
    double dv;

    printf("\n\nFitting exp(2*x) at [-1,+1] by polinomial\n\n");
    printf("Fit type             rms.err max.err    p(0)   dp(0)\n");
    
    //
    // Prepare points
    //
    m = 5;
    a = -1;
    b = +1;
    n = 1000;
    x.setlength(n);
    y.setlength(n);
    w.setlength(n);
    for(i = 0; i <= n-1; i++)
    {
        x(i) = a+(b-a)*i/(n-1);
        y(i) = exp(2*x(i));
        w(i) = 1.0;
    }
    
    //
    // Fitting:
    // a) f(x)=exp(2*x) at [-1,+1]
    // b) by 5th degree polynomial
    // c) without constraints
    //
    polynomialfit(x, y, n, m, info, p, rep);
    barycentricdiff1(p, 0.0, v, dv);
    printf("Unconstrained        %7.4lf %7.4lf %7.4lf %7.4lf\n",
        double(rep.rmserror),
        double(rep.maxerror),
        double(v),
        double(dv));
    
    //
    // Fitting:
    // a) f(x)=exp(2*x) at [-1,+1]
    // b) by 5th degree polynomial
    // c) constrained: p(0)=1
    //
    xc.setlength(1);
    yc.setlength(1);
    dc.setlength(1);
    xc(0) = 0;
    yc(0) = 1;
    dc(0) = 0;
    polynomialfitwc(x, y, w, n, xc, yc, dc, 1, m, info, p, rep);
    barycentricdiff1(p, 0.0, v, dv);
    printf("Constrained, p(0)=1  %7.4lf %7.4lf %7.4lf %7.4lf\n",
        double(rep.rmserror),
        double(rep.maxerror),
        double(v),
        double(dv));
    
    //
    // Fitting:
    // a) f(x)=exp(2*x) at [-1,+1]
    // b) by 5th degree polynomial
    // c) constrained: dp(0)=2
    //
    xc.setlength(1);
    yc.setlength(1);
    dc.setlength(1);
    xc(0) = 0;
    yc(0) = 2;
    dc(0) = 1;
    polynomialfitwc(x, y, w, n, xc, yc, dc, 1, m, info, p, rep);
    barycentricdiff1(p, 0.0, v, dv);
    printf("Constrained, dp(0)=2 %7.4lf %7.4lf %7.4lf %7.4lf\n",
        double(rep.rmserror),
        double(rep.maxerror),
        double(v),
        double(dv));
    
    //
    // Fitting:
    // a) f(x)=exp(2*x) at [-1,+1]
    // b) by 5th degree polynomial
    // c) constrained: p(0)=1, dp(0)=2
    //
    xc.setlength(2);
    yc.setlength(2);
    dc.setlength(2);
    xc(0) = 0;
    yc(0) = 1;
    dc(0) = 0;
    xc(1) = 0;
    yc(1) = 2;
    dc(1) = 1;
    polynomialfitwc(x, y, w, n, xc, yc, dc, 2, m, info, p, rep);
    barycentricdiff1(p, 0.0, v, dv);
    printf("Constrained, both    %7.4lf %7.4lf %7.4lf %7.4lf\n",
        double(rep.rmserror),
        double(rep.maxerror),
        double(v),
        double(dv));
    printf("\n\n");
    return 0;
}

