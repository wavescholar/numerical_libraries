#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "ratint.h"

int main(int argc, char **argv)
{
    int m;
    int n;
    int d;
    ap::real_1d_array x;
    ap::real_1d_array y;
    ap::real_1d_array w;
    ap::real_1d_array xc;
    ap::real_1d_array yc;
    ap::integer_1d_array dc;
    barycentricfitreport rep;
    int info;
    barycentricinterpolant r;
    int i;
    int j;
    double a;
    double b;
    double v;
    double dv;

    printf("\n\nFitting exp(2*x) at [-1,+1] by:\n1. constrained/unconstrained Floater-Hormann functions\n");
    printf("\n");
    printf("Fit type                rms.err max.err    p(0)   dp(0)  DBest\n");
    
    //
    // Prepare points
    //
    m = 5;
    a = -1;
    b = +1;
    n = 10000;
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
    // b) by 5 Floater-Hormann functions
    // c) without constraints
    //
    barycentricfitfloaterhormann(x, y, n, m, info, r, rep);
    barycentricdiff1(r, 0.0, v, dv);
    printf("Unconstrained FH        %7.4lf %7.4lf %7.4lf %7.4lf      %0ld\n",
        double(rep.rmserror),
        double(rep.maxerror),
        double(v),
        double(dv),
        long(rep.dbest));
    
    //
    // Fitting:
    // a) f(x)=exp(2*x) at [-1,+1]
    // b) by 5 Floater-Hormann functions
    // c) constrained: p(0)=1
    //
    xc.setlength(1);
    yc.setlength(1);
    dc.setlength(1);
    xc(0) = 0;
    yc(0) = 1;
    dc(0) = 0;
    barycentricfitfloaterhormannwc(x, y, w, n, xc, yc, dc, 1, m, info, r, rep);
    barycentricdiff1(r, 0.0, v, dv);
    printf("Constrained FH, p(0)=1  %7.4lf %7.4lf %7.4lf %7.4lf      %0ld\n",
        double(rep.rmserror),
        double(rep.maxerror),
        double(v),
        double(dv),
        long(rep.dbest));
    
    //
    // Fitting:
    // a) f(x)=exp(2*x) at [-1,+1]
    // b) by 5 Floater-Hormann functions
    // c) constrained: dp(0)=2
    //
    xc.setlength(1);
    yc.setlength(1);
    dc.setlength(1);
    xc(0) = 0;
    yc(0) = 2;
    dc(0) = 1;
    barycentricfitfloaterhormannwc(x, y, w, n, xc, yc, dc, 1, m, info, r, rep);
    barycentricdiff1(r, 0.0, v, dv);
    printf("Constrained FH, dp(0)=2 %7.4lf %7.4lf %7.4lf %7.4lf      %0ld\n",
        double(rep.rmserror),
        double(rep.maxerror),
        double(v),
        double(dv),
        long(rep.dbest));
    
    //
    // Fitting:
    // a) f(x)=exp(2*x) at [-1,+1]
    // b) by 5 Floater-Hormann functions
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
    barycentricfitfloaterhormannwc(x, y, w, n, xc, yc, dc, 2, m, info, r, rep);
    barycentricdiff1(r, 0.0, v, dv);
    printf("Constrained FH, both    %7.4lf %7.4lf %7.4lf %7.4lf      %0ld\n",
        double(rep.rmserror),
        double(rep.maxerror),
        double(v),
        double(dv),
        long(rep.dbest));
    printf("\n\n");
    return 0;
}

