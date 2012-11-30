#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "spline1d.h"

int main(int argc, char **argv)
{
    ap::real_1d_array x;
    ap::real_1d_array y;
    ap::real_1d_array w;
    ap::real_1d_array xc;
    ap::real_1d_array yc;
    ap::integer_1d_array dc;
    int n;
    int i;
    int info;
    spline1dinterpolant s;
    double t;
    spline1dfitreport rep;

    
    //
    // Fitting by constrained Hermite spline
    //
    printf("FITTING BY CONSTRAINED HERMITE SPLINE\n\n");
    printf("F(x)=sin(x)      function being fitted\n");
    printf("[0, pi]          interval\n");
    printf("M=6              number of basis functions to use\n");
    printf("S(0)=0           first constraint\n");
    printf("S(pi)=0          second constraint\n");
    printf("N=100            number of points to fit\n");
    
    //
    // Create and fit:
    // * X  contains points
    // * Y  contains values
    // * W  contains weights
    // * XC contains constraints locations
    // * YC contains constraints values
    // * DC contains derivative indexes (0 = constrained function value)
    //
    n = 100;
    x.setlength(n);
    y.setlength(n);
    w.setlength(n);
    for(i = 0; i <= n-1; i++)
    {
        x(i) = ap::pi()*i/(n-1);
        y(i) = sin(x(i));
        w(i) = 1;
    }
    xc.setlength(2);
    yc.setlength(2);
    dc.setlength(2);
    xc(0) = 0;
    yc(0) = 0;
    dc(0) = 0;
    xc(0) = ap::pi();
    yc(0) = 0;
    dc(0) = 0;
    spline1dfithermitewc(x, y, w, n, xc, yc, dc, 2, 6, info, s, rep);
    
    //
    // Output results
    //
    if( info>0 )
    {
        printf("\nOK, we have finished\n\n");
        printf("     x   F(x)   S(x)  Error\n");
        t = 0;
        while(ap::fp_less(t,0.999999*ap::pi()))
        {
            printf("%6.3lf %6.3lf %6.3lf %6.3lf\n",
                double(t),
                double(sin(t)),
                double(spline1dcalc(s, t)),
                double(fabs(spline1dcalc(s, t)-sin(t))));
            t = ap::minreal(ap::pi(), t+0.25);
        }
        printf("%6.3lf %6.3lf %6.3lf %6.3lf\n\n",
            double(t),
            double(sin(t)),
            double(spline1dcalc(s, t)),
            double(fabs(spline1dcalc(s, t)-sin(t))));
        printf("rms error is %6.3lf\n",
            double(rep.rmserror));
        printf("max error is %6.3lf\n",
            double(rep.maxerror));
        printf("S(0) = S(pi) = 0 (exactly)\n\n");
    }
    else
    {
        printf("\nSomething wrong, Info=%0ld",
            long(info));
    }
    return 0;
}

