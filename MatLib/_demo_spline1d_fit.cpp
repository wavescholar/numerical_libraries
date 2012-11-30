#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "spline1d.h"

int main(int argc, char **argv)
{
    ap::real_1d_array x;
    ap::real_1d_array y;
    int n;
    int i;
    int info;
    spline1dinterpolant s;
    double t;
    spline1dfitreport rep;

    
    //
    // Fitting by unconstrained natural cubic spline
    //
    printf("FITTING BY UNCONSTRAINED NATURAL CUBIC SPLINE\n\n");
    printf("F(x)=sin(x)      function being fitted\n");
    printf("[0, pi]          interval\n");
    printf("M=4              number of basis functions to use\n");
    printf("N=100            number of points to fit\n");
    
    //
    // Create and fit
    //
    n = 100;
    x.setlength(n);
    y.setlength(n);
    for(i = 0; i <= n-1; i++)
    {
        x(i) = ap::pi()*i/(n-1);
        y(i) = sin(x(i));
    }
    spline1dfitcubic(x, y, n, 4, info, s, rep);
    
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
    }
    else
    {
        printf("\nSomething wrong, Info=%0ld",
            long(info));
    }
    return 0;
}

