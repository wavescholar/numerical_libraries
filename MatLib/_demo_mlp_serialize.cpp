#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "mlpbase.h"

int main(int argc, char **argv)
{
    multilayerperceptron network1;
    multilayerperceptron network2;
    multilayerperceptron network3;
    ap::real_1d_array x;
    ap::real_1d_array y;
    ap::real_1d_array r;
    int rlen;
    double v1;
    double v2;

    
    //
    // Generate two networks filled with small random values.
    // Use MLPSerialize/MLPUnserialize to make network copy.
    //
    mlpcreate0(1, 1, network1);
    mlpcreate0(1, 1, network2);
    mlpserialize(network1, r, rlen);
    mlpunserialize(r, network2);
    
    //
    // Now Network1 and Network2 should be identical.
    // Let's demonstrate it.
    //
    printf("Test serialization/unserialization\n");
    x.setlength(1);
    y.setlength(1);
    x(0) = 2*ap::randomreal()-1;
    mlpprocess(network1, x, y);
    v1 = y(0);
    printf("Network1(X) = %0.2lf\n",
        double(y(0)));
    mlpprocess(network2, x, y);
    v2 = y(0);
    printf("Network2(X) = %0.2lf\n",
        double(y(0)));
    if( ap::fp_eq(v1,v2) )
    {
        printf("Results are equal, OK.\n");
    }
    else
    {
        printf("Results are not equal... Strange...");
    }
    return 0;
}

