#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "mlpbase.h"

int main(int argc, char **argv)
{
    multilayerperceptron net;

    mlpcreate0(2, 1, net);
    mlprandomize(net);
    return 0;
}

