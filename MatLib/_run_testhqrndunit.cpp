#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "testhqrndunit.h"

int main(int argc, char **argv)
{
    unsigned seed;
    if( argc==2 )
        seed = (unsigned)atoi(argv[1]);
    else
    {
        time_t t;
        seed = (unsigned)time(&t);
    }
    srand(seed);
    try
    {
        if(!testhqrndunit_test())
            return 1;
    }
    catch(ap::ap_error e)
    {
        printf("ap::ap_error:\n'%s'\n", e.msg.c_str());
        return 1;
    }
    return 0;
}

