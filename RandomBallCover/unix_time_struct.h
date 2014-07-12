#pragma once
#include <time.h>

#ifndef _INC_KL_UNIX_TIME_STRUCT
#define _INC_KL_UNIX_TIME_STRUCT

//BBCREVISIT 070614- time stuff from Linux - just to get code compiling.
typedef long                    __kernel_time_t;
typedef long                    __kernel_suseconds_t;
struct timeval {
        __kernel_time_t         tv_sec;         /* seconds */
       __kernel_suseconds_t    tv_usec;        /* microseconds */
 };
 inline int gettimeofday(timeval *tp, void* nu =NULL)
 {
	 time_t time_of_day;
	struct tm *tm_buf;
	time_of_day = time( NULL );
	tm_buf=localtime(&time_of_day);


	 tp->tv_sec = tm_buf->tm_sec;
	 tp->tv_usec = clock();

	 return  tp->tv_usec;

 }
 //END BBCREVISIT



#endif