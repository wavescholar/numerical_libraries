/* -------------------------------------------------------------

This file is a component of SDPA
Copyright (C) 2004 SDPA Project

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

------------------------------------------------------------- */
/*-----------------------------------------
  rsdpa_tool.cpp
-----------------------------------------*/

#include "rsdpa_tool.h"
// #include <sys/times.h>
// #include <sys/time.h>
#include <time.h>

// #include <unistd.h>
#ifndef CLK_TCK
#define  CLK_TCK  sysconf(_SC_CLK_TCK)
#endif

// These are constant.
// Do Not Change .
int IZERO =  0;
int IONE  =  1;
int IMONE = -1;
double DZERO =  0.0;
double DONE  =  1.0;
double DMONE = -1.0;

double rTime::rGetUserTime()
{
  return clock();
  #if 0
  struct tms TIME;
  times(&TIME);
  return (double)TIME.tms_utime/(double)CLK_TCK; 
  #endif
}

void rTime::rSetTimeVal(rrealtime& targetVal)
{
  time(&targetVal.ltime);
  _ftime(&targetVal.tstruct);
}

double rTime::rGetRealTime(rrealtime& start, rrealtime& end)
{
  const long int second = end.ltime-start.ltime;
  const long int usecond = end.tstruct.millitm - start.tstruct.millitm;
  return ((double)second) + ((double)usecond)*(1.0e-6);
}

