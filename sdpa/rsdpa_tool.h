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
/*--------------------------------------------------
  rsdpa_tool.h
--------------------------------------------------*/

#ifndef __rsdpa_tool_h__
#define __rsdpa_tool_h__

#include "rsdpa_right.h"

#include <iostream>
#include <time.h>
#include <sys/types.h>
#include <sys/timeb.h>
// #include <sys/time.h>
#include <string>

#if 1
#define rMessage(message) \
cout << message << " :: line " << __LINE__ \
  << " in " << __FILE__ << endl
#else
#define rMessage(message)
#endif

#define rError(message) \
cout << message << " :: line " << __LINE__ \
  << " in " << __FILE__ << endl; \
//exit(false)

#if 0
#define rNewCheck() rMessage("new invoked");
#else
#define rNewCheck() ;
#endif

#define REVERSE_PRIMAL_DUAL 1


// These are constant. Do NOT change
extern int IZERO   ; // =  0;
extern int IONE    ; // =  1;
extern int IMONE   ; // = -1;
extern double DZERO; // =  0.0;
extern double DONE ; // =  1.0;
extern double DMONE; // = -1.0;

struct rrealtime {
  time_t ltime;
  _timeb tstruct;
};

class rTime
{
public:
  static double rGetUserTime();
  static void rSetTimeVal(rrealtime& targetVal);
  static double rGetRealTime(rrealtime& start,
			     rrealtime& end);
};

#if 1 // count time with process time
#define rTimeStart(START__) \
   static clock_t START__; START__ = clock();
#define rTimeEnd(END__) \
   static clock_t END__;   END__ = clock();
#define rTimeCal(START__,END__) (((double)(END__ - START__)) / CLOCKS_PER_SEC);

#else // count time with real time
#define rTimeStart(START__) \
   static rrealtime START__; rTime::rSetTimeVal(START__)
#define rTimeEnd(END__) \
   static rrealtime END__; rTime::rSetTimeVal(END__)
#define rTimeCal(START__,END__) rTime::rGetRealTime(START__,END__)
#endif

#endif // __rsdpa_tool_h__
