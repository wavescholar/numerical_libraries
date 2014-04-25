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
  rsdpa_include.h
--------------------------------------------------*/

#ifndef __rsdpa_include_h__
#define __rsdpa_include_h__

#include "rsdpa_right.h"

// if you use ATLAS, you need to set 0
// otherwise (for example, BLAS in clapack.tgz), set 1
//           and edit Makefile to change LAPACK_LIB

#define NON_ATLAS_SDPA 1

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
//
#include "mkl_lapack.h"
#include "mkl_blas.h"
extern "C" {
#include "f2c.h"
//#if NON_ATLAS_SDPA
//#include "blaswrap.h"
//#endif
//#include "fblaswr.h"
//#include "cblas.h"
//#include "clapack.h"
};



using namespace std;

#define _SUCCESS true
#define FAILURE false

#include "rsdpa_tool.h"

#endif // __rsdpa_include_h__
