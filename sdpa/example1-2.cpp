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
/* The beginning of the ``example1-2.cpp''. */
#include <stdio.h>
#include <stdlib.h>

#include "sdpa-lib.hpp"
#include "sdpa-lib2.hpp"

extern "C" int spda_EX12_ ()
{
	

	SDPA	Problem1;

	 Problem1.CheckMatrix     = false;
        Problem1.Method          = KSH;
	
		/////////////////////////////////problem 2
		strcpy(Problem1.ParameterFileName, "C:/KL/Packages/ConvexOptimization/SDPA_INTEL_BLAS/sdpa/param.sdpa");
		Problem1.ParameterFile = fopen("C:/KL/Packages/ConvexOptimization/SDPA_INTEL_BLAS/sdpa/example2.dat-s", "r");
		strcpy(Problem1.InputFileName, "C:/KL/Packages/ConvexOptimization/SDPA_INTEL_BLAS/sdpa/example2.dat-s");
     	Problem1.InputFile = fopen("C:/KL/Packages/ConvexOptimization/SDPA_INTEL_BLAS/sdpa/example12.dat-s","r");
        strcpy(Problem1.OutputFileName, "C:/KL/Packages/ConvexOptimization/SDPA_INTEL_BLAS/sdpa/example2.result_loc");
        Problem1.OutputFile = fopen("C:/KL/Packages/ConvexOptimization/SDPA_INTEL_BLAS/sdpa/example2.result_loc","w");

	strcpy(Problem1.InitialFileName, "C:/KL/Packages/ConvexOptimization/SDPA_INTEL_BLAS/sdpa/example1.ini");
	Problem1.InitialFile = fopen("C:/KL/Packages/ConvexOptimization/SDPA_INTEL_BLAS/sdpa/example1.ini","r");

	Problem1.InitialPoint = true;		/* If one want to set an initial point */



	Problem1.DisplayInformation = stdout;
	
	SDPA_initialize(Problem1);
	SDPA_Solve(Problem1);

	fclose(Problem1.InputFile);
	fclose(Problem1.OutputFile);
	
	Problem1.Delete();

	return 0;
};	
/* The end of the ``example1-2.cpp''. */
