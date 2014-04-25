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
/* The beginning of the ``example1-1.cpp''. */
#include <stdio.h>
#include <stdlib.h>
extern "C" __declspec(dllexport) int spda_EX12_(void);
extern "C" __declspec(dllexport) int  sdpa_EX21(void);
extern "C" __declspec(dllexport) int sdpa_EX22(void);
extern "C" __declspec(dllexport) int sdpa_EX2r(void);
extern "C" __declspec(dllexport) int sdpa_EX3(void);
extern "C" __declspec(dllexport) int sdpa_EX4(void);
extern "C" __declspec(dllexport) int SDPA_TEST(void);

#include "sdpa-lib.hpp"
#include "sdpa-lib2.hpp"
extern "C"	int spda_EX12_();
extern "C"	int sdpa_EX21();
extern "C"	int sdpa_EX22();
extern "C" int sdpa_EX2r();
extern "C"	int sdpa_EX3();
extern "C"	int sdpa_EX4();
extern "C"	int sdpa_EX5();
extern "C"	int sdpa_EX6();

//int main ()
//{
//	SDPA_TEST();
//}



extern "C" int SDPA_TEST ()
{	int ans22=sdpa_EX22();
	int ans21=sdpa_EX21();
	int ans2r=sdpa_EX2r();	
	int ans4=sdpa_EX4();
int ans5=sdpa_EX5();

//these tests do not work as of 080909
//int ans6=sdpa_EX6();
//	int ans12=spda_EX12_();
//	
//    int ans3=sdpa_EX3();
	/*    if (argc != )
        {
                fprintf(stderr, "%s [Input] [Output] \n", argv[0]);
                exit(EXIT_FAILURE);
        }*/
	//	"Example 1: mDim = 3, nBLOCK = 1, {2}"
 //  3  =  mDIM
 //  1  =  nBOLCK
 //  2  = bLOCKsTRUCT
	//{48, -8, 20}
	//{ {-11,  0}, { 0, 23} }
	//{ { 10,  4}, { 4,  0} }
	//{ {  0,  0}, { 0, -8} }
	//{ {  0, -8}, {-8, -2} }
//
//"Example 1: mDim = 3, nBLOCK = 1, {2}"
//   3  =  mDIM
//   1  =  nBOLCK
//   2  = bLOCKsTRUCT
//	{48, -8, 20}
//	0 1 1 1 -11
//	0 1 2 2 23
//	1 1 1 1 10
//	1 1 1 2 4
//	2 1 2 2 -8
//	3 1 1 2 -8
//	3 1 2 2 -2
        SDPA	Problem1;

        Problem1.Method          = KSH;
	
        strcpy(Problem1.ParameterFileName, "C:/KL/Packages/ConvexOptimization/SDPA_INTEL_BLAS/sdpa/param.sdpa");
     
		Problem1.ParameterFile = fopen("C:/KL/Packages/ConvexOptimization/SDPA_INTEL_BLAS/sdpa/example1.dat-s", "r");
        
		strcpy(Problem1.InputFileName, "C:/KL/Packages/ConvexOptimization/SDPA_INTEL_BLAS/sdpa/example1.dat-s");
     	Problem1.InputFile = fopen("C:/KL/Packages/ConvexOptimization/SDPA_INTEL_BLAS/sdpa/example1.dat-s","r");

        strcpy(Problem1.OutputFileName, "C:/KL/Packages/ConvexOptimization/SDPA_INTEL_BLAS/sdpa/example1.result_loc");
        Problem1.OutputFile = fopen("C:/KL/Packages/ConvexOptimization/SDPA_INTEL_BLAS/sdpa/example1.result_loc","w");

        Problem1.DisplayInformation = stdout;
	
        SDPA_initialize(Problem1);
        SDPA_Solve(Problem1);

        fclose(Problem1.InputFile);
        fclose(Problem1.OutputFile);
	
        Problem1.Delete();



        return 0;
};	



/* The end of the ``example1-1.cpp''. */
