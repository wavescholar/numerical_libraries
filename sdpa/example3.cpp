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
/* The beginning of the ``example3.cpp''. */
#include <stdio.h>
#include <stdlib.h>

#include "sdpa-lib.hpp"
#include "sdpa-lib2.hpp"

void PrintMatrix(int nRow, int nCol, double* element, FILE* fpOut);

extern "C" int sdpa_EX3 ()
{int argc;
char *argv[3];
	double*  	element;

        if (argc != 3)
        {
                fprintf(stderr, "%s [Input] [Output] \n", argv[0]);
                exit(EXIT_FAILURE);
        }

        SDPA	Problem1;

        Problem1.Method          = KSH;
	
        strcpy(Problem1.ParameterFileName, "param.sdpa");
        Problem1.ParameterFile = fopen(Problem1.ParameterFileName, "r");
        strcpy(Problem1.InputFileName, argv[1]);
        Problem1.InputFile = fopen(Problem1.InputFileName,"r");
        strcpy(Problem1.OutputFileName, argv[2]);
        Problem1.OutputFile = fopen(Problem1.OutputFileName,"w");

        Problem1.DisplayInformation = stdout;
	
        SDPA_initialize(Problem1);
        SDPA_Solve(Problem1);

	fprintf(stdout, "yMat = \n");
	// Problem1.printResultYMat(stdout);
	fprintf(stdout, "{\n");
	for (int l=0; l<Problem1.nBLOCK; ++l) {
	  element = Problem1.getResultYMat(l);
	  int nRow = Problem1.bLOCKsTRUCT[l];
	  int nCol = nRow;
	  if (nRow<0) {
	    nCol = -nCol;
	    nRow = 1;
	  }
	  PrintMatrix(nRow,nCol,element,stdout);
	}
	fprintf(stdout, "}\n");

        fclose(Problem1.InputFile);
        fclose(Problem1.OutputFile);
	
        Problem1.Delete();

        return 0;
};	

#define FORMAT "%+8.16e"

//void PrintMatrix(int nRow, int nCol, double* element, FILE* fpOut)
//{
//  fprintf(fpOut,"{");
//  for (int i=1; i<=nRow-1; ++i) {
//    if (i==1) {
//      fprintf(fpOut," ");
//    } else {
//      fprintf(fpOut,"  ");
//    }
//    fprintf(fpOut,"{");
//    for (int j=1; j<=nCol-1; ++j) {
//      fprintf(fpOut, FORMAT",",element[(i-1)+nRow*(j-1)]);
//    }
//    fprintf(fpOut,FORMAT" },\n",element[(i-1)+nRow*(nCol-1)]);
//  }
//  if (nRow>1) {
//    fprintf(fpOut,"  {");
//  }
//  for (int j=1; j<=nCol-1; ++j) {
//    fprintf(fpOut,FORMAT",",element[(nRow-1)+nRow*(j-1)]);
//  }
//  fprintf(fpOut,FORMAT" }",element[(nRow-1)+nRow*(nCol-1)]);
//  if (nRow>1) {
//    fprintf(fpOut,"   }\n");
//  } else {
//    fprintf(fpOut,"\n");
//  }
//}

/* The end of the ``example3.cpp''. */
