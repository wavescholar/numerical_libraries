/* This file is part of the Random Ball Cover (RBC) library.
 * (C) Copyright 2011, Lawrence Cayton [lcayton@tuebingen.mpg.de]
 */

#ifndef BRUTE_H
#define BRUTE_H

#include<stdlib.h>
#include "defs.h"

void brutePar(matrix,matrix,unint*,float*);
void bruteK(matrix,matrix,unint**,float**,unint);
void bruteKHeap(matrix, matrix,unint**,float**, unint);
void bruteMap(matrix,matrix,rep*,unint*,unint*,float*);
void bruteMapK(matrix,matrix,rep*,unint*,unint**,float**,unint);
void bruteList(matrix,matrix,rep*,intList*,unint,unint*,float*);
void bruteListK(matrix,matrix,rep*,intList*,unint,unint**,float**,unint);
void rangeCount(matrix,matrix,float*,unint*);
#endif
