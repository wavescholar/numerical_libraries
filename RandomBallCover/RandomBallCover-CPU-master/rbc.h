/* This file is part of the Random Ball Cover (RBC) library.
 * (C) Copyright 2011, Lawrence Cayton [lcayton@tuebingen.mpg.de]
 */
#ifndef RBC_H
#define RBC_H

#include<stdint.h>
#include "defs.h"

void buildExact(matrix,matrix*,rep*,unint);
void searchExact(matrix,matrix,matrix,rep*,unint*,real*);
void searchExactK(matrix,matrix,matrix,rep*,unint**,real**,unint);
void searchExactManyCores(matrix,matrix,matrix,rep*,unint*,real*);
void searchExactManyCoresK(matrix,matrix,matrix,rep*,unint**,real**,unint);

void buildOneShot(matrix,matrix*,rep*,unint);
void searchOneShot(matrix,matrix,matrix,rep*,unint*);
void searchOneShotK(matrix,matrix,matrix,rep*,unint**,real**,unint);

void pickReps(matrix,matrix*);
void searchStats(matrix,matrix,matrix,rep*,double*);
void reshuffleX(matrix y, matrix x, rep *ri, unint numReps);

void freeRBC(matrix r, rep *ri);
#endif
