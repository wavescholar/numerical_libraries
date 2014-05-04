/* This file is part of the Random Ball Cover (RBC) library.
 * (C) Copyright 2011, Lawrence Cayton [lcayton@tuebingen.mpg.de]
 */

#ifndef UTILS_H
#define UTILS_H

#include "defs.h"
#include<sys/time.h>
#include<stdint.h>
#include<stdlib.h>

void swap(unint*,unint*);
void randPerm(unint,unint*);
unint randBetween(unint,unint);
void printMat(matrix);
void printMatWithIDs(matrix,unint*);
void printCharMat(charMatrix);
void printIntMat(intMatrix);
void printVector(real*,unint);
void copyVector(real*,real*,unint);
void copyRow(matrix *x, matrix *y, unint xi, unint yi);
double timeDiff(struct timeval,struct timeval);
void copyMat(matrix*,matrix*);

void addToList(intList*,unint);
void createList(intList*);
void createSizedList(intList*,unint);
void destroyList(intList*);
void printList(intList*);

void createHeap(heap*,unint);
void destroyHeap(heap*);
void replaceMax(heap*,heapEl);
void heapSort(heap*,unint*,real*);
void reInitHeap(heap*);

void initMat(matrix *x, unint r, unint c);
size_t sizeOfMat(matrix x);
#endif
