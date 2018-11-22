#ifndef UTILS_H
#define UTILS_H

//This file is part of the bregman ball tree package.
//(c) copyright 2008, Lawrence Cayton

#include "bbtree.h"
#include<stdio.h>

#define ERRORTOL 1e-12
#define CLOSEENOUGH 1e-3


//Used in the test program:
#define USEL2 0
#define USEKL 1
#define USEKLD 2
#define USEIS 3


typedef struct sortitem{
  double val;
  int ind;
} sortitem;

void mean(double*,double**,int*,int,int);
int genRand(int);

int compInds(const void*, const void*);
void sortWithInds(double*,int*,int);

void readData(double**,int,int,char*);
void writeDoubs(int,char*,double,...);
void writeTree(treenode*,int,char*);
void writeNode(treenode*,int,FILE*);
treenode* readTree(char*);
treenode* readNode(int,FILE*);

void insert(int*,double*,int,double,int);

#endif

