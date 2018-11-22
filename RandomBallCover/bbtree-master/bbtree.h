#ifndef BBTREE_H
#define BBTREE_H

//This file is part of the bregman ball tree package.
//(c) copyright 2008, Lawrence Cayton

#include "bregdiv.h"

typedef struct treenode{
  double* mu;
  double* mup;
  double divToMu; //used in search procedure
  double R;
  int* inds; /* indices stored in this node */
  int n;
  int isLeaf;
  struct treenode *l;
  struct treenode *r;
} treenode;

treenode* buildTree(double**,int,int,bregdiv,int);
treenode* setupNode(int,int);
treenode* rBuild(double**,int*,int,int,bregdiv,int);
void findSplit(int*, int*, int*,int*, double**, int*,int,int,bregdiv);
void findSplitMaxVar(int*, int*, int*,int*, double**, int*,int,int,bregdiv);
void findSplitKmeans(int*, int*, int*,int*, double**, int*,int,int,bregdiv);
int calcFarthest(double **, int*, int, int, int,bregdiv);
void calcDvector(double *,double**,int*,int,int,int,bregdiv);
void calcDvectorNonInd(double *,double **,int*,int,int,double*,bregdiv);
void deleteTree(treenode*);
#endif
