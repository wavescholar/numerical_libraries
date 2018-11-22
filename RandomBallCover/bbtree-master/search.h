#ifndef SEARCH_H
#define SEARCH_H

//This file is part of the bregman ball tree package.
//(c) copyright 2008, Lawrence Cayton

void multisearch(treenode*,double**,double**,bregdiv,int,int,int,double,int,int);
void rsearch(treenode*,double*,double*,double**,bregdiv,int,double,int);
int needToSearch(treenode*,double*,double*,bregdiv,int,double,double);
int recBinSearch(double,double,treenode*,double*,double*,bregdiv, int,double,double);

#endif
