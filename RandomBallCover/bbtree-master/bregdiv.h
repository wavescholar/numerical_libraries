#ifndef BREGDIV_H
#define BREGDIV_H

//This file is part of the bregman ball tree package.
//(c) copyright 2008, Lawrence Cayton


typedef struct bregdiv{ 
  double (*div)();
  void (*divmat)();
  void (*gradf)();
  void (*gradfstar)();
} bregdiv;

bregdiv l2squared();
double calcl2div(double*,double*,int);
void calcl2gradf(double*,double*,int);
void calcl2gradfstar(double*,double*,int);
void calcl2divmat(double**,double**,double**,int,int,int);

bregdiv kl();
double calckldiv(double*,double*,int);
void calcklgradf(double*,double*,int);
void calcklgradfstar(double*,double*,int);
void calckldivmat(double**,double**,double**,int,int,int);

bregdiv dkl();
double calcdkldiv(double*,double*,int);
void calcdklgradf(double*,double*,int);
void calcdklgradfstar(double*,double*,int);
void calcdkldivmat(double**,double**,double**,int,int,int);

bregdiv is();
double calcisdiv(double*,double*,int);
void calcisgradf(double*,double*,int);
void calcisgradfstar(double*,double*,int);
void calcisdivmat(double**,double**,double**,int,int,int);

#endif
