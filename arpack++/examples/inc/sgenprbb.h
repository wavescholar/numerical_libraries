/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE SGenPrbB.h
   Very simple template class intended to illustrate how to 
   use ARPACK++ to find some few eigenvalues and eigenvectors 
   of symmetric generalized problems in shift and invert,
   buckling and Cayley modes.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef SGENPRBB_H
#define SGENPRBB_H

#include "blas1c.h"
#include "lapackc.h"
#include "smatrixc.h"
#include "smatrixd.h"

template<class T>
class SymGenProblemB {

 private:

  int  n, decsize;
  int  *ipiv;
  T    shift;
  T    *Ad, *Adl, *Adu, *Adu2;

  void FactorDataDeallocate();
  // Eliminates the data structure used on matrix factorization.

  void FactorOP();
  // Factors (A-shift*B).

 public:

  SymMatrixC<T> A;
  SymMatrixD<T> B;

  void MultAv(T* v, T* w) { A.MultMv(v,w); }
  // Matrix-vector product w <- A*v,
  // where the matrix is the 1 dimensional discrete Laplacian on
  // the interval [0,1] with zero Dirichlet boundary conditions.

  void MultBv(T* v, T* w) { B.MultMv(v,w); }
  // Matrix-vector product w <- B*v,
  // where the matrix is the 1 dimensional mass matrix
  // on the interval [0,1].

  void MultOPv(T* v, T* w);
  // Matrix-vector product w <- inv(A-shift*B)*v.

  SymGenProblemB(int nx, T shiftp);
  // Constructor.

  ~SymGenProblemB();
  // Destructor

}; // SymGenProblemB


template<class T>
inline void SymGenProblemB<T>::FactorDataDeallocate()
{

  delete[] Ad;
  delete[] Adl;
  delete[] Adu;
  delete[] Adu2;
  delete[] ipiv;

} // FactorDataDeallocate.


template<class T>
void SymGenProblemB<T>::FactorOP()
{

  int  i, ierr;
  T    h, r1, r2;

  const T one  = 1.0;
  const T two  = 2.0;
  const T four = 4.0;
  const T six  = 6.0;

  if (decsize != n) {
    decsize = n;
    FactorDataDeallocate();
    Ad   = new T[n];
    Adl  = new T[n];
    Adu  = new T[n];
    Adu2 = new T[n];
    ipiv = new int[n];
  }

  h  = one/T(n+1);
  r1 =  two/h - shift*h*four/six;
  r2 = -one/h - shift*h*one/six;

  for (i=0; i<n; i++) {
    Ad[i]  = r1;
    Adl[i] = r2;
  }

  copy(n, Adl, 1, Adu, 1);
  gttrf(n, Adl, Ad, Adu, Adu2, ipiv, ierr);

} // FactorOP.


template<class T>
void SymGenProblemB<T>::MultOPv(T* v, T* w)
{

  int  ierr;
  char *type = "N";

  copy(n, v, 1, w, 1);
  gttrs(type, n, 1, Adl, Ad, Adu, Adu2, ipiv, w, n, ierr);

} // MultOPv.


template<class T>
inline SymGenProblemB<T>::SymGenProblemB(int nx, T shiftp): A(nx), B(nx)
{

  shift   = shiftp;
  decsize = 0;
  Ad      = 0;
  Adl     = 0;
  Adu     = 0;
  Adu2    = 0;
  ipiv    = 0;
  n       = A.ncols(); 
  FactorOP();

} // Constructor.


template<class T>
inline SymGenProblemB<T>::~SymGenProblemB()
{

  FactorDataDeallocate();

} // Destructor.


#endif // SGENPRBB_H

