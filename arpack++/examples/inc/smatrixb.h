/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE SMatrixB.h
   Class template for the one dimensional discrete Laplacian on
   the interval [0,1] with zero Dirichlet boundary conditions.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef SMATRIXB_H
#define SMATRIXB_H

#include "matprod.h"
#include "blas1c.h"
#include "lapackc.h"

template<class T>
class SymMatrixB: public MatrixWithProduct<T> {

 private:

  T    shift;
  T    *Ad, *Adl, *Adu, *Adu2;
  int  *ipiv;
  int  decsize;

  void FactorDataDeallocate();

 public:

  void FactorOP();

  void MultMv(T* v, T* w);

  void MultOPv(T* v, T* w);

  SymMatrixB(int nv);

  SymMatrixB(int nv, T shiftv);

  virtual ~SymMatrixB();

}; // SymMatrixB.


template<class T>
inline void SymMatrixB<T>::FactorDataDeallocate()
// Eliminates the data structure used on matrix factorization.

{

  delete[] Ad;
  delete[] Adl;
  delete[] Adu;
  delete[] Adu2;
  delete[] ipiv;

} // FactorDataDeallocate.


template<class T>
void SymMatrixB<T>::FactorOP()
/*
  Factors (M-shift*I).
*/

{

  int  i, ierr;
  T    h2;

  const T one = 1.0;
  const T two = 2.0;

  if (decsize != ncols()) {
    decsize = ncols();
    FactorDataDeallocate();
    Ad   = new T[ncols()];
    Adl  = new T[ncols()];
    Adu  = new T[ncols()];
    Adu2 = new T[ncols()];
    ipiv = new int[ncols()];
  }

  h2 = T((ncols()+1)*(ncols()+1));

  for (i=0; i<ncols(); i++) {
    Ad[i]  = two*h2 - shift;
    Adl[i] = -one*h2;
  }

  copy(ncols(), Adl, 1, Adu, 1);
  gttrf(ncols(), Adl, Ad, Adu, Adu2, ipiv, ierr);

} // FactorOP.


template<class T>
void SymMatrixB<T>::MultMv(T* v, T* w)
/*
  Matrix-vector multiplication w <- M*v.
*/

{

  int  j;
  T    h2;

  const T two = 2.0;

  w[0] =  two*v[0] - v[1];
  for (j=1; j<ncols()-1; j++) {
    w[j] = - v[j-1] + two*v[j] - v[j+1];
  }
  w[ncols()-1] = - v[ncols()-2] + two*v[ncols()-1];

  // Scaling the vector w by (1 / h^2).

  h2 = T((ncols()+1)*(ncols()+1));
  scal(ncols(), h2, w, 1L);

  return;

} //  MultMv.


template<class T>
void SymMatrixB<T>::MultOPv(T* v, T* w)
/*
  Matrix-vector multiplication w <- inv(M-shift*I)*v.
*/

{

  int  ierr;
  char *type = "N";

  copy(ncols(), v, 1, w, 1);
  gttrs(type, ncols(), 1, Adl, Ad, Adu, Adu2, ipiv, w, ncols(), ierr);

} // MultOPv


template<class T>
inline SymMatrixB<T>::SymMatrixB(int nval): MatrixWithProduct<T>(nval)
// Constructor

{

  decsize = 0;
  Ad      = 0;
  Adl     = 0;
  Adu     = 0;
  Adu2    = 0;
  ipiv    = 0;
  shift   = 0.0;

} // Constructor.


template<class T>
inline SymMatrixB<T>::
SymMatrixB(int nv, T shiftv): MatrixWithProduct<T>(nv)
// Constructor with shift.

{

  decsize = 0;
  Ad      = 0;
  Adl     = 0;
  Adu     = 0;
  Adu2    = 0;
  ipiv    = 0;
  shift   = shiftv;
  FactorOP();

} // Constructor with shift.


template<class T>
inline SymMatrixB<T>::~SymMatrixB()
// Destructor

{

  FactorDataDeallocate();

} // Destructor.


#endif // SMATRIXB_H

