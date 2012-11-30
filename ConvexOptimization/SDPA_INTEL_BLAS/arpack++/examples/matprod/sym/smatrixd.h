/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE SMatrixD.h
   Class template for the 1-dimensional mass matrix
   on the interval [0,1].

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef SMATRIXD_H
#define SMATRIXD_H

#include "matprod.h"
#include "blas1c.h"
#include "lapackc.h"

template<class T>
class SymMatrixD: public MatrixWithProduct<T> {

 private:

  T    *Ad, *Adl, *Adu, *Adu2;
  int  *ipiv;
  int  decsize;

  void FactorDataDeallocate();

 public:

  void FactorM();

  void SolveM(T* v);

  void MultMv(T* v, T* w);

  SymMatrixD(int nv);

  virtual ~SymMatrixD();

}; // SymMatrixD.


template<class T>
inline void SymMatrixD<T>::FactorDataDeallocate()
// Eliminates the data structure used on matrix factorization.

{

  delete[] Ad;
  delete[] Adl;
  delete[] Adu;
  delete[] Adu2;
  delete[] ipiv;

} // FactorDataDeallocate.


template<class T>
void SymMatrixD<T>::FactorM()
// Factors M.

{

  int  i, ierr;
  T    h, r1, r2;

  const T one  = 1.0;
  const T four = 4.0;
  const T six  = 6.0;

  if (decsize != ncols()) {
    decsize = ncols();
    FactorDataDeallocate();
    Ad   = new T[ncols()];
    Adl  = new T[ncols()];
    Adu  = new T[ncols()];
    Adu2 = new T[ncols()];
    ipiv = new int[ncols()];
  }

  h  = one/T(ncols()+1);
  r2 = h/six;
  r1 = r2*four;

  for (i=0; i<ncols(); i++) {
    Ad[i]  = r1;
    Adl[i] = r2;
  }

  copy(ncols(), Adl, 1, Adu, 1);
  gttrf(ncols(), Adl, Ad, Adu, Adu2, ipiv, ierr);

} // FactorM.


template<class T>
inline void SymMatrixD<T>::SolveM(T* v)
// Solves M*w = v. v is overwritten with vector w.

{

  int  ierr;
  char *type = "N";

  gttrs(type, ncols(), 1, Adl, Ad, Adu, Adu2, ipiv, v, ncols(), ierr);

} // SolveM.

template<class T>
void SymMatrixD<T>::MultMv(T* v, T* w)
//  Performs w <- M*v.

{

  int  j;
  T    h;

  const T one  = 1.0;
  const T four = 4.0;
  const T six  = 6.0;

  w[0] = four*v[0] + v[1];
  for (j=1; j<ncols()-1; j++) {
    w[j] = v[j-1] + four*v[j] + v[j+1];
  }
  w[ncols()-1] = v[ncols()-2] + four*v[ncols()-1];

  // Scaling the vector w by h.

  h = one / (T(ncols()+1)*six);
  scal(ncols(), h, w, 1L);

  return;

} //  MultMv.


template<class T>
inline SymMatrixD<T>:: SymMatrixD(int nval): MatrixWithProduct<T>(nval)
// Constructor.

{

  decsize = 0;
  Ad      = 0;
  Adl     = 0;
  Adu     = 0;
  Adu2    = 0;
  ipiv    = 0;

} // Constructor.


template<class T>
inline SymMatrixD<T>::~SymMatrixD()
// Destructor.

{

  FactorDataDeallocate();

} // Destructor.


#endif // SMATRIXD_H

