#ifndef EXPOKIT_H
#define EXPOKIT_H

typedef std::complex<double> cplx;

namespace expokit
{

typedef void(*matvecfunc_double)(void *data, double *in, double *out);
typedef void(*matvecfunc_cplx)(void *data, cplx *in, cplx *out);

extern "C"
{
#include "f2c.h"

//int dgexpv(integer *n, integer *m, doublereal *t, 
//	doublereal *v, doublereal *w, doublereal *tol, doublereal *anorm, 
//	doublereal *wsp, integer *lwsp, integer *iwsp, integer *liwsp, S_fp 
//	matvec, void* matvecdata, integer *itrace, integer *iflag)
int dgexpv(int *n, int *m, double *t, 
	double *v, double *w, double *tol, double *anorm, 
	double *wsp, int *lwsp, int *iwsp, int *liwsp, matvecfunc_double 
	matvec, void* matvecdata, int *itrace, int *iflag);

//int dgpadm_(integer *ideg, integer *m, doublereal *t, 
//	doublereal *h__, integer *ldh, doublereal *wsp, integer *lwsp, 
//	integer *ipiv, integer *iexph, integer *ns, integer *iflag)
int dgpadm_(int *ideg, int *m, double *t, 
	double *h__, int *ldh, double *wsp, int *lwsp, 
	int *ipiv, int *iexph, int *ns, int *iflag);

//int dgpadm(int *ideg, int *m, double *t, 
//	double *h__, int *ldh, double *wsp, int *lwsp, 
//	int *ipiv, int *iexph, int *ns, int *iflag);

int dsexpv(int *n, int *m, double *t, 
	double *v, double *w, double *tol, double *anorm, 
	double *wsp, int *lwsp, int *iwsp, int *liwsp, matvecfunc_double 
	matvec, void *matvecdata, int *itrace, int *iflag);

int zgexpv(int *n, int *m, double *t, 
	cplx *v, cplx *w, double *tol, double *
	anorm, cplx *wsp, int *lwsp, int *iwsp, int *
	liwsp, matvecfunc_cplx matvec, void *matvecdata, int *itrace, int *iflag);

int zhexpv(int *n, int *m, double *t, 
	cplx *v, cplx *w, double *tol, double *
	anorm, cplx *wsp, int *lwsp, int *iwsp, int *
	liwsp, matvecfunc_cplx matvec, void *matvecdata, int *itrace, int *iflag);

}

};

#endif

