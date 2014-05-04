#ifndef __kl_fast_gauss_transform__
#define __kl_fast_gauss_transform__

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <limits.h>
#include "KCenterClustering.h"
#include "GaussTransform.h"
#include "ImprovedFastGaussTransformChooseTruncationNumber.h"
#include "ImprovedFastGaussTransformChooseParameters.h"
#include "ImprovedFastGaussTransform2.h"

#include "kl_matrix.h"

#include <crtdbg.h>
#include <iostream>

class klFastGaussTransform
{
	klMatrix<double> _X;
	klMatrix<double> _Y;
public:

	klFastGaussTransform(klMatrix<double> x,klMatrix<double> y) : _X(x), _Y(y)
	{
	}

	klVector<double> operator()(double h)
	{
		if(_X.getColumns() !=_Y.getColumns())
			throw " klFastGaussTransform error , incompatible dimensions";

		int dim = _X.getColumns();
		int N = _X.getRows();
		int M = _Y.getRows();;
		double* X  = _X.getMemory();
		double* Y = _Y.getMemory();

		double q = N;
		double epsilon = 1.0e-3f; //Error bound
		double KLimit = ceil(0.2f * sqrt((double)dim)* 100/h);

		ImprovedFastGaussTransformChooseParameters IFGTP(dim,h,epsilon,KLimit);
		_CrtCheckMemory( );

		int K = IFGTP.K; //number of clusters
		int p_max = IFGTP.p_max; // maximum truncation number
		double r = IFGTP.r; //cutoff radius


		//BBC 050414
		int *pClusterIndex = new int[K];
		for(int i=0;i<N;i++)
			*(pClusterIndex + i) =0;

		/*int *pClusterIndex = new int[N];
		for(int i=0;i<N;i++)
			*(pClusterIndex + i) =0;*/

		KCenterClustering KCC(dim,N,X,pClusterIndex,K);
		_CrtCheckMemory( );
		double *pClusterCenters = new double[K*dim];
		double *pClusterRadii = new double[K];

		KCC.Cluster();
		KCC.ComputeClusterCenters(K, pClusterCenters,pClusterIndex,pClusterRadii);	
		double rx = KCC.MaxClusterRadius; // maximum radius of the clusters
		_CrtCheckMemory( );
		// Initially the truncation number was chosen based on an estimate
		//of the maximum cluster radius. But now since we have already run
		// the clustering algorithm we know the actual maximum cluster radius.
		ImprovedFastGaussTransformChooseTruncationNumber IFGTCTN(dim,h,epsilon,rx);
		_CrtCheckMemory( );
		p_max = IFGTCTN.p_max;

		double *pWeights = new double[N];
		for(int i=0;i<N;i++)
			*(pWeights+i) = 1.0f/q;

		klVector<double> results(M);
		double *pGaussTransform =  results.getMemory();

		ImprovedFastGaussTransform IFGT(dim,N,M,X,h,pWeights,Y,double(p_max),double(K),pClusterIndex,pClusterCenters,pClusterRadii,r,epsilon,pGaussTransform);
		_CrtCheckMemory( );
		IFGT.Evaluate();

		delete pClusterIndex;
		delete pClusterCenters;
		delete pClusterRadii;
		delete pWeights;

		return results;
	}
};

#endif 