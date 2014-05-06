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
#include "kl_stat.h"
#include "kl_random_number_generator.h"
#include "kl_time_series.h"
#include "kl_multivariate_random_variable.h"
#include "kl_sample_population.h"
#include "kl_principal_components.h"
#include "kl_large_dataset.h"
#include "kl_regression.h"
#include "kl_multiclass_svm.h"
#include "kl_wavelet.h"
#include "kl_ML_helper_fns.h"
#include "kl_bregman_divergence.h"
#include "kl_util.h"
#include "kl_unit_tests.h"
#include "kl_matrix_facorizations.h"
#include "kl_img_pca.h"
#include "kl_matlab_dependent_unit_tests.h"
#include "kl_matlab_iface.h"
#include "kl_arpack.h"
#include "kl_fast_gauss_transform.h"
#include <crtdbg.h>
#include <iostream>


//March 01, 2011.
//BBC Sadly, I could not get this code to work reliably.  After 32+ hours, I belive it is time to give up.
//I did not learn anything by debugging the code, nor did I learn anything wrapping up an interface.  
//The time would have been better spent with a reading programme centered on the multipole methods, and doing 
//some upfron prototyping in Matlab.  This effort needs to be revitalized in the future and approached from that
//direction. 

void testFastGaussTransform(const char* fileName)
{
	ios_base::openmode wMode = ios_base::app;
	ofstream _tex(fileName);
	time_t time_of_day;
	struct tm *tm_buf;
	time_of_day = time( NULL );
	tm_buf=localtime(&time_of_day);
	startLatexDoc("Linear Regression","KL Software Libraries",asctime(tm_buf),_tex,""); 

	char * arg = new char[256];

	unsigned int i;
	unsigned int j;

	unsigned int kernelDim=64;
	unsigned int featureDim=4;

	unsigned int testPoints = 10;

	klVector<double> meanVector1(featureDim);  meanVector1 = 0;
	
	klMatrix<double> covariance = IdentityMatrix<double>(featureDim);// klGenerateRandomSymmetricPositiveDefiniteMatrix<double>(featureDim,fastSeed() );
	
	double sigma = mean(covariance.diag().getMemory(),featureDim);

	klNormalMultiVariate<double> X1(meanVector1,covariance );

	klSamplePopulation<double> X(kernelDim,featureDim);
	klSamplePopulation<double> Y(testPoints,featureDim);

	klVector<double> one(featureDim); one = 1.0;
	for(j=0;j<kernelDim;j++)
	{
		X.setRow(j,one);//X.setRow(j,X1());

	}
	klVector<double> y(featureDim); y=1.0;
	for(j=0;j<testPoints;j++)
	{	y = y+  0.1;
		Y.setRow(j,y);
	}
	LatexPrintMatrix(X,"X",_tex);
	LatexPrintMatrix(Y,"Y",_tex);

	klFastGaussTransform klfg(X,Y);
	double kernelBandwidth =sigma;
	klVector<double> dist =klfg(kernelBandwidth);

	klMatrix<double> dm(dist.getMemory(),dist.getRows(),1);
	LatexPrintMatrix(dm,"D",_tex);

	endLatexDoc(_tex);
	flushall();
	_tex.close();
	delete arg;
}


klVector<double> klFastGaussTransformTest()
{
	unsigned int i;
	unsigned int j;
	unsigned int numTrainingPoints=64;
	unsigned int numFeatureDimensions=5;
	unsigned int numTestPoints = 64;
	double q = numTrainingPoints;//Total Weight of Gaussians

	klVector<double> meanVector1(numFeatureDimensions);  meanVector1 = 0.0;
	klVector<double> meanVector2(numFeatureDimensions);  meanVector2 = 0.5;
	klMatrix<double> covarianceMatrix1 =klGenerateRandomSymmetricPositiveDefiniteMatrix<double>(numFeatureDimensions,fastSeed() );
	klMatrix<double> covarianceMatrix2 =klGenerateRandomSymmetricPositiveDefiniteMatrix<double>(numFeatureDimensions,fastSeed() );

	klNormalMultiVariate<double> X1(meanVector1,covarianceMatrix1 );
	klNormalMultiVariate<double> X2(meanVector2, covarianceMatrix2);

	klSamplePopulation<double> x(numTrainingPoints,numFeatureDimensions);
	klSamplePopulation<double> y(numTrainingPoints,numFeatureDimensions);
	double* X  =x.getMemory();
	double* Y = y.getMemory();

	klVector<double> one (3);one = 1.0;
	klVector<double> oneplus(3); oneplus = 1.1;

	for(j=0;j<numTrainingPoints;j++)
	{
		x.setRow(j,X1());

	}
	for(j=0;j<numTrainingPoints;j++)
	{
		y.setRow(j,X2());
	}

	double h = 1.0; //Kernel Bandwidth
	double epsilon = 1.0e-3f; //Error bound
	
	double KLimit = ceil(0.2f * sqrt((double)numFeatureDimensions)* 100/h);

	ImprovedFastGaussTransformChooseParameters IFGTP(numFeatureDimensions,h,epsilon,KLimit);
	
	//_CrtCheckMemory( );

	int K = IFGTP.K; //number of clusters
	int p_max = IFGTP.p_max; // maximum truncation number
	double r = IFGTP.r; //cutoff radius

	int *pClusterIndex = new int[numTrainingPoints];
	for(int i=0;i<numTrainingPoints;i++)
		*(pClusterIndex + i) =0;

	KCenterClustering KCC(numFeatureDimensions,numTrainingPoints,X,pClusterIndex,K);

	double *pClusterCenters = new double[K*numFeatureDimensions];
	double *pClusterRadii = new double[K];

	KCC.Cluster();
	KCC.ComputeClusterCenters(K, pClusterCenters,pClusterIndex,pClusterRadii);	
	double rx = KCC.MaxClusterRadius; // maximum radius of the clusters

	// Initially the truncation number was chosen based on an estimate
	//of the maximum cluster radius. But now since we have already run
	// the clustering algorithm we know the actual maximum cluster radius.
	ImprovedFastGaussTransformChooseTruncationNumber IFGTCTN(numFeatureDimensions,h,epsilon,rx);
	p_max = IFGTCTN.p_max;
	
	double *pWeights = new double[numTrainingPoints];
	for(int i=0;i<numTrainingPoints;i++)
		*(pWeights+i) = 1.0f/q;

	klVector<double> results(numTestPoints);
	double *pGaussTransform =  results.getMemory();
			
	ImprovedFastGaussTransform IFGT(numFeatureDimensions,numTrainingPoints,numTestPoints,X,h,pWeights,Y,double(p_max),double(K),pClusterIndex,pClusterCenters,pClusterRadii,r,epsilon,pGaussTransform);
	IFGT.Evaluate();
	
	for(int i=0;i<numTestPoints;i++)
		std::cout<<*(pGaussTransform +i)<<std::endl;

	return results;
}





void klFastGaussTransformTestParams()
{
	unsigned int i;
	unsigned int j;
	unsigned int numTrainingPoints=64;
	unsigned int numFeatureDimensions=5;
	unsigned int numTestPoints = 64;
	double q = numTrainingPoints;//Total Weight of Gaussians

	klVector<double> meanVector1(numFeatureDimensions);  meanVector1 = 0.0;
	klVector<double> meanVector2(numFeatureDimensions);  meanVector2 = 0.5;
	klMatrix<double> covarianceMatrix1 =klGenerateRandomSymmetricPositiveDefiniteMatrix<double>(numFeatureDimensions,fastSeed() );
	klMatrix<double> covarianceMatrix2 =klGenerateRandomSymmetricPositiveDefiniteMatrix<double>(numFeatureDimensions,fastSeed() );

	klNormalMultiVariate<double> X1(meanVector1,covarianceMatrix1 );
	klNormalMultiVariate<double> X2(meanVector2, covarianceMatrix2);

	klSamplePopulation<double> x(numTrainingPoints,numFeatureDimensions);
	klSamplePopulation<double> y(numTrainingPoints,numFeatureDimensions);
	double* X  =x.getMemory();
	double* Y = y.getMemory();

	klVector<double> one (3);one = 1.0;
	klVector<double> oneplus(3); oneplus = 1.1;

	for(j=0;j<numTrainingPoints;j++)
	{
		x.setRow(j,X1());

	}
	for(j=0;j<numTrainingPoints;j++)
	{
		y.setRow(j,X2());
	}

	double h = 1.0; //Kernel Bandwidth
	double epsilon = 1.0e-3f; //Error bound
	
	for(int i =1; i<128;i++)
	{
		numFeatureDimensions =i;
		double KLimit = ceil(0.2f * sqrt((double)numFeatureDimensions)* 100/h);
		ImprovedFastGaussTransformChooseParameters IFGTP(numFeatureDimensions,h,epsilon,KLimit);
		int K = IFGTP.K; //number of clusters
		int p_max = IFGTP.p_max; // maximum truncation number
		double r = IFGTP.r; //cutoff radius

		int *pClusterIndex = new int[numTrainingPoints];
		for(int i=0;i<numTrainingPoints;i++)
			*(pClusterIndex + i) =0;

		KCenterClustering KCC(numFeatureDimensions,numTrainingPoints,X,pClusterIndex,K);

		double *pClusterCenters = new double[K*numFeatureDimensions];
		double *pClusterRadii = new double[K];

		KCC.Cluster();
		KCC.ComputeClusterCenters(K, pClusterCenters,pClusterIndex,pClusterRadii);	
		double rx = KCC.MaxClusterRadius; 
		ImprovedFastGaussTransformChooseTruncationNumber IFGTCTN(numFeatureDimensions,h,epsilon,rx);
		p_max = IFGTCTN.p_max;

		ImprovedFastGaussTransformChooseTruncationNumber IFGTCTN2(i,h,epsilon,rx);
		p_max = IFGTCTN2.p_max;
		cerr<<"n = "<<i<<" K = "<<K<<" p_max = "<<IFGTP.p_max<<" r = "<<r<<" p_max after cluster= "<<IFGTCTN.p_max<<endl;
	}

	
}

