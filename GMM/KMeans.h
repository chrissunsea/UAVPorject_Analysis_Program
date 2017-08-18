/*****************************************************************************/
/**Author: Gaojie Sun                                                       **/
/**Date: 21/07/2017                                                         **/
/**Supervisor: Dongryeol Ryu                                                **/
/**Originazation: The University of Melbourne                               **/
/*****************************************************************************/

/*****************************************************************************/
/**Documentation Description:                                               **/
/**   This file is used to define the class that abstracts the KMeans       **/
/** algorithm. This algorithm is used to classify the data into k clusters. **/
/** In this project, the Gaussian Mixture Model is based on this method to  **/
/** initialise the GMM's initial parameters.                                **/
/*****************************************************************************/
#pragma once
#include <fstream>
#include <stdint.h>

class KMeans
{
	public:
		enum InitMode {
			InitRandom,
			InitManual,
			InitUniform,
		};

		~KMeans();
		KMeans(int dimension = 1, int clusterNum = 1);

		int GetInitMode();
		int GetMaxIterNum();
		double  GetEndError();
		double* GetMean(int i);

		void SetInitMode(int i);
		void SetMaxIterNum(int i);
		void SetEndError(double f);
		void SetMean(int i, const double* u);

		void Init(double* data, uint32_t dataPointNum);
		void Cluster(double* data, uint32_t dataPointNum, int *dataLabels) ;

	private:
		int initMode;
		int maxIterNum;
		int dimension;
		int clusterNum;
		double** means;
		double   endError;

		double GetLabel(const double* x, int* label);
		double CalcDistance(const double* x, const double* u, int dimNum);
};