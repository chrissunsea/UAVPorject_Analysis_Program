/*****************************************************************************/
/**Author: Gaojie Sun                                                       **/
/**Date: 21/07/2017                                                         **/
/**Supervisor: Dongryeol Ryu                                                **/
/**Originazation: The University of Melbourne                               **/
/*****************************************************************************/

/*****************************************************************************/
/**Documentation Description:                                               **/
/**   This file is used to define the class that abstracts the fitting data **/
/** to GMM. This algorithm is used to fit the data into Gaussian Mixture    **/
/** Model with K distributions. In this project, the Gaussian Mixture Model **/
/**is used to classify the crop temperature information and soil temperature**/
/** information.                                                            **/
/*****************************************************************************/
#pragma once
#include <fstream>
#include <stdint.h>

class GMM
{
	public:
		~GMM();
		GMM(int dimNum = 1, int mixNum = 1);

		void Copy(GMM* gmm);

		void SetMaxIterNum(int i);
		void SetEndError(double f);

		int GetDimNum();
		int GetMixNum();
		int GetMaxIterNum();
		double GetEndError();

		double* Mean(int i);
		double Prior(int i);
		double* Variance(int i);

		void setMean(int i,double *val);
		void setPrior(int i,double val);
		void setVariance(int i,double *val);

		double GetProbability(const double* sample);

		void Init(double *data, uint32_t dataPointNum);
		void Train(double *data, uint32_t dataPointNum);

	private:
		int dimension;
		int mixtureNum;
		int maxIterNum;	
		
		double*  priors;
		double** means;
		double** vars;

		double* minVars;
		double endError;

	private:
		// Return the "j"th pdf, p(x|j).
		double GetProbability(const double* x, int j);
		void Allocate();
		void Dispose();
};
