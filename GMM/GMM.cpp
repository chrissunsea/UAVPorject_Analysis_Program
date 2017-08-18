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
#include <string.h>
#include <math.h>
#include <iostream>
#include "GMM.h"
#include "KMeans.h"
using namespace std;
/*****************************************************************************/
/** Method Name: GMM (Constructor)                                          **/
/** Method Input:                                                           **/
/**       (1)int dimension: this means the dimension of one data item.      **/
/**       (2)int mixtureNum: this means how many normal distributions the   **/
/**          GMM will have.                                                 **/
/** Method Output: Generate an object used to implement the GMM algorithm.  **/
/** Method Description:                                                     **/
/**      This constructor will generate an object used to implement the     **/
/** GMM algorithm. In this constructor, the data item's dimension and the   **/
/** distribution number are initialised. In addition, the maximum iteration **/
/** is set to 100 (usually 100 is enough. Matlab is also set default        **/
/** iteration number to 100).                                               **/
/**      Importantly,the this->endError can impact the accuracy of the final**/
/** result. Small value means more accurate but more time-consuming. Big    **/
/** value means less accurate but less time-consuming.                      **/
/*****************************************************************************/
GMM::GMM(int dimension, int mixtureNum) {
	if (dimension <= 0 || mixtureNum <= 0) {
		cout<<"Errors(GMM): The dimension number and mixture number should be greater than 0!"<<endl;
		exit(-1);
	}

	this->dimension  = dimension;
	this->mixtureNum = mixtureNum;

	this->maxIterNum = 100;
	this->endError   = 0.000001; //This value may change for high accuracy

	Allocate();

	for (int i = 0; i < this->mixtureNum; i++) {
		this->priors[i] = 1.0 / this->mixtureNum;
		for (int d = 0; d < this->dimension; d++) {
			this->means[i][d] = 0;
			this->vars[i][d]  = 1;
		}
	}
}
/*****************************************************************************/
/** Method Name: ~KMeans (Deconstructor)                                    **/
/** Method Input: NULL                                                      **/
/** Method Output:Destory the object and claim back the memory allocated for**/
/**               this object.                                              **/
/*****************************************************************************/
GMM::~GMM() {
	Dispose();
}
/*****************************************************************************/
/** Method Name: Allocate                                                   **/
/** Method Input: NULL                                                      **/
/** Method Output:                                                          **/
/**     This method is used to allocate the memory for the object.          **/
/*****************************************************************************/
void GMM::Allocate() {
	this->minVars = new double [this->dimension];
	this->priors  = new double [this->mixtureNum];
	this->means   = new double*[this->mixtureNum];
	this->vars    = new double*[this->mixtureNum];
	
	for (int i = 0; i < this->mixtureNum; i++) {
		this->means[i] = new double[this->dimension];
		this->vars[i]  = new double[this->dimension];
	}
}
/*****************************************************************************/
/** Method Name: Dispose                                                    **/
/** Method Input: NULL                                                      **/
/** Method Output:                                                          **/
/**     This method is used to claim back the memory allocated for the      **/
/** object.                                                                 **/
/*****************************************************************************/
void GMM::Dispose() {
	delete[] this->priors;
	delete[] this->minVars;
	for (int i = 0; i < this->mixtureNum; i++) {
		delete[] this->means[i];
		delete[] this->vars[i];
	}
	delete[] this->means;
	delete[] this->vars;
}
/*****************************************************************************/
/** Method Name: Copy                                                       **/
/** Method Input:                                                           **/
/**    (1) GMM* gmm: The source object.                                     **/
/** Method Description: This method is used to copy a GMM object.           **/
/*****************************************************************************/
void GMM::Copy(GMM* gmm) {

	if (this->mixtureNum != gmm->mixtureNum || this->dimension != gmm->dimension) {
		cout<<"Errors(GMM): Cannot copy gmm since the dimension and the mixture number of the 2 objects are different!"<<endl;
		exit(-1);
	}

	for (int i = 0; i < this->mixtureNum; i++) {
		this->priors[i] = gmm->Prior(i);
		memcpy(this->means[i], gmm->Mean(i), sizeof(double) * this->dimension);
		memcpy(this->vars[i], gmm->Variance(i), sizeof(double) * this->dimension);
	}
	memcpy(this->minVars, gmm->minVars, sizeof(double) * this->dimension);
}
/*****************************************************************************/
/** Method Name: GetProbability                                             **/
/** Method Input:                                                           **/
/**    (1) const double* sample: The data reocrd that needs to compute its  **/
/**        probability in the GMM.                                          **/
/** Method Description: This method is used to compute the probability of   **/
/** the given data. The probability of one data record is a weighted average**/
/** probability of multiple normal distributions.                           **/
/*****************************************************************************/
double GMM::GetProbability(const double* sample) {
	double p = 0;
	for (int i = 0; i < this->mixtureNum; i++) {
		p += this->priors[i] * GetProbability(sample, i);
	}
	return p;
}
/*****************************************************************************/
/** Method Name: GetProbability                                             **/
/** Method Input:                                                           **/
/**    (1) const double* x: The data reocrd that needs to compute its       **/
/**        probability in a certain normal distribution of the GMM.         **/
/**    (2) int j: the ID of the normal distribution of the GMM.             **/
/** Method Description: This method is used to compute the probability of   **/
/** the given data in a certain normal distribution of the GMM. This method **/
/** is called by double GMM::GetProbability(const double* sample).          **/
/*****************************************************************************/
double GMM::GetProbability(const double* x, int j) {
	const double PI=3.1415926535897932;
	double p = 1;
	for (int d = 0; d < this->dimension; d++) {
		p *= 1 / sqrt(2 * PI * this->vars[j][d]);
		p *= exp(-0.5 * (x[d] - this->means[j][d]) * (x[d] - this->means[j][d]) / this->vars[j][d]);
	}
	return p;
}
/*****************************************************************************/
/** Method Name: Train                                                      **/
/** Method Input:                                                           **/
/**    (1) double *data: This is the data set that need to be applied to the**/
/**        GMM. These data is flatted into one dimension. In fact, one data **/
/**        item may contain multiple dimensions (multiple data fileds).These**/
/**        data items can be separated into individual items according to   **/
/**        the inner data this->dimension initialised in the constructor.   **/
/**    (2) uint32_t dataPointNum: This is the number of the data item.      **/
/** Method Output:                                                          **/
/**    (1) (object's data member) double** means: The mu(s) of the GMM.     **/
/**	   (2) (object's data member)double** vars: The Sigma of the GMM.       **/
/** Method Description:                                                     **/
/**    This method, first applying the KMeans on the data to cluster the    **/
/**    data.Then the mu and sigma of each cluster are the initial parameters**/
/**    After that, using the EM algorithm to train the GMM.                 **/
/*****************************************************************************/
void GMM::Train(double *data, uint32_t dataPointNum) {
	//Using kMeans to initialise the GMM's paramemters
	Init(data,dataPointNum);
	// Reestimation
	uint32_t dataIndex = 0;
	int unchanged = 0;
	int iterNum = 0;
	double lastL = 0;
	double currL = 0;
	bool loop = true;

	double*  x = new double[this->dimension];
	double*  next_priors = new double [this->mixtureNum];
	double** next_vars   = new double*[this->mixtureNum];
	double** next_means  = new double*[this->mixtureNum];

	for (int i = 0; i < this->mixtureNum; i++) {
		next_means[i] = new double[this->dimension];
		next_vars[i]  = new double[this->dimension];
	}

	double p = 0;
	double pj = 0;
	while (loop) {
		// Clear buffer for reestimation
		memset(next_priors, 0, sizeof(double) * this->mixtureNum);
		for (int i = 0; i < this->mixtureNum; i++) {
			memset(next_vars[i], 0, sizeof(double) * this->dimension);
			memset(next_means[i], 0, sizeof(double) * this->dimension);
		}

		lastL = currL;
		currL = 0;

		// Predict
		for (int k = 0; k < dataPointNum; k++) {

			dataIndex = k * this->dimension;
			for(int j = 0; j < this->dimension; j++){
				x[j]=data[dataIndex + j];
			}

			p = GetProbability(x);

			for (int j = 0; j < this->mixtureNum; j++) {
				pj = GetProbability(x, j) * this->priors[j] / p;

				next_priors[j] += pj;

				for (int d = 0; d < this->dimension; d++) {
					next_means[j][d] += pj * x[d];
					next_vars[j][d] += pj* x[d] * x[d];
				}
			}

			currL += (p > 1E-20) ? log10(p) : -20;
		}
		currL /= dataPointNum;

		// Reestimation: generate new priors, means and variances.
		for (int j = 0; j < this->mixtureNum; j++) {
			this->priors[j] = next_priors[j] / dataPointNum;
			if (this->priors[j] > 0) {
				for (int d = 0; d < this->dimension; d++) {
					this->means[j][d] = next_means[j][d] / next_priors[j];
					this->vars[j][d] = next_vars[j][d] / next_priors[j] - this->means[j][d] * this->means[j][d];
					if (this->vars[j][d] < this->minVars[d]){
						this->vars[j][d] = this->minVars[d];
					}
				}
			}
		}
		// Terminal conditions
		iterNum++;
		if (fabs(currL - lastL) < this->endError * fabs(lastL)) {
			unchanged++;
		}else {
			unchanged = 0;
		}
		if (iterNum >= this->maxIterNum || unchanged >= 3) {
			loop = false;
		}
		if(iterNum == this->maxIterNum && unchanged != 3){
			cout<<"GMM(Warning): KMean have not converged!"<<endl;
		}
	}
	delete[] next_priors;
	for (int i = 0; i < this->mixtureNum; i++) {
		delete[] next_means[i];
		delete[] next_vars[i];
	}
	delete[] next_means;
	delete[] next_vars;
	delete[] x;
}
/*****************************************************************************/
/** Method Name: Init                                                       **/
/** Method Input:                                                           **/
/**    (1) double *data: This is the data set that need to be applied to the**/
/**        GMM. These data is flatted into one dimension. In fact, one data **/
/**        item may contain multiple dimensions (multiple data fileds).These**/
/**        data items can be separated into individual items according to   **/
/**        the inner data this->dimension initialised in the constructor.   **/
/**    (2) uint32_t dataPointNum: This is the number of the data item.      **/
/** Method Output:                                                          **/
/**    (1) (object's data member) double** means: The mu(s) of the GMM.     **/
/**	   (2) (object's data member)double** vars: The Sigma of the GMM.       **/
/** Method Description:                                                     **/
/**    This method, applying the KMeans on the data to cluster the data.Then**/
/**    the mu and sigma of each cluster are the initial parameters.         **/
/*****************************************************************************/
void GMM::Init(double *data, uint32_t dataPointNum) {

	if (this->dimension <=0 || this->mixtureNum <= 0) {
		cout<<"Errors(GMM): Initialization failed, the dimension number and components number should be greater than 0!"<<endl;
		exit(-1);
	}

	if(dataPointNum < this->mixtureNum){
		cout<<"Errors(GMM): the data number should be greater than components number!"<<endl;
		exit(-1);
	}

	const double MIN_VAR = 1E-10;

	int* dataLabels;
	dataLabels = new int[dataPointNum];
	KMeans* kmeans = new KMeans(this->dimension, this->mixtureNum);
	kmeans->Cluster(data, dataPointNum, dataLabels);

	uint32_t* counts = new uint32_t[this->mixtureNum];
	double* overallMeans = new double[this->dimension];	// Overall mean of training data
	for (int i = 0; i < this->mixtureNum; i++) {
		counts[i] = 0;
		this->priors[i] = 0;
		memcpy(this->means[i], kmeans->GetMean(i), sizeof(double) * this->dimension);
		memset(this->vars[i], 0, sizeof(double) * this->dimension);
	}
	memset(overallMeans, 0, sizeof(double) * this->dimension);
	memset(this->minVars, 0, sizeof(double) * this->dimension);

	int label = -1;
	uint32_t dataIndex = 0;
	double* x = new double[this->dimension];
	
	for (uint32_t i = 0; i < dataPointNum; i++) {

		dataIndex = i * this->dimension;
		for(int j = 0; j < this->dimension; j++) {
			x[j]=data[dataIndex + j];
		}

		label=dataLabels[i];
		// Count each Gaussian
		counts[label]++;
		double* m = kmeans->GetMean(label);
		for (int d = 0; d < this->dimension; d++) {
			this->vars[label][d] += (x[d] - m[d]) * (x[d] - m[d]);
		}
		// Count the overall mean and variance.
		for (int d = 0; d < this->dimension; d++) {
			overallMeans[d] += x[d];
			this->minVars[d] += x[d] * x[d];
		}
	}
	// Compute the overall variance (* 0.01) as the minimum variance.
	for (int d = 0; d < this->dimension; d++) {
		overallMeans[d] /= dataPointNum;
		this->minVars[d] = max(MIN_VAR, 0.01 * (this->minVars[d] / dataPointNum - overallMeans[d] * overallMeans[d]));
	}
	// Initialize each Gaussian.
	for (int i = 0; i < this->mixtureNum; i++) {
		this->priors[i] = 1.0 * counts[i] / dataPointNum;
		if (this->priors[i] > 0) {
			for (int d = 0; d < this->dimension; d++) {
				this->vars[i][d] = this->vars[i][d] / counts[i];
				// A minimum variance for each dimension is required.
				if (this->vars[i][d] < this->minVars[d]) {
					this->vars[i][d] = this->minVars[d];
					cout<<"smaller than min var"<<endl;
				}
			}
		} else {
			memcpy(this->vars[i], this->minVars, sizeof(double) * this->dimension);
			cout << "[WARNING] Gaussian " << i << " of GMM is not used!\n";
		}
	}
	delete kmeans;
	delete[] x;
	delete[] counts;
	delete[] overallMeans;
	delete[] dataLabels;
}

void GMM::SetMaxIterNum(int i) { 
	this->maxIterNum = i; 
}

void GMM::SetEndError(double f) { 
	this->endError = f; 
}

int GMM::GetDimNum() {
	return this->dimension; 
}

int GMM::GetMixNum() {
	return this->mixtureNum; 
}

int GMM::GetMaxIterNum() {
	return this->maxIterNum;
}

double GMM::GetEndError() {
	return this->endError; 
}

double GMM::Prior(int i) {
	return this->priors[i];
}

double* GMM::Mean(int i) {
	return this->means[i];
}

double* GMM::Variance(int i) {
	return this->vars[i];
}

void GMM::setPrior(int i,double val) {
	this->priors[i]=val;
}

void GMM::setMean(int i,double *val) {
	for(int j = 0; j < this->dimension; j++) {
		this->means[i][j] = val[j];
	}
}

void GMM::setVariance(int i,double *val) {
	for(int j = 0; j < this->dimension; j++){
		this->vars[i][j] = val[j];
	}
}
