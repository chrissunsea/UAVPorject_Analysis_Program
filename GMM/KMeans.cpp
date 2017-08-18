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
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include "KMeans.h"
using namespace std;
/*****************************************************************************/
/** Method Name: KMeans (Constructor)                                       **/
/** Method Input:                                                           **/
/**       (1)int dimension: this means the dimension of one data item.      **/
/**       (2)int clusterNum: this means the cluster number that the data    **/
/**          will be classified.                                            **/
/** Method Output:Generate an object used to implement the KMeans algorithm.**/
/** Method Description:                                                     **/
/**      This constructor will generate an object used to implement the     **/
/** KMeans algorithm. In this constructor, the data item's dimension and the**/
/** cluster number are initialised.In addition, the maximum iteration is set**/
/** to 100 (usually 100 is enough. Matlab is also set default iteration     **/
/** number to 100).                                                         **/
/**      Importantly,the this->endError can impact the accuracy of the final**/
/** result. Small value means more accurate but more time-consuming. Big    **/
/** value means less accurate but less time-consuming.                      **/
/*****************************************************************************/
KMeans::KMeans(int dimension, int clusterNum) {
	if(dimension <= 0 || clusterNum <= 0) {
		cout<<"Errors(Using KMeans): dimension or cluster number should be greater than 0!"<<endl;
		exit(-1);
	}
	this->dimension  = dimension;
	this->clusterNum = clusterNum;
	this->means = new double*[clusterNum];
	for (int i = 0; i < clusterNum; i++) {
		this->means[i] = new double[dimension];
		memset(this->means[i], 0, sizeof(double) * dimension);
	}
	this->initMode   = InitRandom;
	this->endError   = 0.001; //This value may change for high accuracy
	this->maxIterNum = 100;
}
/*****************************************************************************/
/** Method Name: ~KMeans (Deconstructor)                                    **/
/** Method Input: NULL                                                      **/
/** Method Output:Destory the object and claim back the memory allocated for**/
/**               this object.                                              **/
/** Method Description:                                                     **/
/**       The memory allocated for each KMeans object only includes the     **/
/** means value of 2-dimension array. Therefore, when destory the object,   **/
/** only the 2-d array needs to be recycled.                                **/
/*****************************************************************************/
KMeans::~KMeans() {
	for (int i = 0; i < this->clusterNum; i++) {
		delete[] this->means[i];
	}
	delete[] this->means;
}
/*****************************************************************************/
/** Method Name: Cluster                                                    **/
/** Method Input:                                                           **/
/**      (1)double *data: The data array containing all the data. When one  **/
/**         data item contains multidimension (one data record contains     **/
/**         many attributes), these data should be flatted into continuous  **/
/**         data in the data array.                                         **/
/**      (2)uint32_t dataPointNum: This attribute means the data record     **/
/**         number. It may be different from the data array's length.       **/
/**      (3)int *dataLabels: This array represents the data record labels   **/
/**         after classification.                                           **/
/**             It is also an empty matrix that has been allocated memory   **/
/**         used to return the result matrix.(This method is not responsible**/
/**        for allocating the memory for the result).                       **/
/**                                                                         **/
/** Method Output:                                                          **/
/**      (1)int *dataLabels: This array represents the data record labels   **/
/**         after classification.                                           **/
/**      (2)double** means (object inner data member): This 2-d array means **/
/**         the mean value of each cluster.                                 **/
/*****************************************************************************/
void KMeans::Cluster(double *data, uint32_t dataPointNum, int *dataLabels) {
	if(this->dimension <= 0 || this->clusterNum <= 0) {
		cout<<"Errors(Using KMeans): dimension or cluster number should be greater than 0!"<<endl;
		exit(-1);
	}
	if(dataPointNum < this->clusterNum) {
		cout<<"Errors(Using KMeans): the data number should be greater than cluster number!"<<endl;
		exit(-1);
	}

	Init(data, dataPointNum);

	int label = -1;
	int unchanged = 0;
	int iterNum  = 0;
	uint32_t index_data = 0;
	double lastCost = 0;
	double currCost = 0;
	bool loop = true;
	uint32_t* counts = new uint32_t[this->clusterNum];
	double* x = new double[this->dimension];
	double** next_means = new double*[this->clusterNum];	// New model for reestimation
	for (int i = 0; i < this->clusterNum; i++) {
		next_means[i] = new double[this->dimension];
	}

	while (loop) {
		memset(counts, 0, sizeof(uint32_t) * this->clusterNum);
		for (int i = 0; i < this->clusterNum; i++) {
			memset(next_means[i], 0, sizeof(double) * this->dimension);
		}
		lastCost = currCost;
		currCost = 0;
		// Classification
		for (uint32_t i = 0; i < dataPointNum; i++) {
			index_data = i*this->dimension;
			for(int j = 0; j < this->dimension; j++) {
				x[j] = data[index_data + j];
			}
			currCost += GetLabel(x, &label);
			counts[label]++;
			for (int d = 0; d < this->dimension; d++) {
				next_means[label][d] += x[d];
			}
		}
		currCost /= dataPointNum;
		// Reestimation
		for (int i = 0; i < this->clusterNum; i++) {
			if (counts[i] > 0) {
				for (int d = 0; d < this->dimension; d++) {
					next_means[i][d] /= counts[i];
				}
				memcpy(this->means[i], next_means[i], sizeof(double) * this->dimension);
			}
		}
		// Terminal conditions
		iterNum++;
		if (fabs(lastCost - currCost) < this->endError * lastCost) {
			unchanged++;
		}else {
			unchanged = 0;
		}
		if (iterNum >= this->maxIterNum || unchanged >= 3) {
			loop = false;
		}
		if(iterNum == this->maxIterNum && unchanged != 3){
			cout<<"KMeans(Warning): KMean have not converged!"<<endl;
		}
	}

	for (uint32_t i = 0; i < dataPointNum; i++) {
		index_data = i*this->dimension;
		for(int j = 0; j < this->dimension; j++) {
			x[j] = data[index_data + j];
		}
		GetLabel(x, &label);
		dataLabels[i] = label;
	}
	delete[] counts;
	delete[] x;
	for (int i = 0; i < this->clusterNum; i++) {
		delete[] next_means[i];
	}
	delete[] next_means;
}
/*****************************************************************************/
/** Method Name: Init                                                       **/
/** Method Input:                                                           **/
/**      (1)double *data: The data array containing all the data. When one  **/
/**         data item contains multidimension (one data record contains     **/
/**         many attributes), these data should be flatted into continuous  **/
/**         data in the data array.                                         **/
/**      (2)uint32_t dataPointNum: This attribute means the data record     **/
/**         number. It may be different from the data array's length.       **/
/** Method Output:                                                          **/
/**      (1)double** means (object inner data member): This 2-d array means **/
/**         the mean value of each cluster.                                 **/
/** Method Description:                                                     **/
/**      This method is used to initialize the object's inner data member,  **/
/** means, it can be initialized in 3 way: Random, Uniform and Mannually.   **/
/*****************************************************************************/
void KMeans::Init(double* data, uint32_t dataPointNum) {
	if(this->dimension <= 0 || this->clusterNum <= 0){
		cout<<"Errors(Using KMeans): dimension or cluster number should be greater than 0!"<<endl;
		exit(-1);
	}
	if(dataPointNum < this->clusterNum){
		cout<<"Errors(Using KMeans): the data number should be greater than cluster number!"<<endl;
		exit(-1);
	}

	uint32_t select = 0;
	uint32_t interval = dataPointNum / this->clusterNum;
	uint32_t index_data = 0;
	double* sample = new double[this->dimension];
	if (this->initMode == InitRandom) {
		srand((unsigned)time(NULL));
		for (int i = 0; i < this->clusterNum; i++) {
			select = i * interval + (interval - 1) * rand() / RAND_MAX;
			index_data = select*this->dimension;
			for(int j = 0; j < this->dimension; j++) {
				sample[j] = data[index_data + j];
			}
			memcpy(this->means[i], sample, sizeof(double) * this->dimension);
		}
	}
	else if (this->initMode == InitUniform) {
		for (int i = 0; i < this->clusterNum; i++) {
			select = i * interval;
			index_data =  select*this->dimension;
			for(int j = 0; j < this->dimension; j++) {
				sample[j] = data[index_data + j];
			}
			memcpy(this->means[i], sample, sizeof(double) * this->dimension);
		}
	}
	delete[] sample;
}
/*****************************************************************************/
/** Method Name: GetLabel                                                   **/
/** Method Input:                                                           **/
/**      (1)const double* sample: This represents one data record.          **/
/**      (2)int* label: This means the result(output) of the data record.   **/
/** Method Output:                                                          **/
/**      (1)int* label: This means the result(output) of the data record.   **/
/**                                                                         **/
/** Method Description:                                                     **/
/**      This method is to get the label of one data record according to the**/
/** distance from this point to the central data points.                    **/
/*****************************************************************************/
double KMeans::GetLabel(const double* sample, int* label) {
	double dist = -1;
	double temp = 0;
	for (int i = 0; i < this->clusterNum; i++) {
		temp = CalcDistance(sample, this->means[i], this->dimension);
		if (temp < dist || dist == -1) {
			dist = temp;
			*label = i;
		}
	}
	return dist;
}
/*****************************************************************************/
/** Method Name: CalcDistance                                               **/
/** Method Input:                                                           **/
/**      (1)const double* x: This represents one data record.               **/
/**      (2)const double* mean: This means one central point position.      **/
/**      (3)int dimension: This means the number of attributes of one data  **/
/**         record.                                                         **/
/** Method Output:                                                          **/
/**      (1)The return value(double):The distance from point 'x' to the     **/
/** central point 'mean'.                                                   **/
/**                                                                         **/
/** Method Description:                                                     **/
/**      This method is to calculate the distance according to given central**/
/** point.                                                                  **/
/*****************************************************************************/
double KMeans::CalcDistance(const double* x, const double* mean, int dimension) {
	double temp = 0;
	for (int d = 0; d < dimension; d++) {
		temp += (x[d] - mean[d]) * (x[d] - mean[d]);
	}
	return sqrt(temp);
}
/*****************************************************************************/
/** Method Name: GetInitMode                                                **/
/** Method Input: NULL                                                      **/
/** Method Output:                                                          **/
/**      The return value(int) representing the initial mode: Random,Uniform**/
/** and Mannually, according to the inner data member: enum Mode.           **/
/*****************************************************************************/
int KMeans::GetInitMode() {
	return this->initMode; 
}
/*****************************************************************************/
/** Method Name: GetMaxIterNum                                              **/
/** Method Input: NULL                                                      **/
/** Method Output:                                                          **/
/**      The return value(int) representing the max iteration number.       **/
/*****************************************************************************/
int KMeans::GetMaxIterNum() { 
	return this->maxIterNum; 
}
/*****************************************************************************/
/** Method Name: GetEndError                                                **/
/** Method Input: NULL                                                      **/
/** Method Output:                                                          **/
/**      The return value(double) representing end error rate.              **/
/*****************************************************************************/
double KMeans::GetEndError() { 
	return this->endError; 
}
/*****************************************************************************/
/** Method Name: GetEndError                                                **/
/** Method Input: NULL                                                      **/
/** Method Output:                                                          **/
/**      The return value(double*) representing the ith cluster's central   **/
/** point.                                                                  **/
/*****************************************************************************/
double* KMeans::GetMean(int i) { 
	return this->means[i]; 
}
/*****************************************************************************/
/** Method Name: SetInitMode                                                **/
/** Method Input: (1) int i: the initial mode.                              **/
/** Method Output:                                                          **/
/**     this->initMode will be set to Random,Uniform or Mannually, according**/
/** to the inner data member: enum Mode.                                    **/
/*****************************************************************************/
void KMeans::SetInitMode(int i) { 
	this->initMode = i; 
}
/*****************************************************************************/
/** Method Name: SetMaxIterNum                                              **/
/** Method Input: (1) int i: the max interation number.                     **/
/** Method Output:                                                          **/
/**     this->maxIterNum will be set to the givern number.                  **/
/*****************************************************************************/
void KMeans::SetMaxIterNum(int i) { 
	this->maxIterNum = i; 
}
/*****************************************************************************/
/** Method Name: SetEndError                                                **/
/** Method Input: (1) double f: the end error rate.                         **/
/** Method Output:                                                          **/
/**     this->endError will be set to the givern number.                    **/
/*****************************************************************************/
void KMeans::SetEndError(double f) { 
	this->endError = f; 
}
/*****************************************************************************/
/** Method Name: SetMean                                                    **/
/** Method Input:                                                           **/
/**      (1) int i: the ith cluster.                                        **/
/**      (2) const double* u: the givern mean number.                       **/
/** Method Output:                                                          **/
/**     this->means will be set to the givern values.                       **/
/*****************************************************************************/
void KMeans::SetMean(int i, const double* u) { 
	memcpy(this->means[i], u, sizeof(double) * this->dimension); 
}
