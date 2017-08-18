#include <iostream>
#include <fstream>
#include "KMeans.h"
#include "GMM.h"
using namespace std;

int main()
{
    double data[] = {0.0, 0.2, 0.4, 0.3, 0.2, 0.4, 0.4, 0.2, 0.4, 0.5, 0.2, 0.4, 5.0, 5.2, 8.4, 6.0, 5.2, 7.4, 4.0, 5.2, 4.4, 10.3, 10.4, 10.5, 10.1, 10.6, 10.7, 11.3, 10.2, 10.9};
    //double data[] = {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,0}; 
    const int size = 30; //Number of samples
    const int dim = 1;   //Dimension of feature
    const int cluster_num = 3; //Cluster number

    KMeans* kmeans = new KMeans(dim,cluster_num);
    int* labels = new int[size];
    kmeans->SetInitMode(KMeans::InitRandom);
	kmeans->Cluster(data,size,labels);

	printf("Clustering result by k-meams:\n");
	for(int i = 0; i < cluster_num; ++i)
	{
	    printf("center %d: %f.\n", i, kmeans->GetMean(i)[0]);
	}
    printf("\n");

    //printf("center1 :%f\ncenter2 :%f\n", kmeans->GetMean(0)[0],kmeans->GetMean(1)[0]);

	delete []labels;
	delete kmeans;

    GMM *gmm = new GMM(dim,cluster_num);
    gmm->Train(data,size); //Training GMM

    printf("\nTest GMM:\n");
    for (int i = 0; i < cluster_num; ++i) {
         printf("Mean %d: %f\n",i, gmm->Mean(i)[0]);
    }
   
	

	delete gmm;

    return 0;
}
