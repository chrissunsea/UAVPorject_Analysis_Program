/*****************************************************************************/
/**Author: Gaojie Sun                                                       **/
/**Date: 21/07/2017                                                         **/
/**Supervisor: Dongryeol Ryu                                                **/
/**Originazation: The University of Melbourne                               **/
/*****************************************************************************/

/*****************************************************************************/
/**Documentation Description:                                               **/
/**   This file is used to define the class that abstracts one of edge      **/
/** detection algorithms, called Canny Edge Detection. The implementation   **/
/** is based on the GaussianBlur class and SobelFilter class.               **/
/*****************************************************************************/
#ifndef _CANNY_EDGE_DETECTION_CPP
#define _CANNY_EDGE_DETECTION_CPP

#include <cmath>
#include "../GMM/GMM.h"
#include "CannyEdgeDetection.hpp"
#include <fstream>
using namespace std;
/*****************************************************************************/
/** Method Name: CannyEdgeDetector (Constructor)                            **/
/** Method Input: NULL                                                      **/
/** Method Output: Generate an object used to implement the canny.          **/
/** Method Description:                                                     **/
/**       This constructor does not need any parameter.It initialises the   **/
/** GaussianBlur object and SobelFilter object.                             **/
/*****************************************************************************/
template <typename SRC_T>
CannyEdgeDetector<SRC_T>::CannyEdgeDetector() {
    this->lowThreshold = 0;
    this->highThreshold = 0;

    const int size = 5;
    const double gausSigma = 1.0;

    this->sobelFilter = new SobelFilter<double>();
    this->gaussianBlur = new GaussianBlur<SRC_T>(size, gausSigma);
}
/*****************************************************************************/
/** Method Name: ~CannyEdgeDetector (Deconstructor)                         **/
/** Method Input: NULL                                                      **/
/** Method Output:Destory the object and claim back the memory allocated for**/
/**               this object.                                              **/
/** Method Description:                                                     **/
/**       The memory allocated for each CannyEdgeDetector object includes   **/
/** the objects of GaussianBlur class and SobelFilter class. Therefore, when**/
/** destory the canny object,these objects also need to be recycled.        **/
/*****************************************************************************/
template <typename SRC_T>
CannyEdgeDetector<SRC_T>::~CannyEdgeDetector() {
    delete this->gaussianBlur;
    delete this->sobelFilter;
    this->gaussianBlur = NULL;
    this->sobelFilter = NULL;
}
/*****************************************************************************/
/** Method Name: DetectEdge                                                 **/
/** Method Input:                                                           **/
/**    (1)SRC_T* srcImageMatrix: The input (gray) image matrix that need to **/
/**       detect the edge.The element type can be any basic type.In addition**/
/**       ,the matrix is a one-dimension array.                             **/
/**    (2)char* resImageMatrix: An empty matrix that has been allocated     **/
/**       memory used to return the result of the edge matrix. (This method **/
/**       is not responsible for allocating the memory for the result). In  **/
/**       addition,edge points are set to 1 and no edge points are sets to 0**/
/**    (3)uint32_t height: The input matrix's height.                       **/
/**    (4)uint32_t width: The input matrix's width.                         **/
/** Method Output:                                                          **/
/**    (1)char* resImageMatrix:The result of the edge matrix.All its element**/
/**       are binary. In addition, edge points are set to 1 and no edge     **/
/**       points are sets to 0.                                             **/
/** Method Description:                                                     **/
/**     This method implements the canny edge detection.It provide a service**/
/** to its caller, which takes the original (gray) image matrix (and its    **/
/** height and width) as input, then return the binary edge matrix (edge    **/
/** points are 1 and no edge points are 0).                                 **/
/**      To be more specific,this method first uses Gaussian blur to smooth **/
/** the original image,and then uses sobel filter to calculate the gradient **/
/** for x and y direction. Next,it also calculates the theta of the gradient**/
/** ( atan2(Gy/Gx) ) and the total Gradient( G=sqrt(Gx*Gx + Gy*Gy) ),then   **/
/** apply them to Non Maximum Spression step to filter the points.After that**/
/** , based on the total gradient G, using Gaussian Mixtrure Model to       **/
/** calculate the high threshold and low threshold. Then, based on these 2  **/
/** thresholds to filter the edge points(If the G of the point is greater   **/
/** than high threshold,it is definitely a edge point. If the G of the point**/
/** is smaller than the low threshold, it is definitely a non edge point.   **/
/** However, if the G of one point is between low threshold and high        **/
/** threshold, it is not clear whether or not it is a edge point). Next,    **/
/** the edge tracking step will tackle these points whose G values are      **/
/** between high threshold and low threshold. The detail method is described**/
/** in method EdgeTrack. Finally, for the rest points, they will be set to  **/
/** no edge points.                                                         **/
/*****************************************************************************/
template <typename SRC_T>
void CannyEdgeDetector<SRC_T>::DetectEdge(SRC_T* srcImageMatrix, char* resImageMatrix, uint32_t height,uint32_t width) {
    
    double* G  = new double[height*width];
    double* Gx = new double[height*width];
    double* Gy = new double[height*width];
    double* theta = new double[height*width];
    double* intermediateRes = new double[height*width];

    this->gaussianBlur->Blur(srcImageMatrix, intermediateRes, height, width);

    this->sobelFilter->Filter(this->sobelFilter->DIRECTION_X, intermediateRes, Gx, height, width);
    this->sobelFilter->Filter(this->sobelFilter->DIRECTION_Y, intermediateRes, Gy, height, width);
    this->sobelFilter->CombineGradient(Gx, Gy, G, height, width, theta);
    
    memset(intermediateRes, 0, height*width*sizeof(double));
    this->NonMaximumSuppression(intermediateRes, theta, Gx, Gy, G, height, width);

    delete[] G;
    delete[] Gx;
    delete[] Gy;
    delete[] theta;

    this->SetDoubleThreshold(intermediateRes, height, width);
    
    list<PointPosition*> strongEdgeList;
    list<PointPosition*> weakEdgeList;

    this->DoubleThresholdFilter(intermediateRes, height, width, &strongEdgeList, &weakEdgeList);

    this->EdgeTrack(intermediateRes, height, width, &strongEdgeList);
    this->CleanUp(intermediateRes, height, width, &weakEdgeList);
    
    strongEdgeList.clear();
    weakEdgeList.clear();

    for (uint32_t i = 0; i < height; ++i) {
        for (uint32_t j = 0; j < width; ++j) {
            resImageMatrix[i*width+j] = intermediateRes[i*width+j] > 0? 1 : 0;
        }
    }

    delete[] intermediateRes;
}
/*****************************************************************************/
/** Method Name: NonMaximumSuppression                                      **/
/** Method Input:                                                           **/
/**    (1) double* resImageMatrix: An empty matrix that has been allocated  **/
/**        memory used to return the result matrix. (This method is not     **/
/**        responsible for allocating the memory for the result).           **/
/**    (2)double* theta: The (Gradient) degree of the tangent [-180,180].   **/
/**                      theta = atan(Gy/Gx).                               **/
/**    (3)double* Gx: The matrix (one-dimension array) for each pixel's     **/ 
/**       gradient of x direction.                                          **/
/**    (4)double* Gy: The matrix (one-dimension array) for each pixel's     **/ 
/**       gradient of y direction.                                          **/
/**    (5)double* G: The combination of the 2 perpendicular direction's     **/
/**       (x and y) gradients.                                              **/
/**    (6)uint32_t height: The input matrix's height.                       **/
/**    (7)uint32_t width: The input matrix's width.                         **/
/** Method Output:                                                          **/
/**    (1)double* resImageMatrix: The matrix is same to the matrix G (the   **/
/**       total gradient of the each pixel), expect that the non local      **/
/**       maximum values are set to 0.                                      **/
/** Method Description:                                                     **/
/**     This method is to filter these non edge points. If the Gradient of  **/
/** one point is not a local maximum point, its Gradient will be set to 0,  **/
/** which means it is not a edge point. Otherwise, its Gradient will not be **/
/** changed. The method used to judge whether one point is a local maximum  **/
/** point is interpolation method (this is more accurate).                  **/
/*****************************************************************************/
template <typename SRC_T>
void CannyEdgeDetector<SRC_T>::NonMaximumSuppression(double* resImageMatrix, double* theta, double* Gx, double* Gy, double* G, uint32_t height, uint32_t width) {

    int gridSize = 3;
    double g1,g2,g3,g4;
    double dWeight, dTmp1, dTmp2;
    double grid[gridSize][gridSize];

    uint32_t rowOffset, dataIndex;
    for (uint32_t i = 1; i < height-1; ++i){
        rowOffset = i*width;
        for (uint32_t j = 1; j < width - 1; ++j){
            grid[0][0] = G[(i-1)*width + j-1];
            grid[0][1] = G[(i-1)*width + j];
            grid[0][2] = G[(i-1)*width + j+1];

            grid[1][0] = G[(i)*width + j-1];
            grid[1][1] = G[(i)*width + j];
            grid[1][2] = G[(i)*width + j+1];

            grid[2][0] = G[(i+1)*width + j-1];
            grid[2][1] = G[(i+1)*width + j];
            grid[2][2] = G[(i+1)*width + j+1];

            dataIndex = rowOffset + j;
            ////////////////////situation 1////////////////////// 
            //|---------------------------------------->X////////
            //|//////       g1                      /////////////  
            //|//////       g2  C   g3              /////////////  
            //|//////               g4              /////////////
            //Y//////////////////////////////////////////////////  
            if( (theta[dataIndex] >= 0 && theta[dataIndex] <= 45) ||
                (theta[dataIndex] >= -180 && theta[dataIndex] < -135) ) {  
                g1 = grid[0][0];
                g2 = grid[1][0];
                g3 = grid[1][2];
                g4 = grid[2][2];

                dWeight = abs(Gy[dataIndex] / G[dataIndex]);
                dTmp1 = (g4-g3)*dWeight + g3;  
                dTmp2 = (g1-g2)*dWeight + g2; 
            }     
            ////////////////////situation 2////////////////////// 
            //|---------------------------------------->X////////
            //|//////       g1  g2                  /////////////  
            //|//////           C                   /////////////  
            //|//////           g3    g4            /////////////
            //Y//////////////////////////////////////////////////  
            else if( (theta[dataIndex] > 45 && theta[dataIndex] <= 90) ||
                     (theta[dataIndex] >= -135 && theta[dataIndex] < -90) ) {  
                g1 = grid[0][0];
                g2 = grid[0][1];
                g3 = grid[2][1];
                g4 = grid[2][2];

                dWeight = abs(Gx[dataIndex] / G[dataIndex]); 
                dTmp1 = (g4-g3)*dWeight + g3;  
                dTmp2 = (g1-g2)*dWeight + g2;
            }
            ////////////////////situation 3////////////////////// 
            //|---------------------------------------->X////////
            //|//////           g1    g2            /////////////  
            //|//////           C                   /////////////  
            //|//////     g3    g4                 //////////////
            //Y//////////////////////////////////////////////////  
            else if( (theta[dataIndex] > 90 && theta[dataIndex] <= 135) ||
                     (theta[dataIndex] >= -90 && theta[dataIndex] < -45) ) { 
                g1 = grid[0][1];
                g2 = grid[0][2];
                g3 = grid[2][0];
                g4 = grid[2][1];

                dWeight = abs(Gx[dataIndex] / G[dataIndex]);
                dTmp1 = (g3-g4)*dWeight + g4;  
                dTmp2 = (g2-g1)*dWeight + g1; 
            } 
            ////////////////////situation 4////////////////////// 
            //|---------------------------------------->X////////
            //|//////                 g1            /////////////  
            //|//////       g3    C   g2            /////////////  
            //|//////       g4                      /////////////  
            //Y////////////////////////////////////////////////// 
            else if( (theta[dataIndex] > 135 && theta[dataIndex] <= 180) ||
                     (theta[dataIndex] >= -45 && theta[dataIndex] < 0) ) { 
                g1 = grid[0][2];
                g2 = grid[1][2];
                g3 = grid[1][0];
                g4 = grid[2][0];

                dWeight = abs(Gx[dataIndex] / G[dataIndex]);
                dTmp1 = (g4-g3)*dWeight + g3;  
                dTmp2 = (g1-g2)*dWeight + g2;
            }

            if( (G[dataIndex] >= dTmp1) && (G[dataIndex] >= dTmp2) ) {
                resImageMatrix[dataIndex] = G[dataIndex];  //potential edge point
            } else {
                resImageMatrix[dataIndex] = 0;  //non edge point
            }
        }
    }
}
/*****************************************************************************/
/** Method Name: SetDoubleThreshold                                         **/
/** Method Input:                                                           **/
/**    (1)double* G: The pixels' gradients matrix.                          **/
/**    (2)uint32_t height: The input matrix's height.                       **/
/**    (3)uint32_t width: The input matrix's width.                         **/
/** Method Output:                                                          **/
/**    (1)this->highThreshold:The object's inner data member--highThreshold **/
/**       will be initialized, according to the input matrix G.             **/
/**    (2)this->lowThreshold:The object's inner data member--lowThreshold   **/
/**       will be initialized, according to the input matrix G.             **/
/** Method Description:                                                     **/
/**    This method aims to find optimal thresholds: low and high. It uses   **/
/** the Gaussian Mixture Model with 2 normal distributions to fit the       **/
/** Gradient Matrix without considering the 0 value. It will take the       **/
/** distribution with the higher mean value from the 2 distributions and the**/
/** high threshold will be generated from this distribution. Then, the low  **/
/** threshould is set to the half of the high threshold.                    **/
/*****************************************************************************/
template <typename SRC_T>
void CannyEdgeDetector<SRC_T>::SetDoubleThreshold(double* G, uint32_t height, uint32_t width) {
    int select = 0;
    int dimension = 1;
    int cluster_num = 2;
    
    double* _G = new double[height*width];
    memset(_G, 0, height*width*sizeof(double));
    uint32_t size = 0;
    for (uint32_t i = 0; i < height*width; ++i) {
        if (G[i] != 0) {
            _G[size++] = G[i];
        }
    }

    GMM *gmm = new GMM(dimension,cluster_num);
    gmm->Train(_G, size); //Training GMM
    select = gmm->Mean(0)[0] > gmm->Mean(1)[0] ? 0:1;
/*************************************************************************/   
/** The paramemter 0.9 and 0.5 need to be adjust to get a better result.**/  
/** Unfortunely, for now, we don't know which parameters are best. In   **/
/** the future, for each image, we may take several paramemters and     **/
/** implement it on the each image, then from those results, we can     **/
/**  a best result.                                                     **/
/*************************************************************************/
    this->highThreshold = gmm->Mean(select)[0] - 0.9*sqrt(gmm->Variance(select)[0]);
    this->lowThreshold = this->highThreshold*0.5;
    delete gmm;
    delete[] _G;
}
/*****************************************************************************/
/** Method Name: DoubleThresholdFilter                                      **/
/** Method Input:                                                           **/
/**   (1)double* G: The pixels' gradients matrix.                           **/
/**   (2)uint32_t height: The input matrix's height.                        **/
/**   (3)uint32_t width: The input matrix's width.                          **/
/**   (4)list<PointPosition*>* strongEdgeList:This is a list used to contain**/
/**     the points' positions(i,j).These points need to meet the requirement**/
/**     that the pixel's gradient should be greater than the highThreshold. **/
/**   (5)list<PointPosition*>* weakEdgeList: This is a list used to contain **/
/**     the points' positions(i,j).These points need to meet the requirement**/
/**     that the pixel's gradient should be greater than the lowThreshold   **/
/**     and smaller than the highThreshold.                                 **/
/** Method output:                                                          **/
/**   (1)list<PointPosition*>* strongEdgeList:The points, which are in this **/
/**     list, are the edge points and their Gradient will be set to 1.      **/
/**   (2)list<PointPosition*>* weakEdgeList:The points, which are in this   **/
/**     list, are the potential edge points and their Gradient will be kept.**/
/**   (3)double* G: the input matrix also is the output matrix. In this     **/
/**      method, the pixel whose gradient is smaller than the low threshold **/
/**      will be set to 0. The pixel whose gradient is greater than the high**/
/**      threshold will be set to 1. In addition, the pixel whose value is  **/
/**      between the low threshold and high threshold will be kept as the   **/
/**      input of the edge track method.                                    **/
/** Method Description:                                                     **/
/**  This method is used to filter the points. If the pixel's gradient is   **/
/**  greater than the highThreshold, then it will be treated as a edge point**/
/**  If the pixel's gradient is smaller than the low threshold, then it will**/
/**  be treated as a non edge point. However, if a pixel's gradient is      **/
/**  between the low and high thresholds, it will be treated as the weak    **/
/**  edge, which means this point is a potential edge points.               **/
/*****************************************************************************/
template <typename SRC_T>
void CannyEdgeDetector<SRC_T>::DoubleThresholdFilter(double* G, uint32_t height, uint32_t width, list<PointPosition*>* strongEdgeList, list<PointPosition*>* weakEdgeList) {
    PointPosition* point = NULL;
    uint32_t rowOffset, dataIndex;
    for (uint32_t i = 0; i < height; ++i) {
        rowOffset = i*width;
        for (uint32_t j = 0; j < width; ++j) {
            dataIndex = rowOffset+j;
            if (G[dataIndex] >= this->highThreshold) {
                G[dataIndex] = 1;
                point = new PointPosition;
                point->x = i;
                point->y = j;
                strongEdgeList->push_back(point);
            }else if (G[dataIndex] < this->lowThreshold) {
                G[dataIndex] = 0;   
            }else {
                point = new PointPosition;
                point->x = i;
                point->y = j;
                weakEdgeList->push_back(point);
            }
        }   
    }
}
/*****************************************************************************/
/** Method Name: EdgeTrack                                                  **/
/** Method Input:                                                           **/
/**   (1)double* G: The pixels' gradients matrix.                           **/
/**   (2)uint32_t height: The input matrix's height.                        **/
/**   (3)uint32_t width: The input matrix's width.                          **/
/**   (4)list<PointPosition*>* strongEdgeList:This is a list used to contain**/
/**     the points' positions(i,j).These points are all the edge points     **/
/**     detected in the DoubleThresholdFilter method.                       **/
/** Method output:                                                          **/
/**   (1)double* G: The pixels's gradients matrix G. In this matrix, the    **/
/** edge points are set to 1 and non edge points are set to 0. In addition, **/
/** some points still in the weak list and their gradient are neither 0 nor **/
/** 1. In this situation, these points are also be the non edge points and  **/
/** they will be cleaned up in the CleanUp method.                          **/
/**                                                                         **/
/** Method Description:                                                     **/
/**   For each point in the strongEdgeList, its neighbour will be checked   **/
/** whether they are greater than the lowThreshold. If its neighbour is     **/
/** greater than that value, then it will be treated as edge points;        **/
/** otherwise, it will be treated as non edge ponts.                        **/
/*****************************************************************************/
template <typename SRC_T>
void CannyEdgeDetector<SRC_T>::EdgeTrack(double* G, uint32_t height, uint32_t width, list<PointPosition*>* strongEdgeList) {
    PointPosition* point = NULL;
    while(strongEdgeList->size() > 0) {
        point = strongEdgeList->front();
        this->FindConnectedWeakEdges(G, point->x, point->y, height, width);
        delete point;
        strongEdgeList->pop_front();
    }
}
/*****************************************************************************/
/** Method Name: FindConnectedWeakEdges                                     **/
/** Method Input:                                                           **/
/**   (1)double* G: The pixels' gradients matrix.                           **/
/**   (2)uint32_t x: the input point's i position.                          **/
/**   (3)uint32_t y: the input point's j position.                          **/
/**   (4)uint32_t height: The input matrix's height.                        **/
/**   (5)uint32_t width: The input matrix's width.                          **/
/** Method output:                                                          **/
/**   (1)double* G: The pixels's gradients matrix G. In this matrix, the    **/
/** edge points are set to 1 and non edge points are set to 0. In addition, **/
/** some points still in the weak list and their gradient are neither 0 nor **/
/** 1. In this situation, these points are also be the non edge points and  **/
/** they will be cleaned up in the CleanUp method.                          **/
/**                                                                         **/
/** Method Description:                                                     **/
/**   This method is to check whether one edge point's neighbours are       **/
/** greater than the lowerThreshold. If its neighbour is greater than       **/
/** lowThreshold, then this neighbour point is set to 1 as a edge point;    **/
/** otherwise, it will be kept its value,which means this neighbour point is**/
/** a non edge point and it will be set to 0 in the CleanUp method.         **/
/*****************************************************************************/
template <typename SRC_T>
void CannyEdgeDetector<SRC_T>::FindConnectedWeakEdges(double* G, uint32_t x, uint32_t y, uint32_t height, uint32_t width) {
    const int LEFT_OFFSET = -3;
    const int RIGHT_OFFSET = 3;
    const int UP_OFFSET = -3;
    const int DOWN_OFFSET = 3;

    int64_t _x = x;
    int64_t _y = y;
    int64_t rowOffset, dataIndex;
    for (int i = UP_OFFSET; i < DOWN_OFFSET; ++i) {
        rowOffset = (_x+i)*width;
        for (int j = LEFT_OFFSET; j < RIGHT_OFFSET; ++j) {
            dataIndex = rowOffset + _y+j;
            if (_x+i >= 0 && _y+j >= 0 && _x+i < height && _y+j < width) {
                if (G[rowOffset + _y+j] > this->lowThreshold) {
                    G[rowOffset + _y+j] = 1;
                    FindConnectedWeakEdges(G, _x+i, _y+j, height, width);
                }
            }
        }
    }
}
/*****************************************************************************/
/** Method Name: CleanUp                                                    **/
/** Method Input:                                                           **/
/**   (1)double* G: The pixels' gradients matrix.                           **/
/**   (2)uint32_t height: The input matrix's height.                        **/
/**   (3)uint32_t width: The input matrix's width.                          **/
/**   (4)list<PointPosition*>* weakEdgeList: This is a list used to contain **/
/**     the points' positions(i,j).These points need to meet the requirement**/
/**     that the pixel's gradient should be greater than the lowThreshold   **/
/**     and smaller than the highThreshold.                                 **/
/** Method Output:                                                          **/
/**   (1)double* G: The pixels's gradients matrix G. In this matrix, the    **/
/** edge points are set to 1 and non edge points are set to 0. In addition, **/
/** some points still in the weak list and their gradient are neither 0 nor **/
/** 1. In this situation, these points are also be the non edge points and  **/
/** they will be cleaned up in the CleanUp method.                          **/
/**                                                                         **/
/** Method Description:                                                     **/
/**     This method is used after the edge track procedure. This is used to **/
/** clean up the weak edge points whose gradient values are greater than the**/
/** lowThreshold and smaller than the highThreshold. After edge track step, **/
/** some weak edge points are set to edge points. However, there may be     **/
/** still some points in weak edge list whose gradient is neither 1 (edge   **/
/** point) nor 0 (non edge points). These points is actually the non edge   **/
/** points since its gradient value is smaller than the highThreshold and   **/
/** all its 8 neighbour are non edge points.                                **/
/*****************************************************************************/
template <typename SRC_T>
void CannyEdgeDetector<SRC_T>::CleanUp(double* G, uint32_t height, uint32_t width, list<PointPosition*>* weakEdgeList) {
    PointPosition* point = NULL;
    while(weakEdgeList->size() > 0) {
        point = weakEdgeList->front();
        if (G[point->x*width + point->y] != 1) {
            G[point->x*width + point->y] = 0;
        }
        delete point;
        weakEdgeList->pop_front();
    }
}

#endif
