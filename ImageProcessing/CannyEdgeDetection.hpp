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

#ifndef _CANNY_EDGE_DETECTION_H
#define _CANNY_EDGE_DETECTION_H

#include <list>
#include <stdint.h>
#include "SobelFilter.hpp"
#include "GaussianSmoothing.hpp"

typedef struct PointPosition {
	uint32_t x;
	uint32_t y;
}PointPosition; //used to record the point position(i,j)
/*****************************************************************************/
/** Class Descriptyion:                                                     **/
/**      This class is used to implement Canny Edge Detection algorithm.    **/
/**                                                                         **/
/** -Data Description:                                                      **/
/**      In this class,it is based on the GaussianBlur class and SobelFilter**/
/**      class. Therefore, the data field contains an object of GaussianBlur**/
/**      class and SobelFilter class. Also, Canny Edge Detection relies on  **/
/**      double thresholds to detect the edge points; therefore, the data   **/
/**      field also contains the lowThreshold and highThreshold.            **/ 
/**                                                                         **/
/** -Operation Description:                                                 **/
/**      Canny Edge Detection mainly contains 5 steps: Gaussian Blur, Sobel **/
/**      Filter, Non Maximum Supression, Edge Tracking and Edge Cleaning.   **/
/*****************************************************************************/	
template <typename SRC_T>
class CannyEdgeDetector{
	private:
		double lowThreshold;
		double highThreshold;

		GaussianBlur<SRC_T>* gaussianBlur;
		SobelFilter<double>* sobelFilter;
		
	public:
		CannyEdgeDetector();
		~CannyEdgeDetector();
		void DetectEdge(SRC_T* srcImageMatrix, char* resImageMatrix, uint32_t height,uint32_t width);
		
	private:
		void SetDoubleThreshold(double* G, uint32_t height, uint32_t width);
		void CleanUp(double* G, uint32_t height, uint32_t width, list<PointPosition*>* weakEdgeList);
		void FindConnectedWeakEdges(double* G, uint32_t x, uint32_t y, uint32_t height, uint32_t width);
		void EdgeTrack(double* G, uint32_t height, uint32_t width, list<PointPosition*>* strongEdgeList);
		void NonMaximumSuppression(double* resImageMatrix, double* theta, double* Gx, double* Gy, double* G, uint32_t height, uint32_t width);
		void DoubleThresholdFilter(double* G, uint32_t height, uint32_t width,list<PointPosition*>* strongEdgeList, list<PointPosition*>* weakEdgeList);
		void KernelConvolution(int direction, double* srcImageMatrix, double* resImageMatrix, uint32_t height, uint32_t width);
};

#include "CannyEdgeDetection.cpp"

#endif