/*****************************************************************************/
/**Author: Gaojie Sun                                                       **/
/**Date: 21/07/2017                                                         **/
/**Supervisor: Dongryeol Ryu                                                **/
/**Originazation: The University of Melbourne                               **/
/*****************************************************************************/

/*****************************************************************************/
/**Documentation Description:                                               **/
/**   This file is used to define the class that abstracts one of the basic **/
/** operations of the image processing--Gaussian Blur or Gaussian Smoothing.**/
/** Most of the other image processing algorithms are based on this,        **/
/** such as canny edge detection.                                           **/
/*****************************************************************************/
#ifndef _GAUSSIANSMOTHING_H
#define _GAUSSIANSMOTHING_H

#include <stdint.h>
/*****************************************************************************/
/** Class Descriptyion:                                                     **/
/**      This class is used to implement Gaussian blur operation based on   **/
/**   the Class ImageFilter(Convolution).                                   **/
/**                                                                         **/
/** -Data Description:                                                      **/
/**      In this class,data filed contains a kernel, the kernel's attributes**/
/** (kernelSize) and the sigma used to generated the element of the gaussian**/
/** .Since the type of the input matrix may vary,this class is designed as a**/
/** template.Therefore,it can support different type of input matrix(int or **/
/** double, etc). In addition, the gaussian kernel type is always be double.**/ 
/**                                                                         **/
/** -Operation Description:                                                 **/
/**     Essentially,the Gaussian Blur is multiplication operation related to**/
/** 2 matrix: input matrix and kernel matrix. Therefore, the class operation**/
/** contains setting kernel matrix and calculating the convolution.         **/
/*****************************************************************************/	
template <typename SRC_T>
class GaussianBlur {
/*****************************************************************************/
/**Data Field:                                                              **/
/**     The class is a template.The inner data of one object contains kernel**/
/** matrix, which is a one dimension array that can simulate the 2 dimension**/
/** array.In addition, the type of the kernel is designed as a template.    **/
/**      It also contains 1 integer attribute to describe the height and    **/
/** width of the kernel since Gaussian kernel is always a square.           **/
/**      A double type data, sigma, used to set the gassian distribution    **/
/** formula to generate the element of the kernel.                          **/
/*****************************************************************************/	
	private:
		int kernelSize;
		double  sigma;
		double* kernel;
/*****************************************************************************/
/**Data Field:                                                              **/
/**     The class is a template.The inner data of one object contains kernel**/
/** matrix, which is a one dimension array that can simulate the 2 dimension**/
/** array.In addition, the type of the kernel is designed as a template.    **/
/**      It also contains 1 integer attribute to describe the height and    **/
/** width of the kernel since Gaussian kernel is always a square.           **/
/**      A double type data, sigma, used to set the gassian distribution    **/
/** formula to generate the element of the kernel.                          **/
/*****************************************************************************/	
	private:
		void SetGaussianKernel();
	public:
		GaussianBlur(const int size, const double sigma);
		~GaussianBlur();
		void Blur(SRC_T* srcImageMatrix, double* resImageMatrix, uint32_t height, uint32_t width);
};
#include "GaussianSmoothing.cpp"

#endif