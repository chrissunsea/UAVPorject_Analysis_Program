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
#ifndef _GAUSSIANSMOTHING_CPP
#define _GAUSSIANSMOTHING_CPP

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "ImageFilter.hpp"
#include "GaussianSmoothing.hpp"
using namespace std;
/*****************************************************************************/
/** Method Name: GaussianBlur (Constructor)                                 **/
/** Method Input:                                                           **/
/**    (1)const int size: The size of the Kernel(2D) of Gaussian Blur.      **/
/**    (2)const double sigma: The Gaussian formula's sigma used to generate **/
/**       the element of the Gaussian Blur kernel's element.                **/
/** Method Output: Generate an object containing the Gaussian Blur Kernel   **/
/**                 used to implement the convolution.                      **/
/** Method Description:                                                     **/
/**       This constructor takes the kernel's size that used to set the     **/
/**       Gaussian kernel's height and width and the sigma of the Gaussian  **/
/**       Distribution formula as the input.                                **/
/**                                                                         **/
/**       Importantly, this class only support the kernel's height and width**/
/**  are all odd.                                                           **/
/*****************************************************************************/
template <typename SRC_T>
GaussianBlur<SRC_T>::GaussianBlur(const int size, const double sigma) {
	if (size <=0 || sigma <= 0) {
		cout<<"Error(GaussianBlur): The kernel size and sigma should be greater than 0"<<endl;		
		exit(-1);
	}
	if (size%2 == 0) {
		cout<<"Error(GaussianBlur): The kernel size should be odd"<<endl;		
		exit(-1);
	}

	this->sigma = sigma;
	this->kernelSize = size;
	this->kernel = new double[this->kernelSize*this->kernelSize];

	this->SetGaussianKernel();
}
/*****************************************************************************/
/** Method Name: ~GaussianBlur (Deconstructor)                              **/
/** Method Input: NULL                                                      **/
/** Method Output:Destory the object and claim back the memory allocated for**/
/**               this object.                                              **/
/** Method Description:                                                     **/
/**       The memory allocated for each GaussianBlur object only includes   **/
/** the Kernel array. Therefore, when destory the object, only the kernel   **/
/** array need to be recycled.                                              **/
/*****************************************************************************/
template <typename SRC_T>
GaussianBlur<SRC_T>::~GaussianBlur() {
	delete[] this->kernel;
	this->kernel = NULL;
}
/*****************************************************************************/
/** Method Name: SetGaussianKernel                                          **/
/** Method Input: NULL                                                      **/
/** Method Output: This method is used to generate the gaussian kernel's    **/
/**                element according to the kernel's size, the sigma and the**/
/**                Gaussian distribution's formula.                         **/
/** Method Description:                                                     **/
/**     According to the position (euclidean metric) of the kernel, the     **/
/**     kernel's element will initialised using Gaussian distribution       **/
/**     formula.                                                            **/
/**                                                                         **/
/**     It is called by the constructor and is a private method cannot be   **/
/**     called by outside.                                                  **/
/*****************************************************************************/
template <typename SRC_T>
void GaussianBlur<SRC_T>::SetGaussianKernel() { 

	if (this->kernelSize <=0 || this->sigma <= 0) {
		cout<<"Error(GaussianBlur): The kernel size and sigma should be greater than 0"<<endl;		
		exit(-1);
	}
	if (this->kernelSize%2 == 0) {
		cout<<"Error(GaussianBlur): The kernel size should be odd"<<endl;		
		exit(-1);
	}

    const double PI=3.141592653589793; 
    double sum=0;
    int rowOffset = 0;
    int center=this->kernelSize/2;

    for(int i = 0; i < this->kernelSize; i++)  {  
        rowOffset = i*this->kernelSize;
        for(int j = 0; j < this->kernelSize; j++) {
            this->kernel[rowOffset + j] = (1/(2*PI*sigma*sigma))*exp(-((i-center)*(i-center)+(j-center)*(j-center))/(2*sigma*sigma));  
            sum += this->kernel[rowOffset + j];
        }  
    }  

    for(int i = 0; i < this->kernelSize; i++) {  
        rowOffset = i*this->kernelSize;
        for(int j = 0; j < this->kernelSize; j++) {  
            this->kernel[rowOffset + j] /= sum;  
        }  
    } 
}
/*****************************************************************************/
/** Method Name: Blur                                                       **/
/** Method Input:                                                           **/
/**     (1)SRC_T* srcImageMatrix: the image matrix that need to be blurred. **/
/**     (2)double* resImageMatrix: An empty matrix that has been allocated  **/
/**       memory used to return the result matrix. (This method is not      **/
/**       responsible for allocating the memory for the result).            **/
/**    (3)uint32_t height: The input matrix's width.                        **/
/**    (4)uint32_t width: The input matrix's width.                         **/
/** Method Output:                                                          **/
/**    (1)double* resImageMatrix: The result of the blur.It always be the   **/
/**               double type in order to handle the general situation.     **/
/** Method Description:                                                     **/
/**    This method is based on an ImageFilter object that can implement     **/
/**    the convolution. Essentially, the gaussian blur is a kind of         **/
/**    convolutions with a special kernel. The mainly use of Gaussian blur  **/
/**    is to reduce the image noise.                                        **/
/*****************************************************************************/
template <typename SRC_T>
void GaussianBlur<SRC_T>::Blur(SRC_T* srcImageMatrix, double* resImageMatrix, uint32_t height, uint32_t width) {
	ImageFilter<double, SRC_T>* imgFilter = new ImageFilter<double, SRC_T>(this->kernel, this->kernelSize, this->kernelSize);
	imgFilter->Convolution(srcImageMatrix, resImageMatrix, height, width);
	delete imgFilter;
}

#endif
