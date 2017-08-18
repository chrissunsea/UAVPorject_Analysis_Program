/*****************************************************************************/
/**Author: Gaojie Sun                                                       **/
/**Date: 21/07/2017                                                         **/
/**Supervisor: Dongryeol Ryu                                                **/
/**Originazation: The University of Melbourne                               **/
/*****************************************************************************/

/*****************************************************************************/
/**Documentation Description:                                               **/
/**   This file is used to define the class that abstracts the Sobel Image  **/
/** Fileter.Based on this fileter object, Sobel Edge Detection can be easily**/
/** implemented and Canny edge detection also can use this object. Therefore**/
/** , it is worh to encapsulate this kind of filter into a individual class.**/
/*****************************************************************************/
#ifndef _SOBELFILTER_CPP
#define _SOBELFILTER_CPP

#include <stdlib.h>
#include <cmath>
#include <cstring>
#include <iostream>
#include "SobelFilter.hpp"
using namespace std;

//X direction kernel
template <typename SRC_T> 
const int SobelFilter<SRC_T>::KERNEL_X[] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
//Y direction kernel
template <typename SRC_T> 
const int SobelFilter<SRC_T>::KERNEL_Y[] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};
//45 degree direction kernel
template <typename SRC_T> 
const int SobelFilter<SRC_T>::KERNEL_45_DEGREE[]  = {0, 1, 2, -1, 0, 1, -2, -1, 0};
//135 degree direction kernel
template <typename SRC_T> 
const int SobelFilter<SRC_T>::KERNEL_135_DEGREE[] = {-2, -1, 0, -1, 0, 1, 0, 1, 2};
/*****************************************************************************/
/** Method Name: SobelFilter (Constructor)                                  **/
/** Method Input: NULL                                                      **/
/** Method Output: Generate an object used to implement the Sobel Filter.   **/
/** Method Description:                                                     **/
/**       This constructor does not need any parameter.It sets the kernel to**/
/** NULL and sets the height and width to 0 by default. Also, a new object  **/
/** of ImageFilter class is constructed to implement the convolution.       **/
/*****************************************************************************/
template <typename SRC_T> 
SobelFilter<SRC_T>::SobelFilter() {
	this->height = 0;
	this->width = 0;
	this->kernel = NULL;
	this->filter = new ImageFilter<int, SRC_T>();
}
/*****************************************************************************/
/** Method Name: ~SobelFilter (Deconstructor)                               **/
/** Method Input: NULL                                                      **/
/** Method Output:Destory the object and claim back the memory allocated for**/
/**               this object.                                              **/
/** Method Description:                                                     **/
/**       The memory allocated for each SobelFilter object only includes    **/
/** the Kernel array (used to mannually set the kernel) and the obkect of   **/
/** the ImageFilter class. Therefore, when destory the object, only the     **/
/** kernel array and the object of the class ImageFilter need to be recycled**/
/*****************************************************************************/
template <typename SRC_T> 
SobelFilter<SRC_T>::~SobelFilter() {
	delete[] this->kernel;
	delete this->filter;
	this->kernel = NULL;
	this->filter = NULL;

}
/*****************************************************************************/
/** Method Name: SobelFilter (Constructor)                                  **/
/** Method Input:                                                           **/
/**    (1)int* kernel: (Mannually set the kernel) The kernel array from     **/
/**          outside is used to initialised the object's inner kernel array.**/
/**    (2)const int kernelHeight: The kernel's height.                      **/
/**    (3)const int kernelWidth: The kernel's width.                        **/
/** Method Output: Generate an object used to implement the Sobel Filter.   **/
/** Method Description:                                                     **/
/**       This constructor takes the kernel array from outside as an input. **/
/**  It uses these data to initialize the object's inner data. Indeed, the  **/
/**  inner data is designed for the user to mannually set the kernel to do  **/
/**  the convolution.                                                       **/
/**                                                                         **/
/**       Importantly, this class only support the kernel's height and width**/
/**  are all odd.                                                           **/
/*****************************************************************************/
template <typename SRC_T> 
SobelFilter<SRC_T>::SobelFilter(int* kernel, int height, int width) {
	if (height <= 0 || width <= 0) {
		cout<<"Error(ImageFilter): kernel's width and height should be greater than 0"<<endl;
		exit(-1);		
	}

	if (height%2 == 0 || width%2 == 0) {
		cout<<"Error(ImageFilter): kernel's width and height should be odd"<<endl;
		exit(-1);
	}

	this->kernel = new int[height * width];
	this->height = height;
	this->width  = width;

	memcpy(this->kernel, kernel, height*width*sizeof(int));
	this->filter = new ImageFilter<int, SRC_T>();
}
/*****************************************************************************/
/** Method Name: Filter                                                     **/
/** Method Input:                                                           **/
/**    (1)int direction: The Sobel kernel's direction. Usually use the enum **/
/**       Direction's member as parameter.                                  **/
/**    (2)SRC_T* srcImageMatrix: The input matrix. It can be any basic type.**/
/**       The input matrax is treated as a one dimension array.             **/
/**    (3)double* resImageMatrix: An empty matrix that has been allocated   **/
/**       memory used to return the result matrix. (This method is not      **/
/**       responsible for allocating the memory for the result).            **/
/**    (4)uint32_t height: The input matrix's width.                        **/
/**    (5)uint32_t width: The input matrix's width.                         **/
/** Method Output:                                                          **/
/**    (1)double* resImageMatrix: The result of the convolution.It always be**/
/**               the double type in order to handle the general situation. **/
/** Method Description:                                                     **/
/**       This method takes the input matrix and convolute it with the      **/
/**       object's inner kernels. The result is an array of double type.    **/
/**                                                                         **/
/**       Importantly, the method uses the replicate approach to handle the **/
/**  boundary.                                                              **/
/*****************************************************************************/
template <typename SRC_T> 
void SobelFilter<SRC_T>::Filter(int direction, SRC_T* srcImageMatrix, double* resImageMatrix, uint32_t height, uint32_t width) {
	switch(direction){
		case DIRECTION_X:
			this->filter->ResetImageFilter(&SobelFilter::KERNEL_X[0], KERNEL_HEIGHT_3, KERNEL_WIDTH_3);break;
		case DIRECTION_Y:
			this->filter->ResetImageFilter(&SobelFilter::KERNEL_Y[0], KERNEL_HEIGHT_3, KERNEL_WIDTH_3);break;
		case DIRECTION_45:
			this->filter->ResetImageFilter(&SobelFilter::KERNEL_45_DEGREE[0], KERNEL_HEIGHT_3, KERNEL_WIDTH_3);break;
		case DIRECTION_135:
			this->filter->ResetImageFilter(&SobelFilter::KERNEL_135_DEGREE[0], KERNEL_HEIGHT_3, KERNEL_WIDTH_3);break;
		case MANUALLY:
			if (this->kernel != NULL) {
				this->filter->ResetImageFilter(this->kernel, this->height, this->width);break;
			}
		default:
			cout<<"Error(Sobel): Wrong direction detection."<<endl;exit(-1);
	}
	this->filter->Convolution(srcImageMatrix, resImageMatrix, height, width);
}
/*****************************************************************************/
/** Method Name: CombineGradient                                            **/
/** Method Input:                                                           **/
/**    (1)double* Gx: One of the 2 perpendicular direction's Gradients (the **/
/**                   result of the convolution of the Sobel Filter).       **/
/**    (2)double* Gy: The other one of the 2 perpendicular direction's      **/
/**                   Gradients (the result of the convolution of the Sobel **/
/**                   Filter).                                              **/
/**    (3)double* G: An empty matrix that has been allocated memory used    **/
/**                  to return the result matrix. (This method is not       **/
/**                  responsible for allocating the memory for the result). **/
/**    (4)uint32_t height: The input matrix's height.                       **/
/**    (5)uint32_t width: The input matrix's width.                         **/
/**    (6)double* theta: the (Gradient) degree of the tangent [-180,180].   **/
/**                      theta = atan(Gy/Gx).                               **/
/** Method Output:                                                          **/
/**    (1)double* G: The result of the combination of the 2 perpendicular   **/
/**                  direction's gradients. It always be the double type.   **/
/** Method Description:                                                     **/
/**       This method takes the 2 perpendicular direction's Gradients (the  **/
/**       results of the convolution of the Sobel Filter) and return the    **/
/**       combination of the 2 gradients and the (Gradient) degree of the   **/
/**       tangent [-180,180].  theta = atan(Gy/Gx).                         **/
/*****************************************************************************/
template <typename SRC_T>
void SobelFilter<SRC_T>::CombineGradient(double* Gx, double* Gy, double* G, uint32_t height, uint32_t width, double* theta) {

	const double PI = 3.141592654;
	uint32_t rowOffset, dataIndex;
	for (uint32_t i = 0; i < height; ++i) {
		rowOffset = i*width;
		for (uint32_t j = 0; j < width; ++j) {
			dataIndex = rowOffset + j;
			G[dataIndex] = sqrt(Gx[dataIndex]*Gx[dataIndex] + Gy[dataIndex]*Gy[dataIndex]);
		}
	}

	if (theta != NULL) {
		for (uint32_t i = 0; i < height; ++i) {
			rowOffset = i*width;
			for (uint32_t j = 0; j < width; ++j) {
				dataIndex = rowOffset + j;
				theta[dataIndex] = atan2(Gy[dataIndex], Gx[dataIndex]) * 180 / PI;
			}
		}		
	}
}

#endif


// int main() {
	
// 	uint16_t matrix[] = {1,2,3,4,5,
// 	                     6,7,8,9,10,
// 	                     11,12,13,14,15,
// 	                     16,17,18,19,20,
// 	                     21,22,23,24,25};
// 	int kernel[] = {-1,-1,-1,1,1,1,1,1,1};
// 	double* resImageMatrix = new double[25];
// 	SobelFilter<uint16_t> filter = SobelFilter<uint16_t>(kernel,3,3);
// 	filter.Filter(SobelFilter<uint16_t>::MANUALLY, matrix, resImageMatrix, 5, 5);

// 	for (int i = 0; i < 5; ++i) {
// 		for (int j = 0; j < 5; ++j) {
// 			cout<<resImageMatrix[i*5+j]<<"   ";
// 		}
// 		cout<<endl;
// 	}
// 	cout<<endl;
// 	filter.Filter(SobelFilter<uint16_t>::DIRECTION_Y, matrix, resImageMatrix,5,5);

// 	for (int i = 0; i < 5; ++i) {
// 		for (int j = 0; j < 5; ++j) {
// 			cout<<resImageMatrix[i*5+j]<<"   ";
// 		}
// 		cout<<endl;
// 	}
// }
