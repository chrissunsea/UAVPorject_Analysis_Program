/*****************************************************************************/
/**Author: Gaojie Sun                                                       **/
/**Date: 21/07/2017                                                         **/
/**Supervisor: Dongryeol Ryu                                                **/
/**Originazation: The University of Melbourne                               **/
/*****************************************************************************/

/*****************************************************************************/
/**Documentation Description:                                               **/
/**   This file is used to define the class that abstracts one of the basic **/
/** operations of the image processing--convolution. Most of the other image**/
/** processing algorithms are based on this operation.                      **/
/*****************************************************************************/
#ifndef _IMAGEFILTER_CPP
#define _IMAGEFILTER_CPP
#include <stdlib.h>
#include <iostream>
#include <cstring>
#include "ImageFilter.hpp"
using namespace std;
/*****************************************************************************/
/** Method Name: ImageFilter (Constructor)                                  **/
/** Method Input: NULL                                                      **/
/** Method Output: Generate an object used to implement the convolution.    **/
/** Method Description:                                                     **/
/**       This constructor does not need any parameter.It sets the kernel to**/
/** NULL and sets the height and width to 0 by default.                     **/
/*****************************************************************************/
template <typename KERNEL_T, typename SRC_T> 
ImageFilter<KERNEL_T, SRC_T>::ImageFilter() {
	this->kernel= NULL;
	this->kernelWidth = 0;
	this->kernelHeight = 0;
}
/*****************************************************************************/
/** Method Name: ImageFilter (Constructor)                                  **/
/** Method Input:                                                           **/
/**    (1)const KERNEL_T* kernel: The kernel array from outside is used to  **/
/**       initialised the object's inner kernel array.                      **/
/**    (2)const int kernelHeight: The kernel's height.                      **/
/**    (3)const int kernelWidth: The kernel's width.                        **/
/** Method Output: Generate an object used to implement the convolution.    **/
/** Method Description:                                                     **/
/**       This constructor takes the kernel array from outside as an input. **/
/**  It uses these data to initialise the object's inner data.              **/
/**                                                                         **/
/**       Importantly, this class only support the kernel's height and width**/
/**  are all odd.                                                           **/
/*****************************************************************************/
template <typename KERNEL_T, typename SRC_T> 
ImageFilter<KERNEL_T, SRC_T>::ImageFilter(const KERNEL_T* kernel, const int kernelHeight, const int kernelWidth) {
	
	if (kernelHeight <= 0 || kernelWidth <= 0) {
		cout<<"Error(ImageFilter): kernel's width and height should be greater than 0";
		exit(-1);		
	}
	if (kernelHeight%2 == 0 || kernelWidth%2 == 0) {
		cout<<"Error(ImageFilter): kernel's width and height should be odd";
		exit(-1);
	}

	this->kernelWidth  = kernelWidth;
	this->kernelHeight = kernelHeight;
	this->kernel = new KERNEL_T[kernelHeight*kernelWidth];
	memcpy(this->kernel, kernel, kernelHeight * kernelWidth * sizeof(KERNEL_T));
}
/*****************************************************************************/
/** Method Name: ~ImageFilter (Deconstructor)                               **/
/** Method Input: NULL                                                      **/
/** Method Output:Destory the object and claim back the memory allocated for**/
/**               this object.                                              **/
/** Method Description:                                                     **/
/**       The memory allocated for each ImageFileter object only includes   **/
/** the Kernel array. Therefore, when destory the object, only the kernel   **/
/** array need to be recycled.                                              **/
/*****************************************************************************/
template <typename KERNEL_T, typename SRC_T> 
ImageFilter<KERNEL_T, SRC_T>::~ImageFilter() {
	delete[] kernel;
	this->kernel =NULL;
}
/*****************************************************************************/
/** Method Name: ResetImageFilter                                           **/
/** Method Input:                                                           **/
/**    (1)const KERNEL_T* kernel: The kernel array from outside is used to  **/
/**       initialised the object's inner kernel array.                      **/
/**    (2)const int kernelHeight: The kernel's height.                      **/
/**    (3)const int kernelWidth: The kernel's width.                        **/
/** Method Output: Reset an object's kernel.                                **/
/** Method Description:                                                     **/
/**       This constructor takes the kernel array from outside as an input. **/
/**  It uses these data to reset the object's inner data.                   **/
/**                                                                         **/
/**       Importantly, this class only support the kernel's height and width**/
/**  are all odd.                                                           **/
/*****************************************************************************/
template <typename KERNEL_T, typename SRC_T> 
void ImageFilter<KERNEL_T, SRC_T>::ResetImageFilter(const KERNEL_T* kernel, const int kernelHeight, const int kernelWidth) {
	
	if (this->kernel != NULL) {
		delete[] this->kernel;
		this->kernel = NULL;
	}
	if (kernelHeight <= 0 || kernelWidth <= 0) {
		cout<<"Error(ResetImageFilter): kernel's width and height should be greater than 0";
		exit(-1);		
	}
	if (kernelHeight%2 == 0 || kernelWidth%2 == 0) {
		cout<<"Error(ResetImageFilter): kernel's width and height should be odd";
		exit(-1);
	}

	this->kernelWidth  = kernelWidth;
	this->kernelHeight = kernelHeight;
	this->kernel = new KERNEL_T[kernelHeight*kernelWidth];
	memcpy(this->kernel, kernel, kernelHeight * kernelWidth * sizeof(KERNEL_T));
}
/*****************************************************************************/
/** Method Name: Convolution                                                **/
/** Method Input:                                                           **/
/**    (1)SRC_T* srcImageMatrix: The input matrix. It can be any basic type.**/
/**       The input matrax is treated as a one dimension array.             **/
/**    (2)double* resImageMatrix: An empty matrix that has been allocated   **/
/**       memory used to return the result matrix. (This method is not      **/
/**       responsible for allocating the memory for the result).            **/
/**    (3)uint32_t height: The input matrix's height.                       **/
/**    (4)uint32_t width: The input matrix's width.                         **/
/** Method Output:                                                          **/
/**    (1)double* resImageMatrix: The result of the convolution.It always be**/
/**               the double type in order to handle the general situation. **/
/** Method Description:                                                     **/
/**       This method takes the input matrix and convolute it with the      **/
/**       object's inner kernel. The result is an array of double type.     **/
/**                                                                         **/
/**       Importantly, the method use the replicate approach to handle the  **/
/**  boundary.                                                              **/
/*****************************************************************************/
template <typename KERNEL_T, typename SRC_T> 
void ImageFilter<KERNEL_T, SRC_T>::Convolution(SRC_T* srcImageMatrix, double* resImageMatrix, uint32_t height, uint32_t width) {
	// Image's height and width are unsigned type.
	// Therefore, it doesn't need to check whether
	// the signed of them are less than 0.
	// It rely on the programmer to ensure the height and
	// width are legal input.
	if (height < this->kernelHeight || width < this->kernelWidth) {
		cout<<"Error(Convolution): image's width and height should be greater than kernel's";
		exit(-1);
	}

	memset(resImageMatrix, 0, width*height*sizeof(double));

	int offsetWidth  = this->kernelWidth / 2;
	int offsetHeight = this->kernelHeight / 2;
	double accumulator;

	uint32_t dataRowOffset = 0;
	uint32_t i, j, pixel_i, pixel_j;
	int kernel_i, kernel_j;
	for (i = offsetHeight; i < height -offsetHeight; ++i) {
		for (j = offsetWidth; j < width - offsetWidth; ++j) {
			accumulator = 0.0;
			for (kernel_i = 0, pixel_i = i - offsetHeight; kernel_i < this->kernelHeight; ++kernel_i, ++pixel_i) { 
				for (kernel_j = 0, pixel_j = j - offsetWidth; kernel_j < this->kernelWidth; ++kernel_j, ++pixel_j) {
					accumulator += this->kernel[kernel_i*this->kernelWidth + kernel_j] * srcImageMatrix[pixel_i*width + pixel_j];
				}
			}
			resImageMatrix[i*width + j] = accumulator;
		}
	}
/*****************************************************************************/
/**                       handle the upper boundary                         **/
/*****************************************************************************/
	int64_t signed_i, signed_j;
	int64_t signed_pixel_i,signed_pixel_j;
	int64_t _i,_j;
	for (signed_i = 0; signed_i < offsetHeight; ++signed_i) {
		for (signed_j = offsetWidth; signed_j < width-offsetWidth; ++signed_j) {
			accumulator = 0.0;
			for (kernel_i = 0, signed_pixel_i = signed_i - offsetHeight; kernel_i < this->kernelHeight; ++kernel_i, ++signed_pixel_i) { 
				for (kernel_j = 0, signed_pixel_j = signed_j - offsetWidth; kernel_j < this->kernelWidth; ++kernel_j, ++signed_pixel_j) {
					_i = 0 > signed_pixel_i? 0:signed_pixel_i;
					accumulator += this->kernel[kernel_i*this->kernelWidth + kernel_j] * srcImageMatrix[_i*width + signed_pixel_j];
				}
			}
			resImageMatrix[signed_i*width + signed_j] = accumulator;
		}
	}
/*****************************************************************************/
/**                       handle the lower boundary                         **/
/*****************************************************************************/
	for (signed_i = height -offsetHeight; signed_i < height; ++signed_i) {
		for (signed_j = offsetWidth; signed_j < width-offsetWidth; ++signed_j) {
			accumulator = 0.0;
			for (kernel_i = 0, signed_pixel_i = signed_i - offsetHeight; kernel_i < this->kernelHeight; ++kernel_i, ++signed_pixel_i) { 
				for (kernel_j = 0, signed_pixel_j = signed_j - offsetWidth; kernel_j < this->kernelWidth; ++kernel_j, ++signed_pixel_j) {
					_i = signed_pixel_i>height-1 ? height-1:signed_pixel_i;
					accumulator += this->kernel[kernel_i*this->kernelWidth + kernel_j] * srcImageMatrix[_i*width + signed_pixel_j];
				}
			}
			resImageMatrix[signed_i*width + signed_j] = accumulator;
		}
	}
/*****************************************************************************/
/**                       handle the left boundary                          **/
/*****************************************************************************/
	for (signed_i = 0; signed_i < height; ++signed_i) {
		for (signed_j = 0; signed_j < offsetWidth; ++signed_j) {
			accumulator = 0.0;
			for (kernel_i = 0, signed_pixel_i = signed_i - offsetHeight; kernel_i < this->kernelHeight; ++kernel_i, ++signed_pixel_i) { 
				for (kernel_j = 0, signed_pixel_j = signed_j - offsetWidth; kernel_j < this->kernelWidth; ++kernel_j, ++signed_pixel_j) {
					if (signed_pixel_i < 0) {
						_i = 0;
					}else if (signed_pixel_i > height-1) {
						_i = height-1;
					}else{
						_i = signed_pixel_i;
					}

					if (signed_pixel_j < 0) {
						_j = 0;
					}else if (signed_pixel_j > width-1) {
						_j = width-1;
					}else{
						_j = signed_pixel_j;
					}
					accumulator += this->kernel[kernel_i*this->kernelWidth + kernel_j] * srcImageMatrix[_i*width + _j];
				}
			}
			resImageMatrix[signed_i*width + signed_j] = accumulator;
		}
	}
/*****************************************************************************/
/**                       handle the right boundary                         **/
/*****************************************************************************/
	for (signed_i = 0; signed_i < height; ++signed_i) {
		for (signed_j = width - offsetWidth; signed_j < width; ++signed_j) {
			accumulator = 0.0;
			for (kernel_i = 0, signed_pixel_i = signed_i - offsetHeight; kernel_i < this->kernelHeight; ++kernel_i, ++signed_pixel_i) { 
				for (kernel_j = 0, signed_pixel_j = signed_j - offsetWidth; kernel_j < this->kernelWidth; ++kernel_j, ++signed_pixel_j) {
					
					if (signed_pixel_i < 0) {
						_i = 0;
					}else if (signed_pixel_i > height-1) {
						_i = height-1;
					}else{
						_i = signed_pixel_i;
					}

					if (signed_pixel_j < 0) {
						_j = 0;
					}else if (signed_pixel_j > width-1) {
						_j = width-1;
					}else{
						_j = signed_pixel_j;
					}
					accumulator += this->kernel[kernel_i*this->kernelWidth + kernel_j] * srcImageMatrix[_i*width + _j];
				}
			}
			resImageMatrix[signed_i*width + signed_j] = accumulator;
		}
	}

}

#endif
