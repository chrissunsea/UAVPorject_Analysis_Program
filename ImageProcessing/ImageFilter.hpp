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
#ifndef _IMAGEFILTER_H
#define _IMAGEFILTER_H

#include <stdint.h>
/*****************************************************************************/
/** Class Descriptyion:                                                     **/
/**      This class is used to implement convolution operation.             **/
/**                                                                         **/
/** -Data Description:                                                      **/
/**      In this class, data filed contains a kernel and the kernel's       **/
/** attributes(height and width).Since the types of the kernel's element and**/
/** the input matrix may vary,this class is designed as a template.Therefore**/
/** ,it can support different type of kernel and the input matrix(int or    **/
/** double, etc).                                                           **/ 
/**                                                                         **/
/** -Operation Description:                                                 **/
/**      Essentially, the convolution is multiplication operation related to**/
/** 2 matrix: input matrix and kernel matrix. Therefore, the class operation**/
/** contains setting kernel matrix and calculating the convolution.         **/
/*****************************************************************************/	
template <typename KERNEL_T, typename SRC_T>
class ImageFilter {
/*****************************************************************************/
/** Operations:                                                             **/
/**     The operations contains setting the kernel and convolution besides  **/
/** the constructor and deconstructor.                                      **/
/*****************************************************************************/	
	public:
		ImageFilter();
		~ImageFilter();
		ImageFilter(const KERNEL_T* kernel, const int kernelHeight, const int kernelWidth);
		void ResetImageFilter(const KERNEL_T* kernel, const int kernelHeight, const int kernelWidth);
		void Convolution(SRC_T* srcImageMatrix, double*  resImageMatrix, uint32_t height, uint32_t width);
/*****************************************************************************/
/**Data Field:                                                              **/
/**     The class is a template.The inner data of one object contains kernel**/
/** matrix, which is a one dimension array that can simulate the 2 dimension**/
/** array.In addition, the type of the kernel is designed as a template.    **/
/**      It also contains 2 integer attributes that describe the height and **/
/** width of the kernel.                                                    **/
/*****************************************************************************/	
	private:
		KERNEL_T* kernel;
		int  kernelHeight;
		int  kernelWidth;
};

#include "ImageFilter.cpp"

#endif