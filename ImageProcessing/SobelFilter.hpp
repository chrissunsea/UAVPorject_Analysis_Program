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
#ifndef _SOBELFILTER_H
#define _SOBELFILTER_H

#include <stdint.h>
#include "ImageFilter.hpp"
/*****************************************************************************/
/** Class Descriptyion:                                                     **/
/**      This class is used to implement Sobel image filter operation       **/
/**      including detect x, y, diagonal and any (mannually set) directions.**/
/**      It is based on the ImageFileter class.                             **/
/**                                                                         **/
/** -Data Description:                                                      **/
/**      The kernels including x, y, diagonal directions are encapsulated   **/
/**      into class by default. These are 4 1-dimension array. In addition, **/
/**      this class also support user-defined kernel to implement other     **/
/**      direction's convolution.                                           **/
/**                                                                         **/
/** -Operation Description:                                                 **/
/**      Essentially, the sobel image filter is a kind of convolution with  **/
/**      special kernels. This class supports to do the x, y and diagonal   **/
/**      direction's convolution and also supports to merge the two         **/
/**      perpendicular direction's convolution value.                       **/
/*****************************************************************************/	
template <typename SRC_T>
class SobelFilter{
	private:
		int  height;
		int  width;
		int* kernel; //mannually set the kernel
		ImageFilter<int, SRC_T>* filter; //implement the convolution

		static const int KERNEL_X[];
		static const int KERNEL_Y[];
		static const int KERNEL_45_DEGREE[];
		static const int KERNEL_135_DEGREE[];
		static const int KERNEL_HEIGHT_3 = 3;
		static const int KERNEL_WIDTH_3 = 3;
		
	public:
		enum Direction {
			DIRECTION_X,
			DIRECTION_Y,
			DIRECTION_45,
			DIRECTION_135,
			MANUALLY
		}; //used by Filter method to specify the filter direction

	public:
		SobelFilter();
		~SobelFilter();
		SobelFilter(int* kernel, int height, int width);
		void Filter(int direction, SRC_T* srcImageMatrix, double* resImageMatrix, uint32_t height, uint32_t width);
		void CombineGradient(double* Gx, double* Gy, double* G, uint32_t height, uint32_t width, double* theta = NULL);
};

#include "SobelFilter.cpp"

#endif