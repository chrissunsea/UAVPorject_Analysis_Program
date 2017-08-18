/*****************************************************************************/
/**Author: Gaojie Sun                                                       **/
/**Date: 21/07/2017                                                         **/
/**Supervisor: Dongryeol Ryu                                                **/
/**Originazation: The University of Melbourne                               **/
/*****************************************************************************/

/*****************************************************************************/
/**Documentation Description:                                               **/
/**   This file is used to define the class that abstracts one of the       **/
/** morphological operations--dilation. This is mainly used to enlarge the  **/
/** edge (edge point needs to be set to 1).                                 **/
/*****************************************************************************/
#ifndef _MORPHOLOGYIMAGEPROCESSOR_H
#define _MORPHOLOGYIMAGEPROCESSOR_H

#include <stdint.h>
/*****************************************************************************/
/** Class Descriptyion:                                                     **/
/**      This class is used to implement dilation operation.                **/
/**                                                                         **/
/** -Data Description:                                                      **/
/**      In this class, data filed contains a kernel and the kernel's       **/
/** attributes(height, width (r) and shape).Since the types of the input    **/
/** matrix that needs to implement the dilation may vary, this class is     **/
/** designed as a template. Therefore,it can support different type of the  **/
/** input matrix (int or double, etc).                                      **/ 
/**                                                                         **/
/** -Operation Description:                                                 **/
/**      The main method for this class is to set the dilation kernel(shape **/
/** and size) and to implement the dilation operation on the input matrix.  **/
/*****************************************************************************/	
template <typename SRC_T>
class MorphologyImageProcessor {
	public:
		MorphologyImageProcessor(int shape, int size);
		MorphologyImageProcessor(int shape, int height, int width);
		~MorphologyImageProcessor();
		void Dilation(SRC_T* srcImageMatrix, SRC_T* resImageMatrix, uint32_t height, uint32_t width);

	public:
		enum Shape{
			SQUARE,
			CIRCLE
		};

	private:
		void InitializeKernel();
		void InitializeKernelSquare();
		void InitializeKernelCircle();

	private:
		int* kernel;
		int height;
		int width;
		int shape;
		int r;
	
};
#include "MorphologyImageProcessor.cpp"

#endif