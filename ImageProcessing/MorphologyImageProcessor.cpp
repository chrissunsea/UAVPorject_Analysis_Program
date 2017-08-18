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
#ifndef _MORPHOLOGYIMAGEPROCESSOR_CPP
#define _MORPHOLOGYIMAGEPROCESSOR_CPP

#include <stdlib.h>
#include <cmath>
#include <string>
#include <iostream>
#include "ImageFilter.hpp"
#include "MorphologyImageProcessor.hpp"
using namespace std;
/*****************************************************************************/
/** Method Name: MorphologyImageProcessor (Constructor)                     **/
/** Method Input:                                                           **/
/**        (1)int shape: this is used to specify the shape of the kernel.   **/
/**           However, for this constructor, it is used to set the shape of **/
/**           circle kernel.                                                **/
/**        (2)int size: this is used to specify the size of the kernel.     **/
/**           For this constructor(circle shape), it is used to set the     **/
/**           radius.                                                       **/
/** Method Output: Generate an object used to implement the dilation.       **/
/*****************************************************************************/
template <typename SRC_T>
MorphologyImageProcessor<SRC_T>::MorphologyImageProcessor(int shape, int size) {
	
	if (shape != MorphologyImageProcessor::CIRCLE) {
		cout<<"Error(MorphologyImageProcessor): The arguments is insuficient!"<<endl;
		exit(-1);
	}
	if (size <= 0) {
		cout<<"Error(MorphologyImageProcessor): The size should be greater than 0!"<<endl;
		exit(-1);
	}

	this->shape = MorphologyImageProcessor::CIRCLE;
	this->height = size*2+1;
	this->width = size*2+1;
	this->r = 0;
	if (shape == MorphologyImageProcessor::CIRCLE) {
		this->r = size;
	}
	
	this->kernel = new int[size*size];
	this->InitializeKernel();
}
/*****************************************************************************/
/** Method Name: MorphologyImageProcessor (Constructor)                     **/
/** Method Input:                                                           **/
/**        (1)int shape: this is used to specify the shape of the kernel.   **/
/**           However, for this constructor, it is used to set the shape of **/
/**           square kernel.                                                **/
/**        (2)int height: this is used to specify the height of the square  **/
/**           kernel.                                                       **/
/**                                                                         **/
/**        (3)int width: this is used to specify the width of the square    **/
/**           kernel.                                                       **/
/** Method Output: Generate an object used to implement the dilation.       **/
/*****************************************************************************/
template <typename SRC_T>
MorphologyImageProcessor<SRC_T>::MorphologyImageProcessor(int shape, int height, int width) {

	if (shape != MorphologyImageProcessor::SQUARE) {
		cout<<"Error(MorphologyImageProcessor): Too many arguments!"<<endl;
		exit(-1);
	}
	if (height <= 0 || width <= 0) {
		cout<<"Error(MorphologyImageProcessor): The kernel's height and width should be greater than 0"<<endl;
		exit(-1);
	}

	this->shape = MorphologyImageProcessor::SQUARE;
	this->height = height;
	this->width = width;
	this->r = 0;
	this->kernel = new int[height*width];
	this->InitializeKernel();
}
/*****************************************************************************/
/** Method Name: ~MorphologyImageProcessor (Deconstructor)                  **/
/** Method Input: NULL                                                      **/
/** Method Output:Destory the object and claim back the memory allocated for**/
/**               this object.                                              **/
/** Method Description:                                                     **/
/**       The memory allocated for each ImageFileter object only includes   **/
/** the Kernel array. Therefore, when destory the object, only the kernel   **/
/** array need to be recycled.                                              **/
/*****************************************************************************/
template <typename SRC_T>
MorphologyImageProcessor<SRC_T>::~MorphologyImageProcessor() {
	delete[] this->kernel;
	this->kernel = NULL;
}
/*****************************************************************************/
/** Method Name: InitializeKernel                                           **/
/** Method Input: NULL                                                      **/
/** Method Output:Initialise the data member - kernel according to the shape**/
/**               of the kernel.                                            **/
/** Method Description:                                                     **/
/**       This method will call the coresponding method to initialise the   **/
/** the Kernel array, according to the shape.                               **/
/*****************************************************************************/
template <typename SRC_T>
void MorphologyImageProcessor<SRC_T>::InitializeKernel() {
	if (this->shape == MorphologyImageProcessor::SQUARE) {
		this->InitializeKernelSquare();
	}else if(this->shape == MorphologyImageProcessor::CIRCLE) {
		this->InitializeKernelCircle();
	}
}
/*****************************************************************************/
/** Method Name: InitializeKernelSquare                                     **/
/** Method Input: NULL                                                      **/
/** Method Output:Initialise the data member - kernel according to the      **/
/**               square shape.                                             **/
/** Method Description:                                                     **/
/**       This method will set the square kernel's elements to 1.            **/
/*****************************************************************************/
template <typename SRC_T>
void MorphologyImageProcessor<SRC_T>::InitializeKernelSquare() {
	for (int i = 0; i < this->height; ++i) {
		for (int j = 0; j < this->width; ++j) {
			this->kernel[i*this->width + j] = 1;
		}
	}
}
/*****************************************************************************/
/** Method Name: InitializeKernelCircle                                     **/
/** Method Input: NULL                                                      **/
/** Method Output:Initialise the data member - kernel according to the      **/
/**               circle shape.                                             **/
/** Method Description:                                                     **/
/**       This method will set some of the square kernel's elements to 1.   **/
/*****************************************************************************/
template <typename SRC_T>
void MorphologyImageProcessor<SRC_T>::InitializeKernelCircle() {
	memset(this->kernel, 0, this->height*this->width*sizeof(int));
	int center_i = this->height / 2;
	int center_j = this->width / 2;
	for (int i = 0; i < this->height; ++i) {
		for (int j = 0; j < this->width; ++j) {
			if (sqrt((i-center_i)*(i-center_i) + (j-center_j)*(j-center_j)) <= this->r) {
				this->kernel[i*this->width+j] = 1;
			}
		}
	}
}
/*****************************************************************************/
/** Method Name: Dilation                                                   **/
/** Method Input:                                                           **/
/**        (1)SRC_T* srcImageMatrix: the input image matrix that need to    **/
/**           implement the dilation operation(edge points are 1 and non    **/
/**           edge points are 0).                                           **/
/**        (2)SRC_T* resImageMatrix: An empty matrix that has been allocated**/
/**           memory used to return the result matrix. (This method is not  **/
/**           responsible for allocating the memory for the result).        **/
/**        (3)uint32_t height: The input matrix's height.                   **/
/**        (4)uint32_t width: The input matrix's width.                     **/
/** Method Output:                                                          **/
/**        (1))SRC_T* resImageMatrix: The result of the dilation.           **/
/**            Edge points are 1 and no edge points are 0.                  **/
/** Method Description:                                                     **/
/**       The output and the input can be the same matrix.                  **/
/*****************************************************************************/
template <typename SRC_T>
void MorphologyImageProcessor<SRC_T>::Dilation(SRC_T* srcImageMatrix, SRC_T* resImageMatrix, uint32_t height, uint32_t width) {
	ImageFilter<int, SRC_T> *imageFilter = new ImageFilter<int, SRC_T>(this->kernel, this->height, this->width);
	double* convRes = new double[height*width];
	imageFilter->Convolution(srcImageMatrix, convRes, height, width);

	for (uint32_t i = 0; i < height; ++i) {
		for (uint32_t j = 0; j < width; ++j) {
			if (convRes[i*width + j] > 0) {
				resImageMatrix[i*width + j] = 1;
			}
		}
	}
	delete[] convRes;
	delete imageFilter;
}
#endif
// int main(int argc, char const *argv[]) {
// 	int matrix[100];
// 	int resMatrix[100];
// 	for (int i = 0; i < 100; ++i) {
// 		if (i/10 > 2 && i/10 < 6) {
// 			matrix[i] = 1;
// 		}else{
// 			matrix[i] = 0;
// 		}
// 	}
// 	cout<<endl<<endl;
// 	for (int i = 0; i < 10; ++i) {
// 		for (int j = 0; j < 10; ++j) {
// 			cout<<matrix[i*10+j];
// 		}
// 		cout<<endl;
// 	}

// 	MorphologyImageProcessor<int> *mor = new MorphologyImageProcessor<int>(MorphologyImageProcessor<int>::SQUARE, 3,3);
// 	mor->Dilation(matrix, resMatrix, 10, 10);

// 	cout<<endl<<endl;
// 	for (int i = 0; i < 10; ++i) {
// 		for (int j = 0; j < 10; ++j) {
// 			cout<<resMatrix[i*10+j];
// 		}
// 		cout<<endl;
// 	}
// 	return 0;
// }

