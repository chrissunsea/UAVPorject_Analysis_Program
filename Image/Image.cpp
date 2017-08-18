/*****************************************************************************/
/**Author: Gaojie Sun                                                       **/
/**Date: 21/07/2017                                                         **/
/**Supervisor: Dongryeol Ryu                                                **/
/**Originazation: The University of Melbourne                               **/
/*****************************************************************************/

/*****************************************************************************/
/**Documentation Description:                                               **/
/**   This file is used to define the class that abstracts the image data   **/
/** and its geoinformation.                                                 **/
/*****************************************************************************/
#ifndef _IMAGE_CPP
#define _IMAGE_CPP

#include <cstring>
#include <iostream>
#include "Image.hpp"
#include "ImageReader.hpp"
using namespace std;
/*****************************************************************************/
/** Method Name: Image (Constructor)                                        **/
/** Method Input:                                                           **/
/**       (1) string pathToFile: The path to the image.                     **/
/** Method Output: Generate an object containing the image matrix and its   **/
/**       related information.                                              **/
/*****************************************************************************/
template<typename SRC_T>
Image<SRC_T>::Image(string pathToFile) {
	this->imageMatrix = NULL;
	//implementation in ./ImageReader.cpp
	DecodeImage(pathToFile, this);
}
/*****************************************************************************/
/** Method Name: ~Image (Deconstructor)                                     **/
/** Method Input: NULL                                                      **/
/**       The memory allocated for each Image object includes image matrix. **/
/** Therefore,when destory the canny object,the memory needs to be recycled.**/
/*****************************************************************************/
template<typename SRC_T>
Image<SRC_T>::~Image() {
	delete[] this->imageMatrix;	
	this->imageMatrix = NULL;
}
/*****************************************************************************/
/** Method Name: SetImage                                                   **/
/** Method Input:                                                           **/
/**     (1) SRC_T* imageMatrix: The image matrix used to set the inner data.**/
/**     (2) uint32_t height: The image's height.                            **/
/**     (3) uint32_t width: The image's width.                              **/
/** Method Output:                                                          **/
/**     (1) (Object's inner data member) SRC_T* imageMatrix:The image matrix**/
/**         used to set the inner data.                                     **/
/**     (2) (Object's inner data member)uint32_t height: The image's height.**/
/**     (3) (Object's inner data member)uint32_t width: The image's width.  **/
/*****************************************************************************/
template<typename SRC_T>
void Image<SRC_T>::SetImage(SRC_T* imageMatrix, uint32_t height, uint32_t width) {
	if(this->imageMatrix != NULL) {
		delete[] this->imageMatrix;
	}
	this->height = height;
	this->width = width;
	this->imageMatrix = new SRC_T[height * width];
	memcpy(this->imageMatrix, imageMatrix, height*width*sizeof(SRC_T));
}
/*****************************************************************************/
/** Method Name: SetConnerCoordinates                                       **/
/** Method Input:                                                           **/
/**      (1) int connerNum: This is to indicate which coordinate that will  **/
/**          be set.                                                        **/
/**      (2) double latitude: The latitude of the coordinate.               **/
/**      (3) double longtitude: The longtitude of the coordinate.           **/
/** Method Output:                                                          **/
/**      (1)this->cornerCoordinates[connerNum]: The corresponding coordinate**/
/** will be initialised.                                                    **/
/*****************************************************************************/
template<typename SRC_T>
void Image<SRC_T>::SetConnerCoordinates(int connerNum, double latitude, double longtitude) {
	this->cornerCoordinates[connerNum].latitude = latitude;
	this->cornerCoordinates[connerNum].longtitude = longtitude;
}
//return the image matrix
template<typename SRC_T>
SRC_T* Image<SRC_T>::GetImageMatrix() {
	return this->imageMatrix;
}
//return the image height
template<typename SRC_T>
uint32_t Image<SRC_T>::GetImageHeight() {
	return this->height;
}
//return the image width
template<typename SRC_T>
uint32_t Image<SRC_T>::GetImageWidth() {
	return this->width;
}
//return the corresponding coordinate
template<typename SRC_T>
Coordinate Image<SRC_T>::GetConnerCoordinates(int connerNum) {
	return this->cornerCoordinates[connerNum];
}

#endif