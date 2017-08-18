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
#ifndef _IMAGE_H
#define _IMAGE_H

#include <stdint.h>
#include <string>
#include "PositionDataStructure.h"
using namespace std;
/*****************************************************************************/
/** Class Descriptyion:                                                     **/
/**      This class is used to abstract the image data.                     **/
/**                                                                         **/
/** -Data Description:                                                      **/
/**      For one image, the data the analysis program needs is the image    **/
/**      matrix, the width and height of the image, and the longtitude and  **/
/**      latitude of the image.                                             **/
/**                                                                         **/
/** -Operation Description:                                                 **/
/**      The operations of this class is mainly the operations that set the **/
/**      object's inner data.                                               **/
/*****************************************************************************/	
template<typename SRC_T>
class Image {
	public:
		~Image();
		Image(string pathToFile);
		void SetImage(SRC_T* imageMatrix, uint32_t height, uint32_t width);
		
		SRC_T* GetImageMatrix();
		uint32_t GetImageHeight();
		uint32_t GetImageWidth();
		void SetConnerCoordinates(int connerNum, double latitude, double longtitude);
		Coordinate GetConnerCoordinates(int connerNum);

	protected:
		uint32_t width;
		uint32_t height;
		SRC_T* imageMatrix;
		//The longtitude and latitude of the image.
		//The sequence of the coordinate is: the lower left,
		//the lower right, the higher right and the higher left.
		Coordinate cornerCoordinates[4];
};

#include "Image.cpp"

#endif