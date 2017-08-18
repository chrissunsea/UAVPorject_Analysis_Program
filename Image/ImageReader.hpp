/*****************************************************************************/
/**Author: Gaojie Sun                                                       **/
/**Date: 21/07/2017                                                         **/
/**Supervisor: Dongryeol Ryu                                                **/
/**Originazation: The University of Melbourne                               **/
/*****************************************************************************/

/*****************************************************************************/
/**Documentation Description:                                               **/
/**   This file is used to define some methods that used to read the image. **/
/** They are C code, not the C++ code.It consider the future work may inclde**/
/** different image formats.                                                **/
/*****************************************************************************/
#ifndef _IMAGE_READER_H
#define _IMAGE_READER_H

#include <stdarg.h>
#include <cstring>
#include <stdint.h>
#include "Image.hpp"

#define TIFF_FORMAT 1
#define DEFAULT    -1
//Get Image's format (tif formate)
int  GetImageFormat(string fileName);
//Get the value type of each pixel (it may be 8 bits, 16 bits or 32 bits).
uint16_t  GetPixelValueType(string pathToImage);
uint16_t  GetPixelValueType_TIFF(string pathToImage);
//Decode Image: This is used to set the image matrix, image height and width 
//of the object Image<SRC_T>* image  
template <typename SRC_T>
void DecodeImage(string pathToImage, Image<SRC_T>* image);
//This is called by DecodeImage used to read tiff format image and set the object
template <typename SRC_T>
void ReadImage_TIFF(string fileName, Image<SRC_T>* image);

void GrFmtSilentTIFFErrorHandler(const char*, const char*, va_list );

#include "ImageReader.cpp"

#endif
