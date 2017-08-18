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
#ifndef _IMAGE_READER_CPP
#define _IMAGE_READER_CPP
#include <stdlib.h>
#include "ImageReader.hpp"
#include "xtiffio.h" //libtiff library
#include "geotiff.h" //libgeotiff library
#include "geo_normalize.h" //libgeotiff library

#include <iostream>
using namespace std;
/*****************************************************************************/
/** Method Name: DecodeImage                                                **/
/** Method Input:                                                           **/
/**       (1) string pathToFile: The path to the image.                     **/
/**       (2) Image<SRC_T>* image: The object used to store the data of the **/
/**           the given image.                                              **/
/** Method Output:                                                          **/
/**       (1) Image<SRC_T>* image: The object containing the data of the    **/
/**           given image.                                                  **/
/** Method Description:                                                     **/
/**        This method provides service that decoding the image and storing **/
/** the image data in the given Image object. Different image format have   **/
/** different decode code. First, it will check the suffix of the image name**/
/** to get the format of the image and then call the corresponding method to**/
/** decode the image.                                                       **/
/*****************************************************************************/
template <typename SRC_T>
void DecodeImage(string pathToImage, Image<SRC_T>* image) {
	int format = GetImageFormat(pathToImage);
	switch(format){
		case TIFF_FORMAT:
			cout<<"tif_format"<<endl;
			ReadImage_TIFF(pathToImage, image);break;
		case DEFAULT:
			cout<<"Error(ImageReader): Unknown image format"<<endl;exit(-1);
	}
}
/*****************************************************************************/
/** Method Name: ReadImage_TIFF                                             **/
/** Method Input:                                                           **/
/**       (1) string pathToFile: The path to the tiff image.                **/
/**       (2) Image<SRC_T>* image: The object used to store the data of the **/
/**           the given image.                                              **/
/** Method Output:                                                          **/
/**       (1) Image<SRC_T>* image: The object containing the data of the    **/
/**           given image.                                                  **/
/** Method Description:                                                     **/
/**        This method provides service that decoding the tiff format image **/
/** and storing the image data in the given Image object.                   **/
/*****************************************************************************/
template <typename SRC_T>
void ReadImage_TIFF(string pathToImage, Image<SRC_T>* image) {
  	
  	/************************************************************/
  	/** The below part is decode the image matrix information. **/
	/************************************************************/

	//These 2 methods provided by libtiff library
	//ignore the warning of unknowned tags
	TIFFSetErrorHandler(GrFmtSilentTIFFErrorHandler);
	TIFFSetWarningHandler(GrFmtSilentTIFFErrorHandler);

	//TIFF and GTIF data type define by libtiff and libgeotiff library.
	TIFF *img=(TIFF*)0;  /* TIFF-level descriptor */
	GTIF *gtif=(GTIF*)0; /* GeoKey-level descriptor */

	img = XTIFFOpen(pathToImage.c_str(), "r"); //libtiff library
	gtif= GTIFNew(img); //libtiff library

	if (!img || !gtif) {
		cout<<"Error(ImageReader): Cannot Open the file: "<<pathToImage<<endl;
		exit(-1);
	}

	uint16_t planarConfiguration;
	uint16_t channel;
	uint32_t height, width;

	//get the tag data from the tiff image
	TIFFGetField(img, TIFFTAG_IMAGEWIDTH,  &width);   //libtiff library
	TIFFGetField(img, TIFFTAG_IMAGELENGTH,  &height); //libtiff library
	TIFFGetField(img, TIFFTAG_PLANARCONFIG,  &planarConfiguration); //libtiff library
	TIFFGetField(img, TIFFTAG_SAMPLESPERPIXEL, &channel); //libtiff library
	
	SRC_T* imageMatrix = new SRC_T[height * width];

	//allocate memory for read buffer
	void* buffer = _TIFFmalloc(TIFFScanlineSize(img)); //libtiff library

	//If the image contains different channel, the data of different channel may be
	//organised into 2 ways: continuous way (one pixel's different channel data stored
	//continuously) and seperate way (one certain channel data stored continuously, but
	//different channel's data stored separately). (See the tiff image format)
	if (planarConfiguration == PLANARCONFIG_CONTIG){
		for (uint32_t i = 0; i < height; ++i){
			TIFFReadScanline(img, buffer, i); //libtiff library
			for (uint64_t j = 0; j < width * channel; j++){
				if(j%channel == 0) {
					imageMatrix[j/channel + i*width] = *(reinterpret_cast<SRC_T*> (buffer) + j);				
				}
			}
		}
	}else if (planarConfiguration == PLANARCONFIG_SEPARATE) {
		for (uint32_t i = 0; i < height; ++i) {
			TIFFReadScanline(img, reinterpret_cast<void*>(&imageMatrix[i* width]), i); //libtiff library
		}
	}else {
		cout<<"Error(ImageReader): Wrong image format: planarConfiguration"<<endl;
		exit(-1);
	}
	image->SetImage(imageMatrix, height, width);
	delete[] imageMatrix;

	cout.setf(ios::showpoint);
  	cout.precision(14);
  	cout.setf(ios::fixed);

  	/************************************************************/
  	/** The below part is decode the longtitude and latitude of**/
  	/** the image.                                             **/
	/************************************************************/
  	
	GTIFDefn defn;
    int connerNum = 4;
    double points[4][2] = { { 0.0, height*1.0 },
    						{ width*1.0, height*1.0 },
    						{ width*1.0, 0.0 },
    						{ 0.0, 0.0 }
                             };
    double x = 0.0,y= 0.0;
    if( GTIFGetDefn( gtif, &defn ) ) { //libtiff library
        for (int i = 0; i < connerNum; ++i) {
        	//libtiff library
        	//Map the points to longtitude and latitude
        	//If the tiff image doesn't contain the geoinformation
        	//then the analysis progtam will quit.
        	if (GTIFImageToPCS( gtif, &points[i][0], &points[i][1] ) && GTIFProj4ToLatLong( &defn, 1, &points[i][0], &points[i][1] )) {
        		image->SetConnerCoordinates(i, points[i][1], points[i][0]);
        	}else {
        		cout<<"Error(ImageReader): There is no geoinformation in this image!"<<endl;
        	}
        }
    } else {
     	cout<<"Error(ImageReader): This image does not contain any geoinformation!"<<endl;
     	exit(-1);

    }

    for (int i = 0; i < connerNum; ++i) {
		cout<<"x: "<<image->GetConnerCoordinates(i).latitude<<" y: "<<image->GetConnerCoordinates(i).longtitude<<endl;
	}
    
	GTIFFree(gtif); //libtiff library
	TIFFClose(img); //libtiff library
}
/*****************************************************************************/
/** Method Name: GetImageFormat                                             **/
/** Method Input:                                                           **/
/**      (1) string pathToImage: the path to the image.                     **/
/** Method Output:                                                          **/
/**      (1)return value: represents the image format.                      **/
/** Method Description:                                                     **/
/**      This method is get the format of the image, according to the image **/
/** name's suffix.                                                          **/
/*****************************************************************************/
int GetImageFormat(string pathToImage){
	string format = string( strrchr(pathToImage.c_str(), '.') + 1 );
	if (format.compare("tif") == 0){
		return TIFF_FORMAT;
	}else{
		return DEFAULT;
	}
}
/*****************************************************************************/
/** Method Name: GetPixelValueType                                          **/
/** Method Input:                                                           **/
/**      (1) string pathToImage: the path to the image.                     **/
/** Method Output:                                                          **/
/**      (1)return value: represents the number of the bit for each pixel.  **/
/** Method Description:                                                     **/
/**      Each pixel may use 8 bits, 16bits or 32 bits.                      **/
/*****************************************************************************/
uint16_t GetPixelValueType(string pathToImage) {
	int format = GetImageFormat(pathToImage);
	switch(format){
		case TIFF_FORMAT:
			cout<<"tif_format"<<endl;
			return GetPixelValueType_TIFF(pathToImage);
		case DEFAULT:
			return -1;
	}
	return -1;
}
/*****************************************************************************/
/** Method Name: GetPixelValueType_TIFF                                     **/
/** Method Input:                                                           **/
/**      (1) string pathToImage: the path to the image.                     **/
/** Method Output:                                                          **/
/**      (1)return value: represents the number of the bit for each pixel.  **/
/** Method Description:                                                     **/
/**      Each pixel may use 8 bits, 16bits or 32 bits. This is to get the   **/
/** pixel value for each pixel in an tiff image.                            **/
/*****************************************************************************/
uint16_t GetPixelValueType_TIFF(string pathToImage) {
	TIFFSetErrorHandler(GrFmtSilentTIFFErrorHandler); //libtiff 
	TIFFSetWarningHandler(GrFmtSilentTIFFErrorHandler); //libtiff

	TIFF* img=TIFFOpen(pathToImage.c_str(), "r"); //libtiff

	if (!img) {
		cout<<"Error(ImageReader): Cannot Open the file: "<<pathToImage<<endl;
		exit(-1);
	}
	uint16_t bitsPerPixel = 0;
	TIFFGetField(img, TIFFTAG_BITSPERSAMPLE, &bitsPerPixel); //libtiff
	TIFFClose(img); //libtiff
	return bitsPerPixel;
}

void GrFmtSilentTIFFErrorHandler(const char*, const char*, va_list ) {
	
}
#endif
