/*****************************************************************************/
/**Author: Gaojie Sun                                                       **/
/**Date: 21/07/2017                                                         **/
/**Supervisor: Dongryeol Ryu                                                **/
/**Originazation: The University of Melbourne                               **/
/*****************************************************************************/

/*****************************************************************************/
/**Documentation Description:                                               **/
/**   This file is used to define the Thermal Image class information that  **/
/** contains the temperature matrix and species object and CWSI matrix.     **/
/** The main operation of this class is to transform the temperature        **/
/** inforamtion to CWSI information.                                        **/
/*****************************************************************************/
#ifndef _THERMAL_IMAGE_CPP
#define _THERMAL_IMAGE_CPP

#include <stdlib.h>
#include <cmath>
#include <time.h>
#include <fstream>
#include <iostream>
#include <tiffio.h>
#include "ThermalImage.hpp"
#include "../GoogleEarth/GoogleEarthImage.hpp"
#include "../ImageProcessing/CannyEdgeDetection.hpp"
#include "../ImageProcessing/MorphologyImageProcessor.hpp"
using namespace std;
/*****************************************************************************/
/** Method Name: ThermalImage (Constructor)                                 **/
/** Method Input:                                                           **/
/**       (1) string pathToFile: The path to the image.                     **/
/**       (2) string speciesFile: The path to the species boundary          **/
/**           information, kml file.                                        **/
/** Method Output: Generate an object containing the thermal image matrix   **/
/**       and its related information.                                      **/
/** Method Description:                                                     **/
/**       This method first is to generate the temperature inforamtion from **/
/** the original image matrix. Then it eliminates the nosie temperature data**/
/** After that, use the denoised data to initialise the species object.     **/
/*****************************************************************************/
template<typename SRC_T>
ThermalImage<SRC_T>::ThermalImage(string pathToFile, string speciesFile) 
: Image<SRC_T>(pathToFile) {
	this->temperatureImageMatrix = new double[this->height * this->width];
	this->CWSI = new double[this->height * this->width];
	memset(this->temperatureImageMatrix, 0, this->height*this->width*sizeof(double));
	this->GenerateTemperatureInfo();

	double* tempEdgeOutImageMatrix = new double[this->height * this->width];
	this->TemperatureNoiseElimination(tempEdgeOutImageMatrix);
	this->species = new Species(speciesFile, pathToFile, tempEdgeOutImageMatrix, this->height, this->width);
	delete[] tempEdgeOutImageMatrix;
	
}
/*****************************************************************************/
/** Method Name: ~ThermalImage (Deconstructor)                              **/
/** Method Input: NULL                                                      **/
/**       The memory allocated for each Species object needs to be recycled,**/
/** when destory the Species object.                                        **/
/*****************************************************************************/
template<typename SRC_T>
ThermalImage<SRC_T>::~ThermalImage() {
	delete[] this->temperatureImageMatrix;
	delete[] this->CWSI;
	delete this->species;
	this->temperatureImageMatrix = NULL;
}
/*****************************************************************************/
/** Method Name: GenerateTemperatureInfo                                    **/
/** Method Input:                                                           **/
/**     (1) (onject's inner data) this->imageMatrix(from the Image class):  **/
/**         the original thermal image matrix.                              **/
/** Method Output:                                                          **/
/**     (1) (onject's inner data) double* temperatureImageMatrix: the       **/
/**         temperatue matrix.                                              **/
/** Method Description:                                                     **/
/**     This method is responsible for transforming the original thermal    **/
/** image data to the temperature data matrix.                              **/
/*****************************************************************************/
template<typename SRC_T>
void ThermalImage<SRC_T>::GenerateTemperatureInfo() {
	if (this->temperatureImageMatrix == NULL) {
		cout<<"Error(ThermalImage): temperatureImageMatrix has not been initialised"<<endl;
		exit(-1);
	}

	uint32_t rowOffset;
	double  MAXVALUE = pow (2, sizeof(SRC_T)*8);
	cout<<"Max pixel value is "<<MAXVALUE<<endl;
	for (uint32_t i = 0; i < this->height; ++i) {
		rowOffset = i*this->width;
		for (uint32_t j = 0; j < this->width; ++j) {
			this->temperatureImageMatrix[rowOffset+j] = this->imageMatrix[rowOffset+j] * 60.0 / MAXVALUE + 20;
		}
	}
}
/*****************************************************************************/
/** Method Name: TemperatureNoiseElimination                                **/
/** Method Input:                                                           **/
/**    (1) double* tempDenoise: the empty matrix that used to contain the   **/
/**        result matrix. This method is not responsible for allocating the **/
/**        memory.                                                          **/
/** Method Output:                                                          **/
/**    (1) double* tempDenoise: the temperature matrix that have been       **/
/**        eliminated the edge temperature information.                     **/
/** Method Description:                                                     **/
/**     This method first uses the Canny Edge Detection to detect the edge. **/
/** Then it widens the edge to ensure the noise data has been eliminated.   **/
/** After that, the by multiplying two matrix: temperature matrix and       **/
/** binarized edge matrix to eliminate the noise data.                      **/
/*****************************************************************************/
template<typename SRC_T>
void ThermalImage<SRC_T>::TemperatureNoiseElimination(double* tempDenoise) {
	char* edge = new char[this->height * this->width];
	this->EdgeDetection(edge);
	this->EdgeDilation(edge);
	this->EdgeTemperatureElimination(edge, tempDenoise);
	delete[] edge;
}
/*****************************************************************************/
/** Method Name: EdgeDetection                                              **/
/** Method Input:                                                           **/
/**    (1) char* resEdge: the empty matrix that used to contain the result  **/
/**        edge information matrix. This method is not responsible for      **/
/**        allocating the memory.                                           **/
/**    (2) (Object's inner data)this->imageMatrix: the original image's     **/
/**        matrix that need to do the edge detection. Actually it is the    **/
/**        inner data from the Image class.                                 **/
/**    (3) (Object's inner data)uint32_t height: The height of the image.   **/
/**    (4) (Object's inner data)uint32_t width: The width of the image.     **/
/** Method Output:                                                          **/
/**    (1) char* resEdge: the edge information matrix where 1 represents the**/
/**        edge point and 0 represents the non edge point.                  **/
/** Method Description:                                                     **/
/**    This method uses the canny edge detection to detect the edge.        **/
/*****************************************************************************/
template<typename SRC_T>
void ThermalImage<SRC_T>::EdgeDetection(char* resEdge) {
	if(this->imageMatrix == NULL) {
		cout<<"Error(Image): There is no image!"<<endl;
		 exit(-1);
	}
	
	CannyEdgeDetector<SRC_T> canny = CannyEdgeDetector<SRC_T>();
	canny.DetectEdge(this->imageMatrix, resEdge, this->height, this->width);
}
/*****************************************************************************/
/** Method Name: EdgeDilation                                               **/
/** Method Input:                                                           **/
/**    (1) char* resEdge: the binarized matrix that indicates the edge      **/
/**        information. The 1 represents the edge points and 0 represents   **/
/**        the non edge points.                                             **/
/** Method Output:                                                          **/
/**    (1) char* resEdge: the widened binarized matrix that indicates the   **/
/**        edge information. The 1 represents the edge points and 0         **/
/**        represents the non edge points.                                  **/
/** Method Description:                                                     **/
/**    This method is used to widen the edge of the image to ensure most of **/
/** the noise is eliminated.                                                **/
/*****************************************************************************/
template<typename SRC_T>
void ThermalImage<SRC_T>::EdgeDilation(char* resEdge) {
	MorphologyImageProcessor<char> *mor = new MorphologyImageProcessor<char>(MorphologyImageProcessor<char>::SQUARE, 5,5);
	mor->Dilation(resEdge, resEdge, this->height, this->width);
	delete mor;
}
/*****************************************************************************/
/** Method Name: EdgeTemperatureElimination                                 **/
/** Method Input:                                                           **/
/**    (1) char* resEdge: the binarized matrix that indicates the edge      **/
/**        information. The 1 represents the edge points and 0 represents   **/
/**        the non edge points.                                             **/
/**    (2) double* tempEdgeOutMatrix: the empty matrix that used to contain **/
/**        the result matrix. This method is not responsible for allocating **/
/**        the memory.                                                      **/
/** Method Output:                                                          **/
/**    (1) double* tempEdgeOutMatrix: the temperature matrix that have been **/
/**        eliminated the edge temperature information.                     **/
/** Method Description:                                                     **/
/**  	By multiplying two matrix: temperature matrix and binarized edge    **/
/** matrix to eliminate the noise data.                                     **/
/*****************************************************************************/
template<typename SRC_T>
void ThermalImage<SRC_T>::EdgeTemperatureElimination(char* resEdge, double* tempEdgeOutMatrix) {
	if (resEdge == NULL || tempEdgeOutMatrix == NULL) {
		cout<<"Error(Image): There is no edge or edge out matrix information"<<endl;
		exit(-1);
	}

	uint32_t dataIndex = 0;
	int mask = 0;
	for (uint32_t i = 0; i < this->height; ++i) {
		dataIndex = i*this->width;
		for (uint32_t j = 0; j < this->width; ++j) {
			mask = resEdge[dataIndex + j] > 0? 0 : 1;
			tempEdgeOutMatrix[dataIndex + j] = this->temperatureImageMatrix[dataIndex + j] * mask;
		}
	}
}
/*****************************************************************************/
/** Method Name: AnalyseCWSI                                                **/
/** Method Input: NULL                                                      **/
/** Method Output:                                                          **/
/**   (1) double* CWSI: the crop water stress index matrix to indicate the  **/
/**       the condition of crop water stress.                               **/
/** Method Description:                                                     **/
/** 	The calculation formula of CWSI is: CWSI = (Tc - Tw)/(Td - Tw). Tc  **/
/** means the pixel's value; Tw means the lowest temperature boundary; Td   **/
/** means the highest temperature boudary.                                  **/
/**   This method first calls another function CalTemperatureBoundary() to  **/
/** calculate the Td and Tw.                                                **/
/**   Then it will apply the CWSI formula.                                  **/
/*****************************************************************************/
template<typename SRC_T>
void ThermalImage<SRC_T>::AnalyseCWSI() {
	this->CalTemperatureBoundary();
	this->AdaptiveCWSI();
}
/*****************************************************************************/
/** Method Name: CalTemperatureBoundary                                     **/
/** Method Input:                                                           **/
/**        (1) (Object's inner data)Species* species: the species object.   **/
/** Method Output:                                                          **/
/**        (1) (Object's inner data)Species* species: the species object, to**/
/**            be more specific, the double* tempWet and double* tempDry    **/
/**            is initialised.                                              **/
/** Method Description:                                                     **/
/** 	The calculation formula of CWSI is: CWSI = (Tc - Tw)/(Td - Tw). Tc  **/
/** means the pixel's value; Tw means the lowest temperature boundary; Td   **/
/** means the highest temperature boudary.This function is used to determine**/
/** the Td and Tw for each species.                                         **/
/*****************************************************************************/
template<typename SRC_T>
void ThermalImage<SRC_T>::CalTemperatureBoundary() {
	this->species->CalTemperatureBoundary();
}
/*****************************************************************************/
/** Method Name: AdaptiveCWSI                                               **/
/** Method Input:                                                           **/
/**    (1) (Object's inner data)ddouble* CWSI: The empty matrix to return   **/
/**        the result. This method is not responsible for allocating the    **/
/**        memory.                                                          **/
/**    (2) (onject's inner data) double* temperatureImageMatrix: the        **/
/**         temperatue matrix.                                              **/
/**    (3) (Object's inner data)uint32_t height: The height of the image.   **/
/**    (4) (Object's inner data)uint32_t width: The width of the image.     **/
/** Method Output:                                                          **/
/**    (1) (Object's inner data)ddouble* CWSI: The CWSI matrix to indicate  **/
/**        the crop water stress index.                                     **/
/**    (2) One CWSI image in the current folder.                            **/
/**    (3) One folder that contains the CWSI image and KML file used to show**/
/**        the result on Google Earth.                                      **/
/** Method Description:                                                     **/
/**    First, the method will calculate the CWSI matrix whose pixels' ranges**/
/**    from the 0 to 1.1. Then the method will generate the CWSI image using**/
/**    the jet color map to map the pixel's value to the corresponding RGB  **/
/**    value.                                                               **/
/*****************************************************************************/
template<typename SRC_T>
void ThermalImage<SRC_T>::AdaptiveCWSI() {
	this->species->AdaptiveCWSI(this->CWSI, this->temperatureImageMatrix, this->height, this->width);

	time_t now = time(0);
    struct tm  tstruct;
    char buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d-%H-%M-%S", &tstruct);
	string fileName = "CWSI_Result_"+string(buf)+".tif";

	this->GenerateCWSIImage(fileName);
}
/*****************************************************************************/
/** Method Name: GenerateCWSIImage                                          **/
/** Method Input:                                                           **/
/**    (1) string fileName: the file name of the CWSI image.                **/
/**    (2) (Object's inner data)ddouble* CWSI: The CWSI matrix to indicate  **/
/**        the crop water stress index.                                     **/
/** Method Output:                                                          **/
/**    (1) One CWSI image in the current folder.                            **/
/** Method Description:                                                     **/
/**    The method will generate the CWSI image using the jet color map to   **/
/**    map the pixel's value to the corresponding RGB value.                **/
/*****************************************************************************/
template<typename SRC_T>
void ThermalImage<SRC_T>::GenerateCWSIImage(string fileName) {
	int samplePerPixel = 3;
	int bytesPerSample = 1;
	char *image=new char[this->width * this->height * samplePerPixel];
	this->ConvertRGBImage(image, samplePerPixel);

    TIFF *output_image;
    // Open the TIFF file
    output_image = TIFFOpen(fileName.c_str(), "w");

    // We need to set some values for basic tags before we can add any data
    TIFFSetField(output_image, TIFFTAG_IMAGEWIDTH, this->width);
    TIFFSetField(output_image, TIFFTAG_IMAGELENGTH, this->height);
    TIFFSetField(output_image, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(output_image, TIFFTAG_SAMPLESPERPIXEL, samplePerPixel);
    TIFFSetField(output_image, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(output_image, TIFFTAG_COMPRESSION, COMPRESSION_DEFLATE);
    TIFFSetField(output_image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);

    // Write the information to the file
    TIFFWriteEncodedStrip(output_image, 0, image, this->width*this->height * 3);

    // Close the file
    TIFFClose(output_image);
	string KMLFileName = fileName.substr(0,strrchr(fileName.c_str(), '.') - fileName.c_str());
    GoogleEarthImage googleEarthImage(KMLFileName);
	googleEarthImage.CopyImageToKMZ(fileName);
	googleEarthImage.GenerateResultKMZ(this->cornerCoordinates, 4);
	delete[] image;
}
/*****************************************************************************/
/** Method Name: ConvertRGBImage                                            **/
/** Method Input:                                                           **/
/**      (1) char* rgbImage: The empty matrix to contain the rgb matrix of  **/
/**          the CWSI image. This method is not responsible for allocating  **/
/**          the memory.                                                    **/
/**      (2) int samplePerPixel: the number of channel of the final image.  **/
/** Method Output:                                                          **/
/**      (1) char* rgbImage: the RGB matrix, for each pixel, the sequence is**/
/**          r,g,b.                                                         **/
/** Medthod Description:                                                    **/
/**       This method map the CWSI matrix to rgb matrix according to the jet**/
/** color map.                                                              **/
/*****************************************************************************/
template<typename SRC_T>
void ThermalImage<SRC_T>::ConvertRGBImage(char* rgbImage, int samplePerPixel) {
	COLOUR rgbColor;
	uint32_t index;
	for (uint32_t i = 0; i < this->height; ++i) {
		for (uint32_t j = 0; j < this->width; ++j) {
			rgbColor = GetColour(this->CWSI[i*this->width + j], 0, 1);
			index = (i*this->width + j) * samplePerPixel;
			rgbImage[index] = rgbColor.r;
			rgbImage[index+1] = rgbColor.g;
			rgbImage[index+2] = rgbColor.b;
		}
	}
}
/*****************************************************************************/
/** Method Name: GetColour                                                  **/
/** Method Input:                                                           **/
/**      (1) double colorData: the value of the data that need to be mapped **/
/**          to the rgb value according to the jet color map.               **/
/**      (2) double colorMin: the min value for the data.                   **/
/**      (3) double colorMin: the max value for the data.                   **/
/** Method Output:                                                          **/
/**      (1) COLOUR c: the corresponding rgb value.                         **/
/*****************************************************************************/
template<typename SRC_T>
COLOUR ThermalImage<SRC_T>::GetColour(double colorData,double colorMin,double colorMax) {
    COLOUR c = {0,0,0}; // black
    if (colorData > colorMax || colorData <= colorMin) { 
       return c;
    }
    
    double diffColorValue = colorMax - colorMin;
    double levelOne = colorMin + (2.0/20)*diffColorValue;   //blue
    double levelTwo = levelOne + (5.0/20)*diffColorValue;;   //blue & green
    double levelThree = levelTwo + (5.0/20)*diffColorValue; //green & red
    double levelFour = levelThree + (5.0/20)*diffColorValue;  //red color
   	double levelFive = colorMax;

    if (colorData >= levelFour) {
    	diffColorValue = levelFive - levelFour;
    	c.g = 0;
    	c.b = 0;
    	c.r = 255 - 128 * ( (colorData - levelFour)/diffColorValue );
    }else if (colorData >= levelThree) {
    	diffColorValue = levelFour - levelThree;
    	c.g = 255 - 255 * ( (colorData - levelThree)/diffColorValue );
    	c.b = 0;
    	c.r = char(255);
    }else if (colorData >= levelTwo) {
    	diffColorValue = levelThree - levelTwo;
    	c.g = char(255);
    	c.b = 255 - 255 * ( (colorData - levelTwo)/diffColorValue );
    	c.r = 255 * ( (colorData - levelTwo)/diffColorValue );
    }else if (colorData > levelOne) {
    	diffColorValue = levelTwo - levelOne;
    	c.g = 255 * ( (colorData - levelTwo)/diffColorValue );
    	c.b = char(255);
    	c.r = 0;
    }else if (colorData > colorMin) {
    	diffColorValue = levelOne - colorMin;
    	c.g = 0;
    	c.b = 255 * ( (colorData - levelTwo)/diffColorValue );
    	c.r = 0;
    }
    return(c);
}

#endif
