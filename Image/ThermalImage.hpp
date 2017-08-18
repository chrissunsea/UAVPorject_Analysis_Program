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
#ifndef _THERMAL_IMAGE_H
#define _THERMAL_IMAGE_H

#include "Image.hpp"
#include "Species.hpp"
#include <string>
#include <stdint.h>

typedef struct {
	char r,g,b;
} COLOUR;
/*****************************************************************************/
/** Class Descriptyion:                                                     **/
/**      This class is used to abstract the thermal image data and operation**/
/**                                                                         **/
/** -Data Description:                                                      **/
/**      For one thermal image, the data the analysis program needs is the  **/
/**      temperature matrix that can be inferred from the pixel's value and **/
/**      the image object's basic information, such as height and width. So **/
/**      it extends the Image class.                                        **/
/**                                                                         **/
/** -Operation Description:                                                 **/
/**      The operations of this class is mainly to transform the temperature**/
/** inforamtion to CWSI information.                                        **/
/*****************************************************************************/	
template<typename SRC_T>
class ThermalImage:public Image<SRC_T>{
	public:
		ThermalImage(string pathToFile, string speciesFile);
		~ThermalImage();
		void AnalyseCWSI();
	private:
		void GenerateTemperatureInfo();
		void CalTemperatureBoundary();
		void TemperatureNoiseElimination(double* tempDenoise);
		void EdgeTemperatureElimination(char* resEdge, double* tempEdgeOutMatrix);
		void EdgeDetection(char* resEdge);
		void EdgeDilation(char* resEdge);
		void AdaptiveCWSI();
		void GenerateCWSIImage(string fileName);
		void ConvertRGBImage(char* rgbImage, int samplePerPixel);
		COLOUR GetColour(double v,double vmin,double vmax);
	private:
		Species* species;
		double* temperatureImageMatrix;
		double* CWSI;
};
#include "ThermalImage.cpp"
#endif