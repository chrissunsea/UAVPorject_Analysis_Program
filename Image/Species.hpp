/*****************************************************************************/
/**Author: Gaojie Sun                                                       **/
/**Date: 21/07/2017                                                         **/
/**Supervisor: Dongryeol Ryu                                                **/
/**Originazation: The University of Melbourne                               **/
/*****************************************************************************/

/*****************************************************************************/
/**Documentation Description:                                               **/
/**   This file is used to define the species information that contains the **/
/** species matrix that describe the species of each pixel and temperature  **/
/** data of each pixel.                                                     **/
/*****************************************************************************/
#ifndef _SPECIES_H
#define _SPECIES_H

#include <stdint.h>
#include <string>
#include "PositionDataStructure.h"
using namespace std;

class Species{
	public:
		~Species();
		Species(string speciesBoundaryFile, string pathToImage, double* temperatureMatrix, uint32_t height, uint32_t width);
		void CalTemperatureBoundary(); //calculate tempWet and tempDry
		void AdaptiveCWSI(double* resMatrix, double* temperatureMatrix, uint32_t height, uint32_t width);
	private:
		void SetSpeciesNum(string speciesBoundaryFile);
		void DecodeCoordinates(string content, int index); //decode the coordinate from the kml file
		void SetPolygonPoints(string pathToImage, string speciesBoundaryFile);
		void SetSpeciesTempData(double* tempData, uint32_t width, uint32_t height);

		void SetSpecisesMatrix();
		bool RayCasting(int speciesID, uint32_t i, uint32_t j);
		int IsIntersectAnt(double x, double y, double X1, double Y1, double X2, double Y2);
	private:
		// the species number in the image
		int speciesNum; 
		//the species matrix. Its size is equal to the original image. 
		//It represents the species of each pixel.
		int* speciesMatrix; 
		//The height of the speciesMatrix
		uint32_t  height;
		//The width of the speciesMatrix
		uint32_t  width;
		//The total pixel point number of each species
		uint32_t* totalPoint;
		//The temperature data of each species
		double**  data;
		//The wet temperature of each species
		double*   tempWet;
		//The dry temperature of each species
		double*   tempDry;
		//The polygon's information
		//The structure is define in ./PositionDataStructure.h
		Polygon*  speciesBoundaryInfo;
};

#endif