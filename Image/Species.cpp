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
#include <stdlib.h>
#include <list>
#include <regex>
#include <string.h>
#include <math.h>
#include <iostream>
#include "xtiffio.h"
#include "geotiff.h"
#include "Species.hpp"
#include "../GMM/GMM.h"
#include "geo_normalize.h"
using namespace std;
/*****************************************************************************/
/** Method Name: Species (Constructor)                                      **/
/** Method Input:                                                           **/
/**       (1) string speciesBoundaryFile: The path to the boundary kml file.**/
/**       (2) string pathToImage: The path to the image file.               **/
/**       (3) double* temperatureMatrix: The temperature matrix of the      **/
/**           original image.                                               **/
/**       (4) uint32_t height: the height of the image.                     **/
/**       (5) uint32_t width: the width of the image.                       **/
/** Method Output: Generate an object species information.                  **/
/*****************************************************************************/
Species::Species(string speciesBoundaryFile, string pathToImage, double* temperatureMatrix, uint32_t height, uint32_t width) {
	this->height = height;
	this->width = width;
	//set the species number form the kml file
	this->SetSpeciesNum(speciesBoundaryFile);
	cout <<"Species Number: "<<this->speciesNum<<endl;

	if (this->speciesNum <= 0) {
		cout<<"There is no species information in the species boundary file!"<<endl;
		exit(-1);
	}

	this->tempWet = new double[this->speciesNum];
	this->tempDry = new double[this->speciesNum];
	this->speciesMatrix = new int [height * width];
	this->totalPoint  = new uint32_t[this->speciesNum];

	memset(this->totalPoint, 0, this->speciesNum * sizeof(uint32_t));
	memset(this->tempWet, 0, this->speciesNum * sizeof(uint32_t));
	memset(this->tempDry, 0, this->speciesNum * sizeof(uint32_t));

	//Decode the species polygon position information
	this->SetPolygonPoints(pathToImage, speciesBoundaryFile);
	//Classify the pixels into the corresponding species
	this->SetSpecisesMatrix();
	cout<<"Species Classfication Done!"<<endl;
	this->SetSpeciesTempData(temperatureMatrix, height, width);
}
/*****************************************************************************/
/** Method Name: ~Species (Deconstructor)                                   **/
/** Method Input: NULL                                                      **/
/**       The memory allocated for each Species object needs to be recycled,**/
/** when destory the Species object.                                        **/
/*****************************************************************************/
Species::~Species() {
	delete[] this->speciesMatrix;
	delete[] this->totalPoint;
	delete[] this->tempWet;
	delete[] this->tempDry;
	for (int i = 0; i < speciesNum; ++i) {	
		delete[] this->data[i];
		delete[] this->speciesBoundaryInfo[i].coordinates;
		delete[] this->speciesBoundaryInfo[i].position;
	}
	delete[] this->data;
	delete[] this->speciesBoundaryInfo;
}
/*****************************************************************************/
/** Method Name: SetSpeciesNum                                              **/
/** Method Input:                                                           **/
/**       (1) string speciesBoundaryFile: The path to the boundary kml file.**/
/** Method Output: Set the inner data filed -- speciesNumber.               **/
/** Method Description:                                                     **/
/**       This method will analyse the kml file and count the number of     **/
/** string "<Polygon>" and set the speciesNum equal to this number.         **/
/*****************************************************************************/
void Species::SetSpeciesNum(string speciesBoundaryFile) {
	this->speciesNum = 0;
	ifstream file(speciesBoundaryFile);
	string lineBuffer;

	while(getline(file, lineBuffer)) {
		if (lineBuffer.find("<Polygon>") != string::npos) {
			this->speciesNum++;
		}
	}
	file.close();
}
/*****************************************************************************/
/** Method Name: SetPolygonPoints                                           **/
/** Method Input:                                                           **/
/**       (1) string pathToImage: The path to the image file.               **/
/**       (2) string speciesBoundaryFile: The path to the boundary kml file.**/
/** Method Output:                                                          **/
/**       (1) this->speciesBoundaryInfo: the species boundary information.  **/
/** Method Description:                                                     **/
/**       This method will analyse the kml file and extract the longtitude  **/
/** and latitude information of each polygon from kml file, then map the    **/
/** longtitude and latitude to the pixel's position in the image matrix.    **/
/*****************************************************************************/
void Species::SetPolygonPoints(string pathToImage, string speciesBoundaryFile) {
	this->speciesBoundaryInfo = new Polygon[this->speciesNum];
	
	TIFF *img=(TIFF*)0;  /* TIFF-level descriptor */
	GTIF *gtif=(GTIF*)0; /* GeoKey-level descriptor */
	GTIFDefn defn;

	img = XTIFFOpen(pathToImage.c_str(), "r"); //libtiff library
	gtif= GTIFNew(img); //libtiff library

	if (!img || !gtif) {
		cout<<"Error(ImageReader): Cannot Open the file: "<<pathToImage<<endl;
		exit(-1);
	}

	cout.setf(ios::showpoint);
  	cout.precision(14);
  	cout.setf(ios::fixed);

	ifstream file(speciesBoundaryFile);
	string lineBuffer;
	int index = -1;
	int coordinateFlag = 0;

	//decode the kml file to extract and analyse the polygon information
	while(getline(file, lineBuffer)) {
		if (coordinateFlag == 1) {
			this->DecodeCoordinates(lineBuffer, index);
 		}
		if (lineBuffer.find("<Polygon>") != string::npos) {
			index++;
		}else if (lineBuffer.find("<coordinates>") != string::npos) {
			coordinateFlag = 1;
		}else {
			coordinateFlag = 0;
		}
	}
	double longtitude, latitude;
	GTIFGetDefn( gtif, &defn );
	for (int i = 0; i < this->speciesNum; ++i) {
		for (int j = 0; j < this->speciesBoundaryInfo[i].pointNumber; ++j) {	
			longtitude = this->speciesBoundaryInfo[i].coordinates[j].longtitude;
			latitude = this->speciesBoundaryInfo[i].coordinates[j].latitude;
			
			//map the longtitude and latitude to the pixel's position in the image
			GTIFProj4FromLatLong(&defn, 1, &longtitude, &latitude);
			GTIFPCSToImage(gtif, &longtitude, &latitude);

			this->speciesBoundaryInfo[i].position[j].i = latitude; //pixel's i position
			this->speciesBoundaryInfo[i].position[j].j = longtitude; //pixel's j position
		}
	}
}
/*****************************************************************************/
/** Method Name: DecodeCoordinates                                          **/
/** Method Input:                                                           **/
/**      (1) string content: the content contains the coordinate(longtitude **/
/**          and latitude) information.                                     **/
/**      (2) int index: the index of the polygon indicating the species     **/
/**          index.                                                         **/
/** Method Output:                                                          **/
/**      (1) this->speciesBoundaryInfo: the species boundary information.   **/
/** Method Description:                                                     **/
/**       This method will use the regular expression to analyse the given  **/
/**string content to set the longtitude and latitude data in the data member**/
/** this->speciesBoundaryInfo.                                              **/
/*****************************************************************************/
void Species::DecodeCoordinates(string content, int index) {
	smatch matchString;
	regex pattern ("[+-]?(0|[1-9]+)\\.[0-9]*");
	list<double> list;
	double d;
	while (regex_search (content, matchString, pattern)){
		d = stod(matchString.str());
		list.push_back(d);
		content = matchString.suffix().str();
	}

	if (list.size() % 2 != 0) {
		cout<<"Error(Species): Decode Species boundary information!"<<endl;
		exit(-1);
	}

	int polygonPointNum = list.size() / 2;

	this->speciesBoundaryInfo[index].pointNumber = polygonPointNum;
	this->speciesBoundaryInfo[index].coordinates = new Coordinate[polygonPointNum];
	this->speciesBoundaryInfo[index].position = new Position[polygonPointNum];

	for (int i = 0; i < polygonPointNum; ++i)	{
		this->speciesBoundaryInfo[index].coordinates[i].longtitude = list.front();
		list.pop_front();
		this->speciesBoundaryInfo[index].coordinates[i].latitude = list.front();
		list.pop_front();
	}
	list.clear();
}
/*****************************************************************************/
/** Method Name: SetSpecisesMatrix                                          **/
/** Method Input: NULL                                                      **/
/** Method Output:                                                          **/
/**      (1) (Object's inner data field)int* speciesMatrix: the matrix      **/
/**          indicate the species of each pixel.                            **/
/** Method Description:                                                     **/
/**   This method uses the RayCast algorithm to check whether the point is  **/
/** within the polygons.                                                    **/
/*****************************************************************************/
void Species::SetSpecisesMatrix() {
	int s = 0;
	for (uint32_t i = 0; i < this->height; ++i) {
		for (uint32_t j = 0; j < this->width; ++j) {
			for (s = 0; s < this->speciesNum; ++s) {
				//use the ray cast algorithm to check whether the point(i,j)
				//is within the polygon s.
				if (this->RayCasting(s, i, j)) {
					this->speciesMatrix[i*this->width + j] = s;
					break;
				}
			}
			if (s == this->speciesNum) {
				this->speciesMatrix[i*this->width + j] = -1;
			} else {
				this->totalPoint[s]++;
			}
		}
	}
}
/*****************************************************************************/
/** Method Name: RayCasting                                                 **/
/** Method Input:                                                           **/
/**      (1) int speciesID: this indicates the index of the polygon(species)**/
/**      (2) uint32_t i: the row number of the pixel.                       **/
/**      (3) uint32_t j: the column number of the pixel.                    **/
/** Method Output:                                                          **/
/**      (1)return value(bool): If the point is within the polygon, it will **/
/** return true; otherwise, it will return false.                           **/
/** Method Description:                                                     **/
/**      This method applies the RayCast Algorithm to judge whether the     **/
/** point is within the polygon. In this algorithm, it determines if a point**/ 
/** lies inside a bounded region by casting infinitely long ray from the    **/
/** given point to see how many edges of the given polygon this ray         **/
/** intersects.If the number of intersections is 0 or an even number, the   **/
/** point lies outside the polygon otherwise inside.                        **/
/*****************************************************************************/
bool Species::RayCasting(int speciesID, uint32_t i, uint32_t j) {
	int count = 0;
	int polygonPointNum = this->speciesBoundaryInfo[speciesID].pointNumber;
	Position* pointer = this->speciesBoundaryInfo[speciesID].position;
    Position* p1;
    Position* p2;
    
    for (int index = 0; index < polygonPointNum-1; index++) {
        p1 = pointer + index;
        p2 = p1 + 1;  
        //judge whether the ray intersects with (p1,p2)     	
       	count += this->IsIntersectAnt(j, i, p1->j, p1->i, p2->j, p2->i);
    }

    return count%2 == 1 ? true : false;
}
/*****************************************************************************/
/** Method Name: IsIntersectAnt                                             **/
/** Method Input:                                                           **/
/**        (1) double x: the x position of the start point of the ray.      **/
/**        (2) double y: the y position of the start point of the ray.      **/
/**        (3) double X1: the x position of p1 of the edge.                 **/
/**        (4) double Y1: the y position of p1 of the edge.                 **/
/**        (5) double X2: the x position of p2 of the edge.                 **/
/**        (6) double Y2: the y position of p2 of the edge.                 **/
/** Method Output:                                                          **/
/**        (1) return value: if point p (x,y) intersects with edge (p1,p2)  **/
/**            then it will return 1; otherwise, it will return 0.          **/
/*****************************************************************************/

int Species::IsIntersectAnt(double x, double y, double X1, double Y1, double X2, double Y2) {

    double minX,maxX,minY,maxY;  
    minX = X1;  
    maxX = X2;  
    if (minX > maxX) {  
        minX = X2;  
        maxX = X1;  
    }

    minY = Y1;  
    maxY = Y2;  
    if (minY > maxY) {  
        minY = Y2;  
        maxY = Y1;  
    }  

    //preprocess to accelerate the judgement
    if (y<minY || y>maxY) {
        return 0;  
    }  
  
    //if the point is on the edge
    if (fabs(maxY - minY) == 0) {  
        return (x >= minX && x <= maxX)? 1:0;  
    }  
  
    //calculate the intersection point 
    double x0 = X1 + (double)(y - Y1)*(X2 - X1)/(Y2 - Y1);  
      
    //if the intersection point is on the right side of point p
    //then return 0
    if (x0 > x) {  
        return 0;  
    }  
    //if the intersection point is same to the start point p
    if (fabs(x-x0) == 0)  
    {  
        return 0;  
    }  
    //if the intersection point goes through the lower end point of the edge
    if (fabs(y-minY) == 0) {  
        return 0;  
    }  
    return 1;  
} 
/*****************************************************************************/
/** Method Name: SetSpeciesTempData                                         **/
/** Method Input:                                                           **/
/**     (1) double* tempData: the temperature matrix that are eliminated    **/
/**         noise data.                                                     **/
/**     (2) uint32_t width: the width of the image.                         **/
/**     (3) uint32_t height: the height of the image.                       **/
/** Method Output:                                                          **/
/**     (1) (Object's inner data)double** data: the temperature data of     **/
/**         each species in one dimension array.                            **/
/**     (2) (Object's inner data)uint32_t* totalPoint: the number of points **/
/**         of each species.                                                **/
/** Method Description:                                                     **/
/**     This method is used to extract the species temperature information  **/
/** that will be used to fit the GMM. Importantly, the temperature          **/
/** data used to fit the GMM does not include the useless data, such as the **/
/** the 0 degree temperature and 20 degree temperatue data(the lowest       **/
/** boundary) of the data.                                                  **/
/*****************************************************************************/
void Species::SetSpeciesTempData(double* tempData, uint32_t width, uint32_t height){
	int speciesCategory;
	uint32_t pointNum = 0;
	this->data = new double*[this->speciesNum];
	for (int i = 0; i < this->speciesNum; ++i) {
		pointNum = this->totalPoint[i];
		this->data[i] = new double[pointNum];
	}
	memset(this->totalPoint, 0, sizeof(uint32_t) * this->speciesNum);
	uint32_t rowOffset = 0;
	for (uint32_t i = 0; i < height; ++i) {
		rowOffset = i * width;
		for (uint32_t j = 0; j < width; ++j) {
			speciesCategory = this->speciesMatrix[rowOffset+j];
			if (speciesCategory == -1) {
				continue;
			}
			if(tempData[rowOffset + j] != 0 &&  tempData[rowOffset + j] != 20) {
				data[speciesCategory][this->totalPoint[speciesCategory]++] = tempData[rowOffset + j];
			}
		}
	}
}
/*****************************************************************************/
/** Method Name: CalTemperatureBoundary                                     **/
/** Method Input:                                                           **/
/**     (1) (Object's inner data)double** data: the temperature data of     **/
/**         each species in one dimension array.                            **/
/** Method Output:                                                          **/
/**    (1)(Object's inner data) double* tempWet: the lowest boundary of the **/
/**       GMM trained by a given species temperature data.                  **/
/**    (2)(Object's inner data) double* tempDry: the highest boundary of the**/
/**       GMM trained by a given species temperature data.                  **/
/** Method Description:                                                     **/
/**     This method uses the given data to fit the GMM and calculate the    **/
/** lowest and highest boundary temperature data.                           **/
/*****************************************************************************/
void Species::CalTemperatureBoundary() {
	int cropIndex = -1;
	GMM* gmm = new GMM(1,2);
	uint32_t totalPoint = 0;
	
	for (uint32_t i = 0; i < this->speciesNum; i++) {
		totalPoint = this->totalPoint[i];
		cout<<"Species"<<i<<": "<<"totalPoint: "<<totalPoint<<endl;
		gmm->Train(this->data[i], totalPoint);

		cropIndex = gmm->Mean(0)[0] > gmm->Mean(1)[0] ? 1:0;

		tempWet[i] = gmm->Mean(cropIndex)[0] - 2.58*sqrt(gmm->Variance(cropIndex)[0]);
		tempDry[i] = gmm->Mean(cropIndex)[0] + 2.58*sqrt(gmm->Variance(cropIndex)[0]);

		cout<<"Wet: "<<tempWet[i] <<" Dry: "<<tempDry[i]<<endl;
	}
	delete gmm;
}
/*****************************************************************************/
/** Method Name: AdaptiveCWSI                                               **/
/** Method Input:                                                           **/
/**    (1) double* resMatrix: The empty matrix to return the result. This   **/
/**        method is not responsible for allocating the memory.             **/
/**    (2) double* temperatureMatrix: The temperature matrix that need to   **/
/**        analyse the CWSI.                                                **/
/**    (3) uint32_t height: The height of the image.                        **/
/**    (4) uint32_t width: The width of the image.                          **/
/**    (5) (Object's inner data field)int*speciesMatrix:the matrix indicates**/
/**        the species of each pixel.                                       **/
/** Method Output:                                                          **/
/**    (1) (Object's inner data field) double* resMatrix: the result CWSI   **/
/**        matrix ranging from 0 to 1.1. Also, 1.1 represent this pixel does**/
/**        not need to be considered.                                       **/
/** Method Description:                                                     **/
/** The calculation formula is: CWSI = (Tc - Tw)/(Td - Tw). Tc means the    **/
/** pixel's value; Tw means the lowest temperature boundary; Td means the   **/
/** highest temperature boudary.                                            **/
/*****************************************************************************/
void Species::AdaptiveCWSI(double* resMatrix, double* temperatureMatrix, uint32_t height, uint32_t width) {
	
	uint32_t speciesCategory = -1;
	uint32_t dataIndex = 0;

	for (uint32_t i = 0; i < height; ++i) {
		dataIndex = i * width;
		for (uint32_t j = 0; j < width; ++j) {
			speciesCategory = this->speciesMatrix[dataIndex + j];

			if(speciesCategory == -1) {
				resMatrix[dataIndex + j] = 1.1;
			}else {
				resMatrix[dataIndex + j] = (temperatureMatrix[dataIndex + j] - this->tempWet[speciesCategory])/(this->tempDry[speciesCategory] - this->tempWet[speciesCategory]);
				if (resMatrix[dataIndex + j] > 1.0) {
					resMatrix[dataIndex + j] = 1.1;
				}
			}
		}
	}
}
