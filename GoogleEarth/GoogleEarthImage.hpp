#ifndef GOOGLE_EARTH_IMAGE_H
#define GOOGLE_EARTH_IMAGE_H

#include <string>
#include "../Image/PositionDataStructure.h"
using namespace std;

class GoogleEarthImage {

	public:
		GoogleEarthImage(string fileName);
		~GoogleEarthImage();
		void CopyImage(string pathToImage);
		void GenerateResultKMZ(Coordinate* CornerCoords, int pointNumber);
		bool CopyImageToKMZ(string pathToImage);
	private:
		string fileName;
		string imageName;
		string imageHerf;
};

#endif