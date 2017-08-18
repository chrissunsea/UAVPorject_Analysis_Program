#include "GoogleEarthImage.hpp"
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <iomanip>

GoogleEarthImage::GoogleEarthImage(string fileName) {
    this->fileName = fileName;
    //unix create directory
	  const int dir_err = mkdir(fileName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	  if (-1 == dir_err) {
    	cout<<"Error creating directory!"<<endl;
    	exit(-1);
	  }
}

GoogleEarthImage::~GoogleEarthImage() {

}

bool GoogleEarthImage::CopyImageToKMZ(string pathToImage) {
  this->imageHerf = this->fileName+string("/")+this->fileName+string(".tif");
  this->imageName = this->fileName+string(".tif");
  ifstream src(pathToImage, std::ios::binary);
  ofstream dest(this->imageHerf, std::ios::binary);
  dest << src.rdbuf();
  return src && dest;
}

void GoogleEarthImage::GenerateResultKMZ(Coordinate* cornerCoords, int pointNumber) {
	ofstream myfile;
	string fileName = this->fileName+"/"+this->fileName+".kml";
  myfile.open (fileName);
  	
  myfile.setf(ios::showpoint);
  myfile.precision(14);
  myfile.setf(ios::fixed);
	myfile << string("<?xml version=\"1.0\"")+string(" encoding=\"UTF-8\"?>")<<endl;
	myfile << string("<kml xmlns=\"http:/")+string("/www.opengis.net/kml/2.2\"")+string(" xmlns:gx=\"http:/")+string("/www.google.com/kml/ext/2.2\"")+string(" xmlns:kml=\"http:/")+string("/www.opengis.net/kml/2.2\"")+string(" xmlns:atom=\"http://www.w3.org/2005/Atom\">")<<endl;
	myfile << string("<GroundOverlay>")<<endl;
	myfile << string("<name>image</name>")<<endl;
	myfile << string("<color>a1ffffff</color>")<<endl;
	myfile << string("<Icon>")<<endl;
	myfile << string("<href>")+this->imageName+string("</href>")<<endl;
	myfile << string("<viewBoundScale>0.75</viewBoundScale>")<<endl;
	myfile << string("</Icon>")<<endl;
	myfile << string("<gx:LatLonQuad>")<<endl;
  myfile << string("<coordinates>")<<endl;
  for (int i = 0; i < pointNumber; ++i) {
  	myfile << cornerCoords[i].longtitude;
  	myfile <<",";
  	myfile << cornerCoords[i].latitude;
  	myfile <<" ";
  }
  myfile << string("</coordinates>")<<endl;
  myfile << string("</gx:LatLonQuad>")<<endl;
	myfile << string("</GroundOverlay>")<<endl;
	myfile << string("</kml>")<<endl;
}
