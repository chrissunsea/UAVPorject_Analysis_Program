#include <iostream>
#include <string.h>

#include <stdlib.h>
#include "Image/ThermalImage.hpp"
using namespace std;

#define PARAMETER_IMAGE    0
#define PARAMETER_BOUNDARY 1
#define ERROE_COMMAND     -1

#define ANALYSE_CWSI 1

void PrintPrompt() {
	cout<<"Please input the image path and boundary path:"<<endl;
	cout<<"format:"<<endl<<"\t-i /path/to/image(s) -b /path/to/boundaryfile(s)"<<endl;
	cout<<"It supports batch processing by passing multiple images and multiple files to this program"<<endl;
}

void DecodeParameters(int argc, char **argv, string* images, string* boundaryFiles) {

	int mode = -1;

	if ((argc - 3)%2 != 0) {
		cout<<"The number of image(s) is less(greater) than that of boundary file(s)"<<endl;
		exit(-1);
	}

	int imageNumber = (argc - 3)/2;
	int imageCount = 0;
	int boundaryCount = 0;
	for (int i = 1; i < argc; i++){
		if (strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "-I") == 0){
			mode = PARAMETER_IMAGE;
		}else if (strcmp(argv[i], "-b") == 0 || strcmp(argv[i], "-B") == 0){
			mode = PARAMETER_BOUNDARY;
		}else if (mode == -1) {
			delete[] images;
			delete[] boundaryFiles;
			PrintPrompt();
			exit(-1);
		}else {
			if (mode == PARAMETER_IMAGE && imageCount < imageNumber) {
				images[imageCount] = string(argv[i]);
				imageCount++;
			}else if (mode == PARAMETER_BOUNDARY && boundaryCount < imageNumber) {
				boundaryFiles[boundaryCount] = string(argv[i]);
				boundaryCount++;
			}else if(imageCount > imageNumber || boundaryCount > imageNumber) {
				mode = -1;
			}
		}
	}
}

void ProcessImage_8bits(string pathToFile, string boundaryfile) {
	cout<<"8 bits image"<<endl;
	ThermalImage<uint8_t>* thermalImage = new ThermalImage<uint8_t>(pathToFile, boundaryfile);
	thermalImage->AnalyseCWSI();
	delete thermalImage;
}

void ProcessImage_16bits(string pathToFile, string boundaryfile) {
	cout<<"16 bits image"<<endl;
	ThermalImage<uint16_t>* thermalImage = new ThermalImage<uint16_t>(pathToFile, boundaryfile);
	thermalImage->AnalyseCWSI();
	delete thermalImage;
}

void ProcessImage_32bits(string pathToFile, string boundaryfile) {
	cout<<"32 bits image"<<endl;
	ThermalImage<uint32_t>* thermalImage = new ThermalImage<uint32_t>(pathToFile, boundaryfile);
	thermalImage->AnalyseCWSI();
	delete thermalImage;
}

int main(int argc, char **argv) {

	if (argc <= 1) {
		PrintPrompt();
	}
	
	int imageNumber = (argc - 3)/2;
	
	string* images = new string[imageNumber];
	string* boundaryfiles = new string[imageNumber];
	
	DecodeParameters(argc, argv, images, boundaryfiles);
	int bitsPerPixel = -1;
	for (int i = 0; i < imageNumber; ++i) {

		bitsPerPixel = GetPixelValueType(images[i]);

		if (bitsPerPixel == 8) {
			ProcessImage_8bits(images[i], boundaryfiles[i]);
		}else if(bitsPerPixel == 16){
			ProcessImage_16bits(images[i], boundaryfiles[i]);
		}else if(bitsPerPixel == 32){
			ProcessImage_32bits(images[i], boundaryfiles[i]);
		}
	}
	delete[] images;
	delete[] boundaryfiles;
    return 0;
}

