#include <string>
#include <iostream>
#include <fstream>
#include <sstream>  

#include <tiffio.h>
#include <stdlib.h>

#include "openslide.h"

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

using namespace std;
using namespace boost::filesystem;

// We assume that WSI has been scanned at 40X magnification
double scan_magnification = 40.0;
unsigned int Biased = 5;
double requested_magnification = 40.0;
unsigned int frame_width = 400;
unsigned int frame_height = 400;

unsigned int WSITiling(path OutputDirectory, path WSIFileNameWithPath, double roi_X, double roi_Y, double roi_width, double roi_height);
bool TileExtraction(const char *inputFilename, long x, long y, const char *outputFilename, vector<double> &imageFeatures);
void MyTest(void);

int main(int argc, char* argv[]) {
	MyTest();
	char c;
	cin >> c;
	return 1;
	
	/*if (argc < 2 ) {
		std::cerr << "Usage: " << argv[0] << " WSI_Name Scan_Magnification Req_Magnification Biased Frame_Width Frame_Height ROI_X ROI_Y ROI_Width ROI_Height" << std::endl << std::flush;
		return -1;
	}*/

	double roi_X = 1;
	double roi_Y = 1;
	double roi_width = -1;
	double roi_height = -1;

	//path fileNameWithPath = argv[1];
	path fileNameWithPath("R:/Beck Lab/Atypia_Andy/HighAgreement/1305.tif");

	if (argc >= 3)
		scan_magnification = atof(argv[2]);
	if (argc >= 4)
		requested_magnification = atof(argv[3]);
	if (argc >= 5)
		Biased = atoi(argv[4]);
	if (argc >= 6)
		frame_width = atoi(argv[5]);
	if (argc >= 7)
		frame_height = atoi(argv[6]);
	if (argc >= 8)
		roi_X = atof(argv[7]);
	if (argc >= 9)
		roi_Y = atof(argv[8]);
	if (argc >= 10)
		roi_width = atof(argv[9]);
	if (argc >= 11)
		roi_height = atof(argv[10]);

	string filePath = fileNameWithPath.parent_path().string();
	string fileNameWithExt = fileNameWithPath.leaf().string();
	//string fileExt =  fileNameWithExt.substr(fileNameWithExt.size() - 3, fileNameWithExt.size());
	string fileExt = fileNameWithPath.extension().string();
	boost::algorithm::to_lower(fileExt);
	
	cout << "File Path: " << filePath << "\nFileNameWithExt: " << fileNameWithExt << "\nFileExt: " << fileExt;

	/*Aperio(.svs, .tif)
	Hamamatsu(.vms, .vmu, .ndpi)
	Leica(.scn)
	MIRAX(.mrxs)
	Sakura(.svslide)
	Trestle(.tif)
	Generic tiled TIFF(.tif)*/
	
	if (fileExt == ".svs" || fileExt == ".mrxs" || fileExt == ".tif" || fileExt == ".vms" || fileExt == ".vmu" || fileExt == ".ndpi" || fileExt == ".scn" || fileExt == ".svslide"){

		cout << "\nStart tiling " << fileNameWithExt << " ... ";

		vector<string> nameToken;
		boost::split(nameToken, fileNameWithExt, boost::is_any_of("."));
		if (nameToken.size() == 3){
			path outputDirectory(filePath);
			outputDirectory /= nameToken[0];
			WSITiling(outputDirectory, fileNameWithPath, roi_X, roi_Y, roi_width, roi_height);
		}
		else{
			path outputDirectory (filePath);
			outputDirectory /= string(fileNameWithExt.substr(0, fileNameWithExt.find_last_of("."))) + string("_Tiles");
			WSITiling(outputDirectory, fileNameWithPath, roi_X, roi_Y, roi_width, roi_height);
		}
	}
	else if (fileExt == ".txt"){
			
		cout << "\nStart reading " << fileNameWithExt << " ... ";

		path WSINamesWithPath(filePath), WSINames(filePath);
		WSINamesWithPath /= string("WSINamesWithPath.txt");
		WSINames /= string("WSINames.txt");

		std::ofstream outputFile1(WSINamesWithPath.make_preferred().string().c_str());
		std::ofstream outputFile2(WSINames.make_preferred().string().c_str());

		std::ifstream inputFile(fileNameWithPath.make_preferred().string().c_str());
		
		string line;
		while (std::getline(inputFile, line)){

			istringstream ss(line);
			string id, svsFilename;
			ss >> id >> svsFilename;
			vector<string> nameToken;
			boost::split(nameToken, svsFilename, boost::is_any_of("."));

			if (nameToken.size() == 3){
				path svsFileNameWithPath(filePath), outputDirectory(filePath);
				outputDirectory /= nameToken[0];
				svsFileNameWithPath /= svsFilename;

				unsigned int noOfTiles = WSITiling(outputDirectory, svsFileNameWithPath, roi_X, roi_Y, roi_width, roi_height);
				if (noOfTiles > 0){
					outputFile1 << outputDirectory.make_preferred().string() << endl;
					outputFile2 << nameToken[0] << endl;
				}
			}
			else {
				path svsFileNameWithPath(filePath), outputDirectory(filePath);
				outputDirectory /= svsFilename.substr(0, svsFilename.find_last_of(".")) + string("_Tiles");
				svsFileNameWithPath /= svsFilename;

				unsigned int noOfTiles = WSITiling(outputDirectory, svsFileNameWithPath, roi_X, roi_Y, roi_width, roi_height);
				if (noOfTiles > 0){
					outputFile1 << outputDirectory.make_preferred().string() << endl;
					outputFile2 << svsFilename << endl;
				}
			}
		}
		outputFile1.close();
		outputFile2.close();
	}	//if (fileExt == "txt")
	else {
		fprintf(stderr, "\nNumber of Arguments are: %d", argc );
		fprintf(stderr, "\nUsage: %s inputWSI_withPath Biased requestedMagnification roi_X roi_Y roi_width roi_height frame_width frame_height\n", argv[0]);
		fprintf(stderr, "x and y: Upper left corner, coordinates in pixels\n");
		fprintf(stderr, "width and height: Number of pixels of the requested rectangle\n\n\n");

		exit(1);
	}
	
	return 0;
}

unsigned int WSITiling(path OutputDirectory, path WSIFileNameWithPath, double roi_X, double roi_Y, double roi_width, double roi_height)
{
	// In case bounding box width and height is not mention, then select the width and height of WSI as bounding box
	if (roi_width < 1 || roi_height < 1){
		// Open input file
		openslide_t *wsi = openslide_open(WSIFileNameWithPath.make_preferred().string().c_str());
		if (NULL == wsi) {
			fprintf(stderr, "\nERROR when opening file %s\n", WSIFileNameWithPath.make_preferred().string().c_str());
			exit(1);
		}

		// What is the best layer to use
		// for getting a portion of the WSI at magnification requested by the user
		//int32_t layer = openslide_get_best_layer_for_downsample(wsi, (scan_magnification / requested_magnification));
		int32_t layer = openslide_get_best_layer_for_downsample(wsi,1);
		if (-1 == layer) {
			fprintf(stderr, "\nERROR: Impossible to find an appropriate layer for scan magnification %lfX\n", scan_magnification);
			exit(1);
		}

		// Get dimensions of the available magnification
		int64_t widthLayer, heightLayer;
		openslide_get_layer_dimensions(wsi, layer, &widthLayer, &heightLayer);
		if ((-1 == widthLayer) || (-1 == heightLayer)) {
			fprintf(stderr, "\nERROR: Impossible to get the width and height of layer %d\n", layer);
			exit(1);
		}
		roi_width = widthLayer;
		roi_height = heightLayer;
	}

	string WSIPath = WSIFileNameWithPath.parent_path().string();
	string WSIFileNameWithExt = WSIFileNameWithPath.leaf().string();
	string WSIFileName = WSIFileNameWithExt.substr(0, WSIFileNameWithExt.find_last_of("."));

	path outputDirectory(OutputDirectory);
	if (!exists(outputDirectory.make_preferred().string())) {
		if (!create_directory(outputDirectory)){
			fprintf(stderr, "\nERROR when creating output directory %s\n", outputDirectory.make_preferred().string());
			exit(1);
		}
		cout << "\nDirectory: " << outputDirectory.make_preferred().string() << " is created." << endl;
	}
	else
		cout << "\nDirectory " << outputDirectory.make_preferred().string() << " already exist." << endl;

	unsigned int noOfTiles = 0;

	ostringstream outputFilenames;
	outputFilenames << WSIFileName << ".csv";
	path outputFilenamesWithPath(outputDirectory);
	outputFilenamesWithPath /= outputFilenames.str();
	std::ofstream outputFile(outputFilenamesWithPath.make_preferred().string().c_str());

	cout << "\nWSI Name: " << WSIFileNameWithExt
		<< "\nOutput File (CSV): " << outputFilenamesWithPath.make_preferred().string()
		<< "\nScan Magnification: " << scan_magnification
		<< "\nRequested Magnification: " << requested_magnification
		<< "\nBiased: " << Biased
		<< "\nframe_width: " << frame_width
		<< "\nframe_height: " << frame_height
		<< "\nBoundingBox X: " << roi_X
		<< "\nBoundingBox Y: " << roi_Y
		<< "\nBoundingBox Width: " << roi_width
		<< "\nBoundingBox Height: " << roi_height
		<< "\n\n";

	for (long y = roi_Y; y + frame_height <= roi_height; y += frame_height) {
		for (long x = roi_X; x + frame_width <= roi_width; x += frame_width) {
			vector<double> imageFeatures;
			ostringstream outputFilename;
			outputFilename << x << "_" << y << ".tiff";
			path outputFilenameWithPath(outputDirectory.make_preferred().string() / outputFilename.str());
			if (TileExtraction(WSIFileNameWithPath.make_preferred().string().c_str(), x, y, outputFilenameWithPath.make_preferred().string().c_str(), imageFeatures)){
				outputFile << outputFilenameWithPath.make_preferred().string() << "," << outputFilename.str() << "," << imageFeatures[0] << "," << imageFeatures[1] << "," << imageFeatures[2] << "," << imageFeatures[3] << endl;
				std::cout << "Writing tile: " << outputFilenameWithPath.make_preferred().string() << endl;
				noOfTiles++;
			}
		}
	}
	outputFile.close();
	return noOfTiles;
}

// Second Tile Extraction Code - Work properly for 40X, 20X scaned Magnification slides
bool TileExtraction(const char *inputFilename, long x, long y, const char *outputFilename, vector<double> &imageFeatures)
{
	string tmpString (outputFilename);
	ostringstream tmpString2;
	tmpString2 << tmpString.substr(0, tmpString.size() - 5) << "_tmp.tiff";
	char  tmpFilename[1000];
	sprintf(tmpFilename, "%s", tmpString2.str().c_str());

	// Open input file
	openslide_t *wsi = openslide_open(inputFilename);
	if (NULL == wsi) {
		fprintf(stderr, "\nERROR when opening file %s\n", inputFilename);
		return false;
	}

	// Get the scan magnification layer of WSI
	int32_t layer = openslide_get_best_layer_for_downsample(wsi, 1);
	if (-1 == layer) {
		fprintf(stderr, "\nERROR: Impossible to find an appropriate layer for scan magnification %lfX\n", scan_magnification);
		return false;
	}

	// define the buffer for reading data
	uint32_t *buff = (uint32_t*)malloc(frame_width * frame_height * sizeof(uint32_t));
	if (NULL == buff) {
		fprintf(stderr, "\nERROR: malloc failed for buff, not enough memory?\n");
		return false;
	}

	openslide_read_region(wsi, buff, x, y, layer, frame_width, frame_height);
	openslide_close(wsi);

	char *buffTiff = (char*)malloc(frame_width * frame_height * 3 * sizeof(char));
	if (NULL == buffTiff) {
		fprintf(stderr, "\nERROR: malloc failed for buffTiff, not enough memory?\n");
		return false;
	}

	// Transfer openslide buffer (type uint32_t, format ARGB)
	// into TIFF buffer (type char, format RGB)
	double redMean = 0, greenMean = 0, blueMean = 0, avgRGBPixel = 0;
	double TotalPixel = frame_width*frame_height;
	for (int i = 0; i < TotalPixel; i++) {
		buffTiff[i * 3] = (char)(buff[i] >> 16); // red
		buffTiff[i * 3 + 1] = (char)(buff[i] >> 8); // green
		buffTiff[i * 3 + 2] = (char)(buff[i]); // blue
		
		unsigned char cRed = buffTiff[i * 3];
		unsigned char cGreen = buffTiff[i * 3 + 1];
		unsigned char cBlue = buffTiff[i * 3 + 2];

		unsigned int redPixel = 0, greenPixel = 0, bluePixel = 0;
		redPixel = (unsigned int)cRed;
		greenPixel = (unsigned int)cGreen;
		bluePixel = (unsigned int)cBlue;

		redMean += redPixel;
		greenMean += greenPixel;
		blueMean += bluePixel;
	}
	redMean /= TotalPixel;
	greenMean /= TotalPixel;
	blueMean /= TotalPixel;
	avgRGBPixel = (redMean + greenMean + blueMean)/3;

	if ((redMean - avgRGBPixel) < Biased){
		free(buff);
		free(buffTiff);
		return false;
	}
	imageFeatures.push_back(redMean);
	imageFeatures.push_back(greenMean);
	imageFeatures.push_back(blueMean);
	imageFeatures.push_back(avgRGBPixel);

	cout << "\nX = " << x << "\tY = " << y << "\tWidth = " << frame_width << "\tHeight = " << frame_height 
		 << "\tRequested Magnification= " << requested_magnification << "\tDownsample = " << scan_magnification / requested_magnification;

	// Save region into a TIFF file
	TIFF *tiffFile = TIFFOpen(tmpFilename, "w");
	if (NULL == tiffFile) {
		fprintf(stderr, "\nERROR: Cannot open output TIFF file %s\n", tmpFilename);
		free(buff);
		free(buffTiff);
		return false;
	}

	// Write the TIFF tags to the file
	TIFFSetField(tiffFile, TIFFTAG_IMAGEWIDTH, frame_width);
	TIFFSetField(tiffFile, TIFFTAG_IMAGELENGTH, frame_height);
	TIFFSetField(tiffFile, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	TIFFSetField(tiffFile, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(tiffFile, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
	TIFFSetField(tiffFile, TIFFTAG_BITSPERSAMPLE, 8);
	TIFFSetField(tiffFile, TIFFTAG_SAMPLESPERPIXEL, 3);

	// Actually write the image
	if (0 == TIFFWriteEncodedStrip(tiffFile, 0, buffTiff, frame_width * frame_height * 3)) {
		fprintf(stderr, "\nERROR: Could not write TIFF image\n");
		free(buff);
		free(buffTiff);
		return false;
	}

	TIFFClose(tiffFile);
	free(buff);
	free(buffTiff);

	// Downsample image from available magnification down to magnification requested by the user
	char command[10000], command2[10000];
	if (scan_magnification == requested_magnification) {
#if !defined( _WIN32 )
		sprintf(command, "mv %s %s", tmpFilename.c_str(), outputFilename);
#else
		sprintf(command, "move %s %s", tmpFilename, outputFilename);
#endif
	}
	else {
#if !defined( _WIN32 )
		sprintf(command, "./downsample %s %lf %s %lf", tmpFilename.c_str(), scan_magnification, outputFilename, requested_magnification);
		sprintf(command2, "rm -f %s", tmpFilename);
#else
		sprintf(command, ".\\downsample.exe %s %lf %s %lf", tmpFilename, scan_magnification, outputFilename, requested_magnification);
		sprintf(command2, "del %s", tmpFilename);
#endif
		cout << "\nDownsampling from " << scan_magnification << " down to " << requested_magnification;
	}

	//printf("\n%s\n", command2);
	system(command);
	system(command2);
	return true;
}

// First Tile Extraction code - Work properly for 40X scaned magnification slides
bool TileExtraction_(const char *inputFilename, long x, long y, const char *outputFilename, vector<double> &imageFeatures)
{
	double downsample_size = scan_magnification / requested_magnification;
	string tmpString = outputFilename;
	char tmpFilename[1000];
	//sprintf(tmpFilename, "%s_tmp.tiff", tmpString.substr(0, tmpString.size() - 4) );
	sprintf(tmpFilename, "%s_tmp.tiff", outputFilename);

	// Open input file
	openslide_t *wsi = openslide_open(inputFilename);
	if (NULL == wsi) {
		fprintf(stderr, "\nERROR when opening file %s\n", inputFilename);
		return false;
	}

	// What is the best layer to use
	// for getting a portion of the WSI at magnification requested by the user
	int32_t layer = openslide_get_best_layer_for_downsample(wsi, downsample_size);
	if (-1 == layer) {
		fprintf(stderr, "\nERROR: Impossible to find an appropriate layer for requested magnification %lfX\n", requested_magnification);
		return false;
	}

	// Get dimensions of the available magnification
	int64_t widthLayer, heightLayer;
	openslide_get_layer_dimensions(wsi, layer, &widthLayer, &heightLayer);
	if ((-1 == widthLayer) || (-1 == heightLayer)) {
		fprintf(stderr, "\nERROR: Impossible to get the width and height of layer %d\n", layer);
		return false;
	}

	// Get the magnification of the selected layer
	double downsampleAvailable = openslide_get_layer_downsample(wsi, layer);
	if (-1.0 == downsampleAvailable) {
		fprintf(stderr, "\nERROR: Impossible to get the magnification of layer %d\n", layer);
		return false;
	}

	double magAvailable = scan_magnification / downsampleAvailable;

	cout << "\nMagnification Requested = " << requested_magnification << "\tDownsample = " << downsample_size << "\tDownsample Available = " << downsampleAvailable << "\tLayer = " << layer;
	cout << "\nMagnification at Scan   = " << magAvailable << "\tWSI_Width  = " << widthLayer << "\tWSI_Hieght = " << heightLayer;

	// Compute enlarged region needed at magnification available for selected layer
	// to be able, after downsampling from available magnification down to requested magnification,
	// to get an image of the exact size requested by the user
	int xEnlarged, yEnlarged;
	int64_t widthEnlarged, heightEnlarged;
	if (magAvailable == requested_magnification) {
		xEnlarged = x;
		yEnlarged = y;
		widthEnlarged = frame_width;
		heightEnlarged = frame_height;
	}
	else {
		printf("\nMagnification Available and Magnification Requested are different.");
		double tmp;

		tmp = (double)(x)+((double)(x)* magAvailable / requested_magnification) / 10.0;
		xEnlarged = (int)tmp;
		if ((double)(xEnlarged) < tmp) {
			// Some fractions after floating point have been lost when casting double to int
			// ==> Add 1 to compensate
			xEnlarged++;
		}

		tmp = (double)(y)+((double)(y)* magAvailable / requested_magnification) / 10.0;
		yEnlarged = (int)tmp;
		if ((double)(yEnlarged) < tmp) {
			// Some fractions after floating point have been lost when casting double to int
			// ==> Add 1 to compensate
			yEnlarged++;
		}

		tmp = (double)(frame_width)* magAvailable / requested_magnification;
		widthEnlarged = (int64_t)tmp;
		if ((double)(widthEnlarged) < tmp) {
			// Some fractions after floating point have been lost when casting double to int
			// ==> Add 1 to compensate
			widthEnlarged++;
		}

		tmp = (double)(frame_height)* magAvailable / requested_magnification;
		heightEnlarged = (int64_t)tmp;
		if ((double)(heightEnlarged) < tmp) {
			// Some fractions after floating point have been lost when casting double to int
			// ==> Add 1 to compensate
			heightEnlarged++;
		}
	}

	// Read data from WSI
	uint32_t *buff = (uint32_t*)malloc(widthEnlarged * heightEnlarged * sizeof(uint32_t));
	if (NULL == buff) {
		fprintf(stderr, "\nERROR: malloc failed for buff, not enough memory?\n");
		return false;
		//exit(1);
	}

	openslide_read_region(wsi, buff, xEnlarged, yEnlarged, layer, widthEnlarged, heightEnlarged);

	openslide_close(wsi);

	char *buffTiff = (char*)malloc(widthEnlarged * heightEnlarged * 3 * sizeof(char));
	if (NULL == buffTiff) {
		fprintf(stderr, "\nERROR: malloc failed for buffTiff, not enough memory?\n");
		return false;
		//exit(1);
	}

	// Transfer openslide buffer (type uint32_t, format ARGB)
	// into TIFF buffer (type char, format RGB)
	double redMean = 0, greenMean = 0, blueMean = 0, avgRGBPixel = 0;
	double TotalPixel = widthEnlarged*heightEnlarged;
	for (int i = 0; i < TotalPixel; i++) {
		buffTiff[i * 3] = (char)(buff[i] >> 16); // red
		buffTiff[i * 3 + 1] = (char)(buff[i] >> 8); // green
		buffTiff[i * 3 + 2] = (char)(buff[i]); // blue

		unsigned char cRed = buffTiff[i * 3];
		unsigned char cGreen = buffTiff[i * 3 + 1];
		unsigned char cBlue = buffTiff[i * 3 + 2];

		unsigned int redPixel = 0, greenPixel = 0, bluePixel = 0;
		redPixel = (unsigned int)cRed;
		greenPixel = (unsigned int)cGreen;
		bluePixel = (unsigned int)cBlue;

		redMean += redPixel;
		greenMean += greenPixel;
		blueMean += bluePixel;
	}
	redMean /= TotalPixel;
	greenMean /= TotalPixel;
	blueMean /= TotalPixel;
	avgRGBPixel = (redMean + greenMean + blueMean) / 3;

	if ((redMean - avgRGBPixel) < Biased){
		cout << "\nTile not selected. \nRed Mean: " << redMean << "\tAvgRGB: " << avgRGBPixel << "\tBiased: " << Biased;
		free(buff);
		free(buffTiff);
		return false;
	}
	imageFeatures.push_back(redMean);
	imageFeatures.push_back(greenMean);
	imageFeatures.push_back(blueMean);
	imageFeatures.push_back(avgRGBPixel);

	cout << "\nRequested Region, X = " << x << "\tY = " << y << "\tWidth = " << frame_width << "\tHeight = " << frame_height << "\tMag. Request = " << requested_magnification;
	cout << "\nModified Region, X = " << xEnlarged << "\tY = " << yEnlarged << "\tWidth = " << widthEnlarged << "\tHeight = " << heightEnlarged << "\tMag. Available = " << magAvailable;

	// Save region into a TIFF file
	TIFF *tiffFile = TIFFOpen(tmpFilename, "w");
	if (NULL == tiffFile) {
		fprintf(stderr, "\nERROR: Cannot open output TIFF file %s\n", tmpFilename);
		free(buff);
		free(buffTiff);
		exit(1);
	}

	// Write the TIFF tags to the file
	TIFFSetField(tiffFile, TIFFTAG_IMAGEWIDTH, widthEnlarged);
	TIFFSetField(tiffFile, TIFFTAG_IMAGELENGTH, heightEnlarged);
	TIFFSetField(tiffFile, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	TIFFSetField(tiffFile, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(tiffFile, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
	TIFFSetField(tiffFile, TIFFTAG_BITSPERSAMPLE, 8);
	TIFFSetField(tiffFile, TIFFTAG_SAMPLESPERPIXEL, 3);

	// Actually write the image
	if (0 == TIFFWriteEncodedStrip(tiffFile, 0, buffTiff, widthEnlarged * heightEnlarged * 3)) {
		fprintf(stderr, "\nERROR: Could not write TIFF image\n");
		free(buff);
		free(buffTiff);
		return false;
		//exit(1);
	}

	TIFFClose(tiffFile);
	free(buff);
	free(buffTiff);


	// Downsample image from available magnification down to magnification requested by the user
	char command[10000];
	if (downsample_size == downsampleAvailable) {
#if !defined( _WIN32 )
		sprintf(command, "mv %s %s", tmpFilename, outputFilename);
#else
		sprintf(command, "move %s %s", tmpFilename, outputFilename);
#endif
	}
	else {
#if !defined( _WIN32 )
		sprintf(command, "./downsample %s %lf %s %lf", tmpFilename, magAvailable, outputFilename, requested_magnification);
#else
		sprintf(command, ".\\downsample.exe %s %lf %s %lf", tmpFilename, magAvailable, outputFilename, requested_magnification);
#endif
		printf("\nDOWNSAMPLING from %lfX down to %lfX\n", magAvailable, requested_magnification);
	}

	printf("%s\n", command);
	system(command);
	return true;
}


void MyTest(void){
	path wsiFile("R:/Beck Lab/Atypia_Andy/HighAgreement/1305.tif");

	double downsample_size = (scan_magnification / requested_magnification);

	openslide_t *wsi = openslide_open(wsiFile.make_preferred().string().c_str());
	if (NULL == wsi) {
		fprintf(stderr, "\nERROR when opening file %s\n", wsiFile.make_preferred().string().c_str());
		exit(1);
	}

	// What is the best layer to use
	// for getting a portion of the WSI at magnification requested by the user
	int32_t layer = openslide_get_best_layer_for_downsample(wsi, downsample_size);
	if (-1 == layer) {
		fprintf(stderr, "\nERROR: Impossible to find an appropriate layer for requested magnification %lfX\n", requested_magnification);
		exit(1);
	}

	// Get dimensions of the available magnification
	int64_t widthLayer, heightLayer;
	openslide_get_layer_dimensions(wsi, layer, &widthLayer, &heightLayer);
	if ((-1 == widthLayer) || (-1 == heightLayer)) {
		fprintf(stderr, "\nERROR: Impossible to get the width and height of layer %d\n", layer);
		exit(1);
	}
	
	// Get the magnification of the selected layer
	double downsampleAvailable = openslide_get_layer_downsample(wsi, layer);
	if (-1.0 == downsampleAvailable) {
		fprintf(stderr, "\nERROR: Impossible to get the magnification of layer %d\n", layer);
		exit(1);
	}
	cout << "\nLayer: " << layer << "\tWidth = " << widthLayer << "\tHieght = " << heightLayer << "\tDownSample: " << downsampleAvailable;

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	layer = 0;
	openslide_get_layer_dimensions(wsi, layer, &widthLayer, &heightLayer);
	if ((-1 == widthLayer) || (-1 == heightLayer)) {
		fprintf(stderr, "\nERROR: Impossible to get the width and height of layer %d\n", layer);
	}
	downsampleAvailable = openslide_get_layer_downsample(wsi, layer);
	if (-1.0 == downsampleAvailable) {
		fprintf(stderr, "\nERROR: Impossible to get the magnification of layer %d\n", layer);
	}
	cout << "\nLayer: " << layer << "\tWidth = " << widthLayer << "\tHieght = " << heightLayer << "\tDownSample: " << downsampleAvailable;
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	layer = -1;
	openslide_get_layer_dimensions(wsi, layer, &widthLayer, &heightLayer);
	if ((-1 == widthLayer) || (-1 == heightLayer)) {
		fprintf(stderr, "\nERROR: Impossible to get the width and height of layer %d\n", layer);
	}
	downsampleAvailable = openslide_get_layer_downsample(wsi, layer);
	if (-1.0 == downsampleAvailable) {
		fprintf(stderr, "\nERROR: Impossible to get the magnification of layer %d\n", layer);
	}
	cout << "\nLayer: " << layer << "\tWidth = " << widthLayer << "\tHieght = " << heightLayer << "\tDownSample: " << downsampleAvailable;
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	layer = -2;
	openslide_get_layer_dimensions(wsi, layer, &widthLayer, &heightLayer);
	if ((-1 == widthLayer) || (-1 == heightLayer)) {
		fprintf(stderr, "\nERROR: Impossible to get the width and height of layer %d\n", layer);
	}
	downsampleAvailable = openslide_get_layer_downsample(wsi, layer);
	if (-1.0 == downsampleAvailable) {
		fprintf(stderr, "\nERROR: Impossible to get the magnification of layer %d\n", layer);
	}
	cout << "\nLayer: " << layer << "\tWidth = " << widthLayer << "\tHieght = " << heightLayer << "\tDownSample: " << downsampleAvailable;
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	layer = -3;
	openslide_get_layer_dimensions(wsi, layer, &widthLayer, &heightLayer);
	if ((-1 == widthLayer) || (-1 == heightLayer)) {
		fprintf(stderr, "\nERROR: Impossible to get the width and height of layer %d\n", layer);
	}
	downsampleAvailable = openslide_get_layer_downsample(wsi, layer);
	if (-1.0 == downsampleAvailable) {
		fprintf(stderr, "\nERROR: Impossible to get the magnification of layer %d\n", layer);
	}
	cout << "\nLayer: " << layer << "\tWidth = " << widthLayer << "\tHieght = " << heightLayer << "\tDownSample: " << downsampleAvailable;
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	layer = -4;
	openslide_get_layer_dimensions(wsi, layer, &widthLayer, &heightLayer);
	if ((-1 == widthLayer) || (-1 == heightLayer)) {
		fprintf(stderr, "\nERROR: Impossible to get the width and height of layer %d\n", layer);
	}
	downsampleAvailable = openslide_get_layer_downsample(wsi, layer);
	if (-1.0 == downsampleAvailable) {
		fprintf(stderr, "\nERROR: Impossible to get the magnification of layer %d\n", layer);
	}
	cout << "\nLayer: " << layer << "\tWidth = " << widthLayer << "\tHieght = " << heightLayer << "\tDownSample: " << downsampleAvailable;

	/*double magAvailable = scan_magnification / downsampleAvailable;
	cout << "\nMagnification Requested = " << requested_magnification << "\tDownsample_Size = " << downsample_size << "\tDownsample Available = " << downsampleAvailable << "\tLayer = " << layer;
	cout << "\nMagnification Available = " << magAvailable << "\tWidth = " << widthLayer << "\tHieght = " << heightLayer;
	*/
}