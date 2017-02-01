#ifndef __COLORCHANNELEXTRACTION_H__
#define __COLORCHANNELEXTRACTION_H__

#include "ITKDeclarations.h"

#include <iostream>

CharImagePointer ColorChannelExtraction(RGBImagePointer rgbImage, std::string outputColorSpace, int colorChannel);
int RegressionTestImage (std::string testImageFilename, std::string baselineImageFilename,int reportErrors, bool differences);

#endif // __COLORCHANNELEXTRACTION_H__