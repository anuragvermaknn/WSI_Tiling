/*
 * RGBChannelExtraction.h
 *
 *  Created on: 20 May 2012
 *      Author: humayun
 */
#ifndef __RGBCHANNELEXTRACTION_H__
#define __RGBCHANNELEXTRACTION_H__

#include "ITKDeclarations.h"

class RedChannelPixelAccessor
{
public:
  typedef RGBPixelType InternalType;
  typedef CharPixelType  ExternalType;
 
  static ExternalType Get(const InternalType& input)
  {
    return static_cast<ExternalType>(input.GetRed());
  }
};

class GreenChannelPixelAccessor
{
public:
  typedef RGBPixelType InternalType;
  typedef CharPixelType  ExternalType;
 
  static ExternalType Get(const InternalType& input)
  {
     return static_cast<ExternalType>(input.GetGreen());
  }
};

class BlueChannelPixelAccessor
{
public:
  typedef RGBPixelType InternalType;
  typedef CharPixelType  ExternalType;
 
  static ExternalType Get(const InternalType& input)
  {
     return static_cast<ExternalType>(input.GetBlue());
  }
}; 

class BlueRatioPixelAccessor
{
public:
  typedef RGBPixelType InternalType;
  typedef CharPixelType  ExternalType;
 
  static ExternalType Get(const InternalType& input)
  {
     return static_cast<ExternalType>( ( 100 * input.GetBlue() ) / ( 1 + input.GetRed() + input.GetGreen() ) *
				( 256 / ( 1 + input.GetRed() + input.GetGreen() + input.GetBlue() ) ) );
  }
}; 

CharImagePointer RedChannelExtraction(RGBImagePointer rgbImage);

CharImagePointer BlueChannelExtraction(RGBImagePointer rgbImage);

CharImagePointer GreenChannelExtraction(RGBImagePointer rgbImage);

CharImagePointer BlueRatioExtraction(RGBImagePointer rgbImage);

CharImagePointer BlueRatio( RGBImagePointer rgbImage);

CharImagePointer HematoxylinChannelExtraction( RGBImagePointer rgbImage, int matrixType = 1 );

CharImagePointer EosinChannelExtraction( RGBImagePointer rgbImage, int matrixType = 1 );

#endif // __RGBCHANNELEXTRACTION_H__