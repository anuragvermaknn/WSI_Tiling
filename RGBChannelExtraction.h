#if defined(_MSC_VER)
#pragma warning (disable : 4786)
#endif

#include <itkImageAdaptor.h>
#include "ITKFunctions.h"

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

CharImagePointer RedChannelExtraction(RGBImagePointer rgbImage )
{
 
  typedef itk::ImageAdaptor<RGBImageType,RedChannelPixelAccessor> ImageAdaptorType;
  ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
  adaptor->SetImage(rgbImage);

  typedef itk::RescaleIntensityImageFilter<ImageAdaptorType,CharImageType>  RescalerType;
  RescalerType::Pointer rescaler = RescalerType::New();
  rescaler->SetOutputMinimum(0);
  rescaler->SetOutputMaximum(255);
  rescaler->SetInput(adaptor);
  rescaler->Update();
   
  return rescaler->GetOutput();
 }
CharImagePointer BlueChannelExtraction(RGBImagePointer rgbImage )
{
 
  typedef itk::ImageAdaptor<RGBImageType,BlueChannelPixelAccessor> ImageAdaptorType;
  ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
  adaptor->SetImage(rgbImage);

  typedef itk::RescaleIntensityImageFilter<ImageAdaptorType,CharImageType>  RescalerType;
  RescalerType::Pointer rescaler = RescalerType::New();
  rescaler->SetOutputMinimum(0);
  rescaler->SetOutputMaximum(255);
  rescaler->SetInput(adaptor);
  rescaler->Update();
   
  return rescaler->GetOutput();
 }
CharImagePointer GreenChannelExtraction(RGBImagePointer rgbImage )
{
 
  typedef itk::ImageAdaptor<RGBImageType,GreenChannelPixelAccessor> ImageAdaptorType;
  ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
  adaptor->SetImage(rgbImage);

  typedef itk::RescaleIntensityImageFilter<ImageAdaptorType,CharImageType>  RescalerType;
  RescalerType::Pointer rescaler = RescalerType::New();
  rescaler->SetOutputMinimum(0);
  rescaler->SetOutputMaximum(255);
  rescaler->SetInput(adaptor);
  rescaler->Update();
   
  return rescaler->GetOutput();
 }
CharImagePointer BlueRatioExtraction(RGBImagePointer rgbImage )
{
 
  typedef itk::ImageAdaptor< RGBImageType, BlueRatioPixelAccessor > ImageAdaptorType;
  ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
  adaptor->SetImage(rgbImage);

  typedef itk::RescaleIntensityImageFilter<ImageAdaptorType,CharImageType>  RescalerType;
  RescalerType::Pointer rescaler = RescalerType::New();
  rescaler->SetOutputMinimum(0);
  rescaler->SetOutputMaximum(255);
  rescaler->SetInput(adaptor);
  rescaler->Update();
   
  return rescaler->GetOutput();
 }
CharImagePointer BlueRatio( RGBImagePointer rgbImage )
{
	RGBImageSizeType size = rgbImage->GetLargestPossibleRegion().GetSize();
	CharImagePointer charImage = CharImageType::New();
	CreateImage(charImage, size);
 
	RGBImageRegionType region = rgbImage->GetLargestPossibleRegion();
	RGBImageLinearConstIteratorWithIndexType it (rgbImage, region);
	CharImageLinearIteratorWithIndexType it2(charImage, region);

	for(it.GoToBegin(); !it.IsAtEnd(); it.NextLine())	
	{
		while(!it.IsAtEndOfLine())
		{
			CharImageIndexType index = it.GetIndex();
			RGBPixelType rgbPixel = it.Get();
			CharPixelType charPixel = ( 100 * rgbPixel.GetBlue() ) / ( 1 + rgbPixel.GetRed() + rgbPixel.GetGreen() ) *
				( 256 / ( 1 + rgbPixel.GetRed() + rgbPixel.GetGreen() + rgbPixel.GetBlue() ) ) ;
			charImage->SetPixel( index, charPixel );
			++it;
		}
	}

	CharToCharRescaleIntensityImageFilter::Pointer rescaler = CharToCharRescaleIntensityImageFilter::New();
	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255);
	rescaler->SetInput(charImage);
	rescaler->Update();

	return rescaler->GetOutput();
}
void ColorDeconvolution( RGBImagePointer rgbImage, void *RImage, void *GImage, void *BImage, unsigned int matrixChoice = 1 )
{
	CharImagePointer *rImage = (CharImagePointer *)RImage;
	CharImagePointer *gImage = (CharImagePointer *)GImage;
	CharImagePointer *bImage = (CharImagePointer *)BImage;
	
	CharImagePointer redImage = RedChannelExtraction( rgbImage );
	CharImagePointer greenImage	= GreenChannelExtraction( rgbImage );
	CharImagePointer blueImage	= BlueChannelExtraction( rgbImage );

	double MODx[3];
	double MODy[3];
	double MODz[3];
	switch( matrixChoice )
	{
	case 1:			// Antoine PSL images
		MODx[0]= 0.8249;	// GL H matrix
		MODy[0]= 0.8040;
		MODz[0]= 0.5969;
		MODx[1]= 0.0825;	// GL E matrix 
		MODy[1]= 0.5454;
		MODz[1]= 0.3680;
		MODx[2]= 0.0;		//  Zero matrix 
		MODy[2]= 0.0;
		MODz[2]= 0.0;
		break;
	case 2:			// Antoine NUH images
		MODx[0]= 0.650;		// GL H matrix
		MODy[0]= 0.704;
		MODz[0]= 0.286;
		MODx[1]= 0.072;		// GL E matrix 
		MODy[1]= 0.990;
		MODz[1]= 0.105;
		MODx[2]= 0.0;		//  Zero matrix 
		MODy[2]= 0.0;
		MODz[2]= 0.0;
		break;
	case 3:			// Rutifork matrix
		MODx[0]= 0.644211;	// GL H matrix
		MODy[0]= 0.716556;
		MODz[0]= 0.266844;
		MODx[1]= 0.092789;	// GL E matrix
		MODy[1]= 0.954111;
		MODz[1]= 0.283111;
		MODx[2]= 0.0;		//  Zero matrix 
		MODy[2]= 0.0;
		MODz[2]= 0.0;
		break;
	default:
		MODx[0]= 0.8249;	// GL H matrix
		MODy[0]= 0.8040;
		MODz[0]= 0.5969;
		MODx[1]= 0.0825;	// GL E matrix 
		MODy[1]= 0.5454;
		MODz[1]= 0.3680;
		MODx[2]= 0.0;		//  Zero matrix 
		MODy[2]= 0.0;
		MODz[2]= 0.0;
		break;
	}	

	/**********************************************************************
	*  Convert the stain vectors into the 'q' vector
	***********************************************************************/
	double leng, A, V, C, log255=log(255.0);
	double cosx[3];
	double cosy[3];
	double cosz[3];
	double len[3];
	double q[9];
	
	for (int i=0; i<3; i++){
		/* normalise vector length */
		cosx[i]=0.0;
		cosy[i]=0.0;
		cosz[i]=0.0;
		len[i]=sqrt(MODx[i]*MODx[i] + MODy[i]*MODy[i] + MODz[i]*MODz[i]);
		if (len[i] != 0.0){
			cosx[i]= MODx[i]/len[i];
			cosy[i]= MODy[i]/len[i];
			cosz[i]= MODz[i]/len[i];
		}
	}

	/* translation matrix */
	if (cosx[1]==0.0){ /* 2nd colour is unspecified */
		if (cosy[1]==0.0){
			if (cosz[1]==0.0){
				cosx[1]=cosz[0];
				cosy[1]=cosx[0];
				cosz[1]=cosy[0];
			}	
		}
	}

	if (cosx[2]==0.0){ /* 3rd colour is unspecified */
		if (cosy[2]==0.0){
			if (cosz[2]==0.0){
				if ((cosx[0]*cosx[0] + cosx[1]*cosx[1])> 1)
					cosx[2]=0.0;
				else
					cosx[2]=sqrt(1.0-(cosx[0]*cosx[0])-(cosx[1]*cosx[1]));

				if ((cosy[0]*cosy[0] + cosy[1]*cosy[1])> 1)
					cosy[2]=0.0;
				else
					cosy[2]=sqrt(1.0-(cosy[0]*cosy[0])-(cosy[1]*cosy[1]));

				if ((cosz[0]*cosz[0] + cosz[1]*cosz[1])> 1)
					cosz[2]=0.0;
				else
					cosz[2]=sqrt(1.0-(cosz[0]*cosz[0])-(cosz[1]*cosz[1]));
			}
		}
	}

	leng= sqrt(cosx[2]*cosx[2] + cosy[2]*cosy[2] + cosz[2]*cosz[2]);

	cosx[2]= cosx[2]/leng;
	cosy[2]= cosy[2]/leng;
	cosz[2]= cosz[2]/leng;

	/* matrix inversion */
	A = cosy[1] - cosx[1] * cosy[0] / cosx[0];
	V = cosz[1] - cosx[1] * cosz[0] / cosx[0];
	C = cosz[2] - cosy[2] * V/A + cosx[2] * (V/A * cosy[0] / cosx[0] - cosz[0] / cosx[0]);
	q[2] = (-cosx[2] / cosx[0] - cosx[2] / A * cosx[1] / cosx[0] * cosy[0] / cosx[0] + cosy[2] / A * cosx[1] / cosx[0]) / C;
	q[1] = -q[2] * V / A - cosx[1] / (cosx[0] * A);
	q[0] = 1.0 / cosx[0] - q[1] * cosy[0] / cosx[0] - q[2] * cosz[0] / cosx[0];
	q[5] = (-cosy[2] / A + cosx[2] / A * cosy[0] / cosx[0]) / C;
	q[4] = -q[5] * V / A + 1.0 / A;
	q[3] = -q[4] * cosy[0] / cosx[0] - q[5] * cosz[0] / cosx[0];
	q[8] = 1.0 / C;
	q[7] = -q[8] * V / A;
	q[6] = -q[7] * cosy[0] / cosx[0] - q[8] * cosz[0] / cosx[0];

	
	/************************************************************************
	* Apply the 'q' vector to the original RGB image to make some new stain images
	*************************************************************************/
		
	RGBImageType::SizeType imageSize = rgbImage->GetLargestPossibleRegion().GetSize();

	for (int i=0; i<imageSize[0]; i++)
	{
		for(int j=0; j<imageSize[1]; j++)
		{
			/* log transform the RGB data */
			CharImageIndexType index;
			index[0] = i;
			index[1] = j;
			
			CharImagePixelType R = redImage->GetPixel( index );
			CharImagePixelType G = greenImage->GetPixel( index );
			CharImagePixelType B = blueImage->GetPixel( index );
	
			double Rlog = -((255.0*log(((double)R+1)/255.0))/log255);
			double Glog = -((255.0*log(((double)G+1)/255.0))/log255);
			double Blog = -((255.0*log(((double)B+1)/255.0))/log255);

			for (int k=0; k<3; k++){
				/* rescale to match original paper values */
				double Rscaled = Rlog * q[k*3];
				double Gscaled = Glog * q[k*3+1];
				double Bscaled = Blog * q[k*3+2];

				double output = exp(-((Rscaled + Gscaled + Bscaled) - 255.0) * log255 / 255.0);
				if(output>255) 
					output=255;
				 
				if (k==0) {
					redImage->SetPixel( index, floor(output+.5) );
				} else if (k==1) {
					greenImage->SetPixel( index, floor(output+.5) );
				} else {
					blueImage->SetPixel( index, floor(output+.5) );
				}
			}
		}
	}
	*rImage = redImage;
	*gImage = greenImage;
	*bImage = blueImage;
}
CharImagePointer HematoxylinChannelExtraction( RGBImagePointer rgbImage, int matrixType = 1)
{
	CharImagePointer hImage, eImage, rImage;
	void *HImage = &hImage;		// Hematazlin
	void *EImage = &eImage;		// Eosin
	void *RImage = &rImage;		// Resdi
	
	ColorDeconvolution( rgbImage, HImage, EImage, RImage, matrixType );
	return hImage;
}
CharImagePointer EosinChannelExtraction( RGBImagePointer rgbImage, int matrixType = 1)
{
	CharImagePointer hImage, eImage, rImage;
	void *HImage = &hImage;		// Hematazlin
	void *EImage = &eImage;		// Eosin
	void *RImage = &rImage;		// Resdi
	
	ColorDeconvolution( rgbImage, HImage, EImage, RImage, matrixType );
	return eImage;
}

