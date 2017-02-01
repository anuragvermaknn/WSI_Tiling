/************************************ Basic Filters and Functions **********************************************/
void CreateImage(CharImagePointer inImage, CharImageSizeType size)
{
	CharImageRegionType region;
	CharImageIndexType start;
	start[0] = 0;
	start[1] = 0;

	region.SetSize(size);
	region.SetIndex(start);

	inImage->SetRegions(region);
	inImage->Allocate();

	//// Make a square
	//for(unsigned int r = 0; r < size[0]; r++)
	//{
	//	for(unsigned int c = 0; c < size[1]; c++)
	//	{
	//		CharImageIndexType pixelIndex;
	//		pixelIndex[0] = r;
	//		pixelIndex[1] = c;
	//		inImage->SetPixel(pixelIndex, 0);		// 255 = White image, 0 = Black image
	//	}
	//}
}

#include "itkCastImageFilter.h"
typedef itk::CastImageFilter< CharImageType, FloatImageType >			CharToFloatCastImageFilterType;
typedef itk::CastImageFilter< FloatImageType, CharImageType >			FloatToCharCastImageFilterType;
typedef itk::CastImageFilter< RGBImageType, RGBImageType >				RGBToRGBCastImageFilterType;
template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP CastFilter( TInImage *inImage )
{
	typedef itk::CastImageFilter< TInImage, TOutImage> CastImageFilterType;
	typedef typename CastImageFilterType::Pointer	CastImageFilterPointer;
	CastImageFilterPointer cast = CastImageFilterType::New();
	cast->SetInput( inImage );
	cast->Update();
	return cast->GetOutput();
}

#include "itkSubtractImageFilter.h"
typedef itk::SubtractImageFilter< CharImageType, CharImageType >		CharSubtractImageFilterType;
typedef itk::SubtractImageFilter< FloatImageType, FloatImageType >		FloatSubtractImageFilterType;
template <class TImage, class TOutImageP>
TOutImageP SubtractFilter( TImage * inImage1, TImage * inImage2 )
{
	typedef itk::SubtractImageFilter< TImage, TImage >	SubtractImageFilterType;
	typedef typename SubtractImageFilterType::Pointer SubtractImageFilterPointer;
	SubtractImageFilterPointer subtract = SubtractImageFilterType::New();
	subtract->SetInput1( inImage1 );
	subtract->SetInput2( inImage2 );
	subtract->Update();
	return subtract->GetOutput();
}

#include "itkImageToHistogramFilter.h"
template <class TImage>
void ComputeHistogram(TImage * inImage, std::string histoFileNameCString )
{
	typedef itk::Statistics::ImageToHistogramFilter< TImage >				HistogramFilterType;
	typedef typename HistogramFilterType::Pointer							HistogramFilterPointer;
	typedef typename HistogramFilterType::HistogramMeasurementVectorType	HistogramMeasurementVectorType;
	typedef typename HistogramFilterType::HistogramSizeType					HistogramSizeType;
	typedef typename HistogramFilterType::HistogramType						HistogramType;
	typedef typename HistogramType::ConstIterator							HistogramConstIterator;

	HistogramSizeType histogramSize( 3 );
	histogramSize[0] = 2;  // number of bins for the Red   channel
	histogramSize[1] = 2;  // number of bins for the Green channel
	histogramSize[2] = 2;  // number of bins for the Blue  channel
 
	HistogramFilterPointer filter = HistogramFilterType::New();
	filter->SetInput( inImage );
	filter->SetAutoMinimumMaximum( true );
	filter->SetHistogramSize( histogramSize );
	filter->SetMarginalScale( 10 ); // Required (could this be set in the filter?)
	filter->Update();
 
	const HistogramType * histogram = filter->GetOutput();
 
	HistogramConstIterator histogramIterator = histogram->Begin();
 
	while( histogramIterator  != histogram->End() )
	{
	std::cout << "Index = " << histogram->GetIndex(histogramIterator.GetMeasurementVector())
			<< " Histogram cell center = " << histogramIterator.GetMeasurementVector() 
			<< " Frequency = " << histogramIterator.GetFrequency() << std::endl;
 
	++histogramIterator ;
	}
 
	HistogramMeasurementVectorType mv(3);
	mv[0] = 255;
	mv[1] = 0;
	mv[2] = 0;
	std::cout << "Frequency = " << histogram->GetFrequency(histogram->GetIndex(mv)) << std::endl;
}

#include "itkMaskedImageToHistogramFilter.h"
template <class TImage>
void ComputeHistogram(TImage * inImage, CharImagePointer mask, std::string histoFileNameCString )
{
	typedef itk::Statistics::MaskedImageToHistogramFilter< TImage, CharImageType >	MaskedHistogramFilterType;
	typedef typename MaskedHistogramFilterType::Pointer								MaskedHistogramFilterPointer;
	typedef typename MaskedHistogramFilterType::HistogramMeasurementVectorType		MaskedHistogramMeasurementVectorType;
	typedef typename MaskedHistogramFilterType::HistogramSizeType					MaskedHistogramSizeType;
	typedef typename MaskedHistogramFilterType::HistogramType						MaskedHistogramType;
	typedef typename MaskedHistogramType::ConstIterator								MaskedHistogramConstIterator;

	MaskedHistogramSizeType histogramSize( 3 );
	histogramSize[0] = 2;  // number of bins for the Red   channel
	histogramSize[1] = 2;  // number of bins for the Green channel
	histogramSize[2] = 2;  // number of bins for the Blue  channel
 
	MaskedHistogramFilterPointer filter = MaskedHistogramFilterType::New();
	filter->SetInput( inImage );
	filter->SetMaskImage( mask );
	filter->SetAutoMinimumMaximum( true );
	filter->SetHistogramSize( histogramSize );
	filter->SetMarginalScale( 10 ); // Required (could this be set in the filter?)
	filter->Update();
 
	const MaskedHistogramType * histogram = filter->GetOutput();
 
	MaskedHistogramConstIterator histogramIterator = histogram->Begin();
 
	while( histogramIterator  != histogram->End() )
	{
	std::cout << "Index = " << histogram->GetIndex(histogramIterator.GetMeasurementVector())
			<< " Histogram cell center = " << histogramIterator.GetMeasurementVector() 
			<< " Frequency = " << histogramIterator.GetFrequency() << std::endl;
 
	++histogramIterator ;
	}
 
	MaskedHistogramMeasurementVectorType mv(3);
	mv[0] = 255;
	mv[1] = 0;
	mv[2] = 0;
	std::cout << "Frequency = " << histogram->GetFrequency(histogram->GetIndex(mv)) << std::endl;
}

#include "itkShrinkImageFilter.h"
typedef itk::ShrinkImageFilter< CharImageType, CharImageType >		CharShrinkImageFilterType;
typedef itk::ShrinkImageFilter< FloatImageType, FloatImageType >	FloatShrinkImageFilterType;
template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP ShrinkingImage( TInImage * inImage )
{
	typedef itk::ShrinkImageFilter< TInImage, TOutImage > ShrinkImageFilterType;
	typedef typename ShrinkImageFilterType::Pointer	ShrinkImageFilterPointer;
	ShrinkImageFilterPointer shrinkFilter = ShrinkImageFilterType::New();
	shrinkFilter->SetInput( inImage );
	shrinkFilter->SetShrinkFactor( 0, 2);
	shrinkFilter->SetShrinkFactor( 1, 2);
	shrinkFilter->Update();

	std::cout << "\n\nShrink ..." 
		<< "\nOld Size: " << inImage->GetLargestPossibleRegion().GetSize() 
		<< "\tSpacing: " << inImage->GetSpacing()
		<< "\nNew Size: " << shrinkFilter->GetOutput()->GetLargestPossibleRegion().GetSize() 
		<< "\tSpacing: " << shrinkFilter->GetOutput()->GetSpacing();

	return shrinkFilter->GetOutput();
}

#include "itkInvertIntensityImageFilter.h"
typedef itk::InvertIntensityImageFilter< CharImageType >			CharInvertIntensityImageFilterType;
typedef itk::InvertIntensityImageFilter< FloatImageType >			FloatInvertIntensityImageFilterType;
template <class TImage, class TOutImageP>
TOutImageP InvertIntensityFilter( TImage* inImage, int forePixel = 255 )
{
	typedef itk::InvertIntensityImageFilter< TImage > InvertIntensityImageFilterType;
	typedef typename InvertIntensityImageFilterType::Pointer InvertIntensityImageFilterPointer;
	InvertIntensityImageFilterPointer filter = InvertIntensityImageFilterType::New();
	filter->SetInput( inImage );
	filter->SetMaximum( forePixel );
	filter->Update();
	return filter->GetOutput();
}

/**************************** Image Enhancmenet, Smoothing and Denoisification ************************************/

#include "itkMeanImageFilter.h"
typedef itk::MeanImageFilter< CharImageType, CharImageType >		CharMeanImageFilterType;
typedef itk::MeanImageFilter< FloatImageType, FloatImageType >		FloatMeanImageFilterType;
template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP MeanFilter( TInImage* inImage, double radius = 2 )
{
	typedef itk::MeanImageFilter< TInImage, TOutImage >	MeanImageFilterType;
	typedef typename MeanImageFilterType::Pointer	MeanImageFilterPointer;
	CharImageSizeType meanRadius;
	meanRadius.Fill( radius );
	MeanImageFilterPointer smoothing = MeanImageFilterType::New();
	smoothing->SetRadius( meanRadius );
	smoothing->SetInput( inImage );
	smoothing->Update();
	return smoothing->GetOutput();
}

#include "itkMedianImageFilter.h"
typedef itk::MedianImageFilter< CharImageType, CharImageType >		CharMedianImageFilterType;
typedef itk::MedianImageFilter< FloatImageType, FloatImageType >	FloatMedianImageFilterType;
template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP MedianFilter( TInImage* inImage, double radius = 2 )
{
	typedef itk::MedianImageFilter< TInImage, TOutImage >	MedianImageFilterType;
	typedef typename MedianImageFilterType::Pointer MedianImageFilterPointer;
	CharImageSizeType medianRadius;
	medianRadius.Fill( radius );
	MedianImageFilterPointer smoothing = MedianImageFilterType::New();
	smoothing->SetRadius( medianRadius );
	smoothing->SetInput( inImage );
	smoothing->Update();
	return smoothing->GetOutput();
}

#include "itkBinaryMedianImageFilter.h"
typedef itk::BinaryMedianImageFilter< CharImageType, CharImageType > BinaryMedianImageFilterType;
CharImagePointer BinaryMedianFilter( CharImagePointer inImage, double radius = 2 )
{
	CharImageSizeType medianRadius;
	medianRadius.Fill( radius );
	BinaryMedianImageFilterType::Pointer smoothing = BinaryMedianImageFilterType::New();
	smoothing->SetRadius( medianRadius );
	smoothing->SetInput( inImage );
	smoothing->Update();
	return smoothing->GetOutput();
}

#include "itkBilateralImageFilter.h"
typedef itk::BilateralImageFilter< CharImageType, CharImageType >	CharBilateralImageFilterType;
typedef itk::BilateralImageFilter< FloatImageType, FloatImageType >	FloatBilateralImageFilterType;
template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP BilateralFilter(TInImage* inImage, double domainSigma = 1, double rangeSigma = 1 )
{
	typedef itk::BilateralImageFilter< TInImage, TOutImage > BilateralImageFilterType;
	typedef typename BilateralImageFilterType::Pointer BilateralImageFilterPointer;
	BilateralImageFilterPointer smoothing = BilateralImageFilterType::New();
	smoothing->SetInput( inImage );
	smoothing->SetDomainSigma( domainSigma );
	smoothing->SetRangeSigma( rangeSigma );
	smoothing->Update();
	return smoothing->GetOutput();
}

#include "itkCurvatureFlowImageFilter.h"
typedef itk::CurvatureFlowImageFilter< CharImageType, CharImageType >	CharCurvatureFlowImageFilterType;
typedef itk::CurvatureFlowImageFilter< FloatImageType, FloatImageType >	FloatCurvatureFlowImageFilterType;
template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP CurvatureFlowFilter(TInImage* inImage )
{
	typedef itk::CurvatureFlowImageFilter< TInImage, TOutImage > CurvatureFlowImageFilterType;
	typedef typename CurvatureFlowImageFilterType::Pointer CurvatureFlowImageFilterPointer;
	CurvatureFlowImageFilterPointer smoothing = CurvatureFlowImageFilterType::New();
	smoothing->SetInput( inImage );
	smoothing->SetNumberOfIterations( 5 );
	smoothing->SetTimeStep( 0.125 );
	smoothing->Update();
	return smoothing->GetOutput();
}

#include "itkDiscreteGaussianImageFilter.h"
typedef itk::DiscreteGaussianImageFilter< CharImageType, CharImageType >	CharDiscreteGaussianImageFilterType;
template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP DiscreateGaussianFilter(TInImage* inImage, double variance = 1 )
{
	typedef itk::DiscreteGaussianImageFilter< TInImage, TOutImage >	DiscreteGaussianImageFilterType;
	typedef typename DiscreteGaussianImageFilterType::Pointer DiscreteGaussianImageFilterPointer;
	DiscreteGaussianImageFilterPointer smoothing = DiscreteGaussianImageFilterType::New();
	smoothing->SetInput( inImage );
	smoothing->SetVariance( variance );
	smoothing->Update();
	return smoothing->GetOutput();
}

#include "itkSmoothingRecursiveGaussianImageFilter.h"
typedef itk::SmoothingRecursiveGaussianImageFilter< CharImageType, CharImageType >	CharSmoothingRecursiveGaussianImageFilterType;
template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP SmoothingRecursiveGaussianFilter(TInImage* inImage, double sigma = 1 )
{
	typedef itk::SmoothingRecursiveGaussianImageFilter< TInImage, TOutImage >	SmoothingRecursiveGaussianImageFilterType;
	typedef typename SmoothingRecursiveGaussianImageFilterType::Pointer SmoothingRecursiveGaussianImageFilterPointer;
	SmoothingRecursiveGaussianImageFilterPointer smoothing = SmoothingRecursiveGaussianImageFilterType::New();
	smoothing->SetInput( inImage );
	smoothing->SetSigma( sigma );
	smoothing->Update();
	return smoothing->GetOutput();
}

#include "itkAnisotropicDiffusionImageFilter.h"
typedef itk::AnisotropicDiffusionImageFilter< CharImageType, CharImageType >			CharAnisotropicDiffusionImageFilterType;
typedef itk::AnisotropicDiffusionImageFilter< FloatImageType, FloatImageType >			FloatAnisotropicDiffusionImageFilterType;

#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
typedef itk::CurvatureAnisotropicDiffusionImageFilter< CharImageType, CharImageType >	CharCurvatureAnisotropicDiffusionImageFilterType;
typedef itk::CurvatureAnisotropicDiffusionImageFilter< FloatImageType, FloatImageType >	FloatCurvatureAnisotropicDiffusionImageFilterType;

template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP CurvatureAnisotropicDiffusionFilter( TInImage* inImage, double conductance = 9.0 )
{
	typedef itk::CurvatureAnisotropicDiffusionImageFilter< TInImage, TOutImage >	CurvatureAnisotropicDiffusionImageFilterType;
	typedef typename CurvatureAnisotropicDiffusionImageFilterType::Pointer CurvatureAnisotropicDiffusionImageFilterPointer;
	CurvatureAnisotropicDiffusionImageFilterPointer smoothing = CurvatureAnisotropicDiffusionImageFilterType::New();
	smoothing->SetTimeStep( 0.125 );
	smoothing->SetNumberOfIterations(  5 );
	smoothing->SetConductanceParameter( conductance );
	smoothing->SetInput( inImage );
	smoothing->Update();
	return smoothing->GetOutput();
}

/************************ Image Gradient, Edge Detection, Derivatives, Speed, Height ***********************************/
#include "itkRecursiveGaussianImageFilter.h"
typedef itk::RecursiveGaussianImageFilter< CharImageType, CharImageType >	CharRecursiveGaussianImageFilterType;
template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP RecursiveGaussianFilter(TInImage* inImage, double direction = 0 )
{
	typedef itk::RecursiveGaussianImageFilter< TInImage, TOutImage > RecursiveGaussianImageFilterType;
	typedef typename RecursiveGaussianImageFilterType::Pointer RecursiveGaussianImageFilterPointer;
	RecursiveGaussianImageFilterPointer smoothing = RecursiveGaussianImageFilterType::New();
	smoothing->SetInput( inImage );
	smoothing->SetDirection( direction );		// 0 = x-axis 
	smoothing->SetSecondOrder();
	smoothing->Update();
	return smoothing->GetOutput();
}

#include "itkGradientAnisotropicDiffusionImageFilter.h"
typedef itk::GradientAnisotropicDiffusionImageFilter< CharImageType, CharImageType >	CharGradientAnisotropicDiffusionImageFilterType;
template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP GradientAnisotropicDiffusionFilter(TInImage* inImage, double conductance = 1.5 )
{
	typedef itk::GradientAnisotropicDiffusionImageFilter< TInImage, TOutImage >	GradientAnisotropicDiffusionImageFilterType;
	typedef typename GradientAnisotropicDiffusionImageFilterType::Pointer GradientAnisotropicDiffusionImageFilterPointer;
	GradientAnisotropicDiffusionImageFilterPointer smoothing = GradientAnisotropicDiffusionImageFilterType::New();
	smoothing->SetInput( inImage );
	smoothing->SetConductanceParameter( conductance );
	smoothing->SetNumberOfIterations( 5 );
	smoothing->SetTimeStep( 0.125 );
	smoothing->Update();
	return smoothing->GetOutput();
}

#include "itkVectorGradientAnisotropicDiffusionImageFilter.h"
typedef itk::VectorGradientAnisotropicDiffusionImageFilter< RGBImageType, RGBImageType > VectorGradientAnisotropicDiffusionImageFilterType;

#include "itkLaplacianSharpeningImageFilter.h"
typedef itk::LaplacianSharpeningImageFilter< FloatImageType, FloatImageType >		FloatLaplacianSharpeningImageFilterType;
template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP LaplacianSharpeningFilter( TInImage* inImage )
{
	typedef itk::LaplacianSharpeningImageFilter< TInImage, TOutImage >		LaplacianSharpeningImageFilterType;
	typedef typename LaplacianSharpeningImageFilterType::Pointer LaplacianSharpeningImageFilterPointer;
	LaplacianSharpeningImageFilterPointer sharpening = LaplacianSharpeningImageFilterType::New();
	sharpening->SetInput( inImage );
	sharpening->Update();
	return sharpening->GetOutput();
}

#include "itkLaplacianRecursiveGaussianImageFilter.h"
typedef itk::LaplacianRecursiveGaussianImageFilter< CharImageType, CharImageType >	CharLaplacianRecursiveGaussianImageFilterType;
template <class TInImage, class TOutImage, class TOutImageP>
TOutImageP LaplacianRecursiveGaussianFilter( TInImage* inImage )
{
	typedef itk::LaplacianRecursiveGaussianImageFilter< TInImage, TOutImage >	LaplacianRecursiveGaussianImageFilterType;
	typedef typename LaplacianRecursiveGaussianImageFilterType::Pointer LaplacianRecursiveGaussianImageFilterPointer;
	LaplacianRecursiveGaussianImageFilterPointer LoGFilter = LaplacianRecursiveGaussianImageFilterType::New();
	LoGFilter->SetInput( inImage );
	LoGFilter->Update();
	return LoGFilter->GetOutput();
}

// Speed
#include "itkSigmoidImageFilter.h"
typedef itk::SigmoidImageFilter< FloatImageType, FloatImageType >					FloatSigmoidImageFilterType;

#include "itkGradientMagnitudeImageFilter.h"
typedef itk::GradientMagnitudeImageFilter< FloatImageType, FloatImageType >			FloatGradientMagnitudeImageFilterType;

#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
typedef itk::GradientMagnitudeRecursiveGaussianImageFilter< CharImageType, CharImageType >	CharGradientMagnitudeRecursiveGaussianImageFilterType;
typedef itk::GradientMagnitudeRecursiveGaussianImageFilter< FloatImageType, FloatImageType > FloatGradientMagnitudeRecursiveGaussianImageFilterType;

/***************************** Segmentation **********************************/
#include "itkScalarConnectedComponentImageFilter.h"
typedef itk::ScalarConnectedComponentImageFilter< CharImageType, CharImageType >	ScalarConnectedComponentImageFilterType;
CharImagePointer ScalarConnectedComponentFilter( CharImagePointer inImage )
{
	ScalarConnectedComponentImageFilterType::Pointer segment = ScalarConnectedComponentImageFilterType::New();
	segment->SetInput( inImage );
	try	{	segment->Update();	}
	catch (itk::ExceptionObject& excp)	{		std::cerr << "\n ScalarConnectedComponent Err 1:  Exception caught " << "\n" << excp << std::endl;	}
	return segment->GetOutput();
}

#include "itkConfidenceConnectedImageFilter.h"
typedef itk::ConfidenceConnectedImageFilter< FloatImageType, FloatImageType >	ConfidenceConnectedImageFilterType;
CharImagePointer ConfidenceConnectedSegmentation( FloatImagePointer inImage, CharImageIndexType centroid, int rad = 3 )
{
	ConfidenceConnectedImageFilterType::Pointer confidenceConnected = ConfidenceConnectedImageFilterType::New();
	confidenceConnected->SetInput( CurvatureFlowFilter<FloatImageType, FloatImageType, FloatImagePointer>( inImage ) );
	confidenceConnected->SetMultiplier( 2.5 );
	//confidenceConnected->SetTimeStep( 0.125 );
	confidenceConnected->SetNumberOfIterations( 5 );
	confidenceConnected->SetReplaceValue( 255 );
	confidenceConnected->SetInitialNeighborhoodRadius( rad );
	confidenceConnected->AddSeed( centroid );

	try	{	confidenceConnected->Update();	}
	catch (itk::ExceptionObject& excp)	{		std::cerr << "\n ConfidenceConnectedSegmentation Err 1:  Exception caught " << "\n" << excp << std::endl;	}

	return CastFilter<FloatImageType, CharImageType, CharImagePointer>( confidenceConnected->GetOutput() );
}

#include "itkFastMarchingImageFilter.h"
typedef  itk::FastMarchingImageFilter< CharImageType, CharImageType >			CharFastMarchingFilterType;
typedef  itk::FastMarchingImageFilter< FloatImageType, FloatImageType >			FloatFastMarchingFilterType;
CharImagePointer FastMarchingSegmentation( CharImagePointer inImage, CharImageIndexType centroid )
{
	CharToFloatCastImageFilterType::Pointer charToFloat = CharToFloatCastImageFilterType::New();

	//typedef itk::CurvatureFlowImageFilter< FloatImageType, FloatImageType > CurvatureFlowImageFilterType;
	//CurvatureFlowImageFilterType::Pointer smoothing = CurvatureFlowImageFilterType::New();
	//smoothing->SetNumberOfIterations( 5 );
	//smoothing->SetTimeStep( 0.125 );

	FloatCurvatureAnisotropicDiffusionImageFilterType::Pointer smoothing = FloatCurvatureAnisotropicDiffusionImageFilterType::New();
	smoothing->SetTimeStep( 0.125 );
	smoothing->SetNumberOfIterations(  5 );
	smoothing->SetConductanceParameter( 9.0 );

	FloatGradientMagnitudeRecursiveGaussianImageFilterType::Pointer  gradientMagnitude = FloatGradientMagnitudeRecursiveGaussianImageFilterType::New();
	gradientMagnitude->SetSigma(  1.0  );		// Sigma = 1.0

	FloatSigmoidImageFilterType::Pointer sigmoid = FloatSigmoidImageFilterType::New();
	sigmoid->SetOutputMinimum( 0.0 );
	sigmoid->SetOutputMaximum( 1.0 );
	sigmoid->SetAlpha( -0.5 );					// SigmoidAlpha = -0.5, -0.3
	sigmoid->SetBeta( 3 );						// SigmoidBeta = 3.0, 2.0

	FloatFastMarchingFilterType::Pointer  fastMarching = FloatFastMarchingFilterType::New();
	fastMarching->SetStoppingValue( 100 );		// StopingCreteria = 100

	typedef FloatFastMarchingFilterType::NodeContainer   NodeContainer;
	typedef FloatFastMarchingFilterType::NodeType        NodeType;
	NodeType node;
	const double seedValue = 0.0;
	node.SetValue( seedValue );
	node.SetIndex( centroid );
	NodeContainer::Pointer seeds = NodeContainer::New();
	seeds->Initialize();
	seeds->InsertElement( 0, node );
	fastMarching->SetTrialPoints(  seeds  );
	fastMarching->SetOutputSize( inImage->GetBufferedRegion().GetSize() );

	FloatToCharThresholdingFilterType::Pointer thresholder = FloatToCharThresholdingFilterType::New();
	thresholder->SetLowerThreshold( 0 );
	thresholder->SetUpperThreshold( 100 );		// 8 - TimeThreshold = 100,200
	thresholder->SetOutsideValue(  0  );
	thresholder->SetInsideValue(  255 );

	charToFloat->SetInput( inImage );
	smoothing->SetInput( charToFloat->GetOutput() );
	gradientMagnitude->SetInput( smoothing->GetOutput() );
	fastMarching->SetInput( gradientMagnitude->GetOutput() );
	thresholder->SetInput( fastMarching->GetOutput() );

	try	{	thresholder->Update();	}
	catch( itk::ExceptionObject & excep )	{	std::cerr << "Exception caught !" << std::endl	<< excep << std::endl; }

	return thresholder->GetOutput();
}

#include "itkFiles//itkLevelSetBasedCellSegmentation.h"
typedef itk::LevelSetBasedCellSegmentation< CharImageType, CharImageType >		CharLevelSetBasedCellSegmentationType;
typedef itk::LevelSetBasedCellSegmentation< FloatImageType, FloatImageType >    FloatLevelSetBasedCellSegmentationType;
CharImagePointer LevelSetBasedCellSegmentation( CharImagePointer inImage, std::vector< CharImageIndexType > &centroids )
{
	CharLevelSetBasedCellSegmentationType::Pointer levelSetSegmentor = CharLevelSetBasedCellSegmentationType::New();

	levelSetSegmentor->SetInput( inImage );
	levelSetSegmentor->SetLargestCellRadius( 10.0 ); // in real coordinates
	levelSetSegmentor->SetSeedValue( 3 );
	levelSetSegmentor->SetIterations( 500 );
	levelSetSegmentor->SetPropagationScaling( 2 );
	levelSetSegmentor->SetCurvatureScaling( 1 );
	levelSetSegmentor->SetAdvectionScaling( 4 );
	levelSetSegmentor->SetMaxRMSChange( 0.01 );

	for( unsigned int seedNo = 0; seedNo < centroids.size(); seedNo++ )
		levelSetSegmentor->seeds[seedNo] = centroids[seedNo];

	try		
	{
		levelSetSegmentor->Update();
	}
	catch (itk::ExceptionObject& excp)
	{
		std::cerr << "Err in LevelSetBasedCellSegmentationFilter 1: - Exception caught " << "\n" << excp << std::endl;	
		return NULL;
	}
	CharImagePointer outImage = levelSetSegmentor->GetOutput();
	return outImage;
}

/******************************** *******************************************************/
#include "itkApproximateSignedDistanceMapImageFilter.h"
typedef itk::ApproximateSignedDistanceMapImageFilter< CharImageType, FloatImageType >	ApproximateSignedDistanceMapImageFilterType;
CharImagePointer ApproximateSignedDistanceMapFilter( CharImagePointer inImage )
{
	ApproximateSignedDistanceMapImageFilterType::Pointer distance = ApproximateSignedDistanceMapImageFilterType::New();
	distance->SetInput( inImage );
	distance->SetInsideValue( 255 );
	distance->SetOutsideValue( 0 );

	FloatToCharRescaleIntensityImageFilterType::Pointer rescaling = FloatToCharRescaleIntensityImageFilterType::New();
	rescaling->SetInput( distance->GetOutput() );
	rescaling->SetOutputMinimum( 0 );
	rescaling->SetOutputMaximum( 255 );

	rescaling->Update();

	return rescaling->GetOutput();
}

#include "itkLabelContourImageFilter.h"
typedef itk::LabelContourImageFilter< CharImageType, CharImageType >					CharLabelContourImageFilterType;
CharImagePointer LabelContourFilter( CharImagePointer inImage )
{
	CharLabelContourImageFilterType::Pointer labelContour = CharLabelContourImageFilterType::New();
	labelContour->SetInput( inImage );

	CharToCharRescaleIntensityImageFilter::Pointer rescaling = CharToCharRescaleIntensityImageFilter::New();
	rescaling->SetInput( labelContour->GetOutput() );
	rescaling->SetOutputMinimum( 0 );
	rescaling->SetOutputMaximum( 255 );

	rescaling->Update();

	return rescaling->GetOutput();
}

#include "itkConnectedComponentImageFilter.h"
typedef itk::ConnectedComponentImageFilter< CharImageType, CharImageType >			CharConnectedComponentImageFilterType;
CharImagePointer ConnectedComponentWithLabelContourFilter( CharImagePointer inImage )
{
	// Label the boundary of segment region
	CharConnectedComponentImageFilterType::Pointer segmenting = CharConnectedComponentImageFilterType::New();
	segmenting->SetInput( inImage );
	segmenting->Update();

	CharLabelContourImageFilterType::Pointer labelContour = CharLabelContourImageFilterType::New();
	labelContour->SetInput( segmenting->GetOutput() );

	CharToCharRescaleIntensityImageFilter::Pointer rescaling = CharToCharRescaleIntensityImageFilter::New();
	rescaling->SetInput( labelContour->GetOutput() );
	rescaling->SetOutputMinimum( 0 );
	rescaling->SetOutputMaximum( 255 );

	rescaling->Update();

	return rescaling->GetOutput();
}

//Relabel region
#include "itkCustomColormapFunction.h"
typedef itk::Function::CustomColormapFunction< CharImagePixelType, RGBImagePixelType >	ColormapType;
#include "itkRelabelComponentImageFilter.h"
typedef itk::RelabelComponentImageFilter< CharImageType, CharImageType >			CharRelabelComponentImageFilterType;
#include "itkScalarToRGBColormapImageFilter.h"
typedef itk::ScalarToRGBColormapImageFilter< CharImageType, RGBImageType>			ColormapFilterType;
RGBImagePointer RelabelColormap( CharImagePointer inImage )
{
	CharRelabelComponentImageFilterType::Pointer relabel = CharRelabelComponentImageFilterType::New();
	relabel->SetInput( inImage );
	
	ColormapType::Pointer largeColormap = ColormapType::New();
	// need to generate the color map for RGBColorMap
	//CreateRandomColormap( 255, largeColormap );

	ColormapFilterType::Pointer colorMapFilter = ColormapFilterType::New();

	colorMapFilter->SetInput( relabel->GetOutput() );
	colorMapFilter->SetColormap( largeColormap );
	colorMapFilter->Update();
	return colorMapFilter->GetOutput();
}

#include "itkHConcaveImageFilter.h"			// HConcave = Original - HMinimaImage
typedef itk::HConcaveImageFilter< CharImageType, CharImageType >						CharHConcaveImageFilterType;
CharImagePointer HConcaveFilter( CharImagePointer inImage, int height )
{
	CharHConcaveImageFilterType::Pointer minima = CharHConcaveImageFilterType::New();
	minima->SetInput( inImage );
	minima->SetHeight( height );
	minima->SetFullyConnected( true );
	minima->Update();
	return minima->GetOutput();
}

#include "itkHMaximaImageFilter.h"
typedef itk::HMaximaImageFilter< CharImageType, CharImageType >							CharHMaximaImageFilterType;
CharImagePointer HMaximaFilter( CharImagePointer inImage, int height )
{
	CharHMaximaImageFilterType::Pointer maxima = CharHMaximaImageFilterType::New();
	maxima->SetInput( inImage );
	maxima->SetFullyConnected( true );
	maxima->SetHeight( height );
	maxima->Update();
	return maxima->GetOutput();
}

#include "itkHConvexImageFilter.h"		// HConcave = Original - HMinimaImage
typedef itk::HConvexImageFilter< CharImageType, CharImageType >							CharHConvexImageFilterType;
CharImagePointer HConvexFilter( CharImagePointer inImage, int height )
{
	CharHConvexImageFilterType::Pointer maxima = CharHConvexImageFilterType::New();
	maxima->SetInput( inImage );
	maxima->SetFullyConnected( true );
	maxima->SetHeight( height );
	maxima->Update();
	return maxima->GetOutput();
}

#include "itkScalarToRGBPixelFunctor.h"
#include "itkUnaryFunctorImageFilter.h"
#include "itkWatershedImageFilter.h"
/*
CharImagePointer WatershedSegmentation( CharImagePointer image, float lowerThreshold, float outputScaleLevel )
{
	FloatGradientMagnitudeRecursiveGaussianImageFilterType::Pointer gradienMagnitudeFilter = FloatGradientMagnitudeRecursiveGaussianImageFilterType::New();
	gradienMagnitudeFilter->SetInput( image );
	gradienMagnitudeFilter->SetSigma( 1.0 );

	typedef  itk::WatershedImageFilter< CharImageType > WatershedImageFilterType;
	WatershedImageFilterType::Pointer watershedFilter = WatershedImageFilterType::New();
	watershedFilter->SetInput( gradienMagnitudeFilter->GetOutput() );
	watershedFilter->SetThreshold( lowerThreshold );
	watershedFilter->SetLevel( outputScaleLevel );

//  Instantiate the filter that will encode the label image
//  into a color image (random color attribution).
  
	typedef itk::Functor::ScalarToRGBPixelFunctor< unsigned long > ColorMapFunctorType;
	typedef itk::UnaryFunctorImageFilter< LabeledImageType, RGBImageType, ColorMapFunctorType > ColorMapFilterType;
	ColorMapFilterType::Pointer colorMapFilter = ColorMapFilterType::New();
	colorMapFilter->SetInput(  watershedFilter->GetOutput() );
	colorMapFilter->Update();

	return watershedFilter->GetOutput();
	return colorMapFilter->GetOutput();
	
}
*/
#include "itkScalarImageKmeansImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkImageRegionIterator.h"
typedef itk::ScalarImageKmeansImageFilter< CharImageType > KMeansFilterType;
CharImagePointer KMeansImageClassifier( CharImagePointer inImage, double &meanA, double &meanB, double &meanC)
{
 
	KMeansFilterType::Pointer kmeansFilter = KMeansFilterType::New();
 
	kmeansFilter->SetInput(inImage);
	kmeansFilter->SetUseNonContiguousLabels( true );
	kmeansFilter->AddClassWithInitialMean( meanA );
	kmeansFilter->AddClassWithInitialMean( meanB );
	kmeansFilter->AddClassWithInitialMean( meanC );
	kmeansFilter->Update();
 
	KMeansFilterType::ParametersType estimatedMeans = kmeansFilter->GetFinalMeans();
 
	const unsigned int numberOfClasses = estimatedMeans.Size();
 
	for(unsigned int i = 0 ; i < numberOfClasses ; ++i)
	{
		std::cout << "cluster[" << i << "] ";
		std::cout << "    estimated mean : " << estimatedMeans[i] << std::endl;
	}
	meanA = estimatedMeans[0];
	meanB = estimatedMeans[1];
	meanC = estimatedMeans[2];
 
	typedef itk::RelabelComponentImageFilter< CharImageType, CharImageType > RelabelFilterType;
 
	RelabelFilterType::Pointer relabeler = RelabelFilterType::New();
	relabeler->SetInput( kmeansFilter->GetOutput() );

	return relabeler->GetOutput();
 
	/*typedef itk::RescaleIntensityImageFilter< CharImagePointer, CharImagePointer > RescaleFilterType;
	RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
	rescaleFilter->SetInput(relabeler->GetOutput());
	rescaleFilter->SetOutputMinimum(0);
	rescaleFilter->SetOutputMaximum(255);
 
	typedef std::vector< unsigned long > SizesType;
 
	const SizesType &  sizes = relabeler->GetSizeOfObjectsInPixels();
 
	SizesType::const_iterator sizeItr = sizes.begin();
	SizesType::const_iterator sizeEnd = sizes.end();
 
	std::cout << "Number of pixels per class " << std::endl;
	unsigned int kclass = 0;
	while( sizeItr != sizeEnd )
	{
		std::cout << "Class " << kclass << " = " << *sizeItr << std::endl;
		++kclass;
		++sizeItr;
	}
	return rescaleFilter->GetOutput();*/
}

#include "itkKdTree.h"
#include "itkKdTreeBasedKmeansEstimator.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkImageToListSampleAdaptor.h"
void KdTreeBasedKmeansEstimator( CharImagePointer inImage, double &meanA, double &meanB, double &meanC )
{
	typedef itk::Statistics::ImageToListSampleAdaptor< CharImageType >   AdaptorType;
	AdaptorType::Pointer adaptor = AdaptorType::New();
	adaptor->SetImage(  inImage );

	// Define the Measurement vector type from the AdaptorType
	typedef AdaptorType::MeasurementVectorType  MeasurementVectorType;
	// Create the K-d tree structure
	typedef itk::Statistics::WeightedCentroidKdTreeGenerator< AdaptorType >	TreeGeneratorType;
	TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();
	treeGenerator->SetSample( adaptor );
	treeGenerator->SetBucketSize( 16 );
	treeGenerator->Update();

	typedef TreeGeneratorType::KdTreeType TreeType;
	typedef itk::Statistics::KdTreeBasedKmeansEstimator<TreeType> EstimatorType;
	EstimatorType::Pointer estimator = EstimatorType::New();

	const unsigned int numberOfClasses = 3;
	EstimatorType::ParametersType initialMeans( numberOfClasses );
	initialMeans[0] = meanA;		// 25
	initialMeans[1] = meanB;		// 125
	initialMeans[2] = meanC;		// 255
	estimator->SetParameters( initialMeans );

	estimator->SetKdTree( treeGenerator->GetOutput() );
	estimator->SetMaximumIteration( 200 );
	estimator->SetCentroidPositionChangesThreshold(0.0);
	estimator->StartOptimization();

	EstimatorType::ParametersType estimatedMeans = estimator->GetParameters();
	for ( unsigned int i = 0 ; i < numberOfClasses ; ++i )
	{
		std::cout << "cluster[" << i << "] " << std::endl;
		std::cout << "    estimated mean : " << estimatedMeans[i] << std::endl;
	}
	meanA = estimatedMeans[0];
	meanB = estimatedMeans[1];
	meanC = estimatedMeans[2];

	//return estimator->GetOutput();
}