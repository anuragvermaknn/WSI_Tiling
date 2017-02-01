#include <boost/filesystem.hpp>
using namespace boost::filesystem;

double string_to_double( const std::string& s )
{
	std::istringstream i(s);
	double x;
	if (!(i >> x))
		return 0;
	return x;
}

#ifndef isNaN
inline  bool isNaN( double x)
{
	return x != x;
}
#endif

void ConvertIndexToVectorData( std::vector< CharImageIndexType > indexes1D, std::vector< std::vector< double > > &data2D)
{
	for( int i=0; i<indexes1D.size(); i++)
	{
		std::vector<double> tmp;
		tmp.push_back( indexes1D[i][0] );
		tmp.push_back( indexes1D[i][1] );
		data2D.push_back( tmp );
	}
}

void TrainingSetSelection( std::vector<std::vector<double> > &t1, std::vector<std::vector<double> > &t2 )
{
	// 1- Select all mitosis
	for( int i = 0; i < t1.size(); i++)
		if( t1[i][ t1[i].size()-1 ] == 0 )				// Mitosis = 0
			t2.push_back( t1[i] );
	// 2- Select randomly nonmitosis
	srand( (unsigned) time( NULL ) );
	for( int i = 0; i < TrainingSetSize; )
	{
		int j =  (int) rand() % t1.size();
		if( t1[j][ t1[j].size()-1 ] == 1 )				// NonMitosis = 1
		{
			t2.push_back( t1[j] );
			i++;
		}
	}
}

void CalculateGTCentroid(std::string inputDir, int totalImages = 35)
{
	for(unsigned int imageNo = 0; imageNo < totalImages; imageNo++)	
	{
		std::vector<std::vector<CharImageIndexType> > annotatedIndexes;
		std::vector< CharImageIndexType > centroidIndexes;
	
		std::stringstream annotatedFile;
		annotatedFile << imageNo << csvExt;
		path annotatedFilePath ( inputDir / annotatedFile.str() );
		Read2DIndexesFromCSVFile( annotatedFilePath.make_preferred().string(), annotatedIndexes);

		for( int i = 0; i < annotatedIndexes.size(); i++ )
		{
			long xSum = 0, ySum = 0;
			for( int j = 0; j < annotatedIndexes[i].size(); j++ )
			{
				CharImageIndexType index = annotatedIndexes[i][j];
				xSum += index[0];
				ySum += index[1];
			}
			CharImageIndexType centroid;
			centroid[0] = xSum / annotatedIndexes[i].size();
			centroid[1] = ySum / annotatedIndexes[i].size();
			centroidIndexes.push_back( centroid );
		}

		std::stringstream centroidFile;
		centroidFile << imageNo << centroidsFile << csvExt;
		path centroidFilePath ( inputDir / centroidFile.str() );
		Write1DIndexesToCSVFile( centroidFilePath.make_preferred().string(), centroidIndexes );
	}
}

void ComputeDetectionTPFNFP( std::vector< CharImageIndexType > &GTCentroids, std::vector< CharImageIndexType > &candidateCentroids, 
			   std::vector<CharImageIndexType> &TPCentroids, std::vector<CharImageIndexType> &FNCentroids, std::vector<CharImageIndexType> &FPCentroids )
{
	for( int i=0; i<candidateCentroids.size(); i++)
	{
		bool centroidFound = false;
		for( int j=0; j<GTCentroids.size() && !centroidFound; j++)
			if( candidateCentroids[i] == GTCentroids[j] )
				centroidFound = true;
		if( centroidFound )
			TPCentroids.push_back( candidateCentroids[i] );
		else
			FPCentroids.push_back( candidateCentroids[i] );
	}
	for( int i=0; i<GTCentroids.size(); i++)
	{
		bool GTFound = false;
		for( int j=0; j<candidateCentroids.size() && !GTFound; j++)
			if( candidateCentroids[j] == GTCentroids[i] )
				GTFound = true;
		if( GTFound )
			FNCentroids.push_back( GTCentroids[i] );
	}
}

void ComputeDetectionFP( std::vector< CharImageIndexType > &GTCentroids, std::vector< CharImageIndexType > &candidateCentroids, std::vector< CharImageIndexType > &FPCentroids)
{
	for( int i=0; i<candidateCentroids.size(); i++)
	{
		bool centroidFound = false;
		for( int j=0; j<GTCentroids.size() && !centroidFound; j++)
			if( candidateCentroids[i] == GTCentroids[j] )
				centroidFound = true;
		if( !centroidFound )
			FPCentroids.push_back( candidateCentroids[i] );
	}
}

void ComputeSegmentationTPFNFP( std::vector<CharImageIndexType> &annotatedIndexes, std::vector<CharImageIndexType> &candidateIndexes, int &TPPixels, int &FNPixels, int &FPPixels)
{
	for( int i=0; i<annotatedIndexes.size(); i++)
	{
		bool annotatedPixelFound = false;
		for( int j=0; j<candidateIndexes.size() && !annotatedPixelFound; j++)
			if( candidateIndexes[j] == annotatedIndexes[i] )
				annotatedPixelFound = true;
		if( annotatedPixelFound )
			TPPixels++;
		else
			FNPixels++;
	}
	for( int i=0; i<candidateIndexes.size(); i++)
	{
		bool candidatePixelFound = false;
		for( int j=0; j<annotatedIndexes.size() && !candidatePixelFound; j++)
			if( candidateIndexes[i] == annotatedIndexes[j] )
				candidatePixelFound = true;
		if(!candidatePixelFound)
			FPPixels++;
	}
}

double ComputePDF( double mean, double variance, double value )
{
	itk::Statistics::GaussianDistribution::Pointer gaussian = itk::Statistics::GaussianDistribution::New();
	gaussian->SetMean(mean);
	gaussian->SetVariance(variance);
	return gaussian->EvaluatePDF(value);
}

double ComputePDFByFormula( double mean, double variance, double value )
{
	double e = exp ( - ( value - mean ) * ( value - mean ) / 2 * variance );
	double sD = sqrt( variance * 2 * 3.14 );
	double pDF = e / sD ;
	return pDF;
}

void ComputeVectorMeanVariance( std::vector<double> &data, double &mean, double &variance )
{
	double sum = 0.0;
	for(int i=0; i<data.size(); i++)
		sum += data[i];
	mean = sum / data.size();

	double temp = 0;
	for(int i=0; i<data.size(); i++)
		temp += (mean-data[i])*(mean-data[i]);
	variance = temp/data.size();
	double standardDeviation = sqrt( variance );
}

void LoadingDataset(std::string directory, unsigned int TotalFiles, std::vector<std::vector<double> > &dataset)
{
	for( int fileNo = 0; fileNo < TotalFiles; fileNo++ )
	{
		std::stringstream featuresFile;
		featuresFile << fileNo << csvExt;	
		path FeaturesPath( directory / featuresFile.str() );
		Read2DDataFromCSVFile( FeaturesPath.make_preferred().string(), dataset );
	}
}

void ComputeHistogramFromCSVFile( std::string directory, int totalImages, int histogramChannel = CandidateDetectionChannel )
{
	std::vector< std::vector< int > > histogramPerMitosis;
	std::vector< int > histogramForAllMitosis( 256,0 );
	int mitosisNumber=0;

	for(unsigned int imageNo = 0; imageNo < totalImages; imageNo++)
	{
		std::vector< std::vector< CharImageIndexType > > annotatedIndexes;

		// Step 1 - Read char image
		std::stringstream imageFile;
		imageFile << imageNo << imageExt;
		path imagePath( directory / imageFile.str() );
		RGBImagePointer rgbImage = RGBFileReading( imagePath.make_preferred().string(), "\nError Reading char image - Exception caught" );

		CharLaplacianRecursiveGaussianImageFilterType::Pointer LoGFilter = CharLaplacianRecursiveGaussianImageFilterType::New();
		CharImagePointer charImage,brImage;

		switch( histogramChannel )
		{
		case 0:
			brImage = BlueRatioExtraction( rgbImage );
			SelectedChannel	= std::string("BR");
			LoGFilter->SetInput( brImage );
			LoGFilter->Update();
			charImage = LoGFilter->GetOutput();
			break;
		case 1: 
			charImage = BlueChannelExtraction( rgbImage );
			SelectedChannel	= std::string("B");
			break;
		case 2:
			charImage = RedChannelExtraction( rgbImage );
			SelectedChannel	= std::string("R");
			break;
		case 3:
			charImage = ColorChannelExtraction( rgbImage, "HSV", 2 );
			SelectedChannel	= std::string("HSV");
			break;
		case 4:
			charImage = ColorChannelExtraction( rgbImage, "Lab", 0 );
			SelectedChannel	= std::string("Lab");
			break;
		case 5:
			charImage = ColorChannelExtraction( rgbImage, "Luv", 0 );
			SelectedChannel	= std::string("Luv");
			break;
		case 6:
			charImage = HematoxylinChannelExtraction( rgbImage );
			SelectedChannel	= std::string("HE");
			//charImage = MedianFilter( HematoxylinChannelExtraction( rgbImage ), 2 );
			//SelectedChannel	= std::string("HE_S");
			break;
		case 7: 
			charImage = GreenChannelExtraction( rgbImage );
			SelectedChannel	= std::string("G");
			break;
		default:
			charImage = HematoxylinChannelExtraction( rgbImage );
			SelectedChannel	= std::string("HE");
		}

		// Step 2 - Read GT (Manullay Annotated) 
		std::stringstream annotatedFile;
		annotatedFile << imageNo << csvExt;
		path annotatedFilePath ( directory / annotatedFile.str() );
		Read2DIndexesFromCSVFile( annotatedFilePath.make_preferred().string(), annotatedIndexes);

		// Step 3 - Calculate Histogram (per mitosis and all mitosis)
		histogramPerMitosis.resize(mitosisNumber+1);
		for(int i=0;i<256;i++)
			histogramPerMitosis[mitosisNumber].push_back(0);

		for( int i = 0; i < annotatedIndexes.size(); i++)
		{
			for ( int j = 0; j < annotatedIndexes[i].size(); j++ )
			{
				int pixelValue = charImage->GetPixel( annotatedIndexes[i][j] );
				histogramForAllMitosis[pixelValue]+=1;
				histogramPerMitosis[mitosisNumber][pixelValue]+=1;
			}
		}
		mitosisNumber++;
	}

	// Writing histogram in CSV file
	std::stringstream histogramFileName;
	histogramFileName << "Histogram_" << SelectedChannel << csvExt;
	path histogramFilePath ( directory / histogramFileName.str() );
	std::ofstream histogramFile( histogramFilePath.make_preferred().string().c_str() );

	// Writing Pixel Values (0 to 255) in CSV file
	for(int i=0; i < histogramForAllMitosis.size(); i++)
		histogramFile << i << ",";
	histogramFile << std::endl;
	// Writing Histogram of all mitosis
	for(int i=0; i < histogramForAllMitosis.size(); i++)
		histogramFile << histogramForAllMitosis[i] << ",";
	histogramFile << std::endl << std::endl << std::endl;

	// Writing Histogram of one mitosis
	for(int i=0; i < histogramPerMitosis.size(); i++)
	{
		histogramFile << std::endl;
		for(int j=0; j < histogramPerMitosis[i].size(); j++)
			histogramFile << histogramPerMitosis[i][j] << ",";
	}
	histogramFile.close();
}

void NormalizeFeatures( std::vector< std::vector< double > > inData2D, std::vector< std::vector< int > > &outData2D, int normalizeRange = 1000 )
{
	for( int r = 0; r < inData2D.size(); r++ )
	{
		std::vector< int > tmp;
		for( int c = 0; c < inData2D[r].size(); c++)
			tmp.push_back( inData2D[r][c] );
		outData2D.push_back( tmp );
	}
	
	// Start from First features then N-1 Features because last column is class label
	for( int c = 0; c < inData2D[0].size()-1; c++)
	{
		std::vector< double > oneFeature;
		// Extract one feature values for all instances
		for( int r = 0; r < inData2D.size(); r++ )
			oneFeature.push_back( inData2D[r][c] );

		double min = *std::min_element(oneFeature.begin(), oneFeature.end());
		double max = *std::max_element(oneFeature.begin(), oneFeature.end());
		
		// Normalize one feature value for all instances
		for( int r = 0; r < inData2D.size(); r++ )
			outData2D[r][c] = normalizeRange * ( oneFeature[r] - min ) / (max - min);
	}
}

void NormalizeFeatures( std::vector< std::vector< double > > inData2D, std::vector< std::vector< double > > &outData2D, int normalizeRange = 1000 )
{
	for( int r = 0; r < inData2D.size(); r++ )
	{
		std::vector< double > tmp;
		for( int c = 0; c < inData2D[r].size(); c++)
			tmp.push_back( inData2D[r][c] );
		outData2D.push_back( tmp );
	}
	
	// Start from First features then N-1 Features because last column is class label
	for( int c = 0; c < inData2D[0].size()-1; c++)
	{
		std::vector< double > oneFeature;
		// Extract one feature values for all instances
		for( int r = 0; r < inData2D.size(); r++ )
			oneFeature.push_back( inData2D[r][c] );

		double min = *std::min_element(oneFeature.begin(), oneFeature.end());
		double max = *std::max_element(oneFeature.begin(), oneFeature.end());
		
		// Normalize one feature value for all instances
		for( int r = 0; r < inData2D.size(); r++ )
			outData2D[r][c] = normalizeRange * ( oneFeature[r] - min ) / (max - min);
	}
}
