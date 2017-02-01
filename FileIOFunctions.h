RGBImagePointer RGBFileReading( std::string fileName, std::string errorMsg )
{
	RGBReaderType::Pointer reader = RGBReaderType::New();
	reader->SetFileName( fileName );
	try {	
		reader->Update();	
	}
	catch ( itk::ExceptionObject& excp )	{	
		std::cerr << errorMsg << excp << std::endl;	
		throw ( std::runtime_error("Access error") );
		return NULL;
	}
	return reader->GetOutput();
}
CharImagePointer CharFileReading( std::string fileName, std::string errorMsg )
{
	CharReaderType::Pointer reader = CharReaderType::New();
	reader->SetFileName( fileName );
	try {	
		reader->Update();	
	}
	catch ( itk::ExceptionObject& excp )	{	
		std::cerr << errorMsg << excp << std::endl;	
		throw ( std::runtime_error("Access error") );
		return NULL;
	}
	return reader->GetOutput();
}
FloatImagePointer FloatFileReading( std::string fileName, std::string errorMsg )
{
	FloatReaderType::Pointer reader = FloatReaderType::New();
	reader->SetFileName( fileName );
	try {	
		reader->Update();	
	}
	catch ( itk::ExceptionObject& excp )	{	
		std::cerr << errorMsg << excp << std::endl;	
		throw ( std::runtime_error("Access error") );
		return NULL;
	}
	return reader->GetOutput();
}
void RGBFileWriting( RGBImagePointer inImage, std::string fileName, std::string errorMsg )
{
	RGBWriterType::Pointer writer = RGBWriterType::New();
	writer->SetFileName( fileName );		
	writer->SetInput( inImage );
	try {	
		writer->Update();	
	}
	catch ( itk::ExceptionObject& excp )	{	
		std::cerr << errorMsg << excp << std::endl;	
		throw ( std::runtime_error("Access error") );
	}
}
void CharFileWriting( CharImagePointer inImage, std::string fileName, std::string errorMsg )
{
	CharWriterType::Pointer writer = CharWriterType::New();
	writer->SetFileName( fileName );		
	writer->SetInput( inImage );
	try {	
		writer->Update();	
	}
	catch ( itk::ExceptionObject& excp )	{	
		std::cerr << errorMsg << excp << std::endl;	
		throw ( std::runtime_error("Access error") );
	}
}
void FloatFileWriting( FloatImagePointer inImage, std::string fileName, std::string errorMsg )
{
	FloatWriterType::Pointer writer = FloatWriterType::New();
	writer->SetFileName( fileName );		
	writer->SetInput( inImage );
	try {	
		writer->Update();	
	}
	catch ( itk::ExceptionObject& excp )	{	
		std::cerr << errorMsg << excp << std::endl;	
		throw ( std::runtime_error("Access error") );
	}
}
void Read1DIndexesFromCSVFile( std::string fileName, std::vector<CharImageIndexType> &data1D )
{
	std::ifstream inCSVFile( fileName.c_str() );
	std::string fileLine;
	int mitosisNumber = 0;
	data1D.resize( 0 );
	while( getline( inCSVFile, fileLine ) )
	{
		std::string lineItems;
		int itemNumber=0;
		CharImageIndexType index;
		std::istringstream lineStream(fileLine);
		while( getline( lineStream, lineItems, ','))
		{
			int val = atoi(lineItems.c_str());
			if(itemNumber%2==0)
				index[0]=val;
			else
			{
				index[1]=val;
				data1D.push_back(index);
			}
			itemNumber++;
		}
		mitosisNumber++;
	}
	inCSVFile.close();
}
void Read2DIndexesFromCSVFile( std::string fileName, std::vector<std::vector<CharImageIndexType> > &data2D )
{
	std::ifstream inCSVFile( fileName.c_str() );
	std::string fileLine;
	int mitosisNumber = 0;
	data2D.resize( 0 );
	while( getline( inCSVFile, fileLine ) )
	{
		data2D.resize( mitosisNumber+1 );
		std::string lineItems;
		int itemNumber=0;
		CharImageIndexType index;
		std::istringstream lineStream(fileLine);
		while( getline( lineStream, lineItems, ','))
		{
			int val = atoi(lineItems.c_str());
			if(itemNumber%2==0)
				index[0]=val;
			else
			{
				index[1]=val;
				data2D[mitosisNumber].push_back(index);
			}
			itemNumber++;
		}
		mitosisNumber++;
	}
	inCSVFile.close();
}
void Read2DDataFromCSVFile(std::string fileName, std::vector< std::vector< double > > &data2D )
{
	std::string line;
	std::ifstream inCSVFile( fileName.c_str() );
	while( getline( inCSVFile, line))
	{
		std::vector<double> tmp;
		std::string lineItems;
		std::istringstream lineStream(line);
		while( getline(lineStream, lineItems, ','))
			tmp.push_back( atof( lineItems.c_str() ) );
		data2D.push_back( tmp );
	}
	inCSVFile.close();
}
void WriteImageToTextFile( CharImageType *inImage, std::string fileName)
{
	std::ofstream file( fileName.c_str());
	CharImageRegionType region = inImage->GetLargestPossibleRegion();
	CharImageLinearConstIteratorWithIndexType it (inImage, region);
	for(it.GoToBegin(); !it.IsAtEnd(); it.NextLine())	
	{
		while(!it.IsAtEndOfLine())
		{
			file << it.Get() << ",";
			++it;
		}
		file << std::endl;
	}
	file.close();
}
void WriteSegmentImageIndexToCSVFile( CharImageType *inImage, std::string fileName)
{
	std::ofstream file( fileName.c_str());
	CharImageRegionType region = inImage->GetLargestPossibleRegion();
	CharImageLinearConstIteratorWithIndexType it (inImage, region);
	for(it.GoToBegin(); !it.IsAtEnd(); it.NextLine())	
	{
		while(!it.IsAtEndOfLine())
		{
			if(it.Get() > 0)
			{
				CharImageIndexType index = it.GetIndex();
				file << index[0] << "," << index[1];
			}
			++it;
		}
	}
	file.close();
}
void Write1DIndexesToCSVFile( std::string fileName, std::vector<CharImageIndexType> &indexes1D )
{
	std::ofstream file(fileName.c_str());
	for( int i=0; i<indexes1D.size(); i++)
	{
		CharImageIndexType index = indexes1D[i];
		file << index[0] << "," << index[1] << std::endl;
	}
	file.close();
}
void Write2DIndexesToCSVFile( std::string fileName, std::vector<std::vector<CharImageIndexType> > indexes2D )
{
	std::ofstream file(fileName.c_str());
	for( int i = 0; i < indexes2D.size(); i++ )
	{
		for( int j = 0; j < indexes2D[i].size(); j++ )
		{
			CharImageIndexType index = indexes2D[i][j];
			file << index[0] << "," << index[1] << ",";
		}
		file << std::endl;
	}
	file.close();
}
void Write2DDataToCSVFile( std::string fileName, std::vector<std::vector<double> > data2D )
{
	std::ofstream file(fileName.c_str());
	for( int i = 0; i < data2D.size(); i++ )
	{
		for( int j = 0; j < data2D[i].size(); j++ )
			file << data2D[i][j] << "," ;
		file << std::endl;
	}
	file.close();
}
void Write3DDataToCSVFile( std::string fileName, std::vector<std::vector<std::vector<double> > > &data3D )
{
	std::ofstream file(fileName.c_str());
	for( int i = 0; i < data3D.size(); i++ )
	{
		for( int j = 0; j < data3D[i].size(); j++ )
			for( int k = 0; k < data3D[i][j].size(); k++ )
				file << data3D[i][j][k] << "," ;
		file << std::endl;
	}
	file.close();
}