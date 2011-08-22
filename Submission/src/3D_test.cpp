#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTwoPointsVesselSegmenterImageFilter.h"

int main (int argc, char* argv[])
{
	typedef float PixelType;
	typedef itk::Image<PixelType, 3> ImageType;
	typedef itk::Image<unsigned short, 3> OutImageType;
	typedef itk::ImageFileReader<ImageType> ReaderType;
	typedef itk::ImageFileWriter<OutImageType> WriterType;
	typedef itk::TwoPointsVesselSegmenterImageFilter<ImageType,OutImageType> FilterType;

	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( argv[1] );
	reader->Update();

	typedef itk::Point<float, ImageType::ImageDimension> PointType;

	
	PointType p1, p2;
	p1[0] = atof(argv[3]); p1[1] = atof(argv[4]); p1[2] = atof(argv[5]);
	p2[0] = atof(argv[6]); p2[1] = atof(argv[7]); p2[2] = atof(argv[8]);


	FilterType::Pointer filter = FilterType::New();
	filter->SetInput(reader->GetOutput());
	filter->SetStartPoint (p1);
	filter->SetEndPoint (p2);

	try
	{
		filter->Update();
	}
	catch(itk::ExceptionObject & excp)
	{
		std::cerr << "Problem with filter" <<std::endl;
		std::cerr << excp << std::endl;
	}


	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(argv[2]);
	writer->SetInput(filter->GetOutput());
	
	try
	{
		writer->Update();
	}
	catch(itk::ExceptionObject & err)
	{
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
