/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $itkTwoPointsVesselSegmenterImageFilter.txx $
  Language:  C++
  Date:      $Date: 2010/06/26 $
  Version:   $Revision: 1.0 $
  Author:    David Pellow <david.pellow@utoronto.ca>


     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#define NDEBUG

namespace itk 
{

template <class TInputImageType, class TOutputImageType>
TwoPointsVesselSegmenterImageFilter<TInputImageType, TOutputImageType>
::TwoPointsVesselSegmenterImageFilter() :
m_FullNeighborsMode (1), m_ActiveContourRadius (5.0), m_NumberOfBackgroundClasses(4)
{
	m_BackgroundMeans.push_back(-200);
	m_BackgroundMeans.push_back(100);
	m_BackgroundMeans.push_back(200);
	m_BackgroundMeans.push_back(300);

	m_InitialProportions.push_back(0.01);
	m_InitialProportions.push_back(0.04);
	m_InitialProportions.push_back(0.005);
	m_InitialProportions.push_back(0.003);
	m_InitialProportions.push_back(0.002);
	m_InitialProportions.push_back(0.94);
} 


template <class TInputImageType, class TOutputImageType>
void TwoPointsVesselSegmenterImageFilter<TInputImageType, TOutputImageType>::
SetStartPoint (const typename TInputImageType::PointType &StartPoint)
{
	for (unsigned int i=0;i<TInputImageType::ImageDimension;++i)
		m_StartPoint[i] = StartPoint[i];
}


template <class TInputImageType, class TOutputImageType>
void TwoPointsVesselSegmenterImageFilter<TInputImageType, TOutputImageType>::
SetEndPoint (const typename TInputImageType::PointType &EndPoint)
{
	for (unsigned int i=0;i<TInputImageType::ImageDimension;++i)
		m_EndPoint[i] = EndPoint[i];
}


template <class TInputImageType, class TOutputImageType>
void TwoPointsVesselSegmenterImageFilter<TInputImageType, TOutputImageType>
::SetInitialMeans(const std::vector<double> & initialMeans)
{
	m_BackgroundMeans.clear();
	for (unsigned i = 0; i < m_NumberOfBackgroundClasses; i++)
		m_BackgroundMeans.push_back(initialMeans[i]);
	m_BackgroundMeans.resize(m_NumberOfBackgroundClasses);
}


template <class TInputImageType, class TOutputImageType>
void TwoPointsVesselSegmenterImageFilter<TInputImageType, TOutputImageType>
::SetInitialProportions(const std::vector<double> & initialProps)
{
	m_InitialProportions.clear();	
	for (unsigned i = 0; i < m_NumberOfBackgroundClasses+2; i++)
		m_InitialProportions.push_back(initialProps[i]);
	m_InitialProportions.resize(m_NumberOfBackgroundClasses+2);
}


template <class TInputImageType, class TOutputImageType>
void
TwoPointsVesselSegmenterImageFilter<TInputImageType, TOutputImageType>::
GenerateData()
{
	//find shortest path	
	typename itkShortestPathImageFilterType::Pointer dijkstra = itkShortestPathImageFilterType::New();
	dijkstra->SetInput (this->GetInput());
	dijkstra->SetStartPoint (m_StartPoint);
	dijkstra->SetEndPoint (m_EndPoint);
	dijkstra->SetFullNeighborsMode (m_FullNeighborsMode);
	dijkstra->Update();
#ifndef NDEBUG
	std::cout << "Completed shortest path filter" << std::endl;

	typedef itk::ImageFileWriter< TInputImageType > WriterType;
	typename WriterType::Pointer writer = WriterType::New();
	writer->SetFileName( "path.vtk" );
	writer->SetInput (dijkstra->GetOutput());

	try
	{
		writer->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << "Problem encountered while writing ";
		std::cerr << " image file : " <<  "path.vtk" << std::endl;
		std::cerr << excp << std::endl;
	}
#endif	

	typename itkCastImageFilterType::Pointer caster = itkCastImageFilterType::New();
	caster->SetInput (dijkstra->GetOutput());
	caster->Update();
	typename TOutputImageType::Pointer outputImage = caster->GetOutput();


	//compute min cut
	typename itkActiveContourMinCutImageFilterType::Pointer minCut = itkActiveContourMinCutImageFilterType::New();
	minCut->SetInput (this->GetInput());
	minCut->SetObjImage (outputImage);
	minCut->SetUncertaintyRadius (m_ActiveContourRadius);
	minCut->SetNumberOfBackgroundClasses(m_NumberOfBackgroundClasses);
	minCut->SetInitialMeans(m_BackgroundMeans);
	minCut->SetInitialProportions(m_InitialProportions);

	try
	{
		minCut->Update();
	}
	catch( itk::ExceptionObject & excp )
   	{
		std::cerr << "Problem encountered while updating minCut Filter" << std::endl;
		std::cerr << excp << std::endl;
   	}

	std::cout << "Done segmentation" << std::endl;
		
	this->GraftOutput(minCut->GetOutput());
}


template <class TInputImageType, class TOutputImageType>
void
TwoPointsVesselSegmenterImageFilter<TInputImageType, TOutputImageType>::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  
  os
    << indent << "to do" << std::endl;
}

} // end namespace itk 
