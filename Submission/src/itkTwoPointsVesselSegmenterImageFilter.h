#ifndef __itkTwoPointsVesselSegmenterImageFilter_h_
#define __itkTwoPointsVesselSegmenterImageFilter_h_

/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $itkTwoPointsVesselSegmenterImageFilter.h $
  Language:  C++
  Date:      $Date: 2010/06/26 $
  Version:   $Revision: 1.0 $
  Author:    David Pellow <david.pellow@utoronto.ca>


     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif


#include "itkImageToImageFilter.h"
#include "itkActiveContourMinCutImageFilter.h"
#include "itkShortestPathImageFilter.h"
#include "itkCastImageFilter.h"

namespace itk {

template <class TInputImageType, class TOutputImageType>
class ITK_EXPORT TwoPointsVesselSegmenterImageFilter :
    public ImageToImageFilter<TInputImageType, TOutputImageType>
{
public:

	typedef TwoPointsVesselSegmenterImageFilter			Self;
	typedef ImageToImageFilter<TInputImageType,TOutputImageType> 	Superclass;
	typedef SmartPointer<Self>                                   	Pointer;
	typedef SmartPointer<const Self>                             	ConstPointer;


	itkNewMacro(Self);

	/** Run-time type information */
	itkTypeMacro(TwoPointsVesselSegmenterImageFilter, ImageToImageFilter);

	/** Display */
	void PrintSelf( std::ostream& os, Indent indent ) const;


	/**Set the shortest path neighbors mode */
	itkSetMacro (FullNeighborsMode, int);


	/**Set the active contour radius */
	itkSetMacro (ActiveContourRadius, double);

	/** Set the number of background classes */
	itkSetMacro (NumberOfBackgroundClasses, unsigned);


	typedef typename TInputImageType::PointType PointType;
	
	/**Set the index of the start and end points */
	void SetStartPoint (const typename TInputImageType::PointType &StartPoint);
	void SetEndPoint (const typename TInputImageType::PointType &EndPoint);

	/** Set the initial class means and proportions */
	void SetInitialMeans (const std::vector<double> & initialMeans);
	void SetInitialProportions(const std::vector<double> & initialProps);

protected:

	TwoPointsVesselSegmenterImageFilter ();


	void GenerateData();

private:

	TwoPointsVesselSegmenterImageFilter(Self&);   // intentionally not implemented
	void operator=(const Self&);          // intentionally not implemented

	typename TInputImageType::PointType m_StartPoint, m_EndPoint;  

	int m_FullNeighborsMode;
		
	double m_ActiveContourRadius;

	unsigned int m_NumberOfBackgroundClasses;
	std::vector<double> m_InitialProportions;
	std::vector<double> m_BackgroundMeans;

	typedef itk::ShortestPathImageFilter <TInputImageType, TInputImageType> itkShortestPathImageFilterType;
	typedef itk::CastImageFilter<TInputImageType, TOutputImageType> itkCastImageFilterType;
	typedef itk::ActiveContourMinCutImageFilter<TInputImageType, TOutputImageType> itkActiveContourMinCutImageFilterType;

};

} // end of namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTwoPointsVesselSegmenterImageFilter.txx"
#endif

#undef ITK_LEAN_AND_MEAN

#endif
