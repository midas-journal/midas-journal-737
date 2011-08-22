#ifndef __itkActiveContourMinCutImageFilter_h_
#define __itkActiveContourMinCutImageFilter_h_

/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $itkActiveContourMinCutImageFilter.h $
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
#include "itkGMMEstimatorImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"


namespace itk {

template <class TInputImageType, class TOutputImageType>
class ITK_EXPORT ActiveContourMinCutImageFilter :
    public ImageToImageFilter<TInputImageType, TOutputImageType>
{
public:

	typedef ActiveContourMinCutImageFilter				Self;
	typedef ImageToImageFilter<TInputImageType,TOutputImageType>	Superclass;
	typedef SmartPointer<Self>                                   	Pointer;
	typedef SmartPointer<const Self>                             	ConstPointer;

	typedef itk::Image <unsigned short, 3> 				InnerImageType;
	typedef itk::Image<float, 3> 					ImageType;	
	typedef typename TOutputImageType::Pointer 			OutImagePointer;

	itkNewMacro(Self);

	/** Run-time type information */
	itkTypeMacro(ActiveContourMinCutImageFilter, ImageToImageFilter);

	/** Display */
	void PrintSelf( std::ostream& os, Indent indent ) const;

	typedef typename TInputImageType::PixelType PixelType;

	/** Set the initial class means and proportions */
	void SetInitialMeans (const std::vector<double> & initialMeans);
	void SetInitialProportions(const std::vector<double> & initialProps);

	/** Set/Get the radius */
	itkSetMacro (UncertaintyRadius, double);	
	itkGetMacro (UncertaintyRadius, double);

	/** Set the object input image */
	itkSetObjectMacro(ObjImage, TOutputImageType);

	/** Set the number of background classes */
	itkSetMacro (NumberOfBackgroundClasses, unsigned);

protected:

	ActiveContourMinCutImageFilter ();

	void GenerateData();
 
private:

	ActiveContourMinCutImageFilter(Self&);   // intentionally not implemented
	void operator=(const Self&);          // intentionally not implemented

	double m_UncertaintyRadius;

	OutImagePointer m_ObjImage;

	unsigned int m_NumberOfBackgroundClasses;
	std::vector<double> m_InitialProportions;
	std::vector<double> m_BackgroundMeans;
	


	float computeEdgeWeight (const typename TInputImageType::PixelType &a, const typename TInputImageType::PixelType &b, 
					const typename TInputImageType::PixelType &sigma);

	typedef itk::GMMEstimatorImageFilter <TOutputImageType, TOutputImageType> itkGMMEstimatorImageFilterType;
	typedef itk::DanielssonDistanceMapImageFilter <TOutputImageType, InnerImageType> itkDanielssonDistanceMapImageFilterType;
	typedef itk::CastImageFilter <TInputImageType, TOutputImageType> itkCastImageFilterType;
	typedef itk::LabelStatisticsImageFilter  <TInputImageType, TOutputImageType> itkLabelStatisticsImageFilterType;

};

} // end of namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkActiveContourMinCutImageFilter.txx"
#endif

#undef ITK_LEAN_AND_MEAN

#endif
