#ifndef __itkShortestPathImageFilter_h_
#define __itkShortestPathImageFilter_h_

/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $itkShortestPathImageFilter.h $
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
#include "itkRescaleIntensityImageFilter.h"
#include "itkMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter.h"
#include "itkCastImageFilter.h"


namespace itk {

template <class TInputImageType, class TOutputImageType>
class ITK_EXPORT ShortestPathImageFilter :
    public ImageToImageFilter<TInputImageType, TOutputImageType>
{
public:

	typedef ShortestPathImageFilter			             	Self;
	typedef ImageToImageFilter<TInputImageType,TOutputImageType>	Superclass;
	typedef SmartPointer<Self>                                 	Pointer;
	typedef SmartPointer<const Self>                            	ConstPointer;



 	itkNewMacro(Self);

	/** Run-time type information */
	itkTypeMacro(ShortestPathImageFilter, ImageToImageFilter);

	/** Display */	
	void PrintSelf( std::ostream& os, Indent indent ) const;

	typedef typename TInputImageType::PixelType PixelType;


	void SetStartPoint (const typename TInputImageType::PointType & StartPoint);
	void SetEndPoint(const typename TInputImageType::PointType & EndPoint);


	/** Set/Get the euclid and neighbors mode */
	itkSetMacro (FullNeighborsMode, int);
	itkGetMacro (FullNeighborsMode, int);
	itkSetMacro (EuclidMode, int);	
	itkGetMacro (EuclidMode, int);

	void SetModeToEuclid ();
	void SetModeToNonEuclid ();

protected:

	ShortestPathImageFilter ();

	void GenerateData();


private:

	ShortestPathImageFilter(Self&);   // intentionally not implemented
	void operator=(const Self&);          // intentionally not implemented

	int m_EuclidMode;
	int m_FullNeighborsMode;

	typename TInputImageType::PointType m_StartPoint, m_EndPoint;  //compelet Set/Get functions

	typedef MultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter<TInputImageType, TOutputImageType>
									itkMultiScaleVesselnessFilterType;
	typedef RescaleIntensityImageFilter<TOutputImageType,TOutputImageType>
									itkRescaleFilterType;
	typedef itk::CastImageFilter <TInputImageType,TOutputImageType> itkCastImageFilterType;

	double computeDijkstraEdgeWeight (const typename itkRescaleFilterType::OutputPixelType &a_vesselness);

	int computeNodeNum (const typename TInputImageType::SizeType &size,
					const typename TInputImageType::IndexType &index);

	const static int EUCLID_SHORT = 0;
	const static int NON_EUCLID_SHORT = 1;

};

} // end of namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkShortestPathImageFilter.txx"
#endif

#undef ITK_LEAN_AND_MEAN

#endif
