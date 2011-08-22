#ifndef __itkGMMEstimatorImageFilter_h_
#define __itkGMMEstimatorImageFilter_h_

/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: CompositeFilterExample.cxx,v $
  Language:  C++
  Date:      $Date: 2005/11/19 16:31:49 $
  Version:   $Revision: 1.7 $
  Author:    Gavin Baker <gavinb@cs.mu.oz.au>

  Copyright (c) 2005 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

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

#include "itkScalarImageToListAdaptor.h"

#include "itkVector.h"
#include "itkListSample.h"
#include "itkGaussianMixtureModelComponent.h"
#include "itkExpectationMaximizationMixtureModelEstimator.h"

#include "vnl/vnl_vector.h"

namespace itk {

template <class TInputImageType, class TOutputImageType>
class ITK_EXPORT GMMEstimatorImageFilter :
    public ImageToImageFilter<TInputImageType, TOutputImageType>
{
public:

  typedef GMMEstimatorImageFilter		                   Self;
  typedef ImageToImageFilter<TInputImageType,TOutputImageType> Superclass;
  typedef SmartPointer<Self>                                   Pointer;
  typedef SmartPointer<const Self>                             ConstPointer;
  


	

  itkNewMacro(Self);

  /** Run-time type information */
  itkTypeMacro(GMMEstimatorImageFilter, ImageToImageFilter);

  /** Display */
  void PrintSelf( std::ostream& os, Indent indent ) const;

//  Software Guide : BeginLatex
//
//  Here we declare an alias (to save typing) for the image's pixel type,
//  which determines the type of the threshold value.  We then use the
//  convenience macros to define the Get and Set methods for this parameter.
//
//  Software Guide : EndLatex

//  Software Guide : BeginCodeSnippet
  typedef typename TInputImageType::PixelType PixelType;

  

  void Compute ();
  void SetInitialParams (typename std::vector< double > &initial_means , typename std::vector< double > &initial_std, typename std::vector< double > &initial_props);
  


  void SetNumberOfClasses (unsigned int numOfClasses);
  unsigned int GetNumberOfClasses ();
  void GetClassParams (unsigned int classNum, double &mean, double &std);



 

//  Software Guide : EndCodeSnippet

protected:

  GMMEstimatorImageFilter ();

//  Software Guide : BeginLatex
//
//  Now we can declare the component filter types, templated over the
//  enclosing image type:
//
//  Software Guide : EndLatex

//  Software Guide : BeginCodeSnippet


 	
//  Software Guide : EndCodeSnippet

  void GenerateData();

  typedef itk::Vector< typename TInputImageType::PixelType, 1 > MeasurementVectorType;
	typedef itk::Statistics::ScalarImageToListAdaptor <TInputImageType> itkScalarImageToListAdaptorType;
	typedef itkScalarImageToListAdaptorType SampleType;
	typedef itk::Array< double > ParametersType;
	typedef itk::Statistics::GaussianMixtureModelComponent< SampleType >  ComponentType;
	typedef itk::Statistics::ExpectationMaximizationMixtureModelEstimator< SampleType > EstimatorType;
	
  typedef typename EstimatorType::ProportionVectorType         ProportionVectorType;


private:

  GMMEstimatorImageFilter(Self&);   // intentionally not implemented
  void operator=(const Self&);          // intentionally not implemented

//  Software Guide : BeginLatex
//
//  The component filters are declared as data members, all using the smart
//  pointer types.
//
//  Software Guide : EndLatex

//  Software Guide : BeginCodeSnippet

  
  typename itkScalarImageToListAdaptorType::Pointer im2list;
  
  typename EstimatorType::Pointer estimator;

  
  
 
  unsigned int m_NumberOfClasses;
  typename std::vector< typename ComponentType::Pointer > components;
  typename std::vector< typename itk::Array <double> > initialParameters;
  typename itk::Array <double> initialProportions;
};

} // end of namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGMMEstimatorImageFilter.txx"
#endif

#undef ITK_LEAN_AND_MEAN

#endif
