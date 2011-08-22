//#define NDEBUG

namespace itk 
{

//  Software Guide : BeginCodeSnippet
template <class TInputImageType, class TOutputImageType>
GMMEstimatorImageFilter<TInputImageType, TOutputImageType>
::GMMEstimatorImageFilter() :
m_NumberOfClasses (2),
initialParameters (m_NumberOfClasses),
initialProportions (m_NumberOfClasses)
{

  im2list = itkScalarImageToListAdaptorType::New(); 
  estimator = EstimatorType::New();
  estimator->SetMaximumIteration( 200 );

} 

template <class TInputImageType, class TOutputImageType>
void
GMMEstimatorImageFilter<TInputImageType, TOutputImageType>:: SetNumberOfClasses (unsigned int numOfClasses)
{
	m_NumberOfClasses = numOfClasses;
	
	typename itk::Array <double> t_initialProportions (m_NumberOfClasses);
	initialProportions = t_initialProportions;
	initialParameters.resize (m_NumberOfClasses);
	
}
template <class TInputImageType, class TOutputImageType>
unsigned int
GMMEstimatorImageFilter<TInputImageType, TOutputImageType>:: GetNumberOfClasses ()
{
	return m_NumberOfClasses;
}
 
template <class TInputImageType, class TOutputImageType>
void
GMMEstimatorImageFilter<TInputImageType, TOutputImageType>:: GetClassParams (unsigned int classNum, double &mean, double &std)
{
	mean = (components[classNum])->GetFullParameters()[0];
	std = sqrt ((components[classNum])->GetFullParameters()[1]);
}


template <class TInputImageType, class TOutputImageType>
void
GMMEstimatorImageFilter<TInputImageType, TOutputImageType>::
GenerateData()
{

 
  
  im2list->SetImage( this->GetInput() );
  for ( unsigned int i = 0 ; i < m_NumberOfClasses ; i++ )
  {
	components.push_back( ComponentType::New() );
	(components[i])->SetSample( im2list );
	(components[i])->SetParameters( initialParameters[i] );
  }
  estimator->SetSample( im2list );
  
  for ( unsigned int i = 0 ; i < m_NumberOfClasses ; i++)
  {
	estimator->AddComponent( (typename ComponentType::Superclass*)
								(components[i]).GetPointer() );
  }
	  
  estimator->Update();
  
 

#ifndef NDEBUG 
  for ( unsigned int i = 0 ; i < m_NumberOfClasses ; i++ )
    {
		std::cout << "Cluster[" << i << "]" << std::endl;
		std::cout << "    Parameters:" << std::endl;
		std::cout << "         " << (components[i])->GetFullParameters()[0] 
                  << "		   " << (components[i])->GetFullParameters()[1] 		
				<< std::endl;
		std::cout << "    Proportion: ";
		std::cout << "         " << (*estimator->GetProportions())[i] << std::endl;
    }
#endif
  
  //this->GraftOutput( this->GetInput() );
}



template <class TInputImageType, class TOutputImageType>
void GMMEstimatorImageFilter<TInputImageType, TOutputImageType>::
SetInitialParams (typename std::vector< double > &initial_means , typename std::vector< double > &initial_std, typename std::vector< double > &initial_props)
{
	for (unsigned int i=0; i< m_NumberOfClasses;++i)
	{
		ParametersType params (2);
		params [0] = initial_means[i];
		params [1] = initial_std[i];
		initialParameters[i] = params;
		initialProportions[i] = initial_props[i];
	}	
	estimator->SetInitialProportions( initialProportions );
	
}
				  
						 
template <class TInputImageType, class TOutputImageType>
void
GMMEstimatorImageFilter<TInputImageType, TOutputImageType>::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  
  os
    << indent << "to do" << std::endl;
}

} /* end namespace itk */
