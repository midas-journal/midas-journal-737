/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $itkShortestPathImageFilter.txx $
  Language:  C++
  Date:      $Date: 2010/06/26 $
  Version:   $Revision: 1.0 $
  Author:    David Pellow <david.pellow@utoronto.ca>


     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

//boost includes
#include <climits>
#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <vector>
#include <cmath>

#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkConstNeighborhoodIterator.h"  
#include "itkShapedNeighborhoodIterator.h"


//#define NDEBUG



using namespace boost;


namespace itk 
{


template <class TInputImageType, class TOutputImageType>
double ShortestPathImageFilter<TInputImageType, TOutputImageType>::
computeDijkstraEdgeWeight (const typename itkRescaleFilterType::OutputPixelType & vesselness)
{
	return (500 - vesselness);
}


template <class TInputImageType, class TOutputImageType>
int ShortestPathImageFilter<TInputImageType, TOutputImageType>::
computeNodeNum (const typename TInputImageType::SizeType &size,
				const typename TInputImageType::IndexType &index)
{
	int temp_val;
	int res=0;
	int dim=TInputImageType::ImageDimension;
	for (int k=0; k<dim;++k){
		temp_val=1;
		for (int l=0;l<k;++l)
			temp_val=temp_val*size[l];
		temp_val=temp_val*index[k];
		res=res+temp_val;
	}
	return res;
}
			
				 
template <class TInputImageType, class TOutputImageType>
void ShortestPathImageFilter<TInputImageType, TOutputImageType>
::SetModeToEuclid ()
{

	m_EuclidMode = ShortestPathImageFilter<TInputImageType, TOutputImageType>::EUCLID_SHORT;
}


template <class TInputImageType, class TOutputImageType>
void ShortestPathImageFilter<TInputImageType, TOutputImageType>
::SetModeToNonEuclid ()
{

	m_EuclidMode = ShortestPathImageFilter<TInputImageType, TOutputImageType>::NON_EUCLID_SHORT;
}


template <class TInputImageType, class TOutputImageType>
ShortestPathImageFilter<TInputImageType, TOutputImageType>
::ShortestPathImageFilter() :
m_EuclidMode (ShortestPathImageFilter<TInputImageType, TOutputImageType>::NON_EUCLID_SHORT),
m_FullNeighborsMode(0)
{


}


template <class TInputImageType, class TOutputImageType>
void ShortestPathImageFilter<TInputImageType, TOutputImageType>::
SetStartPoint (const typename TInputImageType::PointType &StartPoint)
{
	for (unsigned int i=0;i<TInputImageType::ImageDimension;++i)
		m_StartPoint[i] = StartPoint[i];
}

template <class TInputImageType, class TOutputImageType>
void ShortestPathImageFilter<TInputImageType, TOutputImageType>::
SetEndPoint (const typename TInputImageType::PointType &EndPoint)
{
	for (unsigned int i=0;i<TInputImageType::ImageDimension;++i)
		m_EndPoint[i] = EndPoint[i];
}


template <class TInputImageType, class TOutputImageType>
void
ShortestPathImageFilter<TInputImageType, TOutputImageType>::
GenerateData()
{

	typename itkMultiScaleVesselnessFilterType::Pointer vesselness = itkMultiScaleVesselnessFilterType::New();
		                
	vesselness->SetInput(this->GetInput());
	vesselness->SetSigmaMin(0.5);
	vesselness->SetSigmaMax(1);
	vesselness->SetNumberOfSigmaSteps(5);

	int dim=TInputImageType::ImageDimension;

	typename itkRescaleFilterType::Pointer rescale = itkRescaleFilterType::New();
	rescale->SetInput(vesselness->GetOutput());
	rescale->SetOutputMinimum(0);
	rescale->SetOutputMaximum(500);
	rescale->Update();

#ifndef NDEBUG
	typedef itk::ImageFileWriter<TOutputImageType> WriterType;
	typename WriterType::Pointer vesselnessWriter = WriterType::New();
	vesselnessWriter->SetFileName("vesselnessImage.vtk");
	vesselnessWriter->SetInput (rescale->GetOutput());
  
	try
   	 {
		vesselnessWriter->Update();
   	 }
	catch( itk::ExceptionObject & excp )
   	 {
		std::cerr << "Problem encountered while writing ";
		std::cerr << " image file : " << "vesselnessImage.vtk" << std::endl;
		std::cerr << excp << std::endl;
   	 }
#endif

	vesselness = NULL;
	
	//check for n-d
	const typename TInputImageType::SizeType &size = this->GetInput()->GetRequestedRegion().GetSize();
	int nodesNum=1;
	for (int i=0;i<dim;++i)
		nodesNum=nodesNum*size[i];		
	int edgesNum;
        if (m_FullNeighborsMode==1) 
           edgesNum = nodesNum*((int) pow((double) 3.0,dim));
        else
           edgesNum = nodesNum*(2*dim);

	
	// boost code
	typedef adjacency_list < vecS, vecS, undirectedS, no_property, property < edge_weight_t, double > > graph_t;
	typedef graph_traits < graph_t >::vertex_descriptor vertex_descriptor;
	typedef graph_traits < graph_t >::edge_descriptor edge_descriptor;
	typedef std::pair<int, int> Edge;

	
	graph_t g(nodesNum);
	std::vector <Edge> edgeVec (edgesNum);

	property_map<graph_t, edge_weight_t>::type weightmap = get (edge_weight, g);


	typedef itk::ShapedNeighborhoodIterator<TInputImageType> itkShapedNeighborhoodIteratorType;
	typename itkShapedNeighborhoodIteratorType::RadiusType radius;
	radius.Fill(1);
	itkShapedNeighborhoodIteratorType it(radius, this->GetInput(), this->GetInput()->GetRequestedRegion());

	typedef itk::ShapedNeighborhoodIterator<typename  itkRescaleFilterType::OutputImageType> 
								itkShapedNeighborhoodVesselnessIteratorType;
	typename itkShapedNeighborhoodVesselnessIteratorType::RadiusType vesselness_radius;
	vesselness_radius.Fill(1);
	itkShapedNeighborhoodVesselnessIteratorType vesselness_it(vesselness_radius, rescale->GetOutput(), 
								rescale->GetOutput()->GetRequestedRegion());

	
	//build graph
	typename TInputImageType::IndexType startIndex;
	typename TInputImageType::IndexType endIndex;
	this->GetInput()->TransformPhysicalPointToIndex (m_StartPoint, startIndex);
	this->GetInput()->TransformPhysicalPointToIndex (m_EndPoint, endIndex);
	
	int curNodeNum,toNodeNum;
	int startNodeNum=computeNodeNum(size,startIndex);
	int endNodeNum=computeNodeNum(size,endIndex);

	typename TInputImageType::PixelType dstValue;

	for (it.GoToBegin(), vesselness_it.GoToBegin(); !it.IsAtEnd(); ++it, ++vesselness_it)
	{
		
		const typename TInputImageType::IndexType &index = it.GetIndex ();
		curNodeNum=computeNodeNum(size,index);
				
		unsigned int centerIndex = it.GetCenterNeighborhoodIndex();
		
		typename ShapedNeighborhoodIterator<TInputImageType>::OffsetType offset;
		typename ShapedNeighborhoodIterator<typename itkRescaleFilterType::OutputImageType>::OffsetType
												offset_vessel;

		if (m_FullNeighborsMode==1)
			for (unsigned int d=0; d < 2*centerIndex + 1; d++)
			{
				offset = it.GetOffset(d);
				offset_vessel=vesselness_it.GetOffset(d);
				it.ActivateOffset(offset);
				vesselness_it.ActivateOffset(offset_vessel);
			}
			else
			for (int d=0; d<dim; d++)
			{
				offset = it.GetOffset(centerIndex+(int)pow((double)3.0,d));
				offset_vessel=vesselness_it.GetOffset(centerIndex+(int)pow((double)3,d));
				it.ActivateOffset(offset);
				vesselness_it.ActivateOffset(offset_vessel);
				offset = it.GetOffset(centerIndex-(int)pow((double)3.0,d));
				offset_vessel=vesselness_it.GetOffset(centerIndex-(int)pow((double)3,d));
				it.ActivateOffset(offset);
				vesselness_it.ActivateOffset(offset_vessel);
			}
					
					
		it.DeactivateOffset(it.GetOffset(centerIndex));
		vesselness_it.DeactivateOffset(vesselness_it.GetOffset(centerIndex));
    
		typename ShapedNeighborhoodIterator<TInputImageType>::Iterator org_local_it;

		for( org_local_it = it.Begin(); !org_local_it.IsAtEnd(); ++org_local_it )
		{
			unsigned int i = org_local_it.GetNeighborhoodIndex();
			bool IsInBounds;
			dstValue=it.GetPixel(i, IsInBounds);
			typename TInputImageType::IndexType idx_i = it.GetIndex( i );
			if( !IsInBounds || index == idx_i )
			{
				continue;
			}
			toNodeNum=computeNodeNum(size,idx_i);
			
			const typename itkRescaleFilterType::OutputPixelType & toVesselness  = vesselness_it.GetPixel (i);	
	
			edge_descriptor e; bool inserted;                    
			tie(e, inserted) = add_edge(curNodeNum, toNodeNum, g);
			weightmap[e] = computeDijkstraEdgeWeight (toVesselness);
					
		}
	}
	
	
	std::vector<vertex_descriptor> p(num_vertices(g));
	std::vector<int> d(num_vertices(g));
	vertex_descriptor s = vertex(startNodeNum, g);
	typename property_map<graph_t, vertex_index_t>::type indexmap = get(vertex_index, g);
	
#ifndef NDEBUG
	std::cout << "Finished building dijkstra graph" << std::endl;
#endif

	dijkstra_shortest_paths(g, s, &p[0], &d[0], weightmap, indexmap, 
                        std::less<int>(), closed_plus<int>(), 
                        (std::numeric_limits<int>::max)(), 0,
                        default_dijkstra_visitor());

#ifndef NDEBUG
	std::cout << "Done dijkstra computation" << std::endl;
#endif

	//prepare and update res Image
	
	typename itkCastImageFilterType::Pointer caster = itkCastImageFilterType::New();
	caster->SetInput (this->GetInput());
	caster->Update();

	typename TOutputImageType::Pointer outputImage = caster->GetOutput();

	typedef itk::ImageRegionIteratorWithIndex<TOutputImageType> itkImageRegionIteratorWithIndexType;
	itkImageRegionIteratorWithIndexType resIt (outputImage, outputImage->GetRequestedRegion());

	for (resIt.GoToBegin(); !resIt.IsAtEnd(); ++resIt)
	{
		resIt.Value() = 0;
	}

	
	vertex_descriptor c = vertex(endNodeNum, g);
	typename TOutputImageType::IndexType curIndex;
	for (unsigned int i=0;i<TOutputImageType::ImageDimension;++i)
		curIndex [i] = endIndex[i];
	
	resIt.SetIndex (curIndex);
	resIt.Value() = 1;
	while (indexmap [c] != indexmap[s])
	{
		int tmp_c = indexmap[c];
		int temp_val;
		for (int k=dim-1;k>=0;--k){
			temp_val=1;
			for (int l=0;l<k;++l)
				temp_val=temp_val*size[l];
			curIndex[k]=tmp_c/temp_val;
			tmp_c=tmp_c-curIndex[k]*temp_val;
		}
		
		resIt.SetIndex (curIndex);
		resIt.Value() = 1;
                c  = vertex(indexmap[p[c]], g); //compute parent index
	}

	resIt.SetIndex (startIndex);
	resIt.Value() = 1;
	
	//update res image
	this->GraftOutput(outputImage);

} //Generate Data



template <class TInputImageType, class TOutputImageType>
void
ShortestPathImageFilter<TInputImageType, TOutputImageType>::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
  
  os << indent << "to do" << std::endl;
}

} // end namespace itk 
