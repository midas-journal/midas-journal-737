cmake_minimum_required(VERSION 2.6)


PROJECT (vesselnessSegmenter)


INCLUDE(${CMAKE_SOURCE_DIR}/IJMacros.txt)

ENABLE_TESTING()

INCLUDE(Dart)

SET(PROJECT_HOME ${CMAKE_SOURCE_DIR})

SET(RESULTS_PATH ${CMAKE_SOURCE_DIR}/Results)
SET(DATA_PATH ${CMAKE_SOURCE_DIR}/Data)

SET(Boost_USE_STATIC_LIBS ON)
SET(Boost_USE_MULTITHREAD OFF)
FIND_PACKAGE( Boost 1.37.0)

IF(Boost_FOUND)

   INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})

ELSE(Boost_FOUND)

  MESSAGE(FATAL_ERROR

          "Cannot build applications without Boost.  Please set BOOST_INCLUDEDIR.")

ENDIF(Boost_FOUND)




FIND_PACKAGE(ITK)

IF(ITK_FOUND)

   INCLUDE(${ITK_USE_FILE})

ELSE(ITK_FOUND)

   MESSAGE(FATAL_ERROR

           "Cannot build InsightApplications without ITK.  Please set

ITK_DIR.")

ENDIF(ITK_FOUND)



SET(Required_Packages

ITK

)

FOREACH(Package ${Required_Packages})
  
	LOADPACKAGE(${Package})

ENDFOREACH(Package)

SET (ITK_LIBS

     ITKBasicFilters

     ITKAlgorithms

     ITKCommon

     ITKIO

     ITKNumerics

     ITKNrrdIO

     itkpng

     itktiff

     itkgdcm

     itkjpeg8

     itkjpeg12

     itkjpeg16

     ITKMetaIO

     itkzlib

     ITKDICOMParser

     ITKEXPAT

     itkvnl_inst

     itkvnl_algo

     itkvnl

     itkvcl

     itksys

)









SET( vesselnessSegmenter_3D_SRCS


	src/itkTwoPointsVesselSegmenterImageFilter.h
	src/itkTwoPointsVesselSegmenterImageFilter.txx
	src/itkShortestPathImageFilter.txx
	src/itkShortestPathImageFilter.h
	src/itkMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter.txx
	src/itkMultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter.h
	src/itkHessianSmoothed3DToVesselnessMeasureImageFilter.txx
	src/itkHessianSmoothed3DToVesselnessMeasureImageFilter.h
	src/itkGMMEstimatorImageFilter.txx
	src/itkGMMEstimatorImageFilter.h
	src/itkActiveContourMinCutImageFilter.h
	src/itkActiveContourMinCutImageFilter.txx

	src/3D_test.cpp


)


SET( vesselnessSegmenter_3D_LIB

    3D_test
    ${ITK_LIBS}


)

SET (IMAGE_COMPARE_LIB
	ImageCompare
	${ITK_LIBS}


)


ADD_EXECUTABLE(3D_test ${vesselnessSegmenter_3D_SRCS} )

TARGET_LINK_LIBRARIES( ${vesselnessSegmenter_3D_LIB})

ADD_EXECUTABLE(ImageCompare src/ImageCompare.cxx) 
TARGET_LINK_LIBRARIES( ${IMAGE_COMPARE_LIB})

