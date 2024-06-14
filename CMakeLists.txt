cmake_minimum_required(VERSION 3.16.3...3.19.7 FATAL_ERROR)

project(MeshToLabelMap)

if( MeshToLabelMap_BUILD_SLICER_EXTENSION )
  set( EXTENSION_NAME MeshToLabelMap )
  set( EXTENSION_HOMEPAGE "https://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/MeshToLabelMap" )
  set( EXTENSION_CATEGORY "Shape Analysis" )
  set( EXTENSION_CONTRIBUTORS "Francois Budin (UNC), Ipek Oguz (University of Iowa)" )
  set( EXTENSION_DESCRIPTION "This extension computes a label map from a 3D model." )
  set( EXTENSION_ICONURL "https://raw.githubusercontent.com/NIRALUser/MeshToLabelMap/master/icons/MeshToLabelMapIcon128x128.png" )
  set( EXTENSION_SCREENSHOTURLS "https://slicer.org/slicerWiki/images/thumb/5/54/MeshToLabelMapScreenshot.png/800px-MeshToLabelMapScreenshot.png" )
  set( EXTENSION_STATUS "Release" )
  set( EXTENSION_DEPENDS "NA" ) # Specified as a space separated list or 'NA' if any
  set( EXTENSION_BUILD_SUBDIRECTORY . )
  find_package( Git REQUIRED )
  find_package( Slicer REQUIRED )
  include( ${Slicer_USE_FILE} )
  set( Build_Static OFF CACHE BOOL "Static libraries and executables only")
endif()

find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})

#
# ITK
#
set( ALL_IO
ITKIOBMP
ITKIOBioRad
ITKIOCSV
ITKIODCMTK
ITKIOGDCM
ITKIOGE
ITKIOGIPL
ITKIOHDF5
ITKIOIPL
ITKIOJPEG
ITKIOLSM
ITKIOMRC
ITKIOMesh
ITKIOMeta
ITKIONIFTI
ITKIOPNG
ITKIORAW
ITKIOSiemens
ITKIOSpatialObjects
ITKIOStimulate
ITKIOTIFF
ITKIOTransformHDF5
ITKIOTransformMatlab
ITKIOVTK
ITKIOXML
MGHIO
)

find_package(ITK REQUIRED )
list(APPEND LIST_ITK_IO_USED ITKIONRRD ITKIOTransformInsightLegacy)
foreach( io ${ALL_IO} )
  list( FIND ITK_MODULES_ENABLED ${io} position )
  if( ${position} GREATER -1 ) #not found: ${position}==-1
    list( APPEND LIST_ITK_IO_USED ${io} )
  endif()
endforeach()

set(ITKComponents
  ITKCommon
  ITKIOImageBase
  ITKSmoothing
  ${LIST_ITK_IO_USED}
  )

find_package(ITK COMPONENTS ${ITKComponents} REQUIRED )
include(${ITK_USE_FILE})


find_package(VTK COMPONENTS
  vtkCommonComputationalGeometry
  vtkCommonCore
  vtkCommonDataModel
  vtkCommonExecutionModel
  vtkCommonMath
  vtkCommonMisc
  vtkCommonSystem
  vtkCommonTransforms
  vtkFiltersCore
  vtkFiltersGeneral
  vtkIOCore
  vtkIOGeometry
  vtkIOLegacy
  vtkIOXML
  vtkIOXMLParser
  vtkImagingCore
  vtkImagingStencil
)
include (${VTK_USE_FILE})


# This block builds a module as both a shared object and an executable

set(MyModule_SOURCE 
vtkAttributedPolyDataToImage.h
vtkAttributedPolyDataToImage.cxx
vtkPolyDataToImageStencilOBBTree.cxx
vtkPolyDataToImageStencilOBBTree.h
)

option( Build_Static "Static libraries and executables only" ON )
if( Build_Static )
  set( STATIC "EXECUTABLE_ONLY" )
endif()

SEMMacroBuildCLI(
NAME MeshToLabelMap
${STATIC}
TARGET_LIBRARIES ${ITK_LIBRARIES} ${VTK_LIBRARIES}
ADDITIONAL_SRCS ${MyModule_SOURCE}
INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
)

option( BUILD_TESTING "Build Testing" ON )
IF(BUILD_TESTING)
  include(CTest)
  include(ExternalData)
  ADD_SUBDIRECTORY(Testing)
  ADD_SUBDIRECTORY(MeshToLabelMapPipeline)
ENDIF(BUILD_TESTING)

if( MeshToLabelMap_BUILD_SLICER_EXTENSION )
  set( CPACK_INSTALL_CMAKE_PROJECTS "${CPACK_INSTALL_CMAKE_PROJECTS};${CMAKE_BINARY_DIR};${EXTENSION_NAME};ALL;/" )
  include(${Slicer_EXTENSION_GENERATE_CONFIG})
  include( ${Slicer_EXTENSION_CPACK} )
endif()
