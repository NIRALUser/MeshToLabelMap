project(MeshToLabelMap)

cmake_minimum_required(VERSION 2.8.10 )


if( MeshToLabelMap_BUILD_SLICER_EXTENSION )
  set( EXTENSION_NAME MeshToLabelMap )
  set( EXTENSION_HOMEPAGE
"http://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/MeshToLabelMap" )
  set( EXTENSION_CATEGORY "Shape Analysis" )
  set( EXTENSION_CONTRIBUTORS "Francois Budin (UNC), Ipek Oguz (University of
Iowa)" )
  set( EXTENSION_DESCRIPTION "This extension computes a label map from a 3D
model." )
  set( EXTENSION_ICONURL
"https://raw.githubusercontent.com/NIRALUser/MeshToLabelMap/master/ico
ns/MeshToLabelMapIcon128x128.png" )
  set( EXTENSION_SCREENSHOTURLS 
"http://slicer.org/slicerWiki/images/thumb/5/54/MeshToLabelMapScre
enshot.png/800px-MeshToLabelMapScreenshot.png" )
  set( EXTENSION_STATUS "Release" )
  set( EXTENSION_DEPENDS "NA" ) # Specified as a space separated list or 'NA' if any
  set( EXTENSION_BUILD_SUBDIRECTORY . )
  find_package( Subversion REQUIRED )
  find_package( Git REQUIRED )
  find_package( Slicer REQUIRED )
  include( ${Slicer_USE_FILE} )
endif()

find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})
find_package (VTK REQUIRED)
include (${VTK_USE_FILE})
find_package (ITK REQUIRED)
include (${ITK_USE_FILE})


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
ENDIF(BUILD_TESTING)
