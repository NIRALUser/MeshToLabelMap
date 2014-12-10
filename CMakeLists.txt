project(MeshToLabelMap)

cmake_minimum_required(VERSION 2.8.10 )

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

