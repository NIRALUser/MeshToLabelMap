project(MeshToLabelMap)

cmake_minimum_required(VERSION 2.8.10 )

find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})
include(${GenerateCLP_USE_FILE})
find_package (VTK REQUIRED)
include (${VTK_USE_FILE})
find_package (ITK REQUIRED)
include (${ITK_USE_FILE})

INCLUDE_DIRECTORIES( 
${CMAKE_CURRENT_SOURCE_DIR}
 ${CMAKE_CURRENT_BINARY_DIR}
 )


# This block builds a module as both a shared object and an executable

set(MyModule_SOURCE 
main.cxx 
vtkAttributedPolyDataToImage.h
vtkAttributedPolyDataToImage.cxx
vtkPolyDataToImageStencilOBBTree.cxx
vtkPolyDataToImageStencilOBBTree.h
)

set (VTK_LIBRARIES
  vtkRendering
  vtkGraphics
  vtkImaging
  vtkIO
  vtkCommon
  vtkFiltering 
) 



GENERATECLP(MyModule_SOURCE MeshToLabelMap.xml)

add_executable(MeshToLabelMap ${MyModule_SOURCE})
target_link_libraries( MeshToLabelMap ${ITK_LIBRARIES} ${VTK_LIBRARIES})

