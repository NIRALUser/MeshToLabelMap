cmake_minimum_required(VERSION 2.8)
CMAKE_POLICY(VERSION 2.8)

set(CLP MeshToLabelMap )

# Add standard remote object stores to user's
# configuration.
list( APPEND ExternalData_URL_TEMPLATES
"http://slicer.kitware.com/midas3/api/rest?method=midas.bitstream.download&checksum=%(hash)"
)

set( TestingBinarydirectory ${CMAKE_CURRENT_BINARY_DIR}/Temporary )
set( DummyFile ${TestingBinarydirectory}/dummy.nrrd )#dummy file because it is never created. Only used for tests that are designed to fail
set( InputVTK DATA{Data/Input/model.vtk} )
set( referenceImage DATA{Data/Input/reference.nrrd} )

# Should just display the "Usage"
add_test(NAME ${CLP}Help COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}> --help)
# Verify that we cannot give both the spacing information and a reference image. Test fails by design.
ExternalData_Add_Test( ${CLP}Data NAME ${CLP}SpacingAndReference COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}>
 -m ${InputVTK}
 --spacing -1,-2,-3
 -l ${DummyFile}
 -R ${referenceImage}
)
set_tests_properties( ${CLP}SpacingAndReference PROPERTIES WILL_FAIL true)
# Spacing should be positive. Test fails by design.
ExternalData_Add_Test( ${CLP}Data NAME ${CLP}NegativeSpacing COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}>
 -m ${InputVTK}
 --spacing -1,-2,-3
 -l ${DummyFile}
)
set_tests_properties( ${CLP}NegativeSpacing PROPERTIES WILL_FAIL true)
# Smoothing radius should be >= 1. Test fails by design.
ExternalData_Add_Test( ${CLP}Data NAME ${CLP}BadSmoothingRadius COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}>
 -m ${InputVTK}
 --median
 --median_radius 0,.7,2
 -R ${referenceImage}
 -l ${DummyFile}
)
set_tests_properties( ${CLP}BadSmoothingRadius PROPERTIES WILL_FAIL true)
# No input mesh. Test fails by design.
ExternalData_Add_Test( ${CLP}Data NAME ${CLP}NoInputMesh COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}>
 -l ${DummyFile}
)
set_tests_properties( ${CLP}NoInputMesh PROPERTIES WILL_FAIL true)
# No output label map. Test fails by design.
ExternalData_Add_Test( ${CLP}Data NAME ${CLP}NoOutputLabelMap COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}>
 -m ${InputVTK}
)
set_tests_properties( ${CLP}NoOutputLabelMap PROPERTIES WILL_FAIL true)
# Not 3 components for spacing. Test fails by design.
ExternalData_Add_Test( ${CLP}Data NAME ${CLP}SpacingTwoComponents COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}>
 -m ${InputVTK}
 --spacing 1,1
 -l ${DummyFile}
)
set_tests_properties( ${CLP}SpacingTwoComponents PROPERTIES WILL_FAIL true)
# Not 3 components for smoothing radius. Test fails by design.
ExternalData_Add_Test( ${CLP}Data NAME ${CLP}SmoothingRadiusFourComponents COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}>
 -m ${InputVTK}
 --spacing 1,1,1
 --median_radius 1,1,1,1
 -l ${DummyFile}
)
set_tests_properties( ${CLP}SmoothingRadiusFourComponents PROPERTIES WILL_FAIL true)
# Not 3 components for boundary extension. Test fails by design.
ExternalData_Add_Test( ${CLP}Data NAME ${CLP}BoundaryExtensionOneComponents COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}>
 -m ${InputVTK}
 --spacing 1,1,1
 --boundary_extension 1
 -l ${DummyFile}
)
set_tests_properties( ${CLP}BoundaryExtensionOneComponents PROPERTIES WILL_FAIL true)
# Output label value < 1. Test fails by design.
ExternalData_Add_Test( ${CLP}Data NAME ${CLP}LabelValueTooSmall COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}>
 -m ${InputVTK}
 --spacing 1,1,1
 -l ${DummyFile}
 --pixel_value 0
)
set_tests_properties( ${CLP}LabelValueTooSmall PROPERTIES WILL_FAIL true)
# Output label value > 255. Test fails by design.
ExternalData_Add_Test( ${CLP}Data NAME ${CLP}LabelValueTooLarge COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}>
 -m ${InputVTK}
 --spacing 1,1,1
 -l ${DummyFile}
 --pixel_value 256
)
set_tests_properties( ${CLP}LabelValueTooLarge PROPERTIES WILL_FAIL true)
# Only input mesh and output label map given but no spacing or reference image. Test fails by design.
ExternalData_Add_Test( ${CLP}Data NAME ${CLP}NoSpacing COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}>
 -m ${InputVTK}
 -l ${DummyFile}
)
set_tests_properties( ${CLP}NoSpacing PROPERTIES WILL_FAIL true)
# Bad input mesh. Test fails by design.
ExternalData_Add_Test( ${CLP}Data NAME ${CLP}BadInput COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}>
 -m ${DummyFile}
 -l ${DummyFile}
 --spacing 1,1,1
)
set_tests_properties( ${CLP}BadInput PROPERTIES WILL_FAIL true)

if( NOT VTK_DEBUG_LEAKS )
  # Bad output label map. Test fails by design.
  set( BadOutputFile ${TestingBinarydirectory}/dummy )#No extension. Should create a problem when writing output file
  ExternalData_Add_Test( ${CLP}Data NAME ${CLP}BadOutput COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${CLP}>
   -m ${InputVTK}
   -l ${BadOutputFile}
   --spacing 1,1,1
  )
  set_tests_properties( ${CLP}BadOutput PROPERTIES WILL_FAIL true)


  #Build test that will compare baseline and output label map
  add_executable(${CLP}Test ${CLP}Test.cxx)
  target_link_libraries(${CLP}Test ${CLP}Lib )
  set_target_properties(${CLP}Test PROPERTIES LABELS ${CLP})

  #Compute label map using the boundary box computed on the mesh. No reference image
  set( output${CLP}LabelMapWithBB ${TestingBinarydirectory}/LabelMapWithBB.nrrd )
  ExternalData_Add_Test( ${CLP}Data
    NAME ${CLP}TestLabelMapWithBB
    COMMAND $<TARGET_FILE:${CLP}Test>
    --compare DATA{Data/Baseline/baselineLabelMapWithBB.nrrd}
          ${output${CLP}LabelMapWithBB}
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0
    ModuleEntryPoint
     --input_mesh ${InputVTK}
     --output_labelmap ${output${CLP}LabelMapWithBB}
     --pixel_value 2
     --spacing 0.1,0.15,0.2
     --verbose
    )


  #Compute label map using the boundary box computed on the mesh and extend bounding box by 10x15x20 radius. No reference image
  set( output${CLP}LabelMapWithBBExt101520 ${TestingBinarydirectory}/LabelMapWithBBExt101520.nrrd )
  ExternalData_Add_Test( ${CLP}Data
    NAME ${CLP}TestLabelMapWithBBExt101520
    COMMAND $<TARGET_FILE:${CLP}Test>
    --compare DATA{Data/Baseline/baselineLabelMapWithBBExt101520.nrrd}
          ${output${CLP}LabelMapWithBBExt101520}
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0
    ModuleEntryPoint
     --input_mesh ${InputVTK}
     --output_labelmap ${output${CLP}LabelMapWithBBExt101520}
     --pixel_value 2
     --spacing 0.2,0.2,0.3
     --boundary_extension 10,15,20
     --verbose
    )

  #Compute label map using a reference image
  set( output${CLP}LabelMapWithReference ${TestingBinarydirectory}/LabelMapWithReference.nrrd )
  ExternalData_Add_Test( ${CLP}Data
    NAME ${CLP}TestLabelMapWithReference
    COMMAND $<TARGET_FILE:${CLP}Test>
    --compare DATA{Data/Baseline/baselineLabelMapWithReference.nrrd}
          ${output${CLP}LabelMapWithReference}
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0
    ModuleEntryPoint
     --input_mesh ${InputVTK}
     --output_labelmap ${output${CLP}LabelMapWithReference}
     --pixel_value 3
     --reference_volume ${referenceImage}
     --verbose
    )
  #Compute label map using the boundary box computed on the mesh and applies median filter. No reference image
  set( output${CLP}LabelMapWithBBMedianFilter ${TestingBinarydirectory}/LabelMapWithBBMedianFilter.nrrd )
  ExternalData_Add_Test( ${CLP}Data
    NAME ${CLP}TestLabelMapWithBBMedianFilter
    COMMAND $<TARGET_FILE:${CLP}Test>
    --compare DATA{Data/Baseline/baselineLabelMapWithBBMedianFilter.nrrd}
          ${output${CLP}LabelMapWithBBMedianFilter}
    --compareIntensityTolerance 0
    --compareRadiusTolerance 0
    --compareNumberOfPixelsTolerance 0
    ModuleEntryPoint
     --input_mesh ${InputVTK}
     --output_labelmap ${output${CLP}LabelMapWithBBMedianFilter}
     --pixel_value 4
     --spacing 0.5,0.1,0.4
     --median
     --median_radius 3,3,3
     --verbose
    )
endif()

# Add a build target to populate the real data.
ExternalData_Add_Target( ${CLP}Data )
