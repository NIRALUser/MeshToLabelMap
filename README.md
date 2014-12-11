#MeshToLabelMap

##What is it?

MeshToLabelMap converts VTK polydata (*.vtk or *.vtp) into label maps. It uses vtkImageStencil as well as the specific vtk classes included in this project to perform this operation. This tool uses SlicerExecutionModel to generate its user interface. It can be integrated into [3D Slicer](http://www.slicer.org) to have access to a graphical user interface.
To print the help page, use the flag '--help'

Example of usage:

```
./MeshToLabelMap -m model.vtk -R reference_image.nrrd -l my_output_label_map.nrrd --pixel_value 1
./MeshToLabelMap -m model.vtk --spacing 0.2,0.2,.2 -l my_output_label_map.nrrd --pixel_value 2
```

##Compilation

MeshToLabelMap requires:
- [CMake](http:www.cmake.org) 2.8.10 or more recent
- [SlicerExecutionModel](https://github.com/Slicer/SlicerExecutionModel)
- [ITK](http:www.itk.org)
- [VTK](http:www.vtk.org)

##License

See License.txt

##Source code

Find the source code on [GitHub](https://github.com/NIRALUser/MeshToLabelMap)

