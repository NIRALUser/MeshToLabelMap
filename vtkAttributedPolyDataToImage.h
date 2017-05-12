/* 
 * This is a subclass of vtkPolyDataToImageStencilOBBTree, with extended functionality to build a lookup from the 
 * stencil back to the triangles of the mesh.
 *
 * Author:  Ipek Oguz 
 *
 */

#ifndef _vtkAttributedPolyDataToImage_h
#define _vtkAttributedPolyDataToImage_h

#include "vtkPolyDataToImageStencilOBBTree.h"
#include "vtkIdList.h"
#include "vtkIdTypeArray.h"
#include "vtkPoints.h"
#include "vtkCell.h"
#include "vtkImageData.h"
#include "vtkFloatArray.h"
#include "vtkImageStencil.h"
#include "vtkSmartPointer.h"

class vtkAttributedPolyDataToImage : public vtkPolyDataToImageStencilOBBTree
{
  public:

  static vtkAttributedPolyDataToImage *New();
  vtkTypeMacro(vtkAttributedPolyDataToImage, vtkPolyDataToImageStencilOBBTree);
  //void PrintSelf(ostream& os, vtkIndent indent);

  vtkSetMacro(Attributes, vtkFloatArray *) ;

  vtkSmartPointer<vtkImageData> GetAttributeVolume () ;
  vtkSmartPointer<vtkImageData> GetBinaryVolume () ;
  
  
  protected:

    vtkAttributedPolyDataToImage();
    ~vtkAttributedPolyDataToImage();

    virtual int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  
    //void BuildLookup () ;
    void ComputeAttributeVolume () ;


  private:
  vtkAttributedPolyDataToImage (const vtkAttributedPolyDataToImage&);  // Not implemented.
  void operator=(const vtkAttributedPolyDataToImage&);  // Not implemented.

  vtkSmartPointer<vtkImageData> BinaryVolume, AttributeVolume ;
  vtkSmartPointer<vtkFloatArray> Attributes ;

  vtkSmartPointer<vtkIdTypeArray> faceList ;
  vtkSmartPointer<vtkPoints> pointList ;
  vtkSmartPointer<vtkPolyData> mesh ;
  vtkSmartPointer<vtkImageStencil> stencil ;
  typedef vtkImageStencilSource SuperSuperclass ;
  bool ScanConvertPerformed ;
} ;

#endif
