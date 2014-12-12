
/* 
 * This is a subclass of vtkPolyDataToImageStencil, with extended functionality to build a lookup from the 
 * stencil back to the triangles of the mesh.
 *
 * Author:  Ipek Oguz 
 *
 */

#include "vtkAttributedPolyDataToImage.h"

#include "vtkGarbageCollector.h"
#include "vtkImageStencilData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkOBBTree.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <math.h>

vtkStandardNewMacro(vtkAttributedPolyDataToImage);


vtkAttributedPolyDataToImage::vtkAttributedPolyDataToImage()
{
  this->BinaryVolume = NULL ;
  this->AttributeVolume = NULL ;
  this->Attributes = NULL ;
  this->faceList = vtkIdTypeArray::New () ;
  this->pointList = vtkPoints::New () ;
  this->ScanConvertPerformed = false ;
  this->mesh = NULL ;
  this->stencil = NULL ;
}

//----------------------------------------------------------------------------
vtkAttributedPolyDataToImage::~vtkAttributedPolyDataToImage()
{  
  if ( this->faceList )
    this->faceList->Delete () ;
  if ( this->pointList ) 
    this->pointList->Delete () ;
  if ( this->BinaryVolume ) 
    this->BinaryVolume->Delete () ;

  if ( this->AttributeVolume ) 
    this->AttributeVolume->Delete () ;
  if ( this->stencil )
    {
      //this->GetOutput()->Delete() ;
      //if ( this->stencil->GetStencil() )
      //this->stencil->GetStencil()->Delete () ;
    //  this->stencil->Delete () ;
    }
}

//----------------------------------------------------------------------------
//void vtkAttributedPolyDataToImage::PrintSelf(ostream& os, vtkIndent indent)
//{
//  this->Superclass::PrintSelf(os,indent);
//}

//----------------------------------------------------------------------------
int vtkAttributedPolyDataToImage::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
   // the following is largely copied from the vtkPolyDataToImageStencil code

   this->SuperSuperclass::RequestData(request, inputVector, outputVector);

   vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
   vtkInformation *outInfo = outputVector->GetInformationObject(0);

   // need to build the OBB tree
   vtkPolyData *polydata = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
   vtkImageStencilData *data = vtkImageStencilData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

   this->mesh = polydata ;

   if (this->OBBTree == NULL)
    {
    this->OBBTree = vtkOBBTree::New();
    }
  this->OBBTree->SetDataSet(polydata);
  this->OBBTree->SetTolerance(this->Tolerance);
  this->OBBTree->BuildLocator();

  // for keeping track of progress
  unsigned long count = 0;
  int extent[6];
  data->GetExtent(extent);
  unsigned long target = (unsigned long)
    ((extent[5] - extent[4] + 1)*(extent[3] - extent[2] + 1)/50.0);
  target++;

  // if we have no data then return
  if (!polydata->GetNumberOfPoints())
    {
    return 1;
    }
  
  double *spacing = data->GetSpacing();
  double *origin = data->GetOrigin();

  if( this->GetDebug() )
  {
      std::cout << "Spacing: " << spacing[0] << " " << spacing[1] << " " << spacing[2] << std::endl ;
      std::cout << "Origin: " << origin[0] << " " << origin[1] << " " << origin[2] << std::endl ;
      std::cout << "Extent: " << extent[0] << " " << extent[1] << " " << extent[2] << " " << extent[3] << " " << extent[4] << " " << extent[5] << std::endl ;
  }
  vtkOBBTree *tree = this->OBBTree;
  vtkPoints *points = vtkPoints::New();

  double p0[3],p1[3];

  p1[0] = p0[0] = extent[0]*spacing[0] + origin[0];
  p1[1] = p0[1] = extent[2]*spacing[1] + origin[1];
  p0[2] = extent[4]*spacing[2] + origin[2];
  p1[2] = extent[5]*spacing[2] + origin[2];

  int zstate = tree->InsideOrOutside(p0);
  if (zstate == 0)
    {
    zstate = -1;
    }
  int *zlist = 0;
  int zlistlen = 0;
  int zlistidx = 0;
  if (extent[4] < extent[5])
    {
    tree->IntersectWithLine(p0, p1, points, 0);
    vtkTurnPointsIntoList(points, zlist, zlistlen,
                          extent, origin, spacing, 2);
    }

  for (int idZ = extent[4]; idZ <= extent[5]; idZ++)
    {
    if (zlistidx < zlistlen && idZ >= zlist[zlistidx])
      {
      zstate = -zstate;
      zlistidx++;
      }

    p1[0] = p0[0] = extent[0]*spacing[0] + origin[0];
    p0[1] = extent[2]*spacing[1] + origin[1];
    p1[1] = extent[3]*spacing[1] + origin[1];
    p1[2] = p0[2] = idZ*spacing[2] + origin[2];

    int ystate = zstate;
    int *ylist = 0;
    int ylistlen = 0;
    int ylistidx = 0;
    if (extent[2] != extent[3])
      {
      tree->IntersectWithLine(p0, p1, points, 0);
      vtkTurnPointsIntoList(points, ylist, ylistlen,
                            extent, origin, spacing, 1);
      }

    for (int idY = extent[2]; idY <= extent[3]; idY++)
      {
      if (ylistidx < ylistlen && idY >= ylist[ylistidx])
        {
        ystate = -ystate;
        ylistidx++;
        }

      if (count%target == 0) 
        {
        this->UpdateProgress(count/(50.0*target));
        }
      count++;

      p0[1] = p1[1] = idY*spacing[1] + origin[1];
      p0[2] = p1[2] = idZ*spacing[2] + origin[2];
      p0[0] = extent[0]*spacing[0] + origin[0];
      p1[0] = extent[1]*spacing[0] + origin[0];

      int xstate = ystate;
      int *xlist = 0;
      int xlistlen = 0;
      int xlistidx = 0;

      //tree->IntersectWithLine(p0, p1, points, 0);
      vtkIdList *faces = vtkIdList::New () ;
      tree->IntersectWithLine ( p0, p1, points, faces ) ;
      int nFaces = faces->GetNumberOfIds () ;
      int i ;
      for ( i = 0 ; i < nFaces ; i++ )
      {
        this->faceList->InsertNextValue ( faces->GetId ( i ) ) ;
      }
      int nPoints = points->GetNumberOfPoints () ;
      if ( this->GetDebug() && nPoints != nFaces )
      {
          std::cout << "Problem: nPoints != nFaces. nPoints=" << nPoints << "-- nFaces=" << nFaces << std::endl ;
      }
      double p[3] ;
      for ( i = 0 ; i < nPoints ; i++ )
      {
        points->GetPoint ( i, p ) ;
        this->pointList->InsertNextPoint ( p ) ;
      }
      faces->Delete () ;

      vtkTurnPointsIntoList(points, xlist, xlistlen,
                            extent, origin, spacing, 0);

      // now turn 'xlist' into sub-extents:
      int r1 = extent[0];
      int r2 = extent[1];
      for (xlistidx = 0; xlistidx < xlistlen; xlistidx++)
        {
        xstate = -xstate;
        if (xstate < 0)
          { // sub extent starts
          r1 = xlist[xlistidx];
          }
        else
          { // sub extent ends
          r2 = xlist[xlistidx] - 1;
          data->InsertNextExtent(r1, r2, idY, idZ);
          }
        }
      if (xstate < 0)
        { // if inside at end, cap off the sub extent
        data->InsertNextExtent(r1, extent[1], idY, idZ);
        }      

      if (xlist)
        {
        delete [] xlist;
        }

      } // for idY

    if (ylist)
      {
      delete [] ylist;
      }

    } // for idZ

  if (zlist)
    {
    delete [] zlist;
    }
  points->Delete();

  this->ScanConvertPerformed = true ;
  return 1 ;
}

void vtkAttributedPolyDataToImage::ComputeAttributeVolume () 
{
  double spacing[3], origin[3] ;
  int *extent ; 
  vtkImageStencilData *stencil = this->GetOutput () ;
  stencil->GetSpacing ( spacing ) ;
  stencil->GetOrigin ( origin ) ;
  extent= this->GetOutputWholeExtent () ;

  // create empty image
  this->AttributeVolume = vtkImageData::New () ;
  this->AttributeVolume->SetOrigin ( origin ) ;
  this->AttributeVolume->SetSpacing ( spacing ) ;
  this->AttributeVolume->SetDimensions ( extent[1] - extent[0] + 1, extent[3] - extent[2] + 1, extent[5] - extent[4] + 1 ) ;
  this->AttributeVolume->AllocateScalars(VTK_FLOAT, 1);

  for ( int i = 0 ; i <= extent[1] ; i++ )
  {
    for ( int j = 0 ; j <= extent[3] ; j++ )
    {
      for ( int k = 0 ; k <= extent[5] ; k++ )
      {
        this->AttributeVolume->SetScalarComponentFromDouble ( i, j, k, 0, 0.0 ) ;
      }
    }
  }


  int nPts = this->pointList->GetNumberOfPoints () ;
  double p[3]; 
  vtkIdType faceId ;
  double closestPoint[3], dist2, pCoords[3], weights[3] ;
  int subId, result ;

  for ( int i = 0 ; i < nPts ; i++ )
  {
    this->pointList->GetPoint ( i, p ) ;
    faceId = this->faceList->GetValue ( i ) ;

    vtkIdList *facePoints = vtkIdList::New () ;
    this->mesh->GetCellPoints ( faceId, facePoints ) ;
    result = this->mesh->GetCell ( faceId )->EvaluatePosition ( p, closestPoint, subId, pCoords, dist2, weights ) ;
    if( this->GetDebug() )
    {
        if ( result == 0 )
        {
            std::cout << "outside" << std::endl ;
        }
        else if ( result == -1 )
        {
            std::cout << "numerical error" << std::endl ;
        }
    }
    int gridCoords[3] ;
    double attributeValue, vertAttributes[3] ;
    attributeValue = 0 ;
    for ( int j = 0 ; j < 3 ; j++ )
      {
	gridCoords[j] = ceil ( ( p[j] - origin[j] ) / spacing[j] ) ;
	vertAttributes[j] = this->Attributes->GetValue ( facePoints->GetId ( j ) )  ;
	attributeValue += vertAttributes[j] * pCoords[j] ;
      }
    //std::cout << attributeValue << std::endl ;
    if ( !(i % 2) )
      gridCoords[0]-- ; 

    this->AttributeVolume->SetScalarComponentFromDouble ( gridCoords[0], gridCoords[1], gridCoords[2], 0, attributeValue );
    facePoints->Delete () ;
  }
  //std::cout << "Attribute volume computed. " << std::endl ;
}

vtkImageData * vtkAttributedPolyDataToImage::GetBinaryVolume()
{
  if ( !this->ScanConvertPerformed ) 
    return NULL ;

  double spacing[3], origin[3] ;

  int *extent, size[3] ; 
  vtkImageStencilData *temp = this->GetOutput () ;
  temp->GetSpacing ( spacing ) ;
  temp->GetOrigin ( origin ) ;
  extent = this->GetOutputWholeExtent () ;
  size[0] = extent[1] - extent[0] + 1 ;
  size[1] = extent[3] - extent[2] + 1 ;
  size[2] = extent[5] - extent[4] + 1 ;

  // Create an empty image
  vtkImageData *emptyImage = vtkImageData::New () ;
  emptyImage->SetOrigin ( origin ) ;
  emptyImage->SetSpacing ( spacing ) ;
  emptyImage->SetDimensions ( size ) ;
  emptyImage->AllocateScalars(VTK_INT, 1);
  
  // Use the stencil as a cookie cutter
  this->stencil = vtkImageStencil::New () ;
  this->stencil->SetInputData ( emptyImage ) ;
  //this->GetOutput() ;
  
  this->stencil->SetStencilData ( this->GetOutput () ) ;
  
  this->stencil->ReverseStencilOn () ;
  this->stencil->SetBackgroundValue ( 128 ) ;
  this->stencil->Update () ;
  this->BinaryVolume = this->stencil->GetOutput () ;
  //this->stencil->Delete() ;

  emptyImage->Delete () ;

  return this->BinaryVolume ;
}

vtkImageData * vtkAttributedPolyDataToImage::GetAttributeVolume()
{
  if ( ( !this->ScanConvertPerformed ) || ( !this->Attributes ) ) 
    return NULL ;

  this->ComputeAttributeVolume () ;
  return this->AttributeVolume ;
}




