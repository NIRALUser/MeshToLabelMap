#include "MeshToLabelMapCLP.h"
#include <iostream>
#include <vector>
#include <string>

#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include <itkMedianImageFilter.h>

#include <vtkSmartPointer.h>
#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include <vtkXMLPolyDataReader.h>
#include "vtkImageData.h"
#include <vtkCommand.h>
#include "vtkAttributedPolyDataToImage.h"


void WriteITKImage ( itk::Image < unsigned char, 3 >::Pointer image, std::string fileName ) 
{
  typedef itk::Image < unsigned char, 3 > ImageType ;
  typedef itk::ImageFileWriter < ImageType > ImageWriterType ; 
  ImageWriterType::Pointer itkWriter = ImageWriterType::New() ;
  itkWriter->SetFileName ( fileName ) ;
  itkWriter->SetInput ( image ) ;
  itkWriter->UseCompressionOn() ;
  itkWriter->Write() ;
}

itk::Image < unsigned char, 3 >::Pointer VTK2BinaryITK ( vtkImageData *vtkImage ) 
{
  typedef itk::Image < unsigned char, 3 > ImageType ;
  typedef ImageType::Pointer ImagePointer ;
  // convert the vtk image to an itk image
  int extent[ 6 ] ;
  vtkImage->GetExtent( extent ) ;
  ImageType::SizeType size ;
  size[0] = extent[ 1 ] - extent[0] + 1 ;
  size[1] = extent[ 3 ] - extent[2] + 1 ;
  size[2] = extent[ 5 ] - extent[4] + 1 ;
  
  ImagePointer itkImage = ImageType::New () ;
  
  ImageType::RegionType region;
  region.SetSize( size );
  ImageType::IndexType start;
  start[0] = start[1] = start[2] = 0;  
  region.SetIndex( start );
  itkImage->SetRegions( region );
  double origin[ 3 ] ;
  vtkImage->GetOrigin( origin ) ;
  itkImage->SetOrigin ( origin ) ;
  double spacing[3] ;
  vtkImage->GetSpacing( spacing ) ;
  itkImage->SetSpacing ( spacing ) ;
  itkImage->Allocate () ;
  
  ImageType::IndexType index ;
  int pixel ;

  for ( index[0] = 0 ; index[0] < size[0] ; index[0]++ )
  {
    for ( index[1] = 0 ; index[1] < size[1] ; index[1]++ )
    {
      for ( index[2] = 0 ; index[2] < size[2] ; index[2]++ )
      {
        pixel = vtkImage->GetScalarComponentAsFloat ( index[0], index[1], index[2], 0 ) ;
        if ( pixel != 128 )
        {
          pixel = 0 ;
        }
        itkImage->SetPixel ( index, pixel ) ;
      }
    }
  }
  return itkImage ;
}

void ComputeBoundingBox( vtkSmartPointer<vtkPolyData> mesh , double spacing[ 3 ] , int size[ 3 ] , double origin [ 3 ] )
{
   // we'll need to compute the bounding box for the entire population  
   double largestBoundaries[6] ;
   largestBoundaries[0] = largestBoundaries[2] = largestBoundaries[4] = 1000 ;
   largestBoundaries[1] = largestBoundaries[3] = largestBoundaries[5] = -1000 ;
   double bb[6] ;
   mesh->GetBounds ( bb ) ;
   std::cout << "BB: " << bb[0] << " " << bb[1] << " " << bb[2] << " " << bb[3] << " " << bb[4] << " " << bb[5] << std::endl ;
   if ( bb[0] < largestBoundaries[0] )
   {
      largestBoundaries[0] = bb[0] ;
   }
   if ( bb[2] < largestBoundaries[2] )
   {
      largestBoundaries[2] = bb[2] ;
   }
   if ( bb[4] < largestBoundaries[4] )
   {
      largestBoundaries[4] = bb[4] ;
   }
   if ( bb[1] > largestBoundaries[1] )
   {
      largestBoundaries[1] = bb[1] ;
   }
   if ( bb[3] > largestBoundaries[3] )
   {
      largestBoundaries[3] = bb[3] ;
   }
   if ( bb[5] > largestBoundaries[5] )
   {
      largestBoundaries[5] = bb[5] ;
   }

  std::cout << "Largest boundaries: " << largestBoundaries[0] << " " << largestBoundaries[1] << " " << largestBoundaries[2] << " " << largestBoundaries[3] << " " << largestBoundaries[4] << " " << largestBoundaries[5] << std::endl ;

   // compute image size, origin and spacing given the common bounding box
   for ( int i = 0 ; i < 3 ; i++ )
     {
       origin[i] = floor ( largestBoundaries[2*i] ) - 1  ;
       size[i] = ceil ( ( ceil ( largestBoundaries[2*i+1] ) - origin[i] ) / spacing[i] ) + 1 ;
     }

/*if (extraBorder != 0 )
   {
  size[0] += extraBorder ;
        size[1] += extraBorder ;
        size[2] += extraBorder ;
  origin[0] -= extraBorder/2 ;
        origin[1] -= extraBorder/2 ;
        origin[2] -= extraBorder/2 ;
  std::cout << "Extra border added! " << std::endl ; 
   }*/

   std::cout << "Origin: " << origin[0] << " " << origin[1] << " " << origin[2] << std::endl ;
   std::cout << "Size: " << size[0] << " " << size[1] << " " << size[2] << std::endl ;
   std::cout << "Spacing: " << spacing[0] << " " << spacing[1] << " " << spacing[2] << std::endl ;
}


//class ErrorObserver copied from http://www.vtk.org/Wiki/VTK/Examples/Cxx/Utilities/ObserveError
class ErrorObserver : public vtkCommand
{
  public:
    ErrorObserver():
    Error(false),
    Warning(false),
    ErrorMessage(""),
    WarningMessage("") {}
    static ErrorObserver *New()
    {
      return new ErrorObserver;
    }
    bool GetError() const
    {
      return this->Error;
    }
    bool GetWarning() const
    {
      return this->Warning;
    }
    void Clear()
    {
      this->Error = false;
      this->Warning = false;
      this->ErrorMessage = "";
      this->WarningMessage = "";
    }
    virtual void Execute(vtkObject *vtkNotUsed(caller),
    unsigned long event,
    void *calldata)
    {
      switch(event)
      {
        case vtkCommand::ErrorEvent:
          ErrorMessage = static_cast<char *>(calldata);
          this->Error = true;
          break;
        case vtkCommand::WarningEvent:
          WarningMessage = static_cast<char *>(calldata);
          this->Warning = true;
          break;
      }
    }

    std::string GetErrorMessage()
    {
      return ErrorMessage;
    }
    std::string GetWarningMessage()
    {
      return WarningMessage;
    }
  private:
    bool Error;
    bool Warning;
    std::string ErrorMessage;
    std::string WarningMessage;
};


int ReadVTK( std::string input , vtkSmartPointer<vtkPolyData> &polyData )
{
  vtkSmartPointer<ErrorObserver> errorObserver =
  vtkSmartPointer<ErrorObserver>::New();
  if( input.rfind( ".vtk" ) != std::string::npos )
  {
    vtkSmartPointer< vtkPolyDataReader > polyReader = vtkPolyDataReader::New() ;
    polyReader->AddObserver( vtkCommand::ErrorEvent , errorObserver ) ;
    polyReader->SetFileName( input.c_str() ) ;
    polyData = polyReader->GetOutput() ;
    polyReader->Update() ;
  }
  else if( input.rfind( ".vtp" ) != std::string::npos )
  {
    vtkSmartPointer< vtkXMLPolyDataReader > xmlReader = vtkXMLPolyDataReader::New() ;
    xmlReader->SetFileName( input.c_str() ) ;
    xmlReader->AddObserver( vtkCommand::ErrorEvent , errorObserver ) ;
    polyData = xmlReader->GetOutput() ;
    xmlReader->Update() ;
  }
  else
  {
    std::cerr << "Input file format not handled: " << input << " cannot be read" << std::endl ;
    return 1 ;
  }
  if (errorObserver->GetError())
  {
    std::cout << "Caught error opening " << input << std::endl ;
    return 1 ;
  }
  return 0 ;
}

int main ( int argc, char *argv[] )
{
  typedef itk::Image < unsigned char, 3 > ImageType ;
  typedef ImageType::Pointer ImagePointer ;
  //Parses command line arguments
  PARSE_ARGS ;
  if( mesh.empty() || labelMap.empty() )
  {
    std::cerr << "Please specify an input mesh and and output label map" << std::endl ;
    return EXIT_FAILURE ;
  }
  if( spacingVec.size() != 3 || smoothingRadius.size() != 3 )
  {
    std::cerr << "Spacing and smoothingRadius must have 3 values" << std::endl ;
    return EXIT_FAILURE ;
  }
  //Loads mesh
  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New() ;
  if( ReadVTK( mesh , polyData ) )
  {
    return EXIT_FAILURE ;
  }
  double spacing[ 3 ] ;
  double origin[ 3] ;
  int size[ 3 ] ;
  if( reference.empty() )
  {
    //If no reference image given, computes a bounding box around the mesh and uses the spacing given by the user
    for( int i = 0 ; i < 3 ; i++ )
    {
      spacing[ i ] = spacingVec[ i ] ;
    }
    ComputeBoundingBox( polyData , spacing , size , origin ) ;
  }
  else
  {
    // Reads reference volume and get spacing, dimension and size information for the output label map from it
    typedef itk::ImageFileReader < ImageType > ImageReaderType ;
    ImageReaderType::Pointer referenceVolume = ImageReaderType::New() ;
    referenceVolume->SetFileName( reference ) ;
    referenceVolume->UpdateOutputInformation() ;
    for( unsigned int i = 0 ; i < 3 ; i++ )
    {
      size[ i ] = referenceVolume->GetImageIO()->GetDimensions( i ) ;
      origin[ i ] = referenceVolume->GetImageIO()->GetOrigin( i ) ;
      spacing[ i ] = referenceVolume->GetImageIO()->GetSpacing( i ) ;
    }
  }


     
   // Scan-convert the mesh
   vtkSmartPointer<vtkAttributedPolyDataToImage> scanConverter = vtkSmartPointer<vtkAttributedPolyDataToImage>::New () ;
   scanConverter->SetTolerance ( 0.0 ) ;
   scanConverter->SetInput ( polyData ) ;
   scanConverter->SetOutputOrigin ( origin ) ;
   scanConverter->SetOutputSpacing ( spacing ) ;
   scanConverter->SetOutputWholeExtent ( 0, size[0] - 1, 0, size[1] - 1, 0, size[2] - 1 ) ;
   scanConverter->Update () ;
   vtkImageData *vtkBinaryVolume = scanConverter->GetBinaryVolume () ;
   ImagePointer binaryVolume ;
   binaryVolume = VTK2BinaryITK ( vtkBinaryVolume ) ;
   if( smoothing )
   {
     if( verbose )
     {
       std::cout << "Median Filtering" << std::endl ;
     }
     typedef itk::MedianImageFilter<ImageType,ImageType> MedianFilterType ;
     MedianFilterType::InputSizeType radius;
     radius.Fill(2);
     for( int i = 0 ; i < 3 ; i++ )
     {
       radius.Fill( smoothingRadius[i] ) ;
     }
     MedianFilterType::Pointer medianFilter = MedianFilterType::New() ;
     medianFilter->SetInput( binaryVolume ) ;
     medianFilter->SetRadius( radius ) ;
     medianFilter->Update() ;
     binaryVolume = medianFilter->GetOutput() ;
   }
   WriteITKImage( binaryVolume , labelMap ) ;  
   return EXIT_SUCCESS ;
} 

