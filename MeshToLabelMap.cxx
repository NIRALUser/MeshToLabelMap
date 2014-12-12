#include "MeshToLabelMapCLP.h"
#include <iostream>
#include <vector>
#include <string>

#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkMedianImageFilter.h>

#include <vtkSmartPointer.h>
#include <vtkTransform.h>
#include <vtkMatrix4x4.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkImageData.h>
#include <vtkCommand.h>
#include <vtkTransformPolyDataFilter.h>
#include "vtkAttributedPolyDataToImage.h"


int WriteITKImage ( itk::Image < unsigned char, 3 >::Pointer image, std::string fileName )
{
    try
    {
        typedef itk::Image < unsigned char, 3 > ImageType ;
        typedef itk::ImageFileWriter < ImageType > ImageWriterType ;
        ImageWriterType::Pointer itkWriter = ImageWriterType::New() ;
        itkWriter->SetFileName ( fileName ) ;
        itkWriter->SetInput ( image ) ;
        itkWriter->UseCompressionOn() ;
        itkWriter->Write() ;
    }
    catch( itk::ExceptionObject & err )
    {
        std::cerr << "ExceptionObject caught !" << std::endl ;
        std::cerr << err << std::endl ;
        return 1 ;
    }
    return 0 ;
}

itk::Image < unsigned char, 3 >::Pointer VTK2BinaryITK ( vtkImageData *vtkImage , itk::Matrix< double , 3 , 3 > direction , unsigned char value )
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
  
  ImageType::RegionType region ;
  region.SetSize( size ) ;
  ImageType::IndexType start ;
  start[0] = start[1] = start[2] = 0 ;
  region.SetIndex( start ) ;
  itkImage->SetRegions( region ) ;
  double origin[ 3 ] ;
  vtkImage->GetOrigin( origin ) ;
  itkImage->SetOrigin( origin ) ;
  double spacing[ 3 ] ;
  vtkImage->GetSpacing( spacing ) ;
  itkImage->SetSpacing( spacing ) ;
  itkImage->SetDirection( direction ) ;
  itkImage->Allocate() ;
  ImageType::IndexType index ;
  ImageType::IndexType maxIndex ;
//The copy from size to maxIndex is just remove the warning due to the comparison of signed and unsigned
//variables in the scope below
  for( int i = 0 ; i < 3 ; i++ )
  {
    maxIndex[ i ] = size[ i ] ;
  }
  int pixel ;
  
  for ( index[ 0 ] = 0 ; index[ 0 ] < maxIndex[ 0 ] ; index[ 0 ]++ )
  {
    for ( index[ 1 ] = 0 ; index[ 1 ] < maxIndex[ 1 ] ; index[ 1 ]++ )
    {
      for ( index[ 2 ] = 0 ; index[ 2 ] < maxIndex[ 2 ] ; index[ 2 ]++ )
      {
        pixel = vtkImage->GetScalarComponentAsFloat( index[ 0 ], index[ 1 ] , index[ 2 ] , 0 ) ;
        if ( pixel != 128 )
        {
          pixel = 0 ;
        }
        else
        {
            pixel = value ;
        }
        itkImage->SetPixel ( index , pixel ) ;
      }
    }
  }
  return itkImage ;
}

void ComputeBoundingBoxFromReferenceImage( std::string reference ,
                         double spacing[ 3 ] ,
                         int size[ 3 ] ,
                         double origin [ 3 ] ,
                         itk::Matrix< double , 3 , 3 > &direction
                       )
{
    // Reads reference volume and get spacing, dimension and size information for the output label map from it
    typedef itk::Image < unsigned char, 3 > ImageType ;
    typedef itk::ImageFileReader < ImageType > ImageReaderType ;
    ImageReaderType::Pointer referenceVolume = ImageReaderType::New() ;
    referenceVolume->SetFileName( reference ) ;
    referenceVolume->UpdateOutputInformation() ;
    std::vector< std::vector< double > > directionVec ;
    for( unsigned int i = 0 ; i < 3 ; i++ )
    {
      size[ i ] = referenceVolume->GetImageIO()->GetDimensions( i ) ;
      origin[ i ] = referenceVolume->GetImageIO()->GetOrigin( i ) ;
      spacing[ i ] = referenceVolume->GetImageIO()->GetSpacing( i ) ;
      directionVec.push_back( referenceVolume->GetImageIO()->GetDirection( i ) ) ;
      for( int j = 0 ; j < 3 ; j++ )
      {
          direction[ j ][ i ] = directionVec[ i ][ j ] ;
      }
    }
}

void ComputeBoundingBoxFromPolyData( vtkSmartPointer<vtkPolyData> mesh ,
                         std::vector< double > boundaryExtension ,
                         double spacing[ 3 ] ,
                         int size[ 3 ] ,
                         double origin [ 3 ] ,
                         bool verbose
                       )
{
    // we'll need to compute the bounding box for the mesh
    double largestBoundaries[ 6 ] ;
    int largestPossibleImage = 2000 ;
    mesh->GetBounds ( largestBoundaries ) ;
    if( verbose )
    {
        std::cout << "Input Mesh Bounding Box: " << largestBoundaries[0] << " "
                                                 << largestBoundaries[1] << " "
                                                 << largestBoundaries[2] << " "
                                                 << largestBoundaries[3] << " "
                                                 << largestBoundaries[4] << " "
                                                 << largestBoundaries[5] << std::endl ;
    }
    // compute image size and origin given the bounding box, the given spacing and the boundary extension
    for ( int i = 0 ; i < 3 ; i++ )
    {
        origin[i] = largestBoundaries[ 2 * i ] - boundaryExtension[ i ] * spacing[ i ] ;
        size[i] = ceil( ( largestBoundaries[ 2 * i + 1 ] - origin[ i ] ) / spacing[ i ] ) + boundaryExtension[ i ] ;
        if( size[ i ] > largestPossibleImage )
        {
            std::cerr << "Dimension " << i << " : Very large size! (" << size[ i ] <<
                         "). The tool might be very slow or crash due to lack of memory." << std::endl ;
        }
    }
}

//class ErrorObserver copied from http://www.vtk.org/Wiki/VTK/Examples/Cxx/Utilities/ObserveError
class ErrorObserver : public vtkCommand
{
  public:
    ErrorObserver():
    Error( false ) ,
    Warning( false ) ,
    ErrorMessage( "" ) ,
    WarningMessage( "" ) {}
    static ErrorObserver *New()
    {
      return new ErrorObserver ;
    }
    bool GetError() const
    {
      return this->Error ;
    }
    bool GetWarning() const
    {
      return this->Warning ;
    }
    void Clear()
    {
      this->Error = false ;
      this->Warning = false ;
      this->ErrorMessage = "" ;
      this->WarningMessage = "" ;
    }
    virtual void Execute( vtkObject *vtkNotUsed( caller ) ,
                         unsigned long event ,
                         void *calldata
                        )
    {
      switch( event )
      {
        case vtkCommand::ErrorEvent:
          ErrorMessage = static_cast<char *>( calldata ) ;
          this->Error = true ;
          break ;
        case vtkCommand::WarningEvent:
          WarningMessage = static_cast<char *>( calldata ) ;
          this->Warning = true ;
          break ;
      }
    }
    std::string GetErrorMessage()
    {
      return ErrorMessage ;
    }
    std::string GetWarningMessage()
    {
      return WarningMessage ;
    }
  private:
    bool Error ;
    bool Warning ;
    std::string ErrorMessage ;
    std::string WarningMessage ;
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
  if( errorObserver->GetError() )
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
  if( spacingVec.size() != 3 || smoothingRadius.size() != 3 || boundaryExtension.size() != 3 )
  {
    std::cerr << "Spacing, boundaryExtension and smoothingRadius must have 3 values" << std::endl ;
    return EXIT_FAILURE ;
  }
  if( value < 1 || value > 255 )
  {
      std::cerr << "Label map pixel value has to belong to [1;255]" << std::endl ;
      return EXIT_FAILURE ;
  }
  //Check spacing information
  if( ( spacingVec[ 0 ] != -1 || spacingVec[ 1 ] != -1 || spacingVec[ 2 ] != -1 )  && !reference.empty() )
  {
      std::cerr << "A reference image is given." << std::endl ;
      std::cout << "Reference image spacing is used to set output label map spacing." << std::endl ;
      std::cerr << "Do not specify additional spacing information." << std::endl ;
      return EXIT_FAILURE ;
  }
  if( reference.empty() && ( spacingVec[ 0 ] <= 0 || spacingVec[ 1 ] <= 0 || spacingVec[ 2 ] <= 0 ) )
  {
      std::cerr << "Please provide valid spacing information ( > 0 )" << std::endl ;
      return EXIT_FAILURE ;
  }
  if( smoothing && ( smoothingRadius[ 0 ] < 1 || smoothingRadius[ 1 ] < 1 || smoothingRadius[ 2 ] < 1 ) )
  {
      std::cerr << "Please provide valid smoothing radius information ( >= 1 )" << std::endl ;
      return EXIT_FAILURE ;
  }
  //Loads mesh
  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New() ;
  if( ReadVTK( mesh , polyData ) )
  {
    return EXIT_FAILURE ;
  }
  vtkSmartPointer<vtkMatrix4x4> RASMatrix = vtkSmartPointer<vtkMatrix4x4>::New() ;
  RASMatrix->Identity() ;
  RASMatrix->SetElement( 0 , 0 , -1 ) ;
  RASMatrix->SetElement( 1 , 1 , -1 ) ;
  vtkSmartPointer<vtkTransform> transform =
             vtkSmartPointer<vtkTransform>::New() ;
  transform->SetMatrix( RASMatrix ) ;
  vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter =
             vtkSmartPointer<vtkTransformPolyDataFilter>::New() ;
  transformFilter->SetInputData( polyData ) ;
  transformFilter->SetTransform( transform ) ;
  transformFilter->Update() ;
  polyData->ShallowCopy( transformFilter->GetOutput() ) ;
  double spacing[ 3 ] ;
  double origin[ 3] ;
  int size[ 3 ] ;
  itk::Matrix<double,3,3> direction ;
  direction.SetIdentity() ;
  if( reference.empty() )
  {
    //If no reference image given, computes a bounding box around the mesh and uses the spacing given by the user
    for( int i = 0 ; i < 3 ; i++ )
    {
      spacing[ i ] = spacingVec[ i ] ;
    }
    ComputeBoundingBoxFromPolyData( polyData , boundaryExtension , spacing , size , origin , verbose ) ;
  }
  else
  {
    ComputeBoundingBoxFromReferenceImage( reference , spacing , size , origin , direction ) ;
    itk::Matrix<double,3,3> Mt ;
    Mt = direction.GetTranspose() ;
    if( verbose )
    {
        std::cout << "Direction: " << std::endl ;
        std::cout << direction << std::endl ;
        std::cout << "Transpose direction matrix: " << std::endl ;
        std::cout << Mt << std::endl ;
    }
    itk::Vector<double,3> originPoint(origin)  ;
    itk::Vector<double,3> offset ;
    offset = originPoint - Mt * originPoint ;
    vtkSmartPointer<vtkMatrix4x4> OrientationMatrix = vtkSmartPointer<vtkMatrix4x4>::New() ;
    OrientationMatrix->Identity() ;
    for( int i = 0 ; i < 3 ; i++ )
    {
      for( int j = 0 ; j < 3 ; j++ )
      {
          OrientationMatrix->SetElement( i , j , Mt[ i ][ j ] ) ;
      }
      OrientationMatrix->SetElement( i , 3 , offset[ i ] ) ;
    }
    if( verbose )
    {
      OrientationMatrix->Print( std::cout ) ;
    }
    transform->SetMatrix( OrientationMatrix ) ;
    transformFilter->SetInputData( polyData ) ;
    transformFilter->SetTransform( transform ) ;
    transformFilter->Update() ;
    polyData = transformFilter->GetOutput() ;
  }
  if( verbose )
  {
    std::cout << "Origin: " << origin[0] << " " << origin[1] << " " << origin[2] << std::endl ;
    std::cout << "Size: " << size[0] << " " << size[1] << " " << size[2] << std::endl ;
    std::cout << "Spacing: " << spacing[0] << " " << spacing[1] << " " << spacing[2] << std::endl ;
    std::cout << "Direction: " << std::endl ;
    std::cout << direction << std::endl ;
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
  //we can safely can 'value' which is an int to unsigned char because we checked that its value is between 1 and 255 (1 because 0 is the background value)
  binaryVolume = VTK2BinaryITK ( vtkBinaryVolume , direction , (unsigned char)value ) ;
  if( smoothing )
  {
    if( verbose )
    {
      std::cout << "Median Filtering" << std::endl ;
      std::cout << "Radius: " << smoothingRadius[ 0 ] << " " << smoothingRadius[ 1 ] << " " << smoothingRadius[ 2 ] << std::endl ;
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
  if( WriteITKImage( binaryVolume , labelMap ) )
  {
    return EXIT_FAILURE ;
  }
  return EXIT_SUCCESS ;
}

