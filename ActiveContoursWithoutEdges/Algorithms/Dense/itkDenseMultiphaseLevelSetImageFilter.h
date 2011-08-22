#ifndef __itkDenseMultiphaseLevelSetImageFilter_h_
#define __itkDenseMultiphaseLevelSetImageFilter_h_

#include "itkDenseMultiphaseSegmentationLevelSetImageFilter.h"
#include "itkScalarChanAndVeseSharedFunctionData.h"

#ifdef VIZU
  #include "vtkVisualize2DImplicitFunction.h"
  #include "vtkVisualize3DImplicitFunction.h"
#endif

namespace itk
{


template < class TInputImage, class TFeatureImage, class TFunction,
  typename TOutputPixel = typename TInputImage::PointType::CoordRepType >
class ITK_EXPORT DenseMultiphaseLevelSetImageFilter:
public DenseMultiphaseSegmentationLevelSetImageFilter< TInputImage,
  TFeatureImage, TOutputPixel >
{

public:
  typedef DenseMultiphaseLevelSetImageFilter Self;
  typedef DenseMultiphaseSegmentationLevelSetImageFilter< TInputImage,
    TFeatureImage, TOutputPixel > Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( DenseMultiphaseLevelSetImageFilter,
    DenseMultiphaseSegmentationLevelSetImageFilter );

  void PrintSelf( std::ostream& os, Indent indent ) const; 

  /** Repeat definition from Superclass to satisfy Borland compiler quirks */
  itkStaticConstMacro( ImageDimension, unsigned int,
    TInputImage::ImageDimension );

  /** Inherited typedef from the superclass. */
  typedef TFeatureImage FeatureImageType;
  typedef typename FeatureImageType::PixelType FeaturePixelType;

  /** Output image type typedefs */
  typedef Image< TOutputPixel, ImageDimension > OutputImageType;
  typedef typename OutputImageType::ValueType ValueType;
  typedef typename OutputImageType::IndexType IndexType;

  typedef typename Superclass::TimeStepType   TimeStepType;
  typedef typename Superclass::FiniteDifferenceFunctionType
    FiniteDifferenceFunctionType;

  typedef TInputImage InputImageType;
  typedef TFunction FunctionType;

  typedef ScalarChanAndVeseSharedFunctionData< InputImageType,
    FeatureImageType > SharedDataType;
  typedef typename SharedDataType::Pointer SharedDataPointer;

#ifdef VIZU
  typedef vtkVisualize2DImplicitFunction< InputImageType, FeatureImageType >
    vtkVisualize2DImplicitFunctionType;
#endif

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(OutputHasNumericTraitsCheck,
    (Concept::HasNumericTraits<TOutputPixel>) );
  /** End concept checking */
#endif

  virtual FunctionType *GetTypedDifferenceFunction(
    unsigned int functionIndex )
  {
    return static_cast< FunctionType* >
      ( & ( **this->m_DifferenceFunctions[functionIndex] ) );
  }

  void SetFunctionCount( unsigned int n )
  {
    Superclass::SetFunctionCount( n );

    for( unsigned int i = 0; i < this->FunctionCount; i++ )
    {
      this->m_DifferenceFunctions[i] = new typename
        FiniteDifferenceFunctionType::Pointer();
      (*this->m_DifferenceFunctions[i]) = FunctionType::New();
      SetSegmentationFunction( i, GetTypedDifferenceFunction( i ) );
    }
  }

protected:
  DenseMultiphaseLevelSetImageFilter()
  {
    sharedData = SharedDataType::New();
  }

  ~DenseMultiphaseLevelSetImageFilter()
  {
    for ( unsigned int i = 0; i < this->FunctionCount; i++ )
      delete this->m_DifferenceFunctions[i];
  }

#if VIZU
  vtkVisualize2DImplicitFunctionType func;
#endif

  SharedDataPointer sharedData;

  virtual void Initialize();
  virtual void InitializeIteration();
};

} //end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDenseMultiphaseLevelSetImageFilter.txx"
#endif

#endif
