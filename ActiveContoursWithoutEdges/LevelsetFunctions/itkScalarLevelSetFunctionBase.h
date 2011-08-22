#ifndef __itkScalarLevelSetFunctionBase_h_
#define __itkScalarLevelSetFunctionBase_h_

#include "itkRegionBasedLevelSetFunction.h"
#include "itkScalarChanAndVeseSharedFunctionData.h"
//#include "itkStatisticsLabelObject.h"

namespace itk {

template < class TInputImage,
class TFeatureImage,
class TSharedData = ScalarChanAndVeseSharedFunctionData< TInputImage, TFeatureImage > >
class ITK_EXPORT ScalarLevelSetFunctionBase
: public RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >
{
public:
  typedef ScalarLevelSetFunctionBase Self;
  typedef RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData > Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  // itkNewMacro(Self); // Abstract class and so cannot be instantiated

  /** Run-time type information (and related methods) */
  itkTypeMacro( ScalarLevelSetFunctionBase,
    RegionBasedLevelSetFunction );

  itkStaticConstMacro( ImageDimension, unsigned int,
    TFeatureImage::ImageDimension );

  typedef TInputImage InputImageType;
  typedef typename Superclass::InputImageConstPointer InputImageConstPointer;
  typedef typename Superclass::InputImagePointer InputImagePointer;
  typedef typename Superclass::InputPixelType InputPixelType;
  typedef typename Superclass::InputIndexType InputIndexType;
  typedef typename Superclass::InputIndexValueType InputIndexValueType;
  typedef typename Superclass::InputSizeType InputSizeType;
  typedef typename Superclass::InputSizeValueType InputSizeValueType;
  typedef typename Superclass::InputRegionType InputRegionType;
  typedef typename Superclass::InputPointType InputPointType;

  typedef TFeatureImage FeatureImageType;
  typedef typename FeatureImageType::ConstPointer FeatureImageConstPointer;
  typedef typename Superclass::FeaturePixelType FeaturePixelType;
  typedef typename Superclass::FeatureIndexType FeatureIndexType;
  typedef typename Superclass::FeatureOffsetType FeatureOffsetType;

  typedef typename Superclass::ScalarValueType ScalarValueType;
  typedef typename Superclass::NeighborhoodType NeighborhoodType;
  typedef typename Superclass::FloatOffsetType FloatOffsetType;
  typedef typename Superclass::RadiusType RadiusType;
  typedef typename Superclass::TimeStepType TimeStepType;
  typedef typename Superclass::GlobalDataStruct GlobalDataStruct;
  typedef typename Superclass::PixelType PixelType;
  typedef typename Superclass::VectorType VectorType;

  typedef typename Superclass::SharedDataType SharedDataType;
  typedef typename Superclass::SharedDataPointer SharedDataPointer;

  typedef ImageRegionIteratorWithIndex< InputImageType >
    ImageIteratorType;
  typedef ImageRegionConstIteratorWithIndex< InputImageType >
    ConstImageIteratorType;
  typedef ImageRegionConstIterator< FeatureImageType >
    ConstFeatureIteratorType;

  typedef std::list< unsigned int > ListPixelType;
  typedef typename ListPixelType::const_iterator ListPixelConstIterator;
  typedef typename ListPixelType::iterator ListPixelIterator;
  typedef Image< ListPixelType, ImageDimension > ListImageType;

  void UpdatePixel( const unsigned int& idx,
    NeighborhoodIterator<TInputImage> &iterator,
    InputPixelType &newValue,
    bool &status );

protected:
  ScalarLevelSetFunctionBase() : Superclass(){}
  ~ScalarLevelSetFunctionBase(){}

  void computeOverlapParameters( const FeatureIndexType featIndex,
    unsigned int& s, unsigned int& pr );

  void ComputeParameters();

private:
};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkScalarLevelSetFunctionBase.txx"
#endif

#endif

