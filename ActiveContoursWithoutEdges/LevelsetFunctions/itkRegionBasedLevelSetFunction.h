#ifndef __itkRegionBasedLevelSetFunction_h_
#define __itkRegionBasedLevelSetFunction_h_

#include "itkSegmentationLevelSetFunction.h"

#include "itkScalarChanAndVeseSharedFunctionData.h"
#include "itkSparseMultiphaseLevelSetImageFilter.h"
#include "itkDenseMultiphaseLevelSetImageFilter.h"

namespace itk {

template < class TInputImage, // LevelSetImageType
  class TFeatureImage, // FeatureImageType
  class TSharedData >
class ITK_EXPORT RegionBasedLevelSetFunction: public
SegmentationLevelSetFunction< TInputImage, TFeatureImage >
{
public:
  /** Standard class typedefs. */
  typedef RegionBasedLevelSetFunction Self;
  typedef SegmentationLevelSetFunction< TInputImage, TFeatureImage >
    Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /** Run-time type information (and related methods) */
  itkTypeMacro( RegionBasedLevelSetFunction, SegmentationLevelSetFunction );

  typedef TInputImage InputImageType;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;
  typedef typename InputImageType::Pointer InputImagePointer;
  typedef typename InputImageType::PixelType InputPixelType;
  typedef typename InputImageType::IndexType InputIndexType;
  typedef typename InputImageType::IndexValueType InputIndexValueType;
  typedef typename InputImageType::SizeType InputSizeType;
  typedef typename InputImageType::SizeValueType InputSizeValueType;
  typedef typename InputImageType::RegionType InputRegionType;
  typedef typename InputImageType::PointType InputPointType;

  typedef TFeatureImage FeatureImageType;
  typedef typename FeatureImageType::PixelType FeaturePixelType;
  typedef typename FeatureImageType::IndexType FeatureIndexType;
  typedef typename FeatureImageType::OffsetType FeatureOffsetType;

  /** Extract some parameters from the superclass. */
  typedef typename Superclass::ScalarValueType ScalarValueType;
  typedef typename Superclass::NeighborhoodType NeighborhoodType;
  typedef typename Superclass::FloatOffsetType FloatOffsetType;
  typedef typename Superclass::RadiusType RadiusType;
  typedef typename Superclass::TimeStepType TimeStepType;
  typedef typename Superclass::GlobalDataStruct GlobalDataStruct;
  typedef typename Superclass::PixelType PixelType;
  typedef typename Superclass::VectorType VectorType;
  typedef typename Superclass::ImageType ImageType;

  typedef TSharedData SharedDataType;
  typedef typename SharedDataType::Pointer SharedDataPointer;

  typedef ImageRegionIteratorWithIndex< InputImageType >
    ImageIteratorType;
  typedef ImageRegionConstIteratorWithIndex< InputImageType >
    ConstImageIteratorType;
  typedef ImageRegionConstIterator< FeatureImageType >
    ConstFeatureIteratorType;

  virtual void Initialize(const RadiusType &r)
  {
    Superclass::Initialize(r);

    ScalarValueType null_value = NumericTraits<ScalarValueType>::Zero;
    this->SetAdvectionWeight( null_value );
    this->SetPropagationWeight( null_value );
    this->SetCurvatureWeight( null_value );
  }

  /* This structure is derived from LevelSetFunction and stores intermediate
  values for computing time step sizes */
  struct ACGlobalDataStruct : public GlobalDataStruct
  {
    ScalarValueType m_MaxGlobalChange;
  };

  void SetSharedData( SharedDataPointer sharedDataIn )
  {
    sharedData = sharedDataIn;
  }

  void UpdateSharedData( bool forceUpdate );

  void *GetGlobalDataPointer() const
  {
    ACGlobalDataStruct *ans = new ACGlobalDataStruct();

    ScalarValueType null_value = NumericTraits<ScalarValueType>::Zero;

    ans->m_MaxAdvectionChange   = null_value;
    ans->m_MaxPropagationChange = null_value;
    ans->m_MaxCurvatureChange   = null_value;
	  ans->m_MaxGlobalChange      = null_value;
    return ans;
  }

  TimeStepType ComputeGlobalTimeStep(void *GlobalData) const;

  /** Compute the equation value. */
  PixelType ComputeUpdate(const NeighborhoodType &neighborhood,
    void *globalData, const FloatOffsetType& = FloatOffsetType(0.0));

  void SetInitialImage(InputImageType *f)
  {
    m_InitialImage = f;
  }

  void SetMu( const double& mu);
  double GetMu() const;

  void SetNu( const double& nu);

  double GetNu() const;

  void SetLambda1( const double& lambda1 );
  double GetLambda1() const;

  void SetLambda2( const double& lambda2 );
  double GetLambda2() const;

  void SetGamma(const double& gamma);
  double GetGamma() const;

  void SetEta(const double& eta);
  double GetEta() const;

  void SetEpsilon(const double& e);
  double GetEpsilon() const;

  void SetTau(const double& t);
  double GetTau() const;

  void SetVolume( const double& v);
  double GetVolume() const;

  ScalarValueType GetAdvectionWeight() const
  {
    return itk::NumericTraits< ScalarValueType >::Zero;
  }

  ScalarValueType GetPropagationWeight() const
  {
    return itk::NumericTraits< ScalarValueType >::Zero;
  }

  void SetFunctionId( const unsigned int& functionId );

protected:

  RegionBasedLevelSetFunction();
  virtual ~RegionBasedLevelSetFunction() {}

  /** The image whose features will be used to create a speed image */
  InputImageConstPointer m_InitialImage;

  bool updatedC;
  bool updatedH;

  SharedDataPointer sharedData;

  /* Area regularizer term in CV formulation, what about lambda1 and lambda2?*/
  double m_Nu;
  double m_Lambda1;
  double m_Lambda2;
  double m_Gamma;
	double m_Tau;
	double m_Volume;
  double m_Epsilon; // to be used in Heaviside function
  double m_inv_epsilon;
  double m_2_over_pi;
  unsigned int m_FunctionId;

  double calculateH( const InputPixelType& z);
  double calculatedH( const InputPixelType& z);
  void ComputeHImage();
  double computeRegularizationTerms( /* fill with adequate parameters */);

  double computeGlobalTerm(
    const InputPixelType& imagePixel,
    const InputIndexType& inputIndex );

  virtual double computeInternalTerm(const FeaturePixelType& iValue,
    const FeatureIndexType& iIdx, const unsigned int& fId ) = 0;

  virtual double computeExternalTerm(const FeaturePixelType& iValue,
    const FeatureIndexType& iIdx, const unsigned int& pr ) = 0;

  virtual void computeOverlapParameters( const FeatureIndexType featIndex,
    unsigned int& s, unsigned int& pr ) = 0;

  virtual double computeOverlapTerm( const unsigned int& s );

  virtual void ComputeParameters() = 0;

  virtual void SpecialProcessing(){}

private:
  RegionBasedLevelSetFunction(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRegionBasedLevelSetFunction.txx"
#endif

#endif
