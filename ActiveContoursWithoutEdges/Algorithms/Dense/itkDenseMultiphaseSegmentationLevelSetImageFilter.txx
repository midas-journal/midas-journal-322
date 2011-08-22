#ifndef __itkDenseMultiphaseSegmentationLevelSetImageFilter_txx_
#define __itkDenseMultiphaseSegmentationLevelSetImageFilter_txx_

#include "itkDenseMultiphaseSegmentationLevelSetImageFilter.h"

namespace itk {

template < class TInputImage, class TFeatureImage, class TOutputPixelType >
DenseMultiphaseSegmentationLevelSetImageFilter< TInputImage, TFeatureImage,
  TOutputPixelType >
::DenseMultiphaseSegmentationLevelSetImageFilter()
{
  m_AutoGenerateSpeedAdvection = true;

  // Provide some reasonable defaults to prevent infinite looping.
  this->SetMaximumRMSError( 0.02 );
  this->SetNumberOfIterations( 1000 );
  m_ReverseExpansionDirection = false;

  this->SetNumberOfOutputs( 1 );
  AddOutput( OutputImageType::New() );
  this->SetNumberOfRequiredInputs( 1 );
}

template < class TInputImage, class TFeatureImage, class TOutputPixelType >
void
DenseMultiphaseSegmentationLevelSetImageFilter< TInputImage, TFeatureImage,
TOutputPixelType >
::GenerateSpeedImage()
{
	for( unsigned int i = 0; i < this->FunctionCount; i++ )
		m_SegmentationFunctions[i]->AllocateSpeedImage();
	for( unsigned int i = 0; i < this->FunctionCount; i++ )
		m_SegmentationFunctions[i]->CalculateSpeedImage();
}

template <class TInputImage, class TFeatureImage, class TOutputPixelType >
void
DenseMultiphaseSegmentationLevelSetImageFilter<TInputImage, TFeatureImage,
TOutputPixelType >
::GenerateAdvectionImage()
{
  for( unsigned int i = 0; i < this->FunctionCount; i++ )
    m_SegmentationFunctions[i]->AllocateAdvectionImage();
  for( unsigned int i = 0; i < this->FunctionCount; i++ )
    m_SegmentationFunctions[i]->CalculateAdvectionImage();
}

template < class TInputImage, class TFeatureImage, class TOutputPixelType >
void
DenseMultiphaseSegmentationLevelSetImageFilter< TInputImage, TFeatureImage,
  TOutputPixelType >
::GenerateData()
{
  for( unsigned int i = 0; i < this->FunctionCount; i++ )
  {
    if (m_SegmentationFunctions[i] == 0)
      itkExceptionMacro("No finite difference function was specified.");

    // A positive speed value causes surface expansion, the opposite of the
    // default.  Flip the sign of the propagation and advection weights.
    if (m_ReverseExpansionDirection == true)
      m_SegmentationFunctions[i]->ReverseExpansionDirection();
  }

  // Allocate the images from which speeds will be sampled.
  if ( this->GetState() == Superclass::UNINITIALIZED &&
    m_AutoGenerateSpeedAdvection == true )
  {
    this->GenerateSpeedImage();
    this->GenerateAdvectionImage();
  }

  // Start the solver
  Superclass::GenerateData();

  // Reset all the signs of the weights.
  if (m_ReverseExpansionDirection == true)
    for( unsigned int i = 0; i < this->FunctionCount; i++ )
      m_SegmentationFunctions[i]->ReverseExpansionDirection();
}

template < class TInputImage, class TFeatureImage, class TOutputPixelType >
void
DenseMultiphaseSegmentationLevelSetImageFilter< TInputImage, TFeatureImage,
  TOutputPixelType >
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Name Of Class: " << GetNameOfClass( ) << std::endl;
  os << indent << "m_ReverseExpansionDirection = " <<
    GetReverseExpansionDirection() << std::endl;
  os << indent << "m_AutoGenerateSpeedAdvection = " <<
    GetAutoGenerateSpeedAdvection() << std::endl;
}

} // end namespace itk

#endif
