#ifndef __itkMultiphaseFiniteDifferenceImageFilter_txx_
#define __itkMultiphaseFiniteDifferenceImageFilter_txx_

#include "itkMultiphaseFiniteDifferenceImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkExceptionObject.h"
#include "itkEventObject.h"

namespace itk {

template <class TInputImage, class TOutputImage >
void
MultiphaseFiniteDifferenceImageFilter<TInputImage, TOutputImage >
::GenerateData()
{
  if (this->GetState() == UNINITIALIZED)
  {
    // Set the coefficients for the deriviatives
    double coeffs[ImageDimension];
    if (m_UseImageSpacing)
      for (unsigned int i = 0; i < ImageDimension; i++)
        coeffs[i] = 1.0 / m_LevelSet[0]->GetSpacing()[i];
    else
      for (unsigned int i = 0; i < ImageDimension; i++)
        coeffs[i] = 1.0;

    for( unsigned int i = 0; i < this->FunctionCount; i++ )
      GetDifferenceFunction(i)->SetScaleCoefficients(coeffs);

    // Allocate the output image -- inherited
    this->AllocateOutputs();
//     std::cout << "Allocated output" << std::endl;

    // Copy the input image to the output image.  Algorithms will operate
    // directly on the output image and the update buffer.
    this->CopyInputToOutput();
//     std::cout << "Copy input to output" << std::endl;

    // Perform any other necessary pre-iteration initialization.
    this->Initialize();
//     std::cout << "Initialization complete" << std::endl;

    // Allocate the internal update buffer.  This takes place entirely within
    // the subclass, since this class cannot define an update buffer type.
    this->AllocateUpdateBuffer();
//     std::cout << "Allocated update buffer" << std::endl;

    this->SetStateToInitialized();
    m_ElapsedIterations = 0;
  }

  // Iterative algorithm
  TimeStepType dt;

  while ( ! this->Halt() )
  {
    /* An optional method for precalculating global values, or setting
    up for the next iteration*/
    this->InitializeIteration();
    std::cout << "Initialized for iteration: " << this->m_ElapsedIterations
      << std::endl;

    dt = this->CalculateChange();
//     std::cout << "Calculate change complete" << std::endl;

    this->ApplyUpdate( dt );
    ++m_ElapsedIterations;
//     std::cout << "Applied update: " << this->GetRMSChange() << std::endl;

    // Invoke the iteration event.
    this->InvokeEvent( IterationEvent() );
    if( this->GetAbortGenerateData() )
    {
      this->InvokeEvent( IterationEvent() );
      this->ResetPipeline();
      throw ProcessAborted(__FILE__,__LINE__);
    }
  }

  // Reset the state once execution is completed
  if (m_ManualReinitialization == false)
    this->SetStateToUninitialized();

  // Any further processing of the solution can be done here.
  this->PostProcessOutput();
//   std::cout << "Copy input to output" << std::endl;
}

/**
 *
 */
template <class TInputImage, class TOutputImage >
void
MultiphaseFiniteDifferenceImageFilter< TInputImage,TOutputImage >
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  // copy the output requested region to the input requested region
  Superclass::GenerateInputRequestedRegion();

  // get pointers to the input
  typename Superclass::InputImagePointer  inputPtr  =
    const_cast< TInputImage * >( this->GetInput(0) );

  if ( !inputPtr )
    return;

  // Get the size of the neighborhood on which we are going to operate.  This
  // radius is supplied by the difference function we are using.
  typename FiniteDifferenceFunctionType::RadiusType radius
    = this->GetDifferenceFunction(0)->GetRadius();

  // Try to set up a buffered region that will accommodate our
  // neighborhood operations.  This may not be possible and we
  // need to be careful not to request a region outside the largest
  // possible region, because the pipeline will give us whatever we
  // ask for.

  // get a copy of the input requested region (should equal the output
  // requested region)
  typename TInputImage::RegionType inputRequestedRegion;
  inputRequestedRegion = inputPtr->GetRequestedRegion();

  // pad the input requested region by the operator radius
  inputRequestedRegion.PadByRadius( radius );

  // crop the input requested region at the input's largest possible region
  if ( inputRequestedRegion.Crop( inputPtr->GetLargestPossibleRegion() ) )
  {
    inputPtr->SetRequestedRegion( inputRequestedRegion );
    return;
  }
  else
  {
    // Couldn't crop the region (requested region is outside the largest
    // possible region).  Throw an exception.

    // store what we tried to request (prior to trying to crop)
    inputPtr->SetRequestedRegion( inputRequestedRegion );

    // build an exception
    InvalidRequestedRegionError e(__FILE__, __LINE__);
    e.SetLocation(ITK_LOCATION);
    e.SetDescription("Requested region is (at least partially) outside the "
      "largest possible region.");
    e.SetDataObject(inputPtr);
    throw e;
  }
}

template < class TInputImage, class TOutputImage >
typename MultiphaseFiniteDifferenceImageFilter< TInputImage,
  TOutputImage >::TimeStepType
MultiphaseFiniteDifferenceImageFilter< TInputImage, TOutputImage >
::ResolveTimeStep( const TimeStepType *timeStepList, const bool *valid,
  int size)
{
  TimeStepType min;
  bool flag;
  min = NumericTraits<TimeStepType>::Zero;

  // grab first valid value
  flag = false;
  for (int i = 0; i < size; ++i)
    if (valid[i])
    {
      min = timeStepList[i];
      flag = true;
      break;
    }

  if (!flag)
    // no values!
    throw ExceptionObject(__FILE__, __LINE__);

  // find minimum value
  for (int i = 0; i < size; ++i)
    if ( valid[i] && (timeStepList[i] < min) )
      min = timeStepList[i];

  return min;
}

template <class TInputImage, class TOutputImage >
bool
MultiphaseFiniteDifferenceImageFilter<TInputImage, TOutputImage >
::Halt()
{
  if (m_NumberOfIterations != 0)
    this->UpdateProgress( static_cast<float>( this->GetElapsedIterations() )
      / static_cast<float>( m_NumberOfIterations ) );

  if( (this->GetElapsedIterations() >= m_NumberOfIterations) ||
      ( this->GetMaximumRMSError() > m_RMSChange ) )
    return true;
  else
    return false;
}


template <class TInputImage, class TOutputImage >
void
MultiphaseFiniteDifferenceImageFilter<TInputImage, TOutputImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "ElapsedIterations: " << m_ElapsedIterations << std::endl;
  os << indent << "UseImageSpacing: " << (m_UseImageSpacing ? "On" : "Off") << std::endl;
  os << indent << "State: " << m_State << std::endl;
  os << indent << "MaximumRMSError: " << m_MaximumRMSError << std::endl;
  os << indent << "NumberOfIterations: " << m_NumberOfIterations << std::endl;
  os << indent << "ManualReinitialization: " << m_ManualReinitialization << std::endl;
  os << indent << "RMSChange: " << m_RMSChange << std::endl;
  os << std::endl;

  if ( GetDifferenceFunction(0) )
  {
    os << indent << "DifferenceFunction: " << std::endl;
    for( unsigned int i = 0; i < this->FunctionCount; i++ )
		  GetDifferenceFunction( i )->Print(os,indent.GetNextIndent());
  }
  else
    os << indent << "DifferenceFunction: " << "(None)" << std::endl;
  os << std::endl;
}


}// end namespace itk

#endif
