#ifndef __itkRegionBasedLevelSetFunction_txx_
#define __itkRegionBasedLevelSetFunction_txx_

#include "itkRegionBasedLevelSetFunction.h"

namespace itk
{

template < class TInputImage, class TFeatureImage, class TSharedData >
RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >
::RegionBasedLevelSetFunction()
{
  m_Nu = 0.0;
	m_Lambda1 = 1;
	m_Lambda2 = 1;
	m_Gamma = 0;
	m_Tau = 0;
	m_Volume = 0;
  m_Epsilon = 1.0;
  m_inv_epsilon = 1.0;
  m_2_over_pi = 2.0/vnl_math::pi;
  sharedData = 0;
  updatedC = false;
  updatedH = false;
}

template < class TInputImage, class TFeatureImage, class TSharedData >
void RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >::
SetMu( const double& mu)
{
  this->SetCurvatureWeight(mu);
}


template < class TInputImage, class TFeatureImage, class TSharedData >
double RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >::
GetMu() const
{
  return this->GetCurvatureWeight();
}

template < class TInputImage, class TFeatureImage, class TSharedData >
void RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >::
SetNu( const double& nu)
{
  this->m_Nu = nu;
}

template < class TInputImage, class TFeatureImage, class TSharedData >
double RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >::
GetNu() const
{
  return this->m_Nu;
}

template < class TInputImage, class TFeatureImage, class TSharedData >
void RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >::
SetLambda1( const double& lambda1 )
{
  this->m_Lambda1 = lambda1;
}

template < class TInputImage, class TFeatureImage, class TSharedData >
double RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >::
GetLambda1() const
{
  return this->m_Lambda1;
}

template < class TInputImage, class TFeatureImage, class TSharedData >
void RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >::
SetLambda2( const double& lambda2 )
{
  this->m_Lambda2 = lambda2;
}

template < class TInputImage, class TFeatureImage, class TSharedData >
double RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >::
GetLambda2() const
{
  return this->m_Lambda2;
}

template < class TInputImage, class TFeatureImage, class TSharedData >
void RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >::
SetGamma(const double& gamma)
{
  this->m_Gamma = gamma;
}

template < class TInputImage, class TFeatureImage, class TSharedData >
double RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >::
GetGamma() const
{
  return this->m_Gamma;
}

template < class TInputImage, class TFeatureImage, class TSharedData >
void RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >::
SetEta(const double& eta)
{
  this->SetLaplacianSmoothingWeight(eta);
}

template < class TInputImage, class TFeatureImage, class TSharedData >
double RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >::
GetEta() const
{
  return this->GetLaplacianSmoothingWeight();
}

template < class TInputImage, class TFeatureImage, class TSharedData >
void RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >::
SetEpsilon(const double& e)
{
  assert( e > 0. );
  this->m_Epsilon = e;
  this->m_inv_epsilon = 1. / e;
}

template < class TInputImage, class TFeatureImage, class TSharedData >
double RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >::
GetEpsilon() const
{
  return this->m_Epsilon;
}

template < class TInputImage, class TFeatureImage, class TSharedData >
void RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >::
SetTau(const double& t)
{
  this->m_Tau = t;
}

template < class TInputImage, class TFeatureImage, class TSharedData >
double RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >::
GetTau() const
{
  return this->m_Tau;
}

template < class TInputImage, class TFeatureImage, class TSharedData >
void RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >::
SetVolume( const double& v)
{
  this->m_Volume = v;
}

template < class TInputImage, class TFeatureImage, class TSharedData >
double RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >::
GetVolume() const
{
  return this->m_Volume;
}

template < class TInputImage, class TFeatureImage, class TSharedData >
void RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >::
SetFunctionId( const unsigned int& functionId )
{
  this->m_FunctionId = functionId;
}


/* Computes the Heaviside function and stores it in hVals */
template < class TInputImage, class TFeatureImage, class TSharedData >
void RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >
::ComputeHImage()
{
  // The phi function
  InputImageConstPointer contourImage = this->m_InitialImage;
  InputImagePointer hBuffer = this->sharedData->hVals[this->m_FunctionId];

  // Iterator for the phi function
  ConstImageIteratorType constIt( contourImage,
    contourImage->GetRequestedRegion() );
  ImageIteratorType It( hBuffer, hBuffer->GetRequestedRegion() );

  for( It.GoToBegin(), constIt.GoToBegin(); !constIt.IsAtEnd();
    ++It, ++constIt )
  {
    double hVal = this->calculateH( constIt.Get() );
    It.Set( hVal );
  }
}


template < class TInputImage, class TFeatureImage, class TSharedData >
double RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >
::calculateH( const InputPixelType& z )
{
/* z is the distance from the contour as given by the phi function */

// 	if(z < -m_epsilon)
// 		return 1;
// 	if(z > m_epsilon)
// 		return 0;
//
// 	return 0.5*(1.0 + (-z/m_epsilon) +
//     (1.0/vnl_math::pi)*sin(vnl_math::pi*(-z)/m_epsilon));

/* Opposite convention of inside is negative, therefore use -z */
// 	return 0.5*(1.0 + (2.0/vnl_math::pi) * atan( -z * m_inv_epsilon ) );
  return 0.5*(1.0 + m_2_over_pi * vcl_atan( -z * m_inv_epsilon ) );
}

template < class TInputImage, class TFeatureImage, class TSharedData >
double RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >
::calculatedH( const InputPixelType& z )
{
/* Opposite convention of inside is negative, therefore use -z */
//	if( z < -m_epsilon || z > m_epsilon )
//		return 0;
//  return 0.5/epsilon * ( 1.0 + cos(vnl_math::pi*(-z)/m_epsilon));

// 	return (1.0/vnl_math::pi)*(1.0/(1.0 + (z*z) / (m_epsilon*m_epsilon) ));
  double t = z*m_inv_epsilon;
  return 0.5*m_2_over_pi*(1.0 + t*t);
}

template < class TInputImage, class TFeatureImage, class TSharedData >
void
RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >
::UpdateSharedData( bool forceUpdate )
  {
    if ( forceUpdate )
    {
      updatedC = false;
      updatedH = false;
    }
    if ( !updatedH )
    {
      // update H
      updatedH = true;
      this->ComputeHImage();

      // Must update all H before updating C
      return;
    }
    if ( !updatedC )
    {
      updatedC = true;
      this->ComputeParameters();
    }

    this->SpecialProcessing();

    unsigned int fId = this->m_FunctionId;

    if ( sharedData->cDens[fId] == 0 )
      sharedData->cVals[fId] = 0;
    else
      sharedData->cVals[fId] = sharedData->cNums[fId] / sharedData->cDens[fId];

    if ( fId == 0 )
    {
      if ( sharedData->cBgrndDen == 0 )
        sharedData->cBgrnd = 0;
      else
        sharedData->cBgrnd = sharedData->cBgrndNum / sharedData->cBgrndDen;

      std::cout << "Background: " << std::endl;
//      std::cout << this->sharedData->cBgrndNum << ' ' << this->sharedData->cBgrndDen << ' ';
      std::cout << this->sharedData->cBgrnd << std::endl;
    }
//     std::cout << "Foreground: " << std::endl;
//     std::cout << sharedData->cNums[fId] << ' ' << this->sharedData->cDens[fId]
//       << ' ' << this->sharedData->cVals[fId] << std::endl;
  }

template < class TInputImage, class TFeatureImage, class TSharedData >
typename RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >::TimeStepType
RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >
::ComputeGlobalTimeStep(void *GlobalData) const
{
/* Computing the time-step for stable curve evolution */

  TimeStepType dt;

  ACGlobalDataStruct *d = (ACGlobalDataStruct *)GlobalData;

  if (vnl_math_abs(d->m_MaxCurvatureChange) > 0.0)
  {
    if (d->m_MaxGlobalChange > 0.0)
    {
      dt = vnl_math_min( ( this->m_WaveDT / d->m_MaxGlobalChange ),
      ( this->m_DT / d->m_MaxCurvatureChange ) );
    }
    else
    {
      dt = this->m_DT / d->m_MaxCurvatureChange;
    }
  }
  else
  {
    if (d->m_MaxGlobalChange > 0.0)
    {
//       std::cout << this->m_WaveDT << std::endl;
      dt = this->m_WaveDT / d->m_MaxGlobalChange;
    }
    else
    {
      dt = 0.0;
    }
  }

  // Reset the values
  d->m_MaxAdvectionChange   = NumericTraits<ScalarValueType>::Zero;
  d->m_MaxPropagationChange = NumericTraits<ScalarValueType>::Zero;
  d->m_MaxCurvatureChange   = NumericTraits<ScalarValueType>::Zero;
  d->m_MaxGlobalChange		  = NumericTraits<ScalarValueType>::Zero;

  return dt;
}


template < class TInputImage, class TFeatureImage, class TSharedData >
typename RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >::PixelType
RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >
::ComputeUpdate( const NeighborhoodType &it, void *globalData,
  const FloatOffsetType& offset )
{
  // Access the neighborhood center pixel of phi
  const ScalarValueType inputValue = it.GetCenterPixel();

  ScalarValueType laplacian, laplacian_term, curvature_term, globalTerm;

  // Access the global data structure
  ACGlobalDataStruct *gd = (ACGlobalDataStruct *)globalData;

  // Compute the Hessian matrix and various other derivatives.  Some of these
  // derivatives may be used by overloaded virtual functions.
  gd->m_GradMagSqr = 1.0e-6;
  for( unsigned int i = 0 ; i < ImageDimension; i++)
  {
    const unsigned int positionA =
      static_cast< unsigned int >( this->m_Center + this->m_xStride[i] );
    const unsigned int positionB =
      static_cast< unsigned int >( this->m_Center - this->m_xStride[i] );

    gd->m_dx[i] = 0.5 * ( it.GetPixel( positionA ) - it.GetPixel( positionB ) );
    gd->m_dxy[i][i] =
      it.GetPixel( positionA ) + it.GetPixel( positionB ) - 2.0 * inputValue;
    gd->m_dx_forward[i]  = it.GetPixel( positionA ) - inputValue;
    gd->m_dx_backward[i] = inputValue - it.GetPixel( positionB );
    gd->m_GradMagSqr += gd->m_dx[i] * gd->m_dx[i];

    for( unsigned int j = i+1; j < ImageDimension; j++ )
    {
      const unsigned int positionAa = static_cast<unsigned int>(
        this->m_Center - this->m_xStride[i] - this->m_xStride[j] );
      const unsigned int positionBa = static_cast<unsigned int>(
        this->m_Center - this->m_xStride[i] + this->m_xStride[j] );
      const unsigned int positionCa = static_cast<unsigned int>(
        this->m_Center + this->m_xStride[i] - this->m_xStride[j] );
      const unsigned int positionDa = static_cast<unsigned int>(
        this->m_Center + this->m_xStride[i] + this->m_xStride[j] );

      gd->m_dxy[i][j] = gd->m_dxy[j][i] = 0.25 *(
        it.GetPixel( positionAa ) -
        it.GetPixel( positionBa ) -
        it.GetPixel( positionCa ) +
        it.GetPixel( positionDa ) );
    }
  }

  ScalarValueType dh = calculatedH( inputValue );

  // Computing the curvature term
  // Used to regularized using the length of contour
  // What is CurvatureSpeed?
  if ( ( dh != 0. ) &&
    ( this->m_CurvatureWeight != NumericTraits< ScalarValueType >::Zero ) )
  {
    curvature_term =
      this->m_CurvatureWeight *
      this->CurvatureSpeed(it, offset) *
      this->ComputeCurvatureTerm( it, offset, gd ) *
      gd->m_GradMagSqr * dh;

    gd->m_MaxCurvatureChange =
      vnl_math_max( gd->m_MaxCurvatureChange, vnl_math_abs( curvature_term ) );
  }
  else
    curvature_term = NumericTraits<ScalarValueType>::Zero;

  // Computing the laplacian term
  // Used in maintaining squared distance function
  if( this->m_LaplacianSmoothingWeight != NumericTraits<ScalarValueType>::Zero)
  {
    laplacian = NumericTraits<ScalarValueType>::Zero;

    // Compute the laplacian using the existing second derivative values
    for(unsigned int i = 0;i < ImageDimension; i++)
    {
      laplacian += gd->m_dxy[i][i];
    }

	  // Use the laplacian to maintain signed distance function
    // What is LaplacianSmoothingSpeed ?
    // Why do we have 0.1 * LaplacianSmoothingWeight ?
    // Why do we have to subtract the curvature_term ?
    laplacian_term =
      0.1*( this->m_LaplacianSmoothingWeight *
      LaplacianSmoothingSpeed(it,offset, gd) * laplacian - curvature_term);
  }
  else
    laplacian_term = NumericTraits<ScalarValueType>::Zero;


  // Update value from curvature length and laplacian term
  PixelType updateVal =
    static_cast< PixelType >( curvature_term + laplacian_term );

  /* Compute the globalTerm - rms difference of image with c_0 or c_1*/
  if ( dh != 0. )
    globalTerm = dh * this->computeGlobalTerm( inputValue, it.GetIndex() );
  else
    globalTerm = 0;

  /* Final update value is the local terms of curvature lengths and laplacian
  squared distances - global terms of rms differences of image and piecewise
  constant regions*/
  updateVal = updateVal - globalTerm;

  /* If MaxGlobalChange recorded is lower than the current globalTerm */
  if(vnl_math_abs( gd->m_MaxGlobalChange) < vnl_math_abs( globalTerm ) )
  {
    gd->m_MaxGlobalChange = globalTerm;
  }

  return updateVal;
}

template < class TInputImage, class TFeatureImage, class TSharedData >
double
RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >
::computeOverlapTerm( const unsigned int& s )
{
  return this->m_Gamma * s;
}

/* Computes the fidelity term (eg: (intensity - mean)2 ).
Most of the code is concerned with using the appropriate combination
of Heaviside and dirac delta for each part of the fidelity term.
- the final dH is the dirac delta term corresponding to the current
level set we are updating. */
template < class TInputImage, class TFeatureImage, class TSharedData >
double
RegionBasedLevelSetFunction< TInputImage, TFeatureImage, TSharedData >
::computeGlobalTerm(
const InputPixelType& inputPixel,
const InputIndexType& inputIndex )
{
  unsigned int fId = this->m_FunctionId;
  unsigned int pr = 1; // computes if it belongs to background
  unsigned int s = 0; // accumulates the overlap across all functions

  // Assuming only 1 level set function to be present
  FeatureIndexType featIndex = static_cast< FeatureIndexType >( inputIndex );
  ScalarValueType globalTerm = 0;

  double overlapTerm = 0.;
  // This conditional statement computes the amount of overlap s
  // and the presence of background pr
  if ( this->sharedData->FunctionCount > 1 )
  {
    featIndex = this->sharedData->GetFeatureIndex( fId, inputIndex );
    computeOverlapParameters( featIndex, s, pr );
  }

  const FeaturePixelType featureVal =
    this->m_FeatureImage->GetPixel ( featIndex );

  double inTerm = this->computeInternalTerm( featureVal, featIndex, fId );
  double outTerm = this->m_Lambda2 * this->computeExternalTerm( featureVal, featIndex, pr );
  overlapTerm = this->computeOverlapTerm( s );
  double regularizationTerm = 2 * this->m_Tau * ( sharedData->cDens[fId] - this->m_Volume );
  //regularizationTerm -= this->m_Nu;
  //NOTE: regularizationTerm here MUST take into account the curvature term!!!

  globalTerm = - inTerm + outTerm - overlapTerm - regularizationTerm;

  return globalTerm;
}

} // end namespace

#endif
