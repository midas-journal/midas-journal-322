/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkSegmentationLevelSetImageFilter2.h,v $
  Language:  C++
  Date:      $Date: 2006/04/05 13:59:37 $
  Version:   $Revision: 1.30 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMultiphaseDenseSegmentationLevelSetImageFilterDense_h_
#define __itkMultiphaseDenseSegmentationLevelSetImageFilterDense_h_

#include "itkMultiphaseDenseSegmentationFiniteDifferenceImageFilter.h"
#include "itkSegmentationLevelSetFunction.h"

namespace itk {

/**
 *
 * \class SegmentationLevelSetImageFilter2
 *
 * \brief A base class which defines the API for implementing a special class of
 * image segmentation filters using level set methods.
 *
 * \par OVERVIEW
 * This object defines the framework for a class of segmentation filters which
 * use level set methods.  These filters work by constructing a ``feature
image''
 * onto which the evolving level set locks as it moves.  In the feature image,
 * values that are close to zero are associated with object boundaries.  An
 * original (or preprocessed) image is given to the filter as the feature image
 * and a seed for the level set is given as the input of the filter.  The seed
is
 * converted into a level set embedding which propagates according to the
features
 * calculated from the original image.
 *
 * \par TEMPLATE PARAMETERS
 * There are two required and two optional template parameter for these
 * filters.
 *
 * TInputImage is the image type of the initial model(s) you will input to the
 * filter using SetInput() or SetInitialImage().
 *
 * TFeatureImage is the image type of the image from which the filter will
 * calculate the speed term for segmentation (see INPUTS).
 *
 * TOutputPixel is the data type used for the output image phi, the implicit
 * level set image.  This should really only ever be set as float (default) or
 * double.
 *
 * FunctionCount is the number of level sets which are being evolved by this
filter.
 *
 * \par INPUTS
 * The input to any subclass of this filter is the seed image(s) for the initial
 * level set embedding.  As with other subclasses of the
 * SparseLevelSetImageFilter, the type of the input image is is not important.
 * The (RequestedRegion) size of the seed image must, however, match the
 * (RequestedRegion) size of the feature image.
 *
 *
 * \par A NOTE ON DATA SPACING
 * Input data with anisotropic spacing should be resampled to isotropic voxels
 * for use in this filter.  Interpolation has not been built into the
 * underlying p.d.e. solver.
 *
 * \par
 * Depending on the particular application and filter that you are using, the
 * feature image should be preprocessed with some type of noise reduction
 * filtering.  The feature image input can be of any type, but it will be cast
 * to floating point before calculations are done.
 *
 * \par OUTPUTS
 * The output of any subclass of this filter is a level set embedding as
 * described in SparseFieldLevelSetImageFilter2.  The zero crossings of the
output
 * image give the pixels closest to the level set boundary.  By ITK convention,
 * NEGATIVE values are pixels INSIDE the segmented region and POSITIVE values
are
 * pixels OUTSIDE the segmented region.
 *
 * \par PARAMETERS
 * The MaximumRMSChange parameter is used to determine when the solution has
 * converged.  A lower value will result in a tighter-fitting solution, but
 * will require more computations.  Too low a value could put the solver into
 * an infinite loop unless a reasonable NumberOfIterations parameter is set.
 * Values should always be greater than 0.0 and less than 1.0.
 *
 * \par
 * The NumberOfIterations parameter can be used to halt the solution after a
 * specified number of iterations, overriding the MaximumRMSChange halting
 * criteria.
 *
 * \par
 * The standard convention for ITK level-set segmentation filters is that
 * POSITIVE propagation (speed) and advection terms cause the surface to EXPAND
 * while negative terms cause the surface to CONTRACT.  When the
 * ReverseExpansionDirection parameter is set to TRUE (on), it tells the
 * function object to reverse the standard ITK convention so that NEGATIVE
 * terms cause EXPANSION and positive terms cause CONTRACTION.
 *
 * This parameter can be safely changed as appropriate for a particular
 * application or data set to achieve the desired behavior.
 *
 * \par
 * Unlike the original SegmentationLevelSetImageFilter, this filter does not
provide
 * parameter setting for the level set functions.  You should request a pointer
to
 * functions and set parameters specifically through their interface
 * \par
 *  See LevelSetFunction for more information.*/
template < class TInputImage, class TFeatureImage,
  typename TOutputPixel = typename TInputImage::PointType::CoordRepType >
class ITK_EXPORT DenseMultiphaseSegmentationLevelSetImageFilter
  : public DenseSegmentationFiniteDifferenceImageFilter< TInputImage,
    Image< TOutputPixel,
    ::itk::GetImageDimension< TInputImage>::ImageDimension > >
{
public:
  /** Inherited typedef from the superclass. Needs to be placed befroe
      the next macro. */
  typedef DenseMultiphaseSegmentationLevelSetImageFilter Self;

  /** Repeat definition from Superclass to satisfy Borland compiler
      quirks */
  itkStaticConstMacro( ImageDimension, unsigned int,
    TInputImage::ImageDimension );

  /** Output image type typedefs */
  typedef Image< TOutputPixel, itkGetStaticConstMacro( ImageDimension ) >
    OutputImageType;

  /** Standard class typedefs */
  typedef DenseSegmentationFiniteDifferenceImageFilter<TInputImage,
		OutputImageType > Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  typedef typename Superclass::FiniteDifferenceFunctionType
    FiniteDifferenceFunctionType;

  /** The data type used in numerical computations.  Derived from the output
   *  image type. */
  typedef typename OutputImageType::ValueType ValueType;
  typedef typename OutputImageType::IndexType IndexType;

  /** Inherited typedef from the superclass. */

  typedef typename Superclass::TimeStepType TimeStepType;
  typedef typename Superclass::InputImageType InputImageType;

  /** Local image typedefs */
  typedef TFeatureImage FeatureImageType;

  /** The generic level set function type */
  typedef SegmentationLevelSetFunction< OutputImageType, FeatureImageType >
  SegmentationFunctionType;

  /** The type used for the advection field */
  typedef typename SegmentationFunctionType::VectorImageType VectorImageType;
  typedef typename SegmentationFunctionType::ImageType SpeedImageType;

  /** Run-time type information (and related methods). */
  itkTypeMacro( DenseMultiphaseSegmentationLevelSetImageFilter,
    MultiphaseDenseSegmentationFiniteDifferenceImageFilter );


  /** Set the segmentation function.  In general, this should only be
   *  called by a subclass of this object. It is made public to allow
   *  itk::Command objects access. The method is inline to avoid a
   *  problem with the gcc 2.95 compiler matching the declaration with
   *  the definition. */
  virtual void SetSegmentationFunction( unsigned int functionIndex,
    SegmentationFunctionType *s )
  {
    m_SegmentationFunctions[functionIndex] = s;

    typename SegmentationFunctionType::RadiusType r;
    r.Fill( 1 );

    m_SegmentationFunctions[functionIndex]->Initialize(r);
    this->SetDifferenceFunction( functionIndex, s );
    this->Modified();
  }

  /** Set/Get the maximum number of iterations allowed for the solver.  This
   *  prevents infinite loops if a solution "bounces". */
  void SetMaximumIterations( unsigned int i )
  {
    itkWarningMacro("SetMaximumIterations is deprecated. "
      " Please use SetNumberOfIterations instead.");
    this->SetNumberOfIterations(i);
  }

  unsigned int GetMaximumIterations() const
  {
    itkWarningMacro("GetMaximumIterations is deprecated. "
      " Please use GetNumberOfIterations instead.");
    return this->GetNumberOfIterations();
  }

  /** Set/Get the feature image to be used for speed function of the level set
   *  equation.  Equivalent to calling Set/GetInput(1, ..) */
  virtual void SetFeatureImage(const FeatureImageType *f)
  {
    this->ProcessObject::SetNthInput( 0, const_cast< FeatureImageType * >(f) );
    for( unsigned int i = 0; i < this->FunctionCount; i++)
      m_SegmentationFunctions[i]->SetFeatureImage(f);
  }

  virtual FeatureImageType * GetFeatureImage()
  {
		return (static_cast< FeatureImageType*>(this->ProcessObject::GetInput(0)));
	}

  /** Set/Get the initial level set model.  Equivalent to calling SetInput(..)
   */
  virtual void SetInitialImage(InputImageType *f)
  {
    this->SetNthInput(0, f);
  }

  /** THIS METHOD IS DEPRECATED AND SHOULD NOT BE USED.  This method reverses
   * the speed function direction, effectively changing inside feature values to
   * outside feature values and vice versa. */
  void SetUseNegativeFeaturesOn()
  {
    itkWarningMacro( << "SetUseNegativeFeaturesOn has been deprecated. "
      " Please use ReverseExpansionDirectionOn() instead" );
    this->ReverseExpansionDirectionOn();
  }
  void SetUseNegativeFeaturesOff()
  {
    itkWarningMacro( << "SetUseNegativeFeaturesOff has been deprecated. "
      " Please use ReverseExpansionDirectionOff() instead" );
    this->ReverseExpansionDirectionOff();
  }

  /** THIS METHOD IS DEPRECATED AND SHOULD NOT BE USED. Set/Get the
   * value of the UseNegativeFeatures flag.  This method is
   * deprecated.  Use Set/Get ReverseExpansionDirection instead.*/
  void SetUseNegativeFeatures( bool u )
  {
    itkWarningMacro( << "SetUseNegativeFeatures has been deprecated. "
      " Please use SetReverseExpansionDirection instead" );
    if (u == true)
      {
      this->SetReverseExpansionDirection(false);
      }
    else
      {
      this->SetReverseExpansionDirection(true);
      }
  }
  bool GetUseNegativeFeatures() const
  {
    itkWarningMacro( << "GetUseNegativeFeatures has been deprecated. "
      "Please use GetReverseExpansionDirection() instead" );
    if ( m_ReverseExpansionDirection == false)
      {
      return true;
      }
    else
      {
      return false;
      }
  }

  /** Turn On/Off the flag which determines whether Positive or Negative speed
   * terms will cause surface expansion.  If set to TRUE then negative speed
   * terms will cause the surface to expand and positive speed terms will cause
   * the surface to contract.  If set to FALSE (default) then positive speed
   * terms will
   * cause the surface to expand and negative speed terms will cause the
   * surface to contract.  This method can be safely used to reverse the
   * expansion/contraction as appropriate to a particular application or data
   * set. */
  itkSetMacro(ReverseExpansionDirection, bool);
  itkGetConstMacro(ReverseExpansionDirection, bool);
  // itkBooleanMacro(ReverseExpansionDirection);

  /** Turn On/Off automatic generation of Speed and Advection terms when Update
      is called.  If set to Off, the Speed and Advection images must be set
      explicitly using SetSpeedImage, SetAdvectionImage OR the methods
      GenerateSpeedImage() and GenerateAdvectionImage() should be called prior
      to updating the filter. */
  itkSetMacro(AutoGenerateSpeedAdvection, bool);
  itkGetConstMacro(AutoGenerateSpeedAdvection, bool);
  // itkBooleanMacro(AutoGenerateSpeedAdvection);

  void UseMinimalCurvatureOn()
  {
    this->SetUseMinimalCurvature(true);
  }

  void UseMinimalCurvatureOff()
  {
    this->SetUseMinimalCurvature(false);
  }


  /** Allocate and calculate the speed term image in the SegmentationFunction
      object.  This method is called automatically on filter execution
      unless AutoGenerateSpeedAdvection is set to Off.*/
  void GenerateSpeedImage();

  /** Allocate and calculate the advection term image in the
      SegmentationFunction object  This method is called automatically on
      filter execution unless AutoGenerateSpeedAdvection is set to Off.*/
  void GenerateAdvectionImage();

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(OutputHasNumericTraitsCheck,
                  (Concept::HasNumericTraits<TOutputPixel>));
  /** End concept checking */
#endif

protected:
  DenseMultiphaseSegmentationLevelSetImageFilter();
  virtual ~DenseMultiphaseSegmentationLevelSetImageFilter()
  {
    delete [] m_SegmentationFunctions;
  }

  void SetFunctionCount( unsigned int n )
  {
    Superclass::SetFunctionCount( n );

    m_SegmentationFunctions = new
      SegmentationFunctionType*[ this->FunctionCount ];
  }

  DenseMultiphaseSegmentationLevelSetImageFilter(const Self&);

  virtual void PrintSelf(std::ostream& os, Indent indent) const;

  /** Overrides parent implementation */
  // This function is called at the end of each iteration
  virtual void InitializeIteration()
  {
    Superclass::InitializeIteration();
    // Estimate the progress of the filter
    this->SetProgress( (float) ( (float)this->GetElapsedIterations()
      / (float)this->GetNumberOfIterations()) );
  }

  /** Overridden from ProcessObject to set certain values before starting the
   * finite difference solver and then create an appropriate output */
  void GenerateData();

  /** Flag which sets the inward/outward direction of propagation speed. See
      SetReverseExpansionDirection for more information. */
  bool m_ReverseExpansionDirection;

  /** Flag to indicate whether Speed and Advection images are automatically
   *  generated when running the filter.  Otherwise, a pointer to images must
   *  be explicitly set or GenerateSpeedImage() and/or GenerateAdvectionImage()
   *  called directly before updating the filter */
  bool m_AutoGenerateSpeedAdvection;

private:
  // Define a pointer to multiphase level-set functions
  SegmentationFunctionType **m_SegmentationFunctions;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDenseMultiphaseSegmentationLevelSetImageFilter.txx"
#endif

#endif

