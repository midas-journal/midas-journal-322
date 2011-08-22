/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDenseSegmentationFiniteDifferenceImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2005/08/09 12:15:42 $
  Version:   $Revision: 1.28 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMultiphaseDenseSegmentationFiniteDifferenceImageFilter_txx_
#define __itkMultiphaseDenseSegmentationFiniteDifferenceImageFilter_txx_
#include "itkMultiphaseDenseSegmentationFiniteDifferenceImageFilter.h"

namespace itk
{

  template <class TInputImage, class TOutputImage >
  void
  DenseSegmentationFiniteDifferenceImageFilter<TInputImage, TOutputImage >
  ::CopyInputToOutput()
  {
    OutputImagePointer output = this->GetOutput();
    for ( unsigned int i = 0; i < this->FunctionCount; i++ )
    {
      InputImagePointer input = this->m_LevelSet[i];
      InputPointType origin = input->GetOrigin();
      InputSpacingType spacing = input->GetSpacing();
      InputSizeType size = input->GetLargestPossibleRegion().GetSize();

			OutputIndexType start;
			for ( unsigned int j = 0; j < ImageDimension; j++ )
				start[j] = static_cast<OutputIndexValueType>( origin[j]/spacing[j] );

      OutputRegionType region;
      region.SetSize( size );
			region.SetIndex( start );

      if ( !input || !output )
        itkExceptionMacro ( << "Either input and/or output is NULL." );

      ImageRegionIterator< InputImageType > in ( input,
        input->GetLargestPossibleRegion() );
      ImageRegionIterator< OutputImageType > out ( output, region );

      for ( in.Begin(), out.Begin(); !out.IsAtEnd(); ++in, ++out )
        if ( in.Get() < 0 )
          out.Value() =  static_cast<PixelType> ( this->m_Lookup[i] );
        else
          out.Value() = 0;
    }
  }

  template < class TInputImage, class TOutputImage >
  void
  DenseSegmentationFiniteDifferenceImageFilter<TInputImage, TOutputImage >
  ::AllocateUpdateBuffer()
  {
    for ( unsigned int i = 0; i < this->FunctionCount; i++ )
    {
			InputImagePointer input = this->m_LevelSet[i];
      InputRegionType region = input->GetLargestPossibleRegion();

      m_UpdateBuffers[i]->SetLargestPossibleRegion ( region );
      m_UpdateBuffers[i]->SetRequestedRegion ( region );
      m_UpdateBuffers[i]->SetBufferedRegion ( region );
      m_UpdateBuffers[i]->Allocate();
			m_UpdateBuffers[i]->SetSpacing( input->GetSpacing() );
			m_UpdateBuffers[i]->SetOrigin( input->GetOrigin() );
    }
  }

	template < class TInputImage, class TOutputImage >
  typename
  DenseSegmentationFiniteDifferenceImageFilter< TInputImage,
    TOutputImage >::TimeStepType
  DenseSegmentationFiniteDifferenceImageFilter< TInputImage, TOutputImage >
  ::CalculateChange()
  {
    TimeStepType timeStep = 9999;

    for( unsigned int i = 0; i < this->FunctionCount; i++ )
    {
      OutputImagePointer levelset = this->m_LevelSet[i];
      TimeStepType dt;

      // Get the FiniteDifferenceFunction to use in calculations.
      const typename FiniteDifferenceFunctionType::Pointer df
      = this->GetDifferenceFunction ( i );

      const OutputSizeType radius = df->GetRadius();

      // Break the input into a series of regions.  The first region is free
      // of boundary conditions, the rest with boundary conditions.  We operate
      // on the levelset region because input has been copied to output.
      FaceCalculatorType faceCalculator;
      FaceListType faceList = faceCalculator ( levelset,
        levelset->GetLargestPossibleRegion(), radius );

      void *globalData;

      // Ask the function object for a pointer to a data structure it
      // will use to manage any global values it needs.  We'll pass this
      // back to the function object at each calculation and then
      // again so that the function object can use it to determine a
      // time step for this iteration.
      globalData = df->GetGlobalDataPointer();

      typename FaceListType::iterator fIt;
			for ( fIt = faceList.begin(); fIt != faceList.end(); ++fIt )
			{
				// Process the non-boundary region.
				NeighborhoodIteratorType nD ( radius, levelset, *fIt );
				UpdateIteratorType nU ( m_UpdateBuffers[i], *fIt );

				for ( nD.GoToBegin(), nD.GoToBegin(); !nD.IsAtEnd(); ++nD, ++nU )
					nU.Value() = df->ComputeUpdate ( nD, globalData );
			}

      // Ask the finite difference function to compute the time step for
      // this iteration.  We give it the global data pointer to use, then
      // ask it to free the global data memory.
      dt = df->ComputeGlobalTimeStep ( globalData );
      df->ReleaseGlobalDataPointer ( globalData );

      if ( dt < timeStep )
        timeStep = dt;
    }
    timeStep = 0.08;
    return timeStep;
  }

  template< class TInputImage, class TOutputImage >
  void
  DenseSegmentationFiniteDifferenceImageFilter< TInputImage, TOutputImage >
  ::ApplyUpdate ( TimeStepType dt )
  {
		double rms_change_accumulator = 0, den = 0;
    for( unsigned int i = 0;  i < this->FunctionCount; i++)
    {
      double img_size = 1;
      for( unsigned int j = 0; j < ImageDimension; j++)
        img_size *=
          this->m_LevelSet[i]->GetLargestPossibleRegion().GetSize()[j];
      den += img_size;
    }

    // Updating the output image
    for ( unsigned int i = 0;  i < this->FunctionCount; i++ )
    {
      InputRegionType region = this->m_LevelSet[i]->GetRequestedRegion();

      ImageRegionIterator<UpdateBufferType> u ( m_UpdateBuffers[i], region );
      ImageRegionIterator<OutputImageType>  o ( this->m_LevelSet[i], region );

      for ( u.Begin(), o.Begin(); !u.IsAtEnd(); ++u, ++o )
      {
        u.Value() = o.Value() + static_cast<PixelType> ( u.Value() * dt );
      }

      if ( this->GetElapsedIterations()%reinitCounter == 0 )
      {
        typename ThresholdFilterType::Pointer thresh =
          ThresholdFilterType::New();
        thresh->SetLowerThreshold(
          NumericTraits< OutputPixelType >::NonpositiveMin() );
        thresh->SetUpperThreshold( 0 );
        thresh->SetInsideValue( 1 );
        thresh->SetOutsideValue( 0 );
        thresh->SetInput( m_UpdateBuffers[i] );
        thresh->Update();

        typename MaurerType::Pointer maurer = MaurerType::New();
        maurer->SetInput( thresh->GetOutput() );
        maurer->SetSquaredDistance( 0 );
        maurer->SetUseImageSpacing( 1 );
        maurer->SetInsideIsPositive( 0 );
        maurer->Update();

        ImageRegionIterator< OutputImageType >  it ( maurer->GetOutput(),
          region );
        for ( o.GoToBegin(), it.GoToBegin(); !o.IsAtEnd(); ++o, ++it )
        {
          rms_change_accumulator  += static_cast<double> (
            vnl_math_sqr( o.Value() - it.Value() ) );
          o.Value() = it.Value();
        }
      }
    }

    this->SetRMSChange( vcl_sqrt(rms_change_accumulator / den ) );
  }

  template< class TInputImage, class TOutputImage >
  void
  DenseSegmentationFiniteDifferenceImageFilter< TInputImage, TOutputImage >
  ::PostProcessOutput()
  {
    this->CopyInputToOutput();
  }

}// end namespace itk

#endif
