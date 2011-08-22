#ifndef __itkMultiphaseSparseFieldLevelSetImageFilter2_txx_
#define __itkMultiphaseSparseFieldLevelSetImageFilter2_txx_

namespace itk
{

  template <class TNeighborhoodType>
  SparseFieldCityBlockNeighborList2<TNeighborhoodType>
  ::SparseFieldCityBlockNeighborList2()
  {
    typedef typename NeighborhoodType::ImageType ImageType;
    typename ImageType::Pointer dummy_image = ImageType::New();

    unsigned int i, nCenter;
    int d;
    OffsetType zero_offset;

    for ( i = 0; i < Dimension; ++i )
    {
      m_Radius[i] = 1;
      zero_offset[i] = 0;
    }
    NeighborhoodType it ( m_Radius, dummy_image,
      dummy_image->GetRequestedRegion() );
    nCenter = it.Size() / 2;

    m_Size = 2 * Dimension;
    m_ArrayIndex.reserve ( m_Size );
    m_NeighborhoodOffset.reserve ( m_Size );

    for ( i = 0; i < m_Size; ++i )
    {
      m_NeighborhoodOffset.push_back ( zero_offset );
    }

    for ( d = Dimension - 1, i = 0; d >= 0; --d, ++i )
    {
      m_ArrayIndex.push_back ( nCenter - it.GetStride ( d ) );
      m_NeighborhoodOffset[i][d] = -1;
    }

    for ( d = 0; d < static_cast<int> ( Dimension ); ++d, ++i )
    {
      m_ArrayIndex.push_back ( nCenter + it.GetStride ( d ) );
      m_NeighborhoodOffset[i][d] = 1;
    }

    for ( i = 0; i < Dimension; ++i )
    {
      m_StrideTable[i] = it.GetStride ( i );
    }
  }

  template <class TNeighborhoodType>
  void
  SparseFieldCityBlockNeighborList2<TNeighborhoodType>
  ::Print ( std::ostream &os ) const
  {
    os << "SparseFieldCityBlockNeighborList2: " << std::endl;
    for ( unsigned i = 0; i < this->GetSize(); ++i )
    {
      os << "m_ArrayIndex[" << i << "]: " << m_ArrayIndex[i] << std::endl;
      os << "m_NeighborhoodOffset[" << i << "]: " << m_NeighborhoodOffset[i]
      << std::endl;
    }
  }

  template<class TInputImage, class TOutputImage >
  double MultiphaseSparseFieldLevelSetImageFilter<TInputImage, TOutputImage >
  ::m_ConstantGradientValue = 1.0;

  template<class TInputImage, class TOutputImage >
  ITK_TYPENAME MultiphaseSparseFieldLevelSetImageFilter<TInputImage,
  TOutputImage >::ValueType
  MultiphaseSparseFieldLevelSetImageFilter<TInputImage, TOutputImage >
  ::m_ValueOne = NumericTraits<ITK_TYPENAME
    MultiphaseSparseFieldLevelSetImageFilter< TInputImage, TOutputImage >
    ::ValueType >::One;

  template<class TInputImage, class TOutputImage >
  ITK_TYPENAME MultiphaseSparseFieldLevelSetImageFilter< TInputImage,
  TOutputImage >::ValueType
  MultiphaseSparseFieldLevelSetImageFilter<TInputImage, TOutputImage >
  ::m_ValueZero = NumericTraits<ITK_TYPENAME
    MultiphaseSparseFieldLevelSetImageFilter< TInputImage, TOutputImage>::
    ValueType>::Zero;

  template< class TInputImage, class TOutputImage >
  ITK_TYPENAME MultiphaseSparseFieldLevelSetImageFilter< TInputImage,
  TOutputImage >::StatusType
  MultiphaseSparseFieldLevelSetImageFilter< TInputImage, TOutputImage >
  ::m_StatusNull = NumericTraits< ITK_TYPENAME
    MultiphaseSparseFieldLevelSetImageFilter<TInputImage, TOutputImage >::
    StatusType >::NonpositiveMin();

  template< class TInputImage, class TOutputImage >
  ITK_TYPENAME MultiphaseSparseFieldLevelSetImageFilter< TInputImage,
  TOutputImage >::StatusType
  MultiphaseSparseFieldLevelSetImageFilter< TInputImage, TOutputImage >
  ::m_StatusChanging = -1;

  template< class TInputImage, class TOutputImage >
  ITK_TYPENAME MultiphaseSparseFieldLevelSetImageFilter< TInputImage,
  TOutputImage >::StatusType
  MultiphaseSparseFieldLevelSetImageFilter< TInputImage, TOutputImage >
  ::m_StatusActiveChangingUp = -2;

  template< class TInputImage, class TOutputImage >
  ITK_TYPENAME MultiphaseSparseFieldLevelSetImageFilter< TInputImage,
  TOutputImage >::StatusType
  MultiphaseSparseFieldLevelSetImageFilter< TInputImage, TOutputImage >
  ::m_StatusActiveChangingDown = -3;

  template< class TInputImage, class TOutputImage >
  ITK_TYPENAME MultiphaseSparseFieldLevelSetImageFilter< TInputImage,
  TOutputImage >::StatusType
  MultiphaseSparseFieldLevelSetImageFilter< TInputImage, TOutputImage >
  ::m_StatusBoundaryPixel = -4;

template< class TInputImage, class TOutputImage >
MultiphaseSparseFieldLevelSetImageFilter< TInputImage, TOutputImage >
::MultiphaseSparseFieldLevelSetImageFilter()
{
  currentFunction = 0;
  m_IsoSurfaceValue = m_ValueZero;
  m_NumberOfLayers = ImageDimension;
  this->SetRMSChange ( static_cast< double > ( m_ValueZero ) );
  m_InterpolateSurfaceLocation = true;
  m_BoundsCheckingActive = false;
}

template<class TInputImage, class TOutputImage >
void
MultiphaseSparseFieldLevelSetImageFilter<TInputImage, TOutputImage >
::CopyInputToOutput()
{
  for(unsigned int i = 0; i < this->FunctionCount; i++)
  {
    SparseDataStruct *sparsePtr = sparseData[i];
    OutputImagePointer input = this->m_LevelSet[i];

    sparsePtr->m_ShiftedImage = OutputImageType::New();
    sparsePtr->m_ShiftedImage->SetRegions(
      input->GetRequestedRegion() );
    sparsePtr->m_ShiftedImage->Allocate();
    sparsePtr->m_ShiftedImage->CopyInformation( input );

    ImageRegionIterator<OutputImageType> lIt(
      input, input->GetRequestedRegion() );
    ImageRegionIterator<OutputImageType> sIt(
      sparsePtr->m_ShiftedImage, input->GetRequestedRegion() );
    for ( lIt.GoToBegin(), sIt.GoToBegin(); !lIt.IsAtEnd(); ++sIt, ++lIt )
      sIt.Set( lIt.Get() );

    {
      typename ZeroCrossingFilterType::Pointer
        zeroCrossingFilter = ZeroCrossingFilterType::New();
      zeroCrossingFilter->SetInput( sparsePtr->m_ShiftedImage );
      zeroCrossingFilter->SetBackgroundValue( m_ValueOne );
      zeroCrossingFilter->SetForegroundValue( m_ValueZero );
      zeroCrossingFilter->Update();

      this->m_LevelSet[i] = zeroCrossingFilter->GetOutput();
      this->m_LevelSet[i]->DisconnectPipeline();
    }
  }
}

template< class TInputImage, class TOutputImage >
void
MultiphaseSparseFieldLevelSetImageFilter< TInputImage, TOutputImage >
::InitializeIteration()
{
	Superclass::InitializeIteration();

  // Set the values in the output image for the active layer.
  this->InitializeActiveLayerValues();

  // Initialize layer values using the active layer as seeds.
  this->PropagateAllLayerValues();
}


template< class TInputImage, class TOutputImage >
void
MultiphaseSparseFieldLevelSetImageFilter< TInputImage, TOutputImage >
::Initialize()
{
	for ( unsigned int functionIndex = 0; functionIndex < this->FunctionCount;
		functionIndex++ )
	{
		SparseDataStruct *sparsePtr = sparseData[functionIndex];

		unsigned int i;

		// Allocate the status image.
		sparsePtr->m_StatusImage = StatusImageType::New();
		sparsePtr->m_StatusImage->SetRegions (
			this->m_LevelSet[functionIndex]->GetRequestedRegion() );
		sparsePtr->m_StatusImage->Allocate();
		sparsePtr->m_StatusImage->FillBuffer( m_StatusNull );//NonpositiveMin

		// Initialize the boundary pixels in the status image to
		// m_StatusBoundaryPixel values.  Uses the face calculator to find all of
		// the region faces.
		typedef NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<
			StatusImageType > BFCType;

		BFCType faceCalculator;
		typename BFCType::FaceListType faceList;
		typename BFCType::SizeType sz;
		typename BFCType::FaceListType::iterator fit;

		sz.Fill ( 1 ); // Set the difference function radius here
		faceList = faceCalculator ( sparsePtr->m_StatusImage,
			sparsePtr->m_StatusImage->GetRequestedRegion(), sz );
		fit = faceList.begin();

		for ( ++fit; fit != faceList.end(); ++fit )
		{
			ImageRegionIterator<StatusImageType> statusIt(
			sparsePtr->m_StatusImage, *fit );
			for ( statusIt.GoToBegin(); ! statusIt.IsAtEnd(); ++statusIt )
				statusIt.Set ( m_StatusBoundaryPixel );// -4
		}

		// Erase all existing layer lists.
		for ( i = 0; i < sparsePtr->m_Layers.size(); ++i )
		{
			while ( ! sparsePtr->m_Layers[i]->Empty() )
			{
				sparsePtr->m_LayerNodeStore->Return( sparsePtr->m_Layers[i]->Front());
				sparsePtr->m_Layers[i]->PopFront();
			}
		}

		// Allocate the layers for the sparse field.
		sparsePtr->m_Layers.clear();
		sparsePtr->m_Layers.reserve( 2*m_NumberOfLayers + 1 );

		while ( sparsePtr->m_Layers.size() < ( 2*m_NumberOfLayers+1 ) )
			sparsePtr->m_Layers.push_back( LayerType::New() );

		// Throw an exception if we don't have enough layers.
		if ( sparsePtr->m_Layers.size() < 3 )
		{
			itkExceptionMacro ( << "Not enough layers have been allocated for the"
				"sparse field.  Requires at least one layer." );
		}
	}

	// Construct the active layer and initialize the first layers inside and
	// outside of the active layer.
	this->ConstructActiveLayer();

	for ( unsigned int j = 0; j < this->FunctionCount; j++ )
	{
		SparseDataStruct *sparsePtr = sparseData[j];

		// Construct the rest of the non-active set layers using the first two
		// layers. Inside layers are odd numbers, outside layers are even numbers.
		for ( unsigned int i = 1; i < sparsePtr->m_Layers.size() - 2; ++i )
			this->ConstructLayer( sparsePtr, i, i+2 );
	}

  // Set the values in the output image for the active layer.
  this->InitializeActiveLayerValues();

  // Initialize layer values using the active layer as seeds.
  this->PropagateAllLayerValues();

  // Initialize pixels inside and outside the sparse field layers to positive
  // and negative values, respectively?? This is not necessary for the
  // calculations, but is useful for presenting a more intuitive output to the
  // filter.  See PostProcessOutput method for more information.
  this->InitializeBackgroundPixels();


  for ( unsigned int j = 0; j < this->FunctionCount; j++ )
  {
    SparseDataStruct *sparsePtr = sparseData[j];

    // Delete status image
    sparsePtr->m_ShiftedImage->Delete();
  }
}

template <class TInputImage, class TOutputImage >
void
MultiphaseSparseFieldLevelSetImageFilter< TInputImage, TOutputImage >
::InitializeBackgroundPixels()
{
	for ( unsigned int functionIndex = 0; functionIndex < this->FunctionCount;
		functionIndex++ )
	{
		SparseDataStruct *sparsePtr = sparseData[functionIndex];

		// Assign background pixels OUTSIDE the sparse field layers to a new level
		// set
		// with value greater than the outermost layer.  Assign background pixels
		// INSIDE the sparse field layers to a new level set with value less than
		// the innermost layer.
		const ValueType max_layer = static_cast<ValueType> ( m_NumberOfLayers );

		const ValueType outside_value  = max_layer + m_ConstantGradientValue;
		const ValueType inside_value = - ( max_layer + m_ConstantGradientValue );

		ImageRegionConstIterator<StatusImageType> statusIt (
			sparsePtr->m_StatusImage,
			this->m_LevelSet[functionIndex]->GetRequestedRegion() );

		ImageRegionIterator<OutputImageType> outputIt (
			this->m_LevelSet[functionIndex],
			this->m_LevelSet[functionIndex]->GetRequestedRegion() );

		ImageRegionIterator<OutputImageType> shiftedIt (
			sparsePtr->m_ShiftedImage,
			this->m_LevelSet[functionIndex]->GetRequestedRegion() );

		for ( outputIt = outputIt.Begin(), shiftedIt = shiftedIt.Begin(),
			statusIt = statusIt.Begin(); !outputIt.IsAtEnd();
			++shiftedIt, ++outputIt, ++statusIt )
		{
			if ( statusIt.Get() == m_StatusNull || statusIt.Get() ==
				m_StatusBoundaryPixel )
			{
				if ( shiftedIt.Get() > m_ValueZero )
					outputIt.Set ( outside_value );
				else
					outputIt.Set ( inside_value );
			}
		}
	}
}


template < class TInputImage, class TOutputImage >
void
MultiphaseSparseFieldLevelSetImageFilter< TInputImage, TOutputImage>
::ConstructActiveLayer()
{
	for ( unsigned int functionIndex = 0; functionIndex < this->FunctionCount;
		functionIndex++ )
	{
		SparseDataStruct *sparsePtr = sparseData[functionIndex];

//       We find the active layer by searching for 0's in the zero crossing
//       image (output image).  The first inside and outside layers are also
//       constructed by searching the neighbors of the active layer in the
//       (shifted) input image. Negative neighbors not in the active set are
//       assigned to the inside, positive neighbors are assigned to the outside.
//
//       During construction we also check whether any of the layers of the
//       active set (or the active set itself) is sitting on a boundary pixel
//       location. If this is the case, then we need to do active bounds
//       checking in the solver.

		unsigned int i;
		NeighborhoodIterator< OutputImageType >
			outputIt ( sparsePtr->m_NeighborList.GetRadius(),
			this->m_LevelSet[functionIndex],
			this->m_LevelSet[functionIndex]->GetRequestedRegion() );

		NeighborhoodIterator< StatusImageType >
			statusIt ( sparsePtr->m_NeighborList.GetRadius(),
			sparsePtr->m_StatusImage,
			this->m_LevelSet[functionIndex]->GetRequestedRegion() );

		NeighborhoodIterator< OutputImageType >
			shiftedIt ( sparsePtr->m_NeighborList.GetRadius(),
			sparsePtr->m_ShiftedImage,
			this->m_LevelSet[functionIndex]->GetRequestedRegion() );

		OutputIndexType center_index, offset_index;
		LayerNodeType *node;
		bool bounds_status;
		ValueType value;
		StatusType layer_number;

		OutputIndexType lowerBounds;
		OutputSizeType upperBounds;
		lowerBounds =
			this->m_LevelSet[functionIndex]->GetRequestedRegion().GetIndex();
		upperBounds =
			this->m_LevelSet[functionIndex]->GetRequestedRegion().GetSize();

		for ( outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt )
		{
			if ( outputIt.GetCenterPixel() == m_ValueZero )
			{
				// Grab the neighborhood in the status image.
				center_index = outputIt.GetIndex();
				statusIt.SetLocation ( center_index );

				// Check to see if any of the sparse field touches a boundary.  If so,
				// then activate bounds checking.
				for ( i = 0; i < ImageDimension; i++ )
				{
					if ( ( center_index[i] + static_cast< long > (
						m_NumberOfLayers ) >= static_cast< long>( upperBounds[i] - 1 ) )
						|| center_index[i] - static_cast< long >(
						m_NumberOfLayers ) <= static_cast< long >(lowerBounds[i]) )
					{
						m_BoundsCheckingActive = true;
					}
				}

				// Borrow a node from the store and set its value.
				node = sparsePtr->m_LayerNodeStore->Borrow();
				node->m_Value = center_index;

				// Add the node to the active list and set the status in the status
				// image.
				sparsePtr->m_Layers[0]->PushFront ( node );
				statusIt.SetCenterPixel ( 0 );

				// Grab the neighborhood in the image of shifted input values.
				shiftedIt.SetLocation ( center_index );

				// Search the neighborhood pixels for first inside & outside layer
				// members.  Construct these lists and set status list values.
				for ( i = 0; i < sparsePtr->m_NeighborList.GetSize(); ++i )
				{
					offset_index = center_index
						+ sparsePtr->m_NeighborList.GetNeighborhoodOffset ( i );

					if ( outputIt.GetPixel(
						sparsePtr->m_NeighborList.GetArrayIndex( i ) ) != m_ValueZero )
					{
						value = shiftedIt.GetPixel (
							sparsePtr->m_NeighborList.GetArrayIndex ( i ) );

						if ( value < m_ValueZero ) // Assign to first inside layer.
						{
							layer_number = 1;
						}
						else // Assign to first outside layer
						{
							layer_number = 2;
						}

						statusIt.SetPixel ( sparsePtr->m_NeighborList.GetArrayIndex ( i ),
							layer_number, bounds_status );
						if ( bounds_status == true ) // In bounds.
						{
							node = sparsePtr->m_LayerNodeStore->Borrow();
							node->m_Value = offset_index;
							sparsePtr->m_Layers[layer_number]->PushFront ( node );
						} // else do nothing.
					}
				}
			}
		}
	}
}

template<class TInputImage, class TOutputImage >
void
MultiphaseSparseFieldLevelSetImageFilter< TInputImage, TOutputImage >
::ConstructLayer ( SparseDataStruct *sparsePtr, StatusType from, StatusType to )
{
	LayerNodeType *node;
	bool boundary_status;
	typename LayerType::ConstIterator fromIt;
	NeighborhoodIterator<StatusImageType>
	statusIt ( sparsePtr->m_NeighborList.GetRadius(), sparsePtr->m_StatusImage,
		this->m_LevelSet[sparsePtr->m_Index]->GetRequestedRegion() );

	// For all indices in the "from" layer...
	for ( fromIt = sparsePtr->m_Layers[from]->Begin();
		fromIt != sparsePtr->m_Layers[from]->End();  ++fromIt )
	{
		// Search the neighborhood of this index in the status image for
		// unassigned indicies. Push those indicies onto the "to" layer and
		// assign them values in the status image.  Status pixels outside the
		// boundary will be ignored.
		statusIt.SetLocation ( fromIt->m_Value );
		for ( unsigned int i = 0; i < sparsePtr->m_NeighborList.GetSize(); ++i )
		{
			if ( statusIt.GetPixel (
				sparsePtr->m_NeighborList.GetArrayIndex ( i ) ) == m_StatusNull )
			{
				statusIt.SetPixel ( sparsePtr->m_NeighborList.GetArrayIndex ( i ),
					to, boundary_status );
				if ( boundary_status == true ) // in bounds
				{
					node = sparsePtr->m_LayerNodeStore->Borrow();
					node->m_Value = statusIt.GetIndex()
						+ sparsePtr->m_NeighborList.GetNeighborhoodOffset ( i );
					sparsePtr->m_Layers[to]->PushFront ( node );
				}
			}
		}
	}
}

template <class TInputImage, class TOutputImage >
void
MultiphaseSparseFieldLevelSetImageFilter< TInputImage, TOutputImage >
::AllocateUpdateBuffer()
{
	for ( unsigned int i = 0; i < this->FunctionCount; i++ )
	{
		SparseDataStruct *sparsePtr = sparseData[i];

		// Preallocate the update buffer.  NOTE: There is currently no way to
		// downsize a std::vector. This means that the update buffer will grow
		// dynamically but not shrink.  In newer implementations there may be a
		// squeeze method which can do this.  Alternately, we can implement our
		// own strategy for downsizing.
		sparsePtr->m_UpdateBuffer.clear();
		sparsePtr->m_UpdateBuffer.reserve ( sparsePtr->m_Layers[0]->Size() );
	}
}

template<class TInputImage, class TOutputImage >
void
MultiphaseSparseFieldLevelSetImageFilter<TInputImage, TOutputImage >
::PostProcessOutput()
{
  OutputImagePointer output = this->GetOutput();

  // Set the values in the levelset image for the active layer.
  this->InitializeActiveLayerValues();
  // Initialize layer values using the active layer as seeds.
  this->PropagateAllLayerValues();

  for ( unsigned int i = 0; i < this->FunctionCount; i++ )
  {
    InputImagePointer input = this->m_LevelSet[i];
    InputPointType origin = input->GetOrigin();
    InputSpacingType spacing = input->GetSpacing();

    // Local iterator
    ImageRegionIterator<OutputImageType> outputIt ( this->m_LevelSet[i],
      this->m_LevelSet[i]->GetRequestedRegion() );

    // In the context of the global coordinates
    OutputIndexType start;
    for ( unsigned int j = 0; j < ImageDimension; j++ )
      start[j] = static_cast<OutputIndexValueType>( origin[j]/spacing[j] );

    // Defining sub-region in the global coordinates
    OutputRegionType region;
    region.SetSize( input->GetLargestPossibleRegion().GetSize() );
    region.SetIndex( start );

    if ( !input || !output )
      itkExceptionMacro ( << "Either input and/or output is NULL." );

    ImageRegionIterator< OutputImageType > outIt ( output, region );

    for ( outputIt = outputIt.Begin(), outIt = outIt.Begin(); !outIt.IsAtEnd();
			++outputIt, ++outIt)
      if ( outputIt.Get() < 0 )
        outIt.Value() =  static_cast<OutputPixelType> ( this->m_Lookup[i] );
  }
}

template<class TInputImage, class TOutputImage >
void
MultiphaseSparseFieldLevelSetImageFilter<TInputImage, TOutputImage >
::PrintSelf ( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf ( os, indent );

  os << indent << "m_IsoSurfaceValue: " << m_IsoSurfaceValue << std::endl;
  os << indent << "m_BoundsCheckingActive: " << m_BoundsCheckingActive;

  for ( unsigned int i = 0; i < this->FunctionCount; i++ )
  {
    SparseDataStruct *sparsePtr = sparseData[i];
    os << indent << "m_LayerNodeStore: " << std::endl;
    sparsePtr->m_LayerNodeStore->Print ( os,indent.GetNextIndent() );
    for ( i = 0; i < sparsePtr->m_Layers.size(); i++ )
    {
      os << indent << "m_Layers[" << i << "]: size="
        << sparsePtr->m_Layers[i]->Size() << std::endl;
      os << indent << sparsePtr->m_Layers[i];
    }
    os << indent << "m_UpdateBuffer: size=" <<
      static_cast< unsigned long > ( sparsePtr->m_UpdateBuffer.size() )
      << " capacity = " <<
      static_cast<unsigned long> ( sparsePtr->m_UpdateBuffer.capacity() ) <<
      std::endl;
    }
  }

} // end namespace itk

#endif
