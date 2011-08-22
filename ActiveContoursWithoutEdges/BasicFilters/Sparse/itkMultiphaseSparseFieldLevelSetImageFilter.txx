#ifndef __itkMultiphaseSparseFieldLevelSetImageFilter_txx_
#define __itkMultiphaseSparseFieldLevelSetImageFilter_txx_
#include "itkImageFileWriter.h"

namespace itk
{

template < class TInputImage, class TOutputImage >
typename MultiphaseSparseFieldLevelSetImageFilter< TInputImage,
TOutputImage >::TimeStepType
MultiphaseSparseFieldLevelSetImageFilter< TInputImage, TOutputImage >
::CalculateChange()
{
	TimeStepType minTimeStep = 999999;

  for ( unsigned int i = 0; i < this->FunctionCount; i++ )
  {
		currentFunction = i;

		const typename Superclass::FiniteDifferenceFunctionType::Pointer df =
			this->GetDifferenceFunction ( currentFunction );

		SparseDataStruct *sparsePtr = sparseData[currentFunction];

		typename Superclass::FiniteDifferenceFunctionType::FloatOffsetType offset;
		ValueType norm_grad_phi_squared, dx_forward, dx_backward, forwardValue,
      backwardValue, centerValue;
		const ValueType MIN_NORM      = 1.0e-6;
		void *globalData = df->GetGlobalDataPointer();

		typename LayerType::ConstIterator layerIt;
		NeighborhoodIterator<OutputImageType> outputIt ( df->GetRadius(),
			this->m_LevelSet[currentFunction],
			this->m_LevelSet[currentFunction]->GetRequestedRegion() );
		TimeStepType timeStep;

		if ( m_BoundsCheckingActive == false )
			outputIt.NeedToUseBoundaryConditionOff();

		sparsePtr->m_UpdateBuffer.clear();
		sparsePtr->m_UpdateBuffer.reserve ( sparsePtr->m_Layers[0]->Size() );

		// Calculates the update values for the active layer indicies in this
		// iteration.  Iterates through the active layer index list, applying
		// the level set function to the output image (level set image) at each
		// index.  Update values are stored in the update buffer.
		for ( layerIt = sparsePtr->m_Layers[0]->Begin(); layerIt !=
			sparsePtr->m_Layers[0]->End(); ++layerIt )
		{
			outputIt.SetLocation ( layerIt->m_Value );

			// Calculate the offset to the surface from the center of this
			// neighborhood.  This is used by some level set functions in sampling a
			// speed, advection, or curvature term.
			if ( this->GetInterpolateSurfaceLocation()
				&& ( centerValue = outputIt.GetCenterPixel() ) != 0.0 )
			{
				// Surface is at the zero crossing, so distance to surface is:
				// phi(x) / norm(grad(phi)), where phi(x) is the center of the
				// neighborhood.  The location is therefore
				// (i,j,k) - ( phi(x) * grad(phi(x)) ) / norm(grad(phi))^2
				norm_grad_phi_squared = 0.0;
				for ( unsigned int j = 0; j < ImageDimension; ++j )
				{
					forwardValue  = outputIt.GetNext ( j );
					backwardValue = outputIt.GetPrevious ( j );

					if ( forwardValue * backwardValue >= 0 )
					{
						//  Neighbors are same sign OR at least one neighbor is zero.
						dx_forward  = forwardValue - centerValue;
						dx_backward = centerValue - backwardValue;

						// Pick the larger magnitude derivative.
						if ( ::vnl_math_abs ( dx_forward ) >::vnl_math_abs( dx_backward ) )
							offset[j] = dx_forward;
						else
							offset[j] = dx_backward;
					}
					else //Neighbors are opposite sign, pick the direction of 0 surface.
					{
						if ( forwardValue * centerValue < 0 )
							offset[j] = forwardValue - centerValue;
						else
							offset[j] = centerValue - backwardValue;
					}

					norm_grad_phi_squared += offset[j] * offset[j];
				}

				for ( unsigned int j = 0; j < ImageDimension; ++j )
				{
            // Adding sqrt imagedimension "extends the reach" of the
            // interpolation
            // to surfaces that pass close to the center of cells.  This is a
            // heuristic fudge factor that improves interpolation and reduces
            // "wiggling" at convergence.
            offset[j] = ( offset[j] * centerValue ) * vcl_sqrt( ImageDimension
              + 0.5 )/ ( norm_grad_phi_squared + MIN_NORM );
				}

				sparsePtr->m_UpdateBuffer.push_back ( df->ComputeUpdate ( outputIt,
					globalData, offset ) );
			}
			else // Don't do interpolation
				sparsePtr->m_UpdateBuffer.push_back ( df->ComputeUpdate ( outputIt,
					globalData ) );
		}

		// Ask the finite difference function to compute the time step for
		// this iteration.  We give it the global data pointer to use, then
		// ask it to free the global data memory.
		timeStep = df->ComputeGlobalTimeStep ( globalData );
		df->ReleaseGlobalDataPointer ( globalData );
		if ( timeStep < minTimeStep )
			minTimeStep = timeStep;
	}
	minTimeStep = 0.2;
	return minTimeStep;
}


template<class TInputImage, class TOutputImage >
void
MultiphaseSparseFieldLevelSetImageFilter<TInputImage, TOutputImage >
::ApplyUpdate ( TimeStepType dt )
{
	for ( unsigned int i = 0; i < this->FunctionCount; i++ )
	{
		 currentFunction = i;

		SparseDataStruct *sparsePtr = sparseData[currentFunction];
		unsigned int t;

		StatusType up_to, up_search;
		StatusType down_to, down_search;

		LayerPointerType UpList[2];
		LayerPointerType DownList[2];
		for ( unsigned int j = 0; j < 2; ++j )
		{
			UpList[j]   = LayerType::New();
			DownList[j] = LayerType::New();
		}

		// Process the active layer.  This step will update the values in the
		//active
		// layer as well as the values at indices that *will* become part of the
		// active layer when they are promoted/demoted.  Also records promotions,
		// demotions in the m_StatusLayer for current active layer indices
		// (i.e. those indices which will move inside or outside the active
		// layers).
		this->UpdateActiveLayerValues( dt, UpList[0], DownList[0] );

		// Process the status up/down lists.  This is an iterative process which
		// proceeds outwards from the active layer.  Each iteration generates the
		// list for the next iteration.

		// First process the status lists generated on the active layer.
		this->ProcessStatusList ( UpList[0], UpList[1], 2, 1 );
		this->ProcessStatusList ( DownList[0], DownList[1], 1, 2 );

		down_to = up_to = 0;
		up_search       = 3;
		down_search     = 4;
		unsigned int j = 1;
		unsigned int k = 0;
		while ( down_search < static_cast<StatusType>(sparsePtr->m_Layers.size()))
		{
        this->ProcessStatusList(UpList[j], UpList[k], up_to, up_search );
        this->ProcessStatusList(DownList[j], DownList[k], down_to, down_search);

        if ( up_to == 0 ) up_to += 1;
        else            up_to += 2;
        down_to += 2;

        up_search   += 2;
        down_search += 2;

        // Swap the lists so we can re-use the empty one.
        t = j;
        j = k;
        k = t;
		}

		// Process the outermost inside/outside layers in the sparse field.
		this->ProcessStatusList(UpList[j], UpList[k], up_to, m_StatusNull );
		this->ProcessStatusList(DownList[j], DownList[k], down_to, m_StatusNull);

		// Now we are left with the lists of indicies which must be
		// brought into the outermost layers.  Bring UpList into last inside layer
		// and DownList into last outside layer.
		this->ProcessOutsideList ( UpList[k], static_cast<int> (
			sparsePtr->m_Layers.size() ) -2 );
		this->ProcessOutsideList ( DownList[k], static_cast<int> (
			sparsePtr->m_Layers.size() ) -1 );

		// Finally, we update all of the layer values (excluding the active layer,
		// which has already been updated).
		this->PropagateAllLayerValues();
	}
	currentFunction = 0;
}

template < class TInputImage, class TOutputImage >
void
MultiphaseSparseFieldLevelSetImageFilter<TInputImage, TOutputImage >
::ProcessOutsideList ( LayerType *OutsideList, StatusType ChangeToStatus )
{
    SparseDataStruct *sparsePtr = sparseData[currentFunction];
    LayerNodeType *node;

    // Push each index in the input list into its appropriate status layer
    // (ChangeToStatus) and update the status image value at that index.
    while ( ! OutsideList->Empty() )
    {
      sparsePtr->m_StatusImage->SetPixel ( OutsideList->Front()->m_Value,
                                           ChangeToStatus );
      node = OutsideList->Front();
      OutsideList->PopFront();
      sparsePtr->m_Layers[ChangeToStatus]->PushFront ( node );
    }
}

template <class TInputImage, class TOutputImage >
void
MultiphaseSparseFieldLevelSetImageFilter< TInputImage, TOutputImage >
::ProcessStatusList ( LayerType *InputList, LayerType *OutputList,
                      StatusType ChangeToStatus, StatusType SearchForStatus )
{
    SparseDataStruct *sparsePtr = sparseData[currentFunction];

    unsigned int i;
    bool bounds_status;
    LayerNodeType *node;
    StatusType neighbor_status;
    NeighborhoodIterator<StatusImageType>
    statusIt ( sparsePtr->m_NeighborList.GetRadius(), sparsePtr->m_StatusImage,
               this->m_LevelSet[currentFunction]->GetRequestedRegion() );

    if ( m_BoundsCheckingActive == false )
    {
      statusIt.NeedToUseBoundaryConditionOff();
    }

    // Push each index in the input list into its appropriate status layer
    // (ChangeToStatus) and update the status image value at that index.
    // Also examine the neighbors of the index to determine which need to go
    // onto
    // the output list (search for SearchForStatus).
    while ( ! InputList->Empty() )
    {
      statusIt.SetLocation ( InputList->Front()->m_Value );
      statusIt.SetCenterPixel ( ChangeToStatus );

      node = InputList->Front();  // Must unlink from the input list
      InputList->PopFront();      // _before_ transferring to another list.
      sparsePtr->m_Layers[ChangeToStatus]->PushFront ( node );

      for ( i = 0; i < sparsePtr->m_NeighborList.GetSize(); ++i )
      {
        neighbor_status = statusIt.GetPixel (
          sparsePtr->m_NeighborList.GetArrayIndex ( i ) );

        // Have we bumped up against the boundary?  If so, turn on bounds
        // checking.
        if ( neighbor_status == m_StatusBoundaryPixel )
        {
          m_BoundsCheckingActive = true;
        }

        if ( neighbor_status == SearchForStatus )
        {
          // mark this pixel so we don't add it twice.
          statusIt.SetPixel ( sparsePtr->m_NeighborList.GetArrayIndex ( i ),
                              m_StatusChanging, bounds_status );
          if ( bounds_status == true )
          {
            node = sparsePtr->m_LayerNodeStore->Borrow();
            node->m_Value = statusIt.GetIndex() +
              sparsePtr->m_NeighborList.GetNeighborhoodOffset ( i );
            OutputList->PushFront ( node );
          } // else this index was out of bounds.
        }
      }
    }
}

template <class TInputImage, class TOutputImage >
void
MultiphaseSparseFieldLevelSetImageFilter<TInputImage, TOutputImage >
::UpdateActiveLayerValues ( TimeStepType dt, LayerType *UpList, LayerType
*DownList )
{
	SparseDataStruct *sparsePtr = sparseData[currentFunction];

	// This method scales the update buffer values by the time step and adds
	// them to the active layer pixels.  New values at an index which fall
	// outside of the active layer range trigger that index to be placed on the
	// "up" or "down" status list.  The neighbors of any such index are then
	// assigned new values if they are determined to be part of the active list
	// for the next iteration (i.e. their values will be raised or lowered into
	// the active range).
	const ValueType LOWER_ACTIVE_THRESHOLD = -( m_ConstantGradientValue/2.0 );
	const ValueType UPPER_ACTIVE_THRESHOLD =  m_ConstantGradientValue / 2.0;

	ValueType new_value, temp_value, rms_change_accumulator;
	LayerNodeType *node, *release_node;
	StatusType neighbor_status;
	unsigned int i, idx, counter;
	bool bounds_status, flag;

	typename LayerType::Iterator layerIt;
	typename UpdateBufferType::const_iterator updateIt;

	NeighborhoodIterator< OutputImageType >
		outputIt ( sparsePtr->m_NeighborList.GetRadius(),
    this->m_LevelSet[currentFunction],
    this->m_LevelSet[currentFunction]->GetRequestedRegion() );

	NeighborhoodIterator< StatusImageType >
    statusIt ( sparsePtr->m_NeighborList.GetRadius(),
		sparsePtr->m_StatusImage,
    this->m_LevelSet[currentFunction]->GetRequestedRegion() );

	if ( m_BoundsCheckingActive == false )
	{
		outputIt.NeedToUseBoundaryConditionOff();
		statusIt.NeedToUseBoundaryConditionOff();
	}

	counter = 0;
	rms_change_accumulator = m_ValueZero;
	layerIt = sparsePtr->m_Layers[0]->Begin();
	updateIt = sparsePtr->m_UpdateBuffer.begin();
	while ( layerIt != sparsePtr->m_Layers[0]->End() )
	{
		outputIt.SetLocation ( layerIt->m_Value );
		statusIt.SetLocation ( layerIt->m_Value );

		new_value = this->CalculateUpdateValue ( layerIt->m_Value,
			dt, outputIt.GetCenterPixel(), *updateIt );

      // If this index needs to be moved to another layer, then search its
      // neighborhood for indicies that need to be pulled up/down into the
      // active layer. Set those new active layer values appropriately,
      // checking first to make sure they have not been set by a more
      // influential neighbor.

      //   ...But first make sure any neighbors in the active layer are not
      // moving to a layer in the opposite direction.  This step is necessary
      // to avoid the creation of holes in the active layer.  The fix is simply
      // to not change this value and leave the index in the active set.

      if ( new_value >= UPPER_ACTIVE_THRESHOLD )
      {
        // This index will move UP into a positive (outside) layer.

        // First check for active layer neighbors moving in the opposite
        // direction.
        flag = false;
        for ( i = 0; i < sparsePtr->m_NeighborList.GetSize(); ++i )
        {
          if ( statusIt.GetPixel( sparsePtr->m_NeighborList.GetArrayIndex( i ) )
                  == m_StatusActiveChangingDown )
          {
            flag = true;
            break;
          }
        }
        if ( flag == true )
        {
          ++layerIt;
          ++updateIt;
          continue;
        }

        rms_change_accumulator += vnl_math_sqr (
          new_value-outputIt.GetCenterPixel() );

        // Search the neighborhood for inside indicies.
        temp_value = new_value - m_ConstantGradientValue;
        for ( i = 0; i < sparsePtr->m_NeighborList.GetSize(); ++i )
        {
          idx = sparsePtr->m_NeighborList.GetArrayIndex ( i );
          neighbor_status = statusIt.GetPixel ( idx );
          if ( neighbor_status == 1 )
          {
            // Keep the smallest possible value for the new active node.  This
            // places the new active layer node closest to the zero level-set.
            if ( outputIt.GetPixel ( idx ) < LOWER_ACTIVE_THRESHOLD ||
                ::vnl_math_abs ( temp_value ) < ::vnl_math_abs (
                outputIt.GetPixel ( idx ) ) )
            {
              UpdatePixel ( currentFunction, idx, outputIt, temp_value,
                bounds_status );
            }
          }
        }
        // Push it into the uplist
        node = sparsePtr->m_LayerNodeStore->Borrow();
        node->m_Value = layerIt->m_Value;
        UpList->PushFront ( node );
        statusIt.SetCenterPixel ( m_StatusActiveChangingUp );

        // Now remove this index from the active list.
        release_node = layerIt.GetPointer();
        ++layerIt;
        sparsePtr->m_Layers[0]->Unlink ( release_node );
        sparsePtr->m_LayerNodeStore->Return ( release_node );
      }

      else if ( new_value < LOWER_ACTIVE_THRESHOLD )
      {
        // This index will move DOWN into a negative (inside) layer.

        // First check for active layer neighbors moving in the opposite
        // direction.
        flag = false;
        for ( i = 0; i < sparsePtr->m_NeighborList.GetSize(); ++i )
        {
          if ( statusIt.GetPixel( sparsePtr->m_NeighborList.GetArrayIndex( i ) )
                  == m_StatusActiveChangingUp )
          {
            flag = true;
            break;
          }
        }
        if ( flag == true )
        {
          ++layerIt;
          ++updateIt;
          continue;
        }

        rms_change_accumulator += vnl_math_sqr ( new_value -
          outputIt.GetCenterPixel() );

        // Search the neighborhood for outside indicies.
        temp_value = new_value + m_ConstantGradientValue;
        for ( i = 0; i < sparsePtr->m_NeighborList.GetSize(); ++i )
        {
          idx = sparsePtr->m_NeighborList.GetArrayIndex ( i );
          neighbor_status = statusIt.GetPixel ( idx );
          if ( neighbor_status == 2 )
          {
            // Keep the smallest magnitude value for this active set node.  This
            // places the node closest to the active layer.
            if ( outputIt.GetPixel ( idx ) >= UPPER_ACTIVE_THRESHOLD ||
              ::vnl_math_abs ( temp_value ) < ::vnl_math_abs (
              outputIt.GetPixel ( idx ) ) )
            {
              UpdatePixel ( currentFunction, idx, outputIt, temp_value,
                bounds_status );
            }
          }
        }
        node = sparsePtr->m_LayerNodeStore->Borrow();
        node->m_Value = layerIt->m_Value;
        DownList->PushFront ( node );
        statusIt.SetCenterPixel ( m_StatusActiveChangingDown );

        // Now remove this index from the active list.
        release_node = layerIt.GetPointer();
        ++layerIt;
        sparsePtr->m_Layers[0]->Unlink ( release_node );
        sparsePtr->m_LayerNodeStore->Return ( release_node );
      }
      else
      {
        rms_change_accumulator += vnl_math_sqr ( new_value -
          outputIt.GetCenterPixel() );

        UpdatePixel ( currentFunction, outputIt.Size() >>1, outputIt, new_value,
          bounds_status );

        ++layerIt;
      }
      ++updateIt;
      ++counter;
  }

	// Determine the average change during this iteration.
	if ( counter == 0 )
	{
		this->SetRMSChange ( static_cast<double> ( m_ValueZero ) );
	}
	else
	{
		this->SetRMSChange ( static_cast<double> (
			vcl_sqrt ( ( double ) ( rms_change_accumulator /
			static_cast<ValueType> ( counter ) ) ) ) );
	}
}

template <class TInputImage, class TOutputImage >
void
MultiphaseSparseFieldLevelSetImageFilter< TInputImage, TOutputImage >
::InitializeActiveLayerValues()
  {
    const ValueType CHANGE_FACTOR = m_ConstantGradientValue / 2.0;
    //  const ValueType CHANGE_FACTOR = 0.7;
    const ValueType MIN_NORM      = 1.0e-6;

    for ( unsigned int i = 0; i < this->FunctionCount; i++ )
    {
      SparseDataStruct *sparsePtr = sparseData[i];

      unsigned int center;

      typename LayerType::ConstIterator activeIt;
      ConstNeighborhoodIterator<OutputImageType>
      outputIt ( sparsePtr->m_NeighborList.GetRadius(),
        this->m_LevelSet[i],
        this->m_LevelSet[i]->GetRequestedRegion() );

      sparsePtr->m_UpdateBuffer.clear();
      sparsePtr->m_UpdateBuffer.reserve ( sparsePtr->m_Layers[0]->Size() );

      center = outputIt.Size()/2;
      typename OutputImageType::Pointer output = this->m_LevelSet[i];

      ValueType dx_forward, dx_backward, length, distance;

      // For all indicies in the active layer...
      for ( activeIt = sparsePtr->m_Layers[0]->Begin();
        activeIt != sparsePtr->m_Layers[0]->End(); ++activeIt )
      {
        // Interpolate on the (shifted) input image values at this index to
        // assign an active layer value in the output image.
        outputIt.SetLocation ( activeIt->m_Value );

        length = m_ValueZero;
        for ( unsigned int j = 0; j < ImageDimension; ++j )
        {
          dx_forward = outputIt.GetPixel ( center +
            sparsePtr->m_NeighborList.GetStride( j ) ) -
            outputIt.GetCenterPixel();
          dx_backward = outputIt.GetCenterPixel() -
            outputIt.GetPixel ( center -
            sparsePtr->m_NeighborList.GetStride ( j ) );

          if ( vnl_math_abs ( dx_forward ) > vnl_math_abs ( dx_backward ) )
            length += dx_forward * dx_forward;
          else
            length += dx_backward * dx_backward;
        }
        length = vcl_sqrt ( ( double ) length ) + MIN_NORM;
        distance = outputIt.GetCenterPixel() / length;

        sparsePtr->m_UpdateBuffer.push_back(
          vnl_math_min ( vnl_math_max ( -CHANGE_FACTOR, distance ),
          CHANGE_FACTOR ) );
      }

      for ( activeIt = sparsePtr->m_Layers[0]->Begin();
        activeIt != sparsePtr->m_Layers[0]->End(); ++activeIt )
      {
        output->SetPixel ( activeIt->m_Value ,
          sparsePtr->m_UpdateBuffer.front() );
      }
    }
}


template <class TInputImage, class TOutputImage >
void
MultiphaseSparseFieldLevelSetImageFilter<TInputImage, TOutputImage >
::PropagateAllLayerValues()
{
    for ( unsigned int i = 0; i < this->FunctionCount; i++ )
      PropagateFunctionLayerValues ( i );
}

template <class TInputImage, class TOutputImage >
void
MultiphaseSparseFieldLevelSetImageFilter< TInputImage, TOutputImage >
::PropagateFunctionLayerValues ( unsigned int functionIndex )
{
    SparseDataStruct *sparsePtr = sparseData[functionIndex];

    unsigned int i;

    // Update values in the first inside and first outside layers using the
    // active layer as a seed. Inside layers are odd numbers, outside layers are
    // even numbers.
    this->PropagateLayerValues ( sparsePtr, 0, 1, 3, 1 ); // first inside
    this->PropagateLayerValues ( sparsePtr, 0, 2, 4, 2 ); // first outside

    // Update the rest of the layers.
    for ( i = 1; i < sparsePtr->m_Layers.size() - 2; ++i )
      this->PropagateLayerValues ( sparsePtr, i, i+2, i+4, ( i+2 ) %2 );
}

template <class TInputImage, class TOutputImage >
void
MultiphaseSparseFieldLevelSetImageFilter< TInputImage, TOutputImage >
::PropagateLayerValues ( SparseDataStruct *sparsePtr, StatusType from,
                        StatusType to, StatusType promote, int InOrOut )
{
    unsigned int i;
    ValueType value, value_temp, delta;
    value = NumericTraits<ValueType>::Zero; // warnings
    bool found_neighbor_flag;
    typename LayerType::Iterator toIt;
    LayerNodeType *node;
    StatusType past_end = static_cast<StatusType>(
      sparsePtr->m_Layers.size() )-1;

    // Are we propagating values inward (more negative) or outward (more
    // positive)?
    if ( InOrOut == 1 ) delta = - m_ConstantGradientValue;
    else delta = m_ConstantGradientValue;

    NeighborhoodIterator<OutputImageType>
    outputIt ( sparsePtr->m_NeighborList.GetRadius(),
      this->m_LevelSet[sparsePtr->m_Index],
      this->m_LevelSet[sparsePtr->m_Index]->GetRequestedRegion() );
    NeighborhoodIterator<StatusImageType>
    statusIt ( sparsePtr->m_NeighborList.GetRadius(), sparsePtr->m_StatusImage,
      this->m_LevelSet[sparsePtr->m_Index]->GetRequestedRegion() );

    if ( m_BoundsCheckingActive == false )
    {
      outputIt.NeedToUseBoundaryConditionOff();
      statusIt.NeedToUseBoundaryConditionOff();
    }

    toIt  = sparsePtr->m_Layers[to]->Begin();
    while ( toIt != sparsePtr->m_Layers[to]->End() )
    {
      statusIt.SetLocation ( toIt->m_Value );

      // Is this index marked for deletion? If the status image has
      // been marked with another layer's value, we need to delete this node
      // from the current list then skip to the next iteration.
      if ( statusIt.GetCenterPixel() != to )
      {
        node = toIt.GetPointer();
        ++toIt;
        sparsePtr->m_Layers[to]->Unlink ( node );
        sparsePtr->m_LayerNodeStore->Return ( node );
        continue;
      }

      outputIt.SetLocation ( toIt->m_Value );

      found_neighbor_flag = false;
      for ( i = 0; i < sparsePtr->m_NeighborList.GetSize(); ++i )
      {
        // If this neighbor is in the "from" list, compare its absolute value
        // to to any previous values found in the "from" list.  Keep the value
        // that will cause the next layer to be closest to the zero level set.
        if ( statusIt.GetPixel ( sparsePtr->m_NeighborList.GetArrayIndex ( i ) )
                == from )
        {
          value_temp = outputIt.GetPixel (
                           sparsePtr->m_NeighborList.GetArrayIndex ( i ) );

          if ( found_neighbor_flag == false )
          {
            value = value_temp;
          }
          else
          {
            if ( InOrOut == 1 )
            {
              // Find the largest (least negative) neighbor
              if ( value_temp > value )
              {
                value = value_temp;
              }
            }
            else
            {
              // Find the smallest (least positive) neighbor
              if ( value_temp < value )
              {
                value = value_temp;
              }
            }
          }
          found_neighbor_flag = true;
        }
      }
      if ( found_neighbor_flag == true )
      {
        // Set the new value using the smallest distance
        // found in our "from" neighbors.
        bool bounds_status;
        ValueType updateValue = value + delta;

        UpdatePixel ( sparsePtr->m_Index, outputIt.Size() >>1,
          outputIt, updateValue, bounds_status );

        ++toIt;
      }
      else
      {
        // Did not find any neighbors on the "from" list, then promote this
        // node.  A "promote" value past the end of my sparse field size
        // means delete the node instead.  Change the status value in the
        // status image accordingly.
        node  = toIt.GetPointer();
        ++toIt;
        sparsePtr->m_Layers[to]->Unlink ( node );
        if ( promote > past_end )
        {
          sparsePtr->m_LayerNodeStore->Return ( node );
          statusIt.SetCenterPixel ( m_StatusNull );
        }
        else
        {
          sparsePtr->m_Layers[promote]->PushFront ( node );
          statusIt.SetCenterPixel ( promote );
        }
      }
    }
}

} // end namespace itk

#endif
