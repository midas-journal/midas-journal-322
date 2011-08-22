#ifndef __itkDenseMultiphaseLevelSetImageFilter_txx
#define __itkDenseMultiphaseLevelSetImageFilter_txx

#include "itkDenseMultiphaseLevelSetImageFilter.h"

namespace itk
{
template < class TInputImage, class TFeatureImage, class TFunction,
class TOutputPixelType >
void
DenseMultiphaseLevelSetImageFilter< TInputImage, TFeatureImage,
  TFunction, TOutputPixelType >::
Initialize()
{
  sharedData->SetFunctionCount ( this->FunctionCount );

  if ( this->m_KdTree )
  {
    sharedData->SetKdTree( this->m_KdTree );
  }

  for ( unsigned int i = 0; i < this->FunctionCount; i++ )
  {
    typename FunctionType::Pointer typedPointer =
      GetTypedDifferenceFunction ( i );

    typedPointer->SetFunctionId ( i );

    sharedData->CreateHVals (
      i, this->m_LevelSet[i]->GetSpacing(),
      this->m_LevelSet[i]->GetOrigin(),
      this->m_LevelSet[i]->GetLargestPossibleRegion() );

    // Share the sharedData structure
    typedPointer->SetSharedData ( sharedData );
  }

  sharedData->AllocateListImage(
    this->GetFeatureImage()->GetLargestPossibleRegion(),
    this->GetFeatureImage()->GetSpacing() );
  std::cout << "Allocated list image " << std::endl;

  sharedData->PopulateListImage();
  std::cout << "Populated list image " << std::endl;

  Superclass::Initialize();

  for ( unsigned int i = 0; i < this->FunctionCount; i++ )
  {
    typename FunctionType::Pointer typedPointer =
      GetTypedDifferenceFunction ( i );
    typedPointer->SetInitialImage ( this->m_LevelSet[i] );
    typedPointer->UpdateSharedData ( true );
  }

  for ( unsigned int i = 0; i < this->FunctionCount; i++ )
  {
    typename FunctionType::Pointer typedPointer =
      GetTypedDifferenceFunction ( i );
    typedPointer->UpdateSharedData ( false );
  }

#ifdef VIZU
 func.SetFunctionCount( this->FunctionCount );
 func.SetFeatureImage( this->GetFeatureImage() );
#endif
}

/** Overrides parent implementation */
// This function is called at the end of each iteration
template < class TInputImage, class TFeatureImage, class TFunction,
  class TOutputPixelType >
void
DenseMultiphaseLevelSetImageFilter< TInputImage, TFeatureImage,
  TFunction, TOutputPixelType > ::
InitializeIteration()
{
  Superclass::InitializeIteration();

  for ( unsigned int i = 0; i < this->FunctionCount; i++ )
  {
    typename FunctionType::Pointer typedPointer =
      GetTypedDifferenceFunction ( i );
    typedPointer->UpdateSharedData ( true );
  }

  for ( unsigned int i = 0; i < this->FunctionCount; i++ )
  {
    typename FunctionType::Pointer typedPointer =
      GetTypedDifferenceFunction ( i );
    typedPointer->UpdateSharedData ( false );

    //Visualization routine
#ifdef VIZU
    func.SetInputImage( this->m_LevelSet[i] );
    func.Visualize( this->m_ElapsedIterations, 0, i );
#endif
  }

  // Estimate the progress of the filter
  this->SetProgress ( ( float ) ( ( float ) this->m_ElapsedIterations
    / ( float ) this->m_NumberOfIterations ) );
}

template < class TInputImage, class TFeatureImage, class TFunction,
  class TOutputPixelType >
void
DenseMultiphaseLevelSetImageFilter< TInputImage, TFeatureImage,
  TFunction, TOutputPixelType > ::
PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Class Name: " << GetNameOfClass( ) << std::endl;
}


} /* end namespace itk */

#endif
