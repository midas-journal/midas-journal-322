#ifndef __itkMultiphaseFiniteDifferenceImageFilter_h_
#define __itkMultiphaseFiniteDifferenceImageFilter_h_

#include "itkArray.h"
#include "itkInPlaceImageFilter.h"
#include "itkFiniteDifferenceFunction.h"
#include <vnl/vnl_vector.h>
#include "itkImageRegionIterator.h"

#include "itkVector.h"
#include "itkListSample.h"
#include "itkKdTree.h"
#include "itkKdTreeGenerator.h"

namespace itk {

/**
 * \class FiniteDifferenceImageFilter2
 *
 * \par The Finite Difference Solver Hierarchy
 *
 * This is an alternate version of the ITK finite difference solver (FDS)
framework,
 * supporting the solution of multiple functions, simultaneously.  The
 * FDS framework is a set of classes for creating filters to solve partial
 * differential equations on images using an iterative, finite difference
 * update scheme.
 *
 * \par
 * The high-level algorithm implemented by the framework can be described by
 * the following pseudocode.
 *
 * \code
 *  WHILE NOT convergence:
 *     FOR ALL pixels i
 *		  FOR ALL functions f
 *			min_time_step = min(min_time_step, calculate_change(f, i))
 *		  FOR ALL functions f
 *          update(f, i, time_step)
 * \endcode
 *
 * \par
 * The following equation describes update \f$n+1\f$ at pixel \f$i\f$ on
 * discrete image \f$ u \f$ :
 *
 * \par
 * \f$u_{\mathbf{i}}^{n+1}=u^n_{\mathbf{i}}+\Delta u^n_{\mathbf{i}}\Delta t\f$
 *
 * \par Component objects
 * The FDS hierarchy is comprised of two component object types, variations of
 * which are designed to be plugged together to create filters for different
 * applications.  At the process level are the ``solver'' objects, which are
 * subclasses of FiniteDifferenceImageFilter2.  Solver objects are filters that
 * take image inputs and produce image outputs.  Solver objects require a
 * ``finite difference function'' object to perform the calculation at each
 * image pixel during iteration.  These specialized function objects are
 * subclasses of FiniteDifferenceFunction. FiniteDifferenceFunctions take a
 * neighborhood of pixels as input (in the form of an
 * itk::NeighborhoodIterator) and produce a scalar valued result.
 *
 * \par
 * Filters for different applications are created by defining a function object
 * to handle the numerical calculations and choosing (or creating) a solver
 * object that reflects the requirements and constraints of the application.
 * For example, anisotropic diffusion filters are created by plugging
 * anisotropic diffusion functions into the DenseFiniteDifferenceImageFilter2.
 * The separation between function object and solver object allows us to
 * create, for example, sparse-field, dense-field, and narrow-band
 * implementations of a level-set surface evolution filter can all be
 * constructed by plugging the same function object into three different,
 * specialized solvers.
 *
 * \par Creating new filters in this hierarchy
 * The procedure for creating a filter within the FDS hierarchy is to identify
 * all the virtual methods that need to be defined for your particular
 * application.  In the simplest case, a filter needs only to instantiate a
 * specific function object and define some halting criteria.  For more
 * complicated applications, you may need to define a specialized type of
 * iteration scheme or updating procedure in a higher-level solver object.
 *
 * \par
 * Some simple examples are the specific subclasses of
 * AnisotropicDiffusionImageFilter.  The leaves of the anisotropic diffusion
 * filter tree only define the function object they use for their particular
 * flavor of diffusion.  See CurvatureAnisotropicDiffusionImageFilter and
 * GradientAnisotropicDiffusionImageFilter for details.
 *
 * \par FiniteDifferenceImageFilter2
 * This class defines the generic solver API at the top level of the FDS
 * framework. FiniteDifferenceImageFilter2 is an abstract class that implements
 * the generic, high-level algorithm (described above).
 *
 * \par Inputs and Outputs
 * This filter is an Image to Image filter.  Depending on the specific
 * subclass implementation, finite difference image filters may process a
 * variety of image types.  The input to the filter is the initial
 * value of \f$ u \f$ and the output of the filter is the solution to the
 * p.d.e.
 *
 * \par How to use this class
 * GenerateData() relies on several virtual methods that must be defined by a
 * subclass.  Specifically: \em AllocateUpdateBuffer \em ApplyUpdate
 * \em CalculateChange and \em Halt.  To create a finite difference solver,
 * implement a subclass to define these methods.
 *
 * \par
 * Note that there is no fixed container type for the buffer used to hold
 * the update \f$ \Delta \f$.  The container might be another image, or simply
 * a list of values.  AllocateUpdateBuffer is responsible for creating the
 * \f$ \Delta \f$ container.  CalculateChange populates this buffer and
 * ApplyUpdate adds the buffer values to the output image (solution).  The
 * boolean Halt() (or ThreadedHalt) method returns a true value to stop
iteration.
 *
 * \ingroup ImageFilter
 * \ingroup LevelSetSegmentation
 * \sa DenseFiniteDifferenceImageFilter2 */
 template < class TInputImage, class TOutputImage >
class ITK_EXPORT MultiphaseFiniteDifferenceImageFilter
  : public InPlaceImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef MultiphaseFiniteDifferenceImageFilter Self;
  typedef InPlaceImageFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Run-time type information (and related methods) */
  itkTypeMacro( MultiphaseFiniteDifferenceImageFilter, InPlaceImageFilter );

  /** Input and output image types. */
  typedef TInputImage  InputImageType;
  typedef typename InputImageType::Pointer InputImagePointer;
  typedef typename InputImageType::RegionType InputRegionType;
  typedef typename InputImageType::SizeType InputSizeType;
  typedef typename InputImageType::SpacingType InputSpacingType;
  typedef typename InputImageType::PointType InputPointType;

  typedef TOutputImage OutputImageType;
  typedef typename OutputImageType::Pointer OutputImagePointer;
  typedef typename OutputImageType::PixelType PixelType;
  typedef typename OutputImageType::RegionType OutputRegionType;
  typedef typename OutputImageType::SizeType OutputSizeType;
  typedef typename OutputImageType::SizeValueType OutputSizeValueType;
  typedef typename OutputImageType::IndexType OutputIndexType;
  typedef typename OutputImageType::IndexValueType OutputIndexValueType;

  /** Dimensionality of input and output data is assumed to be the same. */
  itkStaticConstMacro(ImageDimension, unsigned int,
    OutputImageType::ImageDimension);

  /** The value type of the time step.  This is distinct from PixelType
   * because PixelType may often be a vector value, while the TimeStep is
   * a scalar value. */
  typedef FiniteDifferenceFunction<TOutputImage> FiniteDifferenceFunctionType;
  typedef typename FiniteDifferenceFunctionType::Pointer
    FiniteDifferenceFunctionTypePtr;
  typedef typename FiniteDifferenceFunctionType::TimeStepType TimeStepType;

  typedef vnl_vector< unsigned int > VectorType;

  typedef Vector< float, ImageDimension > CentroidVectorType;
  typedef itk::Statistics::ListSample< CentroidVectorType > SampleType;
  typedef itk::Statistics::KdTreeGenerator< SampleType > TreeGeneratorType;
  typedef typename TreeGeneratorType::Pointer TreePointer;
  typedef typename TreeGeneratorType::KdTreeType TreeType;
  typedef typename TreeType::Pointer KdTreePointer;

  typedef enum { UNINITIALIZED = 0, INITIALIZED = 1 } FilterStateType;

  /** This method returns a pointer to a FiniteDifferenceFunction object that
   * will be used by the filter to calculate updates at image pixels.
   * \param functionIndex Index of difference function to return.
   * \returns A FiniteDifferenceObject pointer. */
  virtual const FiniteDifferenceFunctionTypePtr &GetDifferenceFunction(
    unsigned int functionIndex ) const
  {
  	return ( *m_DifferenceFunctions[functionIndex] );
  }

  /** This method sets the pointer to a FiniteDifferenceFunction object that
   * will be used by the filter to calculate updates at image pixels.
   * \param functionIndex Index of difference function to set.
   * \param function Pointer to difference function to set. */
  virtual void SetDifferenceFunction(unsigned int functionIndex,
    FiniteDifferenceFunctionTypePtr function)
  {
    (*m_DifferenceFunctions[functionIndex]) = function;
  }

 /** Set/Get the number of iterations that the filter will run. */
  itkSetMacro(NumberOfIterations, unsigned int);
  itkGetConstReferenceMacro(NumberOfIterations, unsigned int);

  /** Use the image spacing information in calculations. Use this option if you
   *  want derivatives in physical space. Default is UseImageSpacingOff. */
  itkSetMacro(UseImageSpacing,bool);
  itkBooleanMacro(UseImageSpacing);
  itkGetConstReferenceMacro(UseImageSpacing, bool);

  /** Set/Get the maximum error allowed in the solution.  This may not be
      defined for all solvers and its meaning may change with the application.
*/
  itkSetMacro(MaximumRMSError, double);
  itkGetConstReferenceMacro(MaximumRMSError, double);

  /** Set/Get the root mean squared change of the previous iteration. May not
      be used by all solvers. */
  itkSetMacro(RMSChange, double);
  itkGetConstReferenceMacro(RMSChange, double);

  /** Set the state of the filter to INITIALIZED */
  void SetStateToInitialized()
  {
    this->SetState( INITIALIZED );
  }

  /** Set the state of the filter to UNINITIALIZED */
  void SetStateToUninitialized()
  {
    this->SetState( UNINITIALIZED );
  }

  /** Set/Get the state of the filter. */
  itkSetMacro( State, FilterStateType );
  itkGetConstReferenceMacro( State, FilterStateType );

  /** Require the filter to be manually reinitialized (by calling
      SetStateToUninitialized() */
  itkSetMacro( ManualReinitialization, bool );
  itkGetConstReferenceMacro( ManualReinitialization, bool );
  itkBooleanMacro( ManualReinitialization );

  /** Set the number of elapsed iterations of the filter. */
  itkSetMacro( ElapsedIterations, unsigned int );
  /** Get the number of elapsed iterations of the filter. */
  itkGetConstReferenceMacro( ElapsedIterations, unsigned int );

  void SetLevelSet( unsigned int i, InputImagePointer I )
  {
    m_LevelSet[i] = InputImageType::New();
    m_LevelSet[i]->SetRequestedRegion( I->GetRequestedRegion() );
    m_LevelSet[i]->SetBufferedRegion( I->GetBufferedRegion() );
    m_LevelSet[i]->SetLargestPossibleRegion( I->GetLargestPossibleRegion() );
    m_LevelSet[i]->Allocate();
    m_LevelSet[i]->CopyInformation( I );

    ImageRegionIterator< InputImageType > in ( I,
      I->GetLargestPossibleRegion());
    ImageRegionIterator< InputImageType > cp ( m_LevelSet[i],
      I->GetLargestPossibleRegion()  );
    for ( in.Begin(), cp.Begin(); !in.IsAtEnd(); ++in, ++cp )
      cp.Set( in.Get() );
  }

  InputImagePointer GetLevelSet( unsigned int i )
  {
    return m_LevelSet[i];
  }

  void SetLookup ( VectorType lookup )
  {
    m_Lookup = lookup;
  }

  void SetKdTree( KdTreePointer kdtree )
  {
    m_KdTree = kdtree;
  }

protected:
	MultiphaseFiniteDifferenceImageFilter()
  {
    m_UseImageSpacing    = false;
    m_ElapsedIterations  = 0;
    m_NumberOfIterations = NumericTraits<unsigned int>::max();
    m_MaximumRMSError = 0.0;
    m_RMSChange = 100.0;
    m_State = UNINITIALIZED;
    m_ManualReinitialization = false;
    m_KdTree = 0;
    this->InPlaceOff();
  }

  ~MultiphaseFiniteDifferenceImageFilter()
  {
	  delete [] m_DifferenceFunctions;
    delete [] m_LevelSet;
  }

  void SetFunctionCount( unsigned int n )
  {
    FunctionCount = n;
    m_DifferenceFunctions =
      new FiniteDifferenceFunctionTypePtr*[ FunctionCount ];

    // Initialize the images
    m_LevelSet = new InputImagePointer [ FunctionCount ];

    // Initialize the lookup table
    m_Lookup.set_size( FunctionCount );
    for( unsigned int i = 0; i < n; i++ )
      m_Lookup[i] = i+1;
  }

  unsigned int FunctionCount;
  InputImagePointer *m_LevelSet;
  VectorType m_Lookup;
  KdTreePointer m_KdTree;

  unsigned int  m_NumberOfIterations;
  double m_RMSChange;
  double m_MaximumRMSError;
  unsigned int m_ElapsedIterations;
    /** The function that will be used in calculating updates for each pixel. */
  FiniteDifferenceFunctionTypePtr **m_DifferenceFunctions;

  void PrintSelf(std::ostream& os, Indent indent) const;

  /** This method allocates a temporary update container in the subclass. */
  virtual void AllocateUpdateBuffer() = 0;

  /** This method is defined by a subclass to apply changes to the output
   * from an update buffer and a time step value "dt".
   * \param dt Time step value. */
  virtual void ApplyUpdate(TimeStepType dt) = 0;

  /** This method is defined by a subclass to populate an update buffer
   * with changes for the pixels in the output.  It returns a time
   * step value to be used for the update.
   * \returns A time step to use in updating the output with the changes
   * calculated from this method. */
  virtual TimeStepType CalculateChange() = 0;

  /** This method can be defined in subclasses as needed to copy the input
   * to the output. See DenseFiniteDifferenceImageFilter2 for an
   * implementation. */
  virtual void CopyInputToOutput() = 0;

  /** This is the default, high-level algorithm for calculating finite
   * difference solutions.  It calls virtual methods in its subclasses
   * to implement the major steps of the algorithm. */
  virtual void GenerateData();

  /** FiniteDifferenceImageFilter2 needs a larger input requested region than
   * the output requested region.  As such, we need to provide
   * an implementation for GenerateInputRequestedRegion() in order to inform
   * the pipeline execution model.
   *
   * \par
   * The filter will ask for a padded region to perform its neighborhood
   * calculations.  If no such region is available, the boundaries will be
   * handled as described in the FiniteDifferenceFunction defined by the
   * subclass.
   * \sa ProcessObject::GenerateInputRequestedRegion() */
  virtual void GenerateInputRequestedRegion();

  /** This method returns true when the current iterative solution of the
   * equation has met the criteria to stop solving.  Defined by a subclass. */
  virtual bool Halt();

  /** This method is similar to Halt(), and its default implementation in this
   * class is simply to call Halt(). However, this method takes as a parameter
   * a void pointer to the MultiThreader::ThreadInfoStruct structure. If you
   * override this method instead of overriding Halt, you will be able to get
   * the current thread ID and handle the Halt method accordingly. This is
useful
   * if you are doing a lot of processing in Halt that you don't want
parallelized.
   * Notice that ThreadedHalt is only called by the multithreaded filters, so
you
   * still should implement Halt, just in case a non-threaded filter is used.
   */
  virtual bool ThreadedHalt(void *itkNotUsed(threadInfo))
  {
    return this->Halt();
  }

  /** This method is optionally defined by a subclass and is called before
   * the loop of iterations of calculate_change & upate. It does the global
   * initialization, i.e. in the SparseFieldLevelSetImageFilter, initialize
   * the list of layers.
   * */
  virtual void Initialize() { };

  /** This method is optionally defined by a subclass and is called immediately
   * prior to each iterative CalculateChange-ApplyUpdate cycle.  It can be
   * used to set global variables needed for the next iteration (ie. average
   * gradient magnitude of the image in anisotropic diffusion functions), or
   * otherwise prepare for the next iteration.
   * */
  virtual void InitializeIteration()
  {
    for(unsigned int i = 0; i < this->FunctionCount; i++ )
      GetDifferenceFunction( i )->InitializeIteration();
  }

  /** Virtual method for resolving a single time step from a set of time steps
   * returned from processing threads.
   * \return Time step (dt) for the iteration update based on a list
   * of time steps generated from the threaded calculated change method (one
   * for each region processed).
   *
   * \param timeStepList The set of time changes compiled from all the threaded
calls
   * to ThreadedGenerateData.
   * \param valid The set of flags indicating which of "list" elements are
   *  valid
   * \param size The size of "list" and "valid"
   *
   * The default is to return the minimum value in the list. */
  virtual TimeStepType ResolveTimeStep(const TimeStepType* timeStepList,
    const bool* valid,int size );


  /** This method is called after the solution has been generated to allow
   * subclasses to apply some further processing to the output.*/
  virtual void PostProcessOutput() {}


private:
  MultiphaseFiniteDifferenceImageFilter(const Self&);
  //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** Control whether derivatives use spacing of the input image in
      its calculation. */
  bool m_UseImageSpacing;

  /** Indicates whether the filter automatically resets to UNINITIALIZED state
      after completing, or whether filter must be manually reset */
  bool m_ManualReinitialization;


  /** State that the filter is in, i.e. UNINITIALIZED or INITIALIZED */
  FilterStateType m_State;
};

}// end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMultiphaseFiniteDifferenceImageFilter.txx"
#endif

#endif

