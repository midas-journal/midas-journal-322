#include "itkImageFileReader.h"
#include "itkCastImageFilter.h"
#include "itkScalarChanAndVeseLevelSetFunction.h"
#include "itkSparseMultiphaseLevelSetImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkImageFileWriter.h"
#include "anyoption.h"

using namespace itk;

int main(int argc, char**argv)
{
/* 1. CREATE AN OBJECT */
  AnyOption *opt = new AnyOption();

  /* 2. SET PREFERENCES  */
  //opt->noPOSIX(); /* do not check for POSIX style character options */
  //opt->setVerbose(); /* print warnings about unknown options */
  //opt->autoUsagePrint(true); /* print usage for bad options */

  /* 3. SET THE USAGE/HELP   */
  opt->addUsage( "" );
  opt->addUsage( "Usage: " );
  opt->addUsage( "" );
  opt->addUsage( " Feature Image (input image) " );
  opt->addUsage( " Contour Image1 (level set initialization) " );
  opt->addUsage( " Contour Image2 (level set initialization) " );
  opt->addUsage( " Contour Image3 (level set initialization) " );
  opt->addUsage( " Output Image  (output image) " );
  opt->addUsage( " -h   --help                   Prints this help " );
  opt->addUsage( " -i   --iter    10  (default)  Number of Iterations" );
  opt->addUsage( " -r   --rms     0   (default)  rms" );
  opt->addUsage( "      --epsilon 1   (default)  Epsilon parameter in the Heaviside and dirac" );
  opt->addUsage( "      --mu      0   (default)  Length Regularization weight" );
  opt->addUsage( "      --nu      0   (default)  Area Regularization weight" );
  opt->addUsage( "      --l1      1    (default) Inside parameter" );
  opt->addUsage( "      --l2      1    (default) Outside parameter" );
  opt->addUsage( "      --gamma   4000 (default) Overlap penalty weight" );
  opt->addUsage( "      --eta     0 (default)    Laplacian smoothing weight" );
  opt->addUsage( "      --tau     0 (default)    Weight to control volume penalty" );
  opt->addUsage( " -v   --volume  0 (default)    Volume constraint value" );
  opt->addUsage( "" );

  /* 4. SET THE OPTION STRINGS/CHARACTERS */

/* by default all  options  will be checked on the command line
    and from option/resource file */
  /* a flag (takes no argument), supporting long and short form */
  opt->setFlag(  "help", 'h' );
  /* an option (takes an argument), supporting long and short form */
  opt->setOption(  "iter", 'i' );
  opt->setOption(  "rms", 'r' );
  opt->setOption(  "epsilon" );
  opt->setOption(  "mu" );
  opt->setOption(  "nu" );
  opt->setOption(  "l1" );
  opt->setOption(  "l2" );
  opt->setOption(  "gamma" );
  opt->setOption(  "eta" );
  opt->setOption(  "tau" );
  opt->setOption(  "volume", 'v' );

  /* a flag (takes no argument), supporting only short form */
//   opt->setFlag( 'c' );

  /* for options that will be checked only on the command and line not in
  option/resource file */
  /* a flag (takes no argument), supporting long and short form */
//   opt->setCommandFlag(  "zip" , 'z');

  /* for options that will be checked only from the option/resource file */
  /* an option (takes an argument), supporting only long form */
//   opt->setFileOption(  "title" );

  /* 5. PROCESS THE COMMANDLINE AND RESOURCE FILE */

  /* read options from a  option/resource file with ':'
  separated options or flags, one per line */
  opt->processFile( ".options" );
  /* go through the command line and get the options  */
  opt->processCommandArgs( argc, argv );

  if( ! opt->hasOptions())
  { /* print usage if no options */
    opt->printUsage();
    delete opt;
    return EXIT_FAILURE;
  }

  unsigned int nb_iteration = 10;
  double rms = 0.;
  double epsilon = 1.;
  double mu = 0.;
  double nu = 0.;
  double l1 = 1.;
  double l2 = 1.;
  double gamma = 4000;
  double eta = 0.;
  double tau = 0.;
  double volume = 0.;

  /* 6. GET THE VALUES */
  if( opt->getFlag( "help" ) || opt->getFlag( 'h' ) )
  {
    opt->printUsage();
    delete opt;
    return EXIT_FAILURE;
  }
  if( opt->getValue( 'i' ) != NULL  || opt->getValue( "iter" ) != NULL  )
  {
    nb_iteration = atoi( opt->getValue( 'i' ) );
  }
  if( opt->getValue( 'r' ) != NULL  || opt->getValue( "rms" ) != NULL  )
  {
    rms = atof( opt->getValue( "rms" ) );
  }
  if( opt->getValue( "epsilon" ) != NULL  )
  {
    epsilon = atof( opt->getValue( "epsilon" ) );
  }
  if( opt->getValue( "mu" ) != NULL  )
  {
    mu = atof( opt->getValue( "mu" ) );
  }
  if( opt->getValue( "nu" ) != NULL  )
  {
    nu = atof( opt->getValue( "nu" ) );
  }
  if( opt->getValue( "l1" ) != NULL  )
  {
    l1 = atof( opt->getValue( "l1" ) );
  }
  if( opt->getValue( "l2" ) != NULL  )
  {
    l2 = atof( opt->getValue( "l2" ) );
  }
  if( opt->getValue( "gamma" ) != NULL  )
  {
    eta = atof( opt->getValue( "gamma" ) );
  }
  if( opt->getValue( "eta" ) != NULL  )
  {
    eta = atof( opt->getValue( "eta" ) );
  }
  if( opt->getValue( "tau" ) != NULL  )
  {
    tau = atof( opt->getValue( "tau" ) );
  }
  if( opt->getValue( 'v' ) != NULL  || opt->getValue( "volume" ) != NULL  )
  {
    volume = atof( opt->getValue( "volume" ) );
  }

  if( opt->getArgc() != 5 )
  {
    opt->printUsage();
    delete opt;
    return EXIT_FAILURE;
  }

  const unsigned int Dimension = 3;
  typedef float ScalarPixelType;
  typedef itk::Image< ScalarPixelType, Dimension > LevelSetImageType;
  typedef itk::Image< ScalarPixelType, Dimension > FeatureImageType;

  typedef itk::ScalarChanAndVeseLevelSetFunction< LevelSetImageType,
    FeatureImageType > LevelSetFunctionType;
  typedef itk::SparseMultiphaseLevelSetImageFilter< LevelSetImageType,
    FeatureImageType, LevelSetFunctionType, float > MultiLevelSetType;

  typedef itk::ImageFileReader< LevelSetImageType > ContourReaderType;
  typedef itk::ImageFileReader< FeatureImageType > FeatureReaderType;
  typedef itk::ImageFileWriter< LevelSetImageType > WriterType;

  FeatureReaderType::Pointer featureReader = FeatureReaderType::New();
  featureReader->SetFileName( opt->getArgv( 0 ) );
  featureReader->Update();

  ContourReaderType::Pointer contourReader1 = ContourReaderType::New();
  contourReader1->SetFileName( opt->getArgv( 1 ) );
  contourReader1->Update();

  ContourReaderType::Pointer contourReader2 = ContourReaderType::New();
  contourReader2->SetFileName( opt->getArgv( 2 ) );
  contourReader2->Update();

  ContourReaderType::Pointer contourReader3 = ContourReaderType::New();
  contourReader3->SetFileName( opt->getArgv( 3 ) );
  contourReader3->Update();

  FeatureImageType::Pointer featureImage = featureReader->GetOutput();
  LevelSetImageType::Pointer contourImage1 = contourReader1->GetOutput();
  LevelSetImageType::Pointer contourImage2 = contourReader2->GetOutput();
  LevelSetImageType::Pointer contourImage3 = contourReader3->GetOutput();

  MultiLevelSetType::Pointer levelSetFilter = MultiLevelSetType::New();
  levelSetFilter->SetFunctionCount( 3 );
  levelSetFilter->SetFeatureImage( featureImage );
  levelSetFilter->SetLevelSet( 0, contourImage1 );
  levelSetFilter->SetLevelSet( 1, contourImage2 );
  levelSetFilter->SetLevelSet( 2, contourImage3 );
  levelSetFilter->SetNumberOfIterations( nb_iteration );
  levelSetFilter->SetMaximumRMSError( rms );
  levelSetFilter->SetUseImageSpacing( 0 );

  for ( unsigned int i = 0; i < 3; i++ )
  {
  levelSetFilter->GetTypedDifferenceFunction(0)->SetEpsilon( epsilon );
  levelSetFilter->GetTypedDifferenceFunction(0)->SetMu( mu );
  levelSetFilter->GetTypedDifferenceFunction(0)->SetNu( nu );
  levelSetFilter->GetTypedDifferenceFunction(0)->SetLambda1( l1 );
  levelSetFilter->GetTypedDifferenceFunction(0)->SetLambda2( l2 );
  levelSetFilter->GetTypedDifferenceFunction(0)->SetGamma( gamma );
  levelSetFilter->GetTypedDifferenceFunction(0)->SetEta( eta );
  levelSetFilter->GetTypedDifferenceFunction(0)->SetTau( tau );
  levelSetFilter->GetTypedDifferenceFunction(0)->SetVolume( volume );
  }

  levelSetFilter->Update();

#ifndef NDEBUG
  std::cout << levelSetFilter << std::endl;
#endif

  for ( unsigned int i = 0; i < 3; i++ )
  {
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput( levelSetFilter->GetLevelSet( i ) );
    writer->SetFileName( opt->getArgv( 4 ) );

    try
    {
      writer->Update();
    }
    catch( itk::ExceptionObject & excep )
    {
      std::cerr << "Exception caught !" << std::endl;
      std::cerr << excep << std::endl;
      return -1;
    }
  }
  return EXIT_SUCCESS;
}
