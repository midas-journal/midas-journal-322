INCLUDE_DIRECTORIES(
  ${CELLSEGMENTATION_SOURCE_DIR}/ActiveContoursWithoutEdges/Common
  ${CELLSEGMENTATION_SOURCE_DIR}/ActiveContoursWithoutEdges/Algorithms/Dense
  ${CELLSEGMENTATION_SOURCE_DIR}/ActiveContoursWithoutEdges/Algorithms/Sparse
  ${CELLSEGMENTATION_SOURCE_DIR}/ActiveContoursWithoutEdges/BasicFilters/Dense
  ${CELLSEGMENTATION_SOURCE_DIR}/ActiveContoursWithoutEdges/BasicFilters/Sparse
  ${CELLSEGMENTATION_SOURCE_DIR}/ActiveContoursWithoutEdges/LevelsetFunctions
  ${CELLSEGMENTATION_SOURCE_DIR}/Options )

SET( EXAMPLE_LIST
  ScalarSinglePhaseDense2DTest
  ScalarSinglePhaseSparse2DTest
  ScalarMultiPhase2DTest
  ScalarSinglePhase3DTest
  ScalarMultiPhase3DTest
 )

FOREACH( var ${EXAMPLE_LIST} )
  ADD_EXECUTABLE( ${var} ${var}.cxx )
  TARGET_LINK_LIBRARIES( ${var} ITKCommon ITKBasicFilters
  ITKStatistics ITKNumerics ITKIO Option )
ENDFOREACH( var )

ADD_TEST( ScalarSinglePhaseDense2DTest
${CXX_TEST_PATH}/ScalarSinglePhaseDense2DTest
${MY_DATAROOT}/example2D.png ${MY_DATAROOT}/phi.mha
${MY_DATAROOT}/segDense2D.mha )

ADD_TEST( ScalarSinglePhaseSparse2DTest
${CXX_TEST_PATH}/ScalarSinglePhaseSparse2DTest
${MY_DATAROOT}/example2D.png ${MY_DATAROOT}/phi.mha
${MY_DATAROOT}/segSparse2D.mha )

