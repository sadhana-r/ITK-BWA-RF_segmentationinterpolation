cmake_minimum_required(VERSION 3.5.0)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()

# This project is designed to be built outside the Insight source tree.
project(BinaryWeightedAverageFilter)

# Find ITK.
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
if (ITKVtkGlue_LOADED)
  find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
endif()


set(HEADERS
   itkBWAfilter.h
        itkBWAfilter.hxx
        itkComputeInterpolation.h
        itkComputeInterpolation.txx
        itkRFLabelMap.h
        itkRFLabelMap.txx
        itkRandomForest.h
        itkRandomForest.txx
        itkBWAandRFinterpolation.h
        itkBWAandRFinterpolation.txx
./RFLibrary/RFTrain.h
./RFLibrary/RandomForestClassifyImageFilter.h
./RFLibrary/RandomForestClassifyImageFilter.txx
./RFLibrary/RandomForestClassifier.h
./RFLibrary/ImageCollectionToImageFilter.h
./RFLibrary/Library/classification.h
./RFLibrary/Library/classifier.h
./RFLibrary/Library/data.h
./RFLibrary/Library/forest.h
./RFLibrary/Library/imageio.h
./RFLibrary/Library/linearalgebra.h
./RFLibrary/Library/node.h
./RFLibrary/Library/random.h
./RFLibrary/Library/statistics.h
./RFLibrary/Library/trainer.h
./RFLibrary/Library/trainingcontext.h
./RFLibrary/Library/tree.h
./RFLibrary/Library/type.h
./RFLibrary/Library/utility.h
        ./IRISSlicer/IRISSlicer.txx
        ./IRISSlicer/SNAPCommon.h
        ./IRISSlicer/RLEImageRegionConstIterator.h
        ./IRISSlicer/RLEImageConstIterator.h
        ./IRISSlicer/IRISSlicer.h
        ./IRISSlicer/ImageCoordinateTransform.h
        ./IRISSlicer/IRISVectorTypes.h
        ./IRISSlicer/IRISVectorTypes.txx
        ./IRISSlicer/RLEImage.h
        ./IRISSlicer/RLEImage.txx
        ./IRISSlicer/IRISSlicer_RLE.txx
    )

add_executable(${PROJECT_NAME} "BinaryWeightedAverageFilter.cpp" ${HEADERS})

if( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(${PROJECT_NAME} ITKReview ${ITK_LIBRARIES})
else( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(${PROJECT_NAME} ${ITK_LIBRARIES})
endif( "${ITK_VERSION_MAJOR}" LESS 4 )





