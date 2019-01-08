#ifndef ITKBWAANDRFINTERPOLATION_HXX
#define ITKBWAANDRFINTERPOLATION_HXX

#include "itkObjectFactory.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include <math.h>
#include "itkImageFileWriter.h"

namespace itk
{

template< class TInputImage, class TOutputImage>
CombineBWAandRFFilter <TInputImage,TOutputImage>
::CombineBWAandRFFilter ()
{
    this->SetNumberOfRequiredInputs(2);
    this->SetNumberOfRequiredOutputs(2);

    this->SetNthOutput( 0, this->MakeOutput(0) );
    this->SetNthOutput( 1, this->MakeOutput(1) );

}

template< class TInputImage, class TOutputImage>
void
CombineBWAandRFFilter <TInputImage,TOutputImage>
::SetIntensityImage(const TInputImage* image)
{
    this->SetNthInput(0, const_cast<TInputImage*>(image));
}

template< class TInputImage, class TOutputImage>
void
CombineBWAandRFFilter <TInputImage,TOutputImage>
::SetSegmentationImage(const TOutputImage* mask)
{
    this->SetNthInput(1, const_cast<TOutputImage*>(mask));
}

template< class TInputImage, class TOutputImage>
typename TInputImage::ConstPointer CombineBWAandRFFilter <TInputImage,TOutputImage>::GetIntensityImage()
{
    return static_cast< const TInputImage * >
            ( this->ProcessObject::GetInput(0) );
}

template< class TInputImage, class TOutputImage>
typename TOutputImage::ConstPointer CombineBWAandRFFilter <TInputImage,TOutputImage>::GetSegmentationImage()
{
    return static_cast< const TOutputImage * >
            ( this->ProcessObject::GetInput(1) );
}

template< class TInputImage, class TOutputImage>
TOutputImage* CombineBWAandRFFilter <TInputImage, TOutputImage>::GetInterpolation()
{
    return dynamic_cast< TOutputImage * >(this->ProcessObject::GetOutput(0) );
}

template< class TInputImage, class TOutputImage>
TInputImage* CombineBWAandRFFilter <TInputImage, TOutputImage>::GetProbabilityMap()
{
    return dynamic_cast< TInputImage * >(this->ProcessObject::GetOutput(1) );
}

template< typename TInputImage, typename TOutputImage >
DataObject::Pointer CombineBWAandRFFilter <TInputImage, TOutputImage>::MakeOutput(unsigned int idx)
{
    DataObject::Pointer output;

    switch ( idx )
    {
    case 0:
        output = ( TOutputImage::New() ).GetPointer();
        break;
    case 1:
        output = ( TInputImage::New() ).GetPointer();
        break;
    default:
        std::cerr << "No output " << idx << std::endl;
        output = NULL;
        break;
    }
    return output.GetPointer();
}

template<unsigned int VDim>
void ExpandRegion(itk::ImageRegion<VDim> &region, const itk::Index<VDim> &idx)
{
    if(region.GetNumberOfPixels() == 0)
    {
        region.SetIndex(idx);
        for(size_t i = 0; i < VDim; i++)
            region.SetSize(i, 1);
    }
    else {
        for(size_t i = 0; i < VDim; i++)
        {
            if(region.GetIndex(i) > idx[i])
            {
                region.SetSize(i, region.GetSize(i) + (region.GetIndex(i) - idx[i]));
                region.SetIndex(i, idx[i]);
            }
            else if(region.GetIndex(i) + (long) region.GetSize(i) <= idx[i]) {
                region.SetSize(i, 1 + idx[i] - region.GetIndex(i));
            }
        }
    }
}


template< class TInputImage, class TOutputImage>
void
CombineBWAandRFFilter<TInputImage,TOutputImage>
::GenerateData()
{

    typename TInputImage::ConstPointer IntensityImage = this->GetIntensityImage();
    typename TOutputImage::ConstPointer SegmentationImage = this->GetSegmentationImage();

    // Find bounding box
    typename TOutputImage::RegionType bbox;

    // Find the extent of the non-background region of the image
    itk::ImageRegionConstIterator<TOutputImage> it(SegmentationImage, SegmentationImage->GetBufferedRegion());
    for( ; !it.IsAtEnd(); ++it)
        if(it.Value() != 0)
            ExpandRegion(bbox, it.GetIndex());

    // Make sure the bounding box is within the contents of the image
    bbox.Crop(SegmentationImage->GetBufferedRegion());

    std::cout << "Bounding box size is " << bbox.GetSize() << std::endl;

    // Chop off the region
    typedef itk::RegionOfInterestImageFilter<TOutputImage, TOutputImage> TrimFilterInt;
    typename TrimFilterInt::Pointer CropSegmentationImage = TrimFilterInt::New();
    CropSegmentationImage->SetInput(SegmentationImage);
    CropSegmentationImage->SetRegionOfInterest(bbox);
    CropSegmentationImage->Update();

    // Chop off the region
    typedef itk::RegionOfInterestImageFilter<TInputImage, TInputImage> TrimFilterDouble;
    typename TrimFilterDouble::Pointer CropIntensityImage = TrimFilterDouble::New();
    CropIntensityImage->SetInput(IntensityImage);
    CropIntensityImage->SetRegionOfInterest(bbox);
    CropIntensityImage->Update();

    // Determine which ones are the segmented slices
    std::cout << "Determining slicing direction" << std::endl;

    std::vector<int> nEmptySlices;

    // Try first dimension
    itk::ImageSliceConstIteratorWithIndex<TOutputImage> It_dim1( SegmentationImage, SegmentationImage->GetRequestedRegion() );
    It_dim1.SetFirstDirection(1);
    It_dim1.SetSecondDirection(2);

    It_dim1.GoToBegin();
    int counter = 0;
    int value = 0;

    while( !It_dim1.IsAtEnd() )
    {
        while( !It_dim1.IsAtEndOfSlice() )
        {
            while( !It_dim1.IsAtEndOfLine() )
            {
                if (It_dim1.Get() != 0){
                    value = It_dim1.Get();
                }
                ++It_dim1;
            }
            It_dim1.NextLine();
        }
        if( value == 0){
            ++counter; // count the number of empty slices
        }
        It_dim1.NextSlice();
        value = 0;
    }

    nEmptySlices.push_back(counter);

    // Try second dimension
    itk::ImageSliceConstIteratorWithIndex<TOutputImage> It_dim2(SegmentationImage, SegmentationImage->GetRequestedRegion() );
    It_dim2.SetFirstDirection(0);
    It_dim2.SetSecondDirection(2);

    It_dim2.GoToBegin();
    counter = 0;
    value = 0;

    while( !It_dim2.IsAtEnd() )
    {
        while( !It_dim2.IsAtEndOfSlice() )
        {
            while( !It_dim2.IsAtEndOfLine() )
            {
                if (It_dim2.Get() != 0){
                    value = It_dim2.Get();
                }
                ++It_dim2;
            }
            It_dim2.NextLine();
        }
        It_dim2.NextSlice();
        if( value == 0){
            ++ counter;
        }
        value = 0;
    }

    nEmptySlices.push_back(counter);

    // Try third dimension
    itk::ImageSliceConstIteratorWithIndex<TOutputImage> It_dim3(SegmentationImage, SegmentationImage->GetRequestedRegion() );
    It_dim3.SetFirstDirection(0);
    It_dim3.SetSecondDirection(1);

    It_dim3.GoToBegin();
    counter = 0;
    value = 0;

    while( !It_dim3.IsAtEnd() )
    {
        while( !It_dim3.IsAtEndOfSlice() )
        {
            while( !It_dim3.IsAtEndOfLine() )
            {
                if (It_dim3.Get() != 0){
                    value = It_dim3.Get();
                }
                ++It_dim3;
            }
            It_dim3.NextLine();
        }
        It_dim3.NextSlice();
        if( value == 0){
            ++ counter;
        }
        value = 0;
    }

    nEmptySlices.push_back(counter);
    m_SlicingAxis = std::distance(nEmptySlices.begin(), std::max_element(nEmptySlices.begin(), nEmptySlices.end()));

    std::vector<int> dimensions;
    dimensions.push_back(0);
    dimensions.push_back(1);
    dimensions.push_back(2);

    dimensions.erase(dimensions.begin() + m_SlicingAxis);

    //After determining slicing direction, can determine which slices have segmentations
    itk::ImageSliceConstIteratorWithIndex<TOutputImage> It(SegmentationImage,SegmentationImage->GetRequestedRegion() );
    It.SetFirstDirection(dimensions[0]);
    It.SetSecondDirection(dimensions[1]);

    It.GoToBegin();
    counter = 1;
    value = 0;

    while( !It.IsAtEnd() )
    {
        while( !It.IsAtEndOfSlice() )
        {
            while( !It.IsAtEndOfLine() )
            {
                if (It.Get() != 0){
                    value = It.Get();
                }
                ++It;
            }
            It.NextLine();
        }
        It.NextSlice();
        if( value != 0){
            m_SegmentationIndices.push_back(counter);
        }
        ++counter;
        value = 0;
    }
    std::cout << "Slices containing segmentations are " << std::endl;
    std::copy( m_SegmentationIndices.begin(),  m_SegmentationIndices.end(), std::ostream_iterator<int>(std::cout, " "));

    // Create the BWA filter
    std::cout << "\nPerforming BWA interpolation" << std::endl;
    typename BWAFilterType::Pointer BWAfilter = BWAFilterType::New();
    BWAfilter->SetInput(SegmentationImage);
    BWAfilter->SetIntermediateSlicesOnly(m_intermediateslices);
    BWAfilter->SetSlicingAxis(m_SlicingAxis);
    BWAfilter->SetSegmentationIndices(m_SegmentationIndices);
    BWAfilter->SetBoundingBox(bbox);
    BWAfilter->Update();

    typedef itk::ImageFileWriter<TInputImage>ThreeDWriterType;
    typename ThreeDWriterType::Pointer writer1 = ThreeDWriterType::New();
    writer1->SetInput(BWAfilter->GetProbabilityMap());
    writer1->SetFileName("ProbabilityMap_BWA.nii" );
    writer1->Update();

    std::cout << "Generating RF label map" << std::endl;
    typename RFLabelFilterType::Pointer generateRFlabelmap = RFLabelFilterType::New();
    generateRFlabelmap->SetInput(SegmentationImage);
    generateRFlabelmap->SetSlicingAxis(m_SlicingAxis);
    generateRFlabelmap->Update(); // Background is 1, FG is 0

    std::cout << "Performing random forest classification" << std::endl;
    typename RandomForestFilterType::Pointer randomForestClassifier = RandomForestFilterType::New();
    randomForestClassifier->SetInputImage(IntensityImage);
    randomForestClassifier->SetLabelMap(generateRFlabelmap->GetOutput());
    randomForestClassifier->SetBoundingBox(bbox);
    randomForestClassifier->SetIntermediateSlicesOnly(m_intermediateslices);
    randomForestClassifier->SetSegmentationIndices(m_SegmentationIndices);
    randomForestClassifier->SetSlicingAxis(m_SlicingAxis);
    randomForestClassifier->Update();

    typedef itk::ImageFileWriter<TInputImage>ThreeDWriterType;
    typename ThreeDWriterType::Pointer writer = ThreeDWriterType::New();
    writer->SetInput(randomForestClassifier->GetOutput());
    writer->SetFileName("ProbabilityMap_RandomForest.nii" );
    writer->Update();

    // Combine probability maps by taking the average
    typename AddImageFilterType::Pointer addProbabilityMaps = AddImageFilterType::New();
    addProbabilityMaps->SetInput1(BWAfilter->GetProbabilityMap());
    addProbabilityMaps->SetInput2(randomForestClassifier->GetOutput());

    typename MultiplyImageFilterType::Pointer scaleProbabilityMap = MultiplyImageFilterType::New();
    scaleProbabilityMap->SetInput(addProbabilityMaps->GetOutput());
    scaleProbabilityMap->SetConstant(0.5);
    scaleProbabilityMap->Update();

    // Threshold the probability map to obtain the final segmentation
    typename BinaryThresholdFilterType::Pointer thresholdProbabilityMap = BinaryThresholdFilterType::New();
    thresholdProbabilityMap->SetInput(scaleProbabilityMap->GetOutput());
    thresholdProbabilityMap->SetLowerThreshold(0.5);
    thresholdProbabilityMap->SetUpperThreshold(1);
    thresholdProbabilityMap->SetInsideValue(1);
    thresholdProbabilityMap->Update();

    this->GetInterpolation()->Graft(thresholdProbabilityMap->GetOutput());
    this->GetProbabilityMap()->Graft(scaleProbabilityMap->GetOutput());

    std::cout << "Output region size" << this->GetInterpolation()->GetLargestPossibleRegion().GetSize() << std::endl;
    std::cout << "Output pmap region size" << this->GetProbabilityMap()->GetLargestPossibleRegion().GetSize() << std::endl;
}

}// end namespace

#endif // ITKBWAANDRFINTERPOLATION_HXX

