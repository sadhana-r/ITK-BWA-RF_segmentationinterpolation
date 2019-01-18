#ifndef ITKBWAANDRFINTERPOLATION_HXX
#define ITKBWAANDRFINTERPOLATION_HXX

#include <math.h>
#include <itkSize.h>
#include "itkObjectFactory.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkBWAandRFinterpolation.h"
#include "RFLibrary/ImageCollectionToImageFilter.h"


namespace itk
{

template< class ImageScalarType, class ImageVectorType, class TLabelImage>
CombineBWAandRFFilter <ImageScalarType,  ImageVectorType,TLabelImage>
::CombineBWAandRFFilter ()
    : m_Label( 0 ),
      m_LabeledSlices( TLabelImage::ImageDimension ) // initialize with empty sets
{
    this->SetNumberOfRequiredInputs(2);
    this->SetNumberOfRequiredOutputs(2);

    this->SetNthOutput( 0, this->MakeOutput(0) );
    this->SetNthOutput( 1, this->MakeOutput(1) );

}

template< class ImageScalarType, class ImageVectorType, class TLabelImage>
void
CombineBWAandRFFilter <ImageScalarType,  ImageVectorType,TLabelImage>
::AddScalarImage(ImageScalarType *image)
{
    this->SetNthInput(0, const_cast<ImageScalarType*>(image));
}

template< class ImageScalarType, class ImageVectorType, class TLabelImage>
void
CombineBWAandRFFilter <ImageScalarType,  ImageVectorType,TLabelImage>
::AddVectorImage(ImageVectorType *image)
{
    this->SetNthInput(0, const_cast<ImageVectorType*>(image));
}

template< class ImageScalarType, class ImageVectorType, class TLabelImage>
void
CombineBWAandRFFilter <ImageScalarType,  ImageVectorType,TLabelImage>
::SetSegmentationImage(const TLabelImage* mask)
{
    this->SetNthInput(1, const_cast<TLabelImage*>(mask));
}

template< class ImageScalarType, class ImageVectorType, class TLabelImage>
typename TLabelImage::ConstPointer
CombineBWAandRFFilter<ImageScalarType,  ImageVectorType,TLabelImage>
::GetSegmentationImage()
{
    return static_cast< const TLabelImage * >
            ( this->ProcessObject::GetInput(1) );
}

template< class ImageScalarType, class ImageVectorType, class TLabelImage>
TLabelImage * CombineBWAandRFFilter <ImageScalarType,  ImageVectorType,TLabelImage>
::GetInterpolation()
{
    return dynamic_cast< TLabelImage * >(this->ProcessObject::GetOutput(0) );
}

template< class ImageScalarType, class ImageVectorType, class TLabelImage>
ProbabilityType * CombineBWAandRFFilter <ImageScalarType,  ImageVectorType,TLabelImage>
::GetProbabilityMap()
{
    return dynamic_cast< ProbabilityType * >(this->ProcessObject::GetOutput(1) );
}

template< class ImageScalarType, class ImageVectorType, class TLabelImage>
DataObject::Pointer CombineBWAandRFFilter<ImageScalarType,  ImageVectorType,TLabelImage>
::MakeOutput(unsigned int idx)
{
    DataObject::Pointer output;

    switch ( idx )
    {
    case 0:
        output = ( TLabelImage::New() ).GetPointer();
        break;
    case 1:
        output = ( ProbabilityType::New() ).GetPointer();
        break;
    default:
        std::cerr << "No output " << idx << std::endl;
        output = NULL;
        break;
    }
    return output.GetPointer();
}

template< class ImageScalarType, class ImageVectorType, class TLabelImage>
template< typename T2 >
void
CombineBWAandRFFilter <ImageScalarType,  ImageVectorType,TLabelImage>
::ExpandRegion( typename T2::RegionType& region, const typename T2::IndexType& index )
{
    for ( unsigned int a = 0; a < T2::ImageDimension; ++a )
    {
        if ( region.GetIndex( a ) > index[a] )
        {
            region.SetSize( a, region.GetSize( a ) + region.GetIndex( a ) - index[a] );
            region.SetIndex( a, index[a] );
        }
        else if ( region.GetIndex( a ) + (typename T2::IndexValueType)region.GetSize( a ) <= index[a] )
        {
            region.SetSize( a, index[a] - region.GetIndex( a ) + 1 );
        }
        // else it is already within
    }
}

template< class ImageScalarType, class ImageVectorType, class TLabelImage>
void
CombineBWAandRFFilter<ImageScalarType,  ImageVectorType,TLabelImage>
::DetermineSliceOrientations()
{
    m_LabeledSlices.clear();
    m_LabeledSlices.resize( TLabelImage::ImageDimension ); // initialize with empty sets
    m_BoundingBoxes.clear();

    typename TLabelImage::ConstPointer m_Input = this->GetSegmentationImage();
    typename TLabelImage::Pointer m_Output = this->GetInterpolation();

    typename TLabelImage::RegionType region = m_Output->GetRequestedRegion();
    ImageRegionConstIteratorWithIndex< TLabelImage > it( m_Input, region );
    while ( !it.IsAtEnd() )
    {
        typename TLabelImage::IndexType indPrev, indNext;
        const typename TLabelImage::IndexType ind = it.GetIndex();
        const typename TLabelImage::PixelType val = m_Input->GetPixel( ind );
        if ( val != 0 )
        {
            typename TLabelImage::RegionType boundingBox1;
            boundingBox1.SetIndex( ind );
            for ( unsigned int a = 0; a < TLabelImage::ImageDimension; ++a )
            {
                boundingBox1.SetSize( a, 1 );
            } //region of size [1,1,1]
            std::pair< typename BoundingBoxesType::iterator, bool > resBB
                    = m_BoundingBoxes.insert( std::make_pair( val, boundingBox1 ) );
            if ( !resBB.second ) // include this index in existing BB
            {
                ExpandRegion< TLabelImage >( resBB.first->second, ind );
            }

            unsigned int cTrue = 0;
            unsigned int cAdjacent = 0;
            unsigned int axis = 0;
            for ( unsigned int a = 0; a < TLabelImage::ImageDimension; ++a )
            {
                indPrev = ind;
                indPrev[a]--;
                indNext = ind;
                indNext[a]++;
                typename TLabelImage::PixelType prev = 0;
                if ( region.IsInside( indPrev ) )
                {
                    prev = m_Input->GetPixel( indPrev );
                }
                typename TLabelImage::PixelType next = 0;
                if ( region.IsInside( indNext ) )
                {
                    next = m_Input->GetPixel( indNext );
                }
                if ( prev == 0 && next == 0 ) // && - isolated slices only, || - flat edges too
                {
                    axis = a;
                    ++cTrue;
                }
                else if ( prev == val && next == val )
                {
                    ++cAdjacent;
                }
            }
            if ( cTrue == 1 && cAdjacent == TLabelImage::ImageDimension - 1 )
                // slice has empty adjacent space only along one axis
            {
                // if ( m_Axis == -1 || m_Axis == int(axis) )
                // {
                m_LabeledSlices[axis][val].insert( ind[axis] );
                //std::cout << ind << std::endl;
                //}
            }
        }
        ++it;
    }

} // >::DetermineSliceOrientations

template< class ImageScalarType, class ImageVectorType, class TLabelImage>
void
CombineBWAandRFFilter <ImageScalarType,  ImageVectorType,TLabelImage>
::GenerateData()
{

    typename TLabelImage::ConstPointer SegmentationImage = this->GetSegmentationImage();

    //IOPixelType pixelType = this->ProcessObject::GetInput(0)->GetPixelType();
    DataObject *dataobj = this->ProcessObject::GetInput(0);

    typename TLabelImage::Pointer m_Output = this->GetInterpolation();
    this->AllocateOutputs();

    this->DetermineSliceOrientations();

    if ( m_BoundingBoxes.size() == 0)
    {
        ImageAlgorithm::Copy< TLabelImage, TLabelImage >( SegmentationImage.GetPointer(), m_Output.GetPointer(),
                                                          m_Output->GetRequestedRegion(), m_Output->GetRequestedRegion() );
        return; // no contours detected - no segmentations drawn
    }

    typename TLabelImage::RegionType bbox;
    // To dertemine slicing direcion,get the size of m_LabelledSlices along each axis
    for ( unsigned i = 0; i < TLabelImage::ImageDimension; i++ ) // loop through axes
    {
        //std::cout << m_LabeledSlices[i].size() <<std::endl;
        if (m_LabeledSlices[i].size())
        { //if the axis has some labelled slices then
            m_SlicingAxis = i;

            // then determine the indices of the segmented slices
            for ( typename LabeledSlicesType::iterator it = m_LabeledSlices[i].begin();it != m_LabeledSlices[i].end(); ++it )
            {
                //if ( m_Label == 0 || m_Label == it->first ) // label needs to be interpolated
                if (1 == it->first ) //only interpolate label of value 1
                {
                    // Get the bounding box for label == 1
                    for ( typename BoundingBoxesType::iterator iBB = m_BoundingBoxes.begin();iBB != m_BoundingBoxes.end(); ++iBB ){
                        //if (m_Label == 0 || m_Label == iBB->first)
                        {
                            std::cout << "label value " << iBB->first << std::endl;
                            bbox = iBB->second;
                            std::cout << "bbox index " << bbox.GetIndex() << std::endl;
                        }
                    }

                    for ( typename SliceSetType::iterator s = it->second.begin();
                          s != it->second.end();
                          ++s ){
                        m_SegmentationIndices.push_back(*s - bbox.GetIndex()[i]);

                    }

                }
            }

        }
    }

    // Perform contour interpolation
    std::cout << "Binary Weighted Averaging " << std::endl;
    typename BWAFilterType::Pointer BWAfilter = BWAFilterType::New();
    BWAfilter->SetInput(SegmentationImage);
    BWAfilter->SetIntermediateSlicesOnly(m_intermediateslices);
    BWAfilter->SetSlicingAxis(m_SlicingAxis);
    BWAfilter->SetSegmentationIndices(m_SegmentationIndices);
    BWAfilter->SetBoundingBox(bbox);
    BWAfilter->Update();

    typename RFLabelFilterType::Pointer generateRFlabelmap = RFLabelFilterType::New();
    generateRFlabelmap->SetInput(SegmentationImage);
    generateRFlabelmap->SetSlicingAxis(m_SlicingAxis);
    generateRFlabelmap->Update(); // Background is 1, FG is 0

    typename RandomForestClassifierType::Pointer randomForestClassifier = RandomForestClassifierType::New();

    ImageScalarType *image = dynamic_cast<ImageScalarType *>(dataobj);
    if(image)
    {
        randomForestClassifier->AddScalarImage(image);
    }
    else
    {
        ImageVectorType *vecImage = dynamic_cast<ImageVectorType *>(dataobj);
        if(vecImage)
        {
            randomForestClassifier->AddVectorImage(vecImage);
        }
        else
        {
            itkAssertInDebugOrThrowInReleaseMacro(
                        "Wrong input type to ImageCollectionConstRegionIteratorWithIndex");
        }
    }

    std::cout << "Random Forest Classification " << std::endl;
    randomForestClassifier->SetLabelMap(generateRFlabelmap->GetOutput());
    randomForestClassifier->SetBoundingBox(bbox);
    randomForestClassifier->SetIntermediateSlicesOnly(m_intermediateslices);
    randomForestClassifier->SetSegmentationIndices(m_SegmentationIndices);
    randomForestClassifier->SetSlicingAxis(m_SlicingAxis);
    randomForestClassifier->Update();

    // Combine probability maps by taking the average
    typename AddImageFilterType::Pointer addProbabilityMaps = AddImageFilterType::New();
    addProbabilityMaps->SetInput1(BWAfilter->GetProbabilityMap());
    addProbabilityMaps->SetInput2(randomForestClassifier->GetOutput());

    typename MultiplyImageFilterType::Pointer scaleProbabilityMap = MultiplyImageFilterType::New();
    scaleProbabilityMap->SetInput(addProbabilityMaps->GetOutput());
    scaleProbabilityMap->SetConstant(0.5);
    scaleProbabilityMap->Update();

    //Threshold the probability map to obtain the final segmentation
    typename BinaryThresholdFilterType::Pointer thresholdProbabilityMap = BinaryThresholdFilterType::New();
    thresholdProbabilityMap->SetInput(scaleProbabilityMap->GetOutput());
    thresholdProbabilityMap->SetLowerThreshold(0.5);
    thresholdProbabilityMap->SetUpperThreshold(1);
    thresholdProbabilityMap->SetInsideValue(1);
    thresholdProbabilityMap->Update();

    this->GetInterpolation()->Graft(thresholdProbabilityMap->GetOutput());
    this->GetProbabilityMap()->Graft(scaleProbabilityMap->GetOutput());
}

}// end namespace

#endif // ITKBWAANDRFINTERPOLATION_HXX

