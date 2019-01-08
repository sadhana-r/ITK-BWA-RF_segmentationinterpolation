#ifndef itkBWAfilter_txx
#define itkBWAfilter_txx

#include <iterator>
#include <vector>
#include "itkObjectFactory.h"
#include "itkImageSliceConstIteratorWithIndex.h"
#include "IRISSlicer/IRISSlicer.h"
#include "itkImage.h"
#include <cmath>
#include <cstring>
#include <algorithm>
#include <iostream>
#include "itkRegionOfInterestImageFilter.h"

namespace itk
{

template< class TInputImage, class TOutputImage>
DataObject::Pointer BinaryWeightedAveragingFilter <TInputImage, TOutputImage >::MakeOutput(unsigned int idx)
{
    DataObject::Pointer output;

    switch ( idx )
    {
    case 0:
        output = ( TInputImage::New() ).GetPointer();
        break;
    case 1:
        output = ( TOutputImage::New() ).GetPointer();
        break;
    default:
        std::cerr << "No output " << idx << std::endl;
        output = NULL;
        break;
    }
    return output.GetPointer();
}

template< class TInputImage, class TOutputImage >
BinaryWeightedAveragingFilter< TInputImage, TOutputImage >::BinaryWeightedAveragingFilter()
{
    this->SetNumberOfRequiredOutputs(2);
    this->SetNumberOfRequiredInputs(1);
    this->SetNthOutput( 0, this->MakeOutput(0) );
    this->SetNthOutput( 1, this->MakeOutput(1) );
}

template< class TInputImage, class TOutputImage>
TInputImage* BinaryWeightedAveragingFilter<TInputImage, TOutputImage >::GetInterpolation()
{
    return dynamic_cast< TInputImage * >(this->ProcessObject::GetOutput(0) );
}

template< class TInputImage, class TOutputImage>
TOutputImage* BinaryWeightedAveragingFilter<TInputImage, TOutputImage >::GetProbabilityMap()
{
    return dynamic_cast< TOutputImage * >(this->ProcessObject::GetOutput(1) );
}


template< class TInputImage, class TOutputImage >
void BinaryWeightedAveragingFilter< TInputImage, TOutputImage >
::GenerateData()
{
    typename TInputImage::ConstPointer input = this->GetInput();

    // Setup output 1
    typename TInputImage::Pointer interpolation = this->GetInterpolation();
    interpolation->SetBufferedRegion( input->GetLargestPossibleRegion() );
    interpolation->Allocate();
    itk::ImageRegionConstIterator<TInputImage> It_interp(input,input->GetLargestPossibleRegion());
    It_interp.GoToBegin();
    while (!It_interp.IsAtEnd())
    {
        int val = input->GetPixel(It_interp.GetIndex());
        interpolation->SetPixel(It_interp.GetIndex(), val);
        ++It_interp;
    }

    // Setup output 2
    typename TOutputImage::Pointer probabilitymap = this->GetProbabilityMap();
    probabilitymap->SetBufferedRegion( input->GetLargestPossibleRegion() );
    probabilitymap->Allocate();
    itk::ImageRegionConstIterator<TInputImage> It_pmap(input,input->GetLargestPossibleRegion());
    It_pmap.GoToBegin();
    while (!It_pmap.IsAtEnd())
    {
        double val = double(input->GetPixel(It_pmap.GetIndex()));
        probabilitymap->SetPixel(It_pmap.GetIndex(), val);
        ++It_pmap;
    }

    typedef itk::RegionOfInterestImageFilter<TInputImage, TInputImage> TrimFilter;
    typename TrimFilter::Pointer TrimImageFilter = TrimFilter::New();
    TrimImageFilter->SetInput(input);
    TrimImageFilter->SetRegionOfInterest(m_boundingbox);
    TrimImageFilter->Update();

    std::vector<int> dimensions;
    dimensions.push_back(0);
    dimensions.push_back(1);
    dimensions.push_back(2);
    dimensions.erase(dimensions.begin() + m_slicingaxis);
    m_firstdirection = dimensions[0];
    m_seconddirection = dimensions[1];

    //Initialize filters to extract slices amd restack to form a 3D volume

    int totalSlices =  m_SegmentationIndices.size();
    typedef IRISSlicer<TInputImage, T2DintImage, TInputImage> IrisSlicerFilterType;
    typename IrisSlicerFilterType::Pointer Slicer[2];

    itk::FixedArray< unsigned int, 3 > layout;
    layout[0] = 1;
    layout[1] = 1;
    layout[2] = 0; // Set third element to zero to get a 3D volume [1,1,0]
    typename TInputImage::PixelType filler = 0;
    typename TOutputImage::PixelType filler_double = 0;

    typename AndImageFilterType::Pointer Intersect[totalSlices-1];
    typename SignedDistanceMapFilterType::Pointer SignedDistanceMapIntersect[totalSlices - 1];

    for ( unsigned int i = 0; i < totalSlices-1; i++ )
    {
        std::cout<< "Interpolating between " <<  m_SegmentationIndices[i] <<" and " <<  m_SegmentationIndices[i+1] <<std::endl;

        const int numSlices =  m_SegmentationIndices[i+1] -  m_SegmentationIndices[i];
        int intermediate_slice = numSlices/2;

        // Extract the first slice
        Slicer[0] = IrisSlicerFilterType::New();
        Slicer[0]->SetInput(TrimImageFilter->GetOutput());
        Slicer[0]->SetSliceIndex( m_SegmentationIndices[i]-1);
        Slicer[0]->SetSliceDirectionImageAxis(m_slicingaxis); //Needs to generalize
        Slicer[0]->SetLineDirectionImageAxis(m_seconddirection);
        Slicer[0]->SetPixelDirectionImageAxis(m_firstdirection);
        Slicer[0]->SetLineTraverseForward(0);
        Slicer[0]->SetPixelTraverseForward(0);
        Slicer[0]->Update();

        // Extract the second slice
        Slicer[1] = IrisSlicerFilterType::New();
        Slicer[1]->SetInput(TrimImageFilter->GetOutput());
        Slicer[1]->SetSliceIndex( m_SegmentationIndices[i+1]-1);
        Slicer[1]->SetSliceDirectionImageAxis(m_slicingaxis);
        Slicer[1]->SetLineDirectionImageAxis(m_seconddirection);
        Slicer[1]->SetPixelDirectionImageAxis(m_firstdirection);
        Slicer[1]->SetLineTraverseForward(0);
        Slicer[1]->SetPixelTraverseForward(0);
        Slicer[1]->Update();

        // Compute connected components
        //Compute A intersect B
        Intersect[i] = AndImageFilterType::New();
        Intersect[i]->SetInput1(Slicer[0]->GetOutput());
        Intersect[i]->SetInput2(Slicer[1]->GetOutput());
        Intersect[i]->Update();

        SignedDistanceMapIntersect[i] = SignedDistanceMapFilterType::New();
        SignedDistanceMapIntersect[i]->SetInput(Intersect[i]->GetOutput());
        SignedDistanceMapIntersect[i]->UseImageSpacingOff();
        SignedDistanceMapIntersect[i]->SquaredDistanceOff();
        SignedDistanceMapIntersect[i]->InsideIsPositiveOn();
        SignedDistanceMapIntersect[i]->Update();

        typename SubtractImageFilterType::Pointer DifferenceImage[2];
        // Compute A\B
        DifferenceImage[0] = SubtractImageFilterType::New();
        DifferenceImage[0]->SetInput1(Slicer[0]->GetOutput());
        DifferenceImage[0]->SetInput2(Intersect[i]->GetOutput());
        DifferenceImage[0]->Update();
        typename T2DintImage::Pointer AdiffB = DifferenceImage[0]->GetOutput();
        AdiffB->DisconnectPipeline();

        // Compute B\A
        DifferenceImage[1] = SubtractImageFilterType::New();
        DifferenceImage[1]->SetInput1(Slicer[1]->GetOutput());
        DifferenceImage[1]->SetInput2(Intersect[i]->GetOutput());
        DifferenceImage[1]->Update();
        typename T2DintImage::Pointer BdiffA = DifferenceImage[1]->GetOutput();
        BdiffA->DisconnectPipeline();

        //Calculate the connected components  in A \ B
        ConnectedComponentImageFilterType::Pointer ConnectedAdifB = ConnectedComponentImageFilterType::New ();
        ConnectedAdifB->SetInput(AdiffB); // A\B
        ConnectedAdifB->Update();

        int numCCAdB = ConnectedAdifB->GetObjectCount();
        typename BinaryThresholdImageFilterType::Pointer ConnectedComponent_AdB[numCCAdB];
        typename OrImageFilterType::Pointer A_dash[numCCAdB];

        //Initialize filters and varaibles to combine results across connected components

        typename TInputImage::Pointer InterpolatedVolume[2];
        typename TOutputImage::Pointer ProbabilityMap[2];
        typename T2DintImage::Pointer InterpolatedVolume_intermediate[2];
        typename T2DdoubleImage::Pointer ProbabilityMap_intermediate[2];

        typename TilerType::Pointer Tiler_Interp[numCCAdB];
        typename DoubleTilerType::Pointer Tiler_ERF[numCCAdB];

        typename SignedDistanceMapFilterType::Pointer SignedDistanceMapA_dash[numCCAdB];

        //Loop through connected components in A\B

        for( int j = 0; j < numCCAdB; j++ )

        {
            //Isolate single connected component
            ConnectedComponent_AdB[j] = BinaryThresholdImageFilterType::New();
            ConnectedComponent_AdB[j]->SetInput(ConnectedAdifB->GetOutput());
            ConnectedComponent_AdB[j]->SetLowerThreshold(j+1);
            ConnectedComponent_AdB[j]->SetUpperThreshold(j+1);
            ConnectedComponent_AdB[j]->SetInsideValue(1);
            ConnectedComponent_AdB[j]->SetOutsideValue(0);
            ConnectedComponent_AdB[j]->Update();

            // A' = (A \intersect B) \u C
            A_dash[j] = OrImageFilterType::New();
            A_dash[j]->SetInput1(ConnectedComponent_AdB[j]->GetOutput());
            A_dash[j]->SetInput2(Intersect[i]->GetOutput());
            A_dash[j]->Update();

            SignedDistanceMapA_dash[j] = SignedDistanceMapFilterType::New();
            SignedDistanceMapA_dash[j]->SetInput(A_dash[j]->GetOutput());
            SignedDistanceMapA_dash[j]->UseImageSpacingOff();
            SignedDistanceMapA_dash[j]->SquaredDistanceOff();
            SignedDistanceMapA_dash[j]->InsideIsPositiveOn();
            SignedDistanceMapA_dash[j]->Update();

            // Calculate the reparameterization function g

            // First calculate the area of difference set A\B
            double AreaAB = 0;
            itk::ImageRegionConstIterator<T2DintImage> It(ConnectedComponent_AdB[j]->GetOutput(),ConnectedComponent_AdB[j]->GetOutput()->GetLargestPossibleRegion());
            It.GoToBegin();
            while (!It.IsAtEnd())
            {
                if (It.Get()!=0) //Or whatever value is your white pixel
                {
                    AreaAB+=1;
                }
                ++It;
            }

            // Compute the reparameterization using an iterator

            // Loop through different values of the parameter t and compute the vaule of g
            const int numElements = 60;
            double gAB [numElements];

            for (int t = 0;t < numElements; t++){

                double I_diffB = 0;
                double B_diffI = 0;
                double AreaIB = 0;
                double t1 = t*(double(1)/(numElements-1));

                itk::ImageRegionConstIterator<T2DdoubleImage> It_sdtA(SignedDistanceMapA_dash[j]->GetOutput(),SignedDistanceMapA_dash[j]->GetOutput()->GetLargestPossibleRegion());
                itk::ImageRegionConstIterator<T2DdoubleImage> It_sdtIntersect(SignedDistanceMapIntersect[i]->GetOutput(),SignedDistanceMapIntersect[i]->GetOutput()->GetLargestPossibleRegion());
                itk::ImageRegionConstIterator<T2DintImage> ItB(Intersect[i]->GetOutput(),Intersect[i]->GetOutput()->GetLargestPossibleRegion());

                It_sdtA.GoToBegin();
                It_sdtIntersect.GoToBegin();
                ItB.GoToBegin();

                while (!It_sdtA.IsAtEnd())
                {
                    double x1 = It_sdtA.Get();
                    double x2 = It_sdtIntersect.Get();
                    double interp = (1-t1)*x1 + t1*x2; // interpolate between the two pixels
                    int B = ItB.Get();

                    if (interp >= 0 & B == 0)
                    {
                        I_diffB+=1;
                    }
                    if (interp < 0 & B == 1)
                    {
                        B_diffI+=1;
                    }

                    ++It_sdtA;\
                    ++It_sdtIntersect;
                    ++ItB;
                }

                AreaIB = I_diffB + B_diffI;
                double g = AreaIB/AreaAB;
                gAB[t] = 1 - g;

            } // End of loop through different paramter values - finished computing gAB

            //Given the reparametrization function, can now interpolate between the different t-values

            // Loop through averaging paramter for interpolation

            if (m_intermediateslices == false){

                Tiler_Interp[j] = TilerType::New();
                Tiler_Interp[j]->SetLayout(layout);
                Tiler_ERF[j] = DoubleTilerType::New();
                Tiler_ERF[j]->SetLayout(layout);
                Tiler_Interp[j]->SetDefaultPixelValue( filler );
                Tiler_ERF[j]->SetDefaultPixelValue( filler_double );

                typename InterpolationFilterType:: Pointer ComputeInterpolation;

                for (int t = 1;t < numSlices; t++){

                    double t1 =  t*(double(1)/numSlices);

                    // Find the reparametrization using inverse of gAB
                    double tdash_max, tdash_min, tdash;
                    for (int n = 0; n < numElements-1; n++){
                        double  gmin_low = gAB[n] - t1;
                        double gmin_high = gAB[n+1] - t1;

                        if ( (gmin_low <= 0) && (gmin_high >= 0)){
                            tdash_min =  n*(double(1)/(numElements-1));
                            tdash_max = (n+1)*(double(1)/(numElements-1));
                            double dist_low = 0 - gmin_low;
                            double dist_high = gmin_high - 0;
                            tdash = dist_high/(dist_low + dist_high)*tdash_min +  dist_low/(dist_low + dist_high)*tdash_max;
                        }
                    }

                    ComputeInterpolation = InterpolationFilterType::New();
                    ComputeInterpolation->SetInputImage1(SignedDistanceMapA_dash[j]->GetOutput());
                    ComputeInterpolation->SetInputImage2(SignedDistanceMapIntersect[i]->GetOutput());
                    ComputeInterpolation->SetConstant(tdash);
                    ComputeInterpolation->Update();

                    T2DintImage::Pointer interp = ComputeInterpolation->GetInterpolation();
                    interp->DisconnectPipeline();
                    T2DdoubleImage::Pointer pmap = ComputeInterpolation->GetProbabilityMap();
                    pmap->DisconnectPipeline();

                    Tiler_Interp[j]->SetInput(t-1,interp);
                    Tiler_ERF[j]->SetInput(t-1,pmap);

                } // End of loop through interpolation using different t values

                Tiler_Interp[j]->Update();
                Tiler_ERF[j]->Update();

                if(j == 0){

                    InterpolatedVolume[0] = Tiler_Interp[j]->GetOutput();
                    ProbabilityMap[0] = Tiler_ERF[j]->GetOutput();
                    InterpolatedVolume[0]->DisconnectPipeline();
                    ProbabilityMap[0]->DisconnectPipeline();
                }
                else{

                    typename ThreeDOrImageFilterType::Pointer CombineComponentsFilterAdB = ThreeDOrImageFilterType::New();
                    CombineComponentsFilterAdB->SetInput1(InterpolatedVolume[0]);
                    CombineComponentsFilterAdB->SetInput2(Tiler_Interp[j]->GetOutput());
                    CombineComponentsFilterAdB->Update();
                    InterpolatedVolume[0] = CombineComponentsFilterAdB->GetOutput();
                    InterpolatedVolume[0]->DisconnectPipeline();

                    typename MaximumImageFilterType::Pointer MaximumImageFilterAdB = MaximumImageFilterType::New ();
                    MaximumImageFilterAdB->SetInput(0, ProbabilityMap[0]);
                    MaximumImageFilterAdB->SetInput(1, Tiler_ERF[j]->GetOutput());
                    MaximumImageFilterAdB ->Update();
                    ProbabilityMap[0] = MaximumImageFilterAdB->GetOutput();
                    ProbabilityMap[0]->DisconnectPipeline();
                }
            }

            if(m_intermediateslices == true){

                double t1 = double(intermediate_slice)/numSlices;
                typename InterpolationFilterType:: Pointer ComputeInterpolation;

                // Find the reparametrization using inverse of gAB
                double tdash_max, tdash_min, tdash;
                for (int n = 0; n < numElements-1; n++){
                    double  gmin_low = gAB[n] - t1;
                    double gmin_high = gAB[n+1] - t1;

                    if ( (gmin_low <= 0) && (gmin_high >= 0)){
                        tdash_min =  n*(double(1)/(numElements-1));
                        tdash_max = (n+1)*(double(1)/(numElements-1));
                        double dist_low = 0 - gmin_low;
                        double dist_high = gmin_high - 0;
                        tdash = dist_high/(dist_low + dist_high)*tdash_min +  dist_low/(dist_low + dist_high)*tdash_max;
                    }
                    ComputeInterpolation = InterpolationFilterType::New();
                    ComputeInterpolation->SetInputImage1(SignedDistanceMapA_dash[j]->GetOutput());
                    ComputeInterpolation->SetInputImage2(SignedDistanceMapIntersect[i]->GetOutput());
                    ComputeInterpolation->SetConstant(tdash);
                    ComputeInterpolation->Update();
                }

                if(j == 0){
                    InterpolatedVolume_intermediate[0] = ComputeInterpolation->GetInterpolation();
                    ProbabilityMap_intermediate[0] = ComputeInterpolation->GetProbabilityMap();
                    InterpolatedVolume_intermediate[0]->DisconnectPipeline();
                    ProbabilityMap_intermediate[0]->DisconnectPipeline();
                }
                else{
                    typename OrImageFilterType::Pointer CombineComponentsFilterAdB = OrImageFilterType::New();
                    CombineComponentsFilterAdB->SetInput1(InterpolatedVolume_intermediate[0]);
                    CombineComponentsFilterAdB->SetInput2(ComputeInterpolation->GetInterpolation());
                    CombineComponentsFilterAdB->Update();
                    InterpolatedVolume_intermediate[0] = CombineComponentsFilterAdB->GetOutput();

                    T2DMaximumImageFilterType::Pointer MaximumImageFilterAdB = T2DMaximumImageFilterType::New ();
                    MaximumImageFilterAdB->SetInput(0, ProbabilityMap_intermediate[0]);
                    MaximumImageFilterAdB->SetInput(1, ComputeInterpolation->GetProbabilityMap());
                    MaximumImageFilterAdB->Update();
                    ProbabilityMap_intermediate[0] = MaximumImageFilterAdB->GetOutput();

                    InterpolatedVolume_intermediate[0]->DisconnectPipeline();
                    ProbabilityMap_intermediate[0]->DisconnectPipeline();
                }
            }

        } // End of loop through connected components in A\B

        //  **********************************************

        // Initialize filters and varaibled need to combine results across connected components
        ConnectedComponentImageFilterType::Pointer ConnectedBdifA = ConnectedComponentImageFilterType::New ();
        ConnectedBdifA->SetInput(BdiffA); // B\A
        ConnectedBdifA->Update();

        int numCCBdA = ConnectedBdifA->GetObjectCount();
        typename BinaryThresholdImageFilterType::Pointer ConnectedComponent_BdA[numCCBdA];
        typename OrImageFilterType::Pointer B_dash[numCCBdA];
        typename SignedDistanceMapFilterType::Pointer SignedDistanceMapB_dash[numCCBdA];

        typename TilerType::Pointer TilerBdA_Interp[numCCBdA];
        typename DoubleTilerType::Pointer TilerBdA_ERF[numCCBdA];

        //  Calculate the connected components  in B\A
        for( int j = 0; j < numCCBdA; j++ )
        {
            // Isolate single connected component
            ConnectedComponent_BdA[j] = BinaryThresholdImageFilterType::New();
            ConnectedComponent_BdA[j]->SetInput(ConnectedBdifA->GetOutput());
            ConnectedComponent_BdA[j]->SetLowerThreshold(j+1);
            ConnectedComponent_BdA[j]->SetUpperThreshold(j+1);
            ConnectedComponent_BdA[j]->SetInsideValue(1);
            ConnectedComponent_BdA[j]->SetOutsideValue(0);
            ConnectedComponent_BdA[j]->Update();

            //B' = (A \intersect B) \u C
            B_dash[j] = OrImageFilterType::New();
            B_dash[j]->SetInput1(ConnectedComponent_BdA[j]->GetOutput());
            B_dash[j]->SetInput2(Intersect[i]->GetOutput());
            B_dash[j]->Update();

            SignedDistanceMapB_dash[j] = SignedDistanceMapFilterType::New();
            SignedDistanceMapB_dash[j]->SetInput(B_dash[j]->GetOutput());
            SignedDistanceMapB_dash[j]->UseImageSpacingOff();
            SignedDistanceMapB_dash[j]->SquaredDistanceOff();
            SignedDistanceMapB_dash[j]->InsideIsPositiveOn();
            SignedDistanceMapB_dash[j]->Update();

            //  Calculate the reparameterization function g

            // Calculate the area of the sym diff between A' and B' - which is the area of C

            // Calcualte the area of the difference set
            double AreaAB = 0;
            itk::ImageRegionConstIterator< T2DintImage > It(ConnectedComponent_BdA[j]->GetOutput(),ConnectedComponent_BdA[j]->GetOutput()->GetLargestPossibleRegion());
            It.GoToBegin();
            while (!It.IsAtEnd())
            {
                if (It.Get()!=0) //Or whatever value is your white pixel
                {
                    AreaAB+=1;
                }
                ++It;
            }

            // Loop through different values of the parameter t and compute the vaule of g
            const int numElements = 60;
            double gAB [numElements];

            for (int t = 0;t < numElements; t++){

                double t1 =  t*(double(1)/(numElements-1));
                double I_diffB = 0;
                double B_diffI = 0;
                double AreaIB = 0;

                itk::ImageRegionConstIterator<T2DdoubleImage> It_sdtIntersect(SignedDistanceMapIntersect[i]->GetOutput(),SignedDistanceMapIntersect[i]->GetOutput()->GetLargestPossibleRegion());
                itk::ImageRegionConstIterator<T2DdoubleImage> It_sdtB(SignedDistanceMapB_dash[j]->GetOutput(),SignedDistanceMapB_dash[j]->GetOutput()->GetLargestPossibleRegion());
                itk::ImageRegionConstIterator<T2DintImage> ItB(B_dash[j]->GetOutput(),B_dash[j]->GetOutput()->GetLargestPossibleRegion());

                It_sdtIntersect.GoToBegin();
                It_sdtB.GoToBegin();
                ItB.GoToBegin();

                while (!It_sdtIntersect.IsAtEnd())
                {
                    double x1 = It_sdtIntersect.Get();
                    double x2 = It_sdtB.Get();

                    double interp = (1-t1)*x1 + t1*x2; // interpolate between the two pixels
                    int B = ItB.Get();

                    if (interp >= 0 & B == 0)
                    {
                        I_diffB+=1;
                    }
                    if (interp < 0 & B == 1)
                    {
                        B_diffI+=1;
                    }

                    ++It_sdtB;\
                    ++It_sdtIntersect;
                    ++ItB;
                }
                AreaIB = I_diffB + B_diffI;
                double g = AreaIB/AreaAB;
                gAB[t] = 1 - g;

            } // End of loop through different paramter values - finished computing gAB

            // Given the reparametrization function, can now interpolate between the different t-values

            // Loop through averaging paramter for interpolation

            if ( m_intermediateslices == false){

                TilerBdA_Interp[j] = TilerType::New();
                TilerBdA_Interp[j]->SetLayout(layout);
                TilerBdA_ERF[j] = DoubleTilerType::New();
                TilerBdA_ERF[j]->SetLayout(layout);
                TilerBdA_ERF[j]->SetDefaultPixelValue( filler_double );
                TilerBdA_Interp[j]->SetDefaultPixelValue( filler );

                typename InterpolationFilterType:: Pointer ComputeInterpolation2;
                T2DintImage::Pointer interp2;
                T2DdoubleImage::Pointer pmap2;

                for (int t = 1;t < numSlices; t++){

                    double t1 = t*(double(1)/numSlices);

                    // Find the reparametrization using inverse of gAB

                    double tdash_max, tdash_min, tdash;
                    for (int n = 0; n < numElements-1; n++){
                        double  gmin_low = gAB[n] - t1;
                        double gmin_high = gAB[n+1] - t1;

                        if ( (gmin_low <= 0) && (gmin_high >= 0)){

                            tdash_min = n*(double(1)/numElements);
                            tdash_max = (n+1)*(double(1)/numElements);
                            double dist_low = 0 - gmin_low;
                            double dist_high = gmin_high - 0;
                            tdash = dist_high/(dist_low + dist_high)*tdash_min +  dist_low/(dist_low + dist_high)*tdash_max;
                        }
                    }

                    ComputeInterpolation2 = InterpolationFilterType::New();
                    ComputeInterpolation2->SetInputImage1(SignedDistanceMapIntersect[i]->GetOutput());
                    ComputeInterpolation2->SetInputImage2(SignedDistanceMapB_dash[j]->GetOutput());
                    ComputeInterpolation2->SetConstant(tdash);
                    ComputeInterpolation2->Update();
                    interp2 = ComputeInterpolation2->GetInterpolation();
                    interp2->DisconnectPipeline();
                    pmap2 = ComputeInterpolation2->GetProbabilityMap();
                    pmap2->DisconnectPipeline();

                    TilerBdA_Interp[j]->SetInput(t-1,interp2);

                    TilerBdA_ERF[j]->SetInput(t-1,pmap2);

                } // Loop through interpolation using different t values

                TilerBdA_Interp[j]->Update();
                TilerBdA_ERF[j]->Update();

                if(j == 0){
                    InterpolatedVolume[1] = TilerBdA_Interp[j]->GetOutput();
                    ProbabilityMap[1] = TilerBdA_ERF[j]->GetOutput();
                    InterpolatedVolume[1]->DisconnectPipeline();
                    ProbabilityMap[1]->DisconnectPipeline();
                }
                else{
                    typename ThreeDOrImageFilterType::Pointer CombineComponentsFilterBdA = ThreeDOrImageFilterType::New();
                    CombineComponentsFilterBdA->SetInput1(InterpolatedVolume[1]);
                    CombineComponentsFilterBdA->SetInput2(TilerBdA_Interp[j]->GetOutput());
                    CombineComponentsFilterBdA->Update();
                    InterpolatedVolume[1] = CombineComponentsFilterBdA->GetOutput();

                    typename MaximumImageFilterType::Pointer MaximumImageFilterBdA = MaximumImageFilterType::New();
                    MaximumImageFilterBdA->SetInput(0, ProbabilityMap[1]);
                    MaximumImageFilterBdA->SetInput(1, TilerBdA_ERF[j]->GetOutput());
                    MaximumImageFilterBdA->Update();
                    ProbabilityMap[1] = MaximumImageFilterBdA->GetOutput();

                    InterpolatedVolume[1]->DisconnectPipeline();
                    ProbabilityMap[1]->DisconnectPipeline();
                }
            }

            if ( m_intermediateslices == true){ // Interpolate only middle slices

                typename InterpolationFilterType:: Pointer ComputeInterpolation2;
                double t1 = double(intermediate_slice)/numSlices;

                // Find the reparametrization using inverse of gAB
                double tdash_max, tdash_min, tdash;
                for (int n = 0; n < numElements-1; n++){
                    double  gmin_low = gAB[n] - t1;
                    double gmin_high = gAB[n+1] - t1;

                    if ( (gmin_low <= 0) && (gmin_high >= 0)){
                        tdash_min =  n*(double(1)/numElements);
                        tdash_max = (n+1)*(double(1)/numElements);
                        double dist_low = 0 - gmin_low;
                        double dist_high = gmin_high - 0;
                        tdash = dist_high/(dist_low + dist_high)*tdash_min +  dist_low/(dist_low + dist_high)*tdash_max;
                    }
                }
                ComputeInterpolation2 = InterpolationFilterType::New();
                ComputeInterpolation2->SetInputImage1(SignedDistanceMapIntersect[i]->GetOutput());
                ComputeInterpolation2->SetInputImage2(SignedDistanceMapB_dash[j]->GetOutput());
                ComputeInterpolation2->SetConstant(tdash);
                ComputeInterpolation2->Update();

                if(j == 0){

                    InterpolatedVolume_intermediate[1] = ComputeInterpolation2->GetInterpolation();
                    InterpolatedVolume_intermediate[1]->DisconnectPipeline();
                    ProbabilityMap_intermediate[1] = ComputeInterpolation2->GetProbabilityMap();
                    ProbabilityMap_intermediate[1]->DisconnectPipeline();
                }
                else{

                    typename OrImageFilterType::Pointer CombineComponentsFilterBdA = OrImageFilterType::New();
                    CombineComponentsFilterBdA->SetInput1(InterpolatedVolume_intermediate[1]);
                    CombineComponentsFilterBdA->SetInput2(ComputeInterpolation2->GetInterpolation());
                    CombineComponentsFilterBdA->Update();
                    InterpolatedVolume_intermediate[1] = CombineComponentsFilterBdA->GetOutput();
                    InterpolatedVolume_intermediate[1]->DisconnectPipeline();

                    typename T2DMaximumImageFilterType::Pointer MaximumImageFilterBdA = T2DMaximumImageFilterType::New ();
                    MaximumImageFilterBdA->SetInput(0, ProbabilityMap_intermediate[1]);
                    MaximumImageFilterBdA->SetInput(1, ComputeInterpolation2->GetProbabilityMap());
                    MaximumImageFilterBdA->Update();
                    ProbabilityMap_intermediate[1] = MaximumImageFilterBdA->GetOutput();
                    ProbabilityMap_intermediate[1]->DisconnectPipeline();
                }
            }

        } // End of loop through connected components in B\A

        itk::ImageSliceConstIteratorWithIndex<TInputImage> It_destination( input, input->GetRequestedRegion() );
        It_destination.SetFirstDirection(m_firstdirection);
        It_destination.SetSecondDirection(m_seconddirection);
        It_destination.GoToBegin();

        if ( m_intermediateslices == false){

            typename ThreeDOrImageFilterType::Pointer CombineInterpolations = ThreeDOrImageFilterType::New();
            CombineInterpolations->SetInput1(InterpolatedVolume[0]);
            CombineInterpolations->SetInput2(InterpolatedVolume[1]);
            CombineInterpolations->Update();
            typename TInputImage::Pointer interpolated_segment = CombineInterpolations->GetOutput();
            interpolated_segment->DisconnectPipeline();

            typename MaximumImageFilterType::Pointer CombineProbabilityMaps = MaximumImageFilterType::New ();
            CombineProbabilityMaps->SetInput(0, ProbabilityMap[0]);
            CombineProbabilityMaps->SetInput(1, ProbabilityMap[1]);
            CombineProbabilityMaps->Update();
            typename TOutputImage::Pointer pmap_segment = CombineProbabilityMaps->GetOutput();
            pmap_segment->DisconnectPipeline();

            typename TInputImage::SizeType sourceDim = interpolated_segment->GetLargestPossibleRegion().GetSize();

//            std::cout << "Copying interpolated region to original image space" << std::endl;
            itk::ImageSliceConstIteratorWithIndex<TInputImage> It_source( interpolated_segment, interpolated_segment->GetLargestPossibleRegion() );
            It_source.SetFirstDirection(0);
            It_source.SetSecondDirection(1);
            It_source.GoToBegin();

            int pixel_val;

            typename TInputImage::IndexType idx;
            typename TInputImage::IndexType source_idx;

            while( !It_source.IsAtEnd() )
            {
                while( !It_source.IsAtEndOfSlice() )
                {

                    while( !It_source.IsAtEndOfLine() )
                    {

                        idx  = It_destination.GetIndex();
                        source_idx = It_source.GetIndex();
                        idx[0] = idx[0] +  m_SegmentationIndices[i];

                        source_idx[0] = sourceDim[0] - source_idx[0] -1;
                        source_idx[1] = sourceDim[1] - source_idx[1] - 1;

                        pixel_val = interpolated_segment->GetPixel(source_idx);
                        interpolation->SetPixel(idx, pixel_val);

                        ++It_source;
                        ++It_destination;
                    }
                    It_source.NextLine();
                    It_destination.NextLine();
                }
                It_source.NextSlice();
                It_destination.NextSlice();
            }

            // Copy probability map
//            std::cout << "Copying probability map to original image space" << std::endl;
            itk::ImageSliceConstIteratorWithIndex<TOutputImage> It_sourceprob( pmap_segment, pmap_segment->GetLargestPossibleRegion() );
            It_sourceprob.SetFirstDirection(0);
            It_sourceprob.SetSecondDirection(1);

            itk::ImageSliceConstIteratorWithIndex<TInputImage> It_destinationprob( input, input->GetRequestedRegion() );
            It_destinationprob.SetFirstDirection(m_firstdirection);
            It_destinationprob.SetSecondDirection(m_seconddirection);

            It_destinationprob.GoToBegin();
            It_sourceprob.GoToBegin();
            typename TOutputImage::IndexType idx_prob;
            typename TOutputImage::IndexType source_idx_prob;
            double prob_val;

            while( !It_sourceprob.IsAtEnd() )
            {
                while( !It_sourceprob.IsAtEndOfSlice() )
                {
                    while( !It_sourceprob.IsAtEndOfLine() )
                    {

                        idx_prob  = It_destinationprob.GetIndex();
                        source_idx_prob = It_sourceprob.GetIndex();
                        idx_prob[0] = idx_prob[0] +  m_SegmentationIndices[i];
                        source_idx_prob[0] = sourceDim[0] - source_idx_prob[0] -1;
                        source_idx_prob[1] = sourceDim[1] - source_idx_prob[1] -1;

                        prob_val = pmap_segment->GetPixel(source_idx_prob);
                        probabilitymap->SetPixel(idx_prob, prob_val);

                        ++It_sourceprob;
                        ++It_destinationprob;
                    }

                    It_sourceprob.NextLine();
                    It_destinationprob.NextLine();
                }
                It_sourceprob.NextSlice();
                It_destinationprob.NextSlice();
            }
        }

        if ( m_intermediateslices == true){

            typename OrImageFilterType::Pointer CombineInterpolations = OrImageFilterType::New();
            CombineInterpolations->SetInput1(InterpolatedVolume_intermediate[0]);
            CombineInterpolations->SetInput2(InterpolatedVolume_intermediate[1]);
            CombineInterpolations->Update();
            typename T2DintImage::Pointer interpolated_segment = CombineInterpolations->GetOutput();
            interpolated_segment->DisconnectPipeline();

            typename T2DMaximumImageFilterType::Pointer CombineProbabilityMaps = T2DMaximumImageFilterType::New ();
            CombineProbabilityMaps->SetInput(0, ProbabilityMap_intermediate[0]);
            CombineProbabilityMaps->SetInput(1, ProbabilityMap_intermediate[1]);
            CombineProbabilityMaps->Update();
            typename T2DdoubleImage::Pointer pmap_segment = CombineProbabilityMaps->GetOutput();
            pmap_segment->DisconnectPipeline();

            typename T2DintImage::SizeType sourceDim = interpolated_segment->GetLargestPossibleRegion().GetSize();
//            std::cout << sourceDim << std::endl;

//            std::cout << "Copying interpolated region to original image space" << std::endl;
            itk::ImageSliceConstIteratorWithIndex<T2DintImage> It_source( interpolated_segment, interpolated_segment->GetLargestPossibleRegion() );
            It_source.SetFirstDirection(0);
            It_source.SetSecondDirection(1);
            It_source.GoToBegin();

            int pixel_val;
            typename TInputImage::IndexType idx;
            typename T2DintImage::IndexType source_idx;

            while( !It_source.IsAtEnd() )
            {
                while( !It_source.IsAtEndOfSlice() )
                {

                    while( !It_source.IsAtEndOfLine() )
                    {

                        idx  = It_destination.GetIndex();
                        source_idx = It_source.GetIndex();
                        idx[0] = idx[0] +  m_SegmentationIndices[i] + intermediate_slice -1;

                        source_idx[0] = sourceDim[0] - source_idx[0] -1;
                        source_idx[1] = sourceDim[1] - source_idx[1] - 1;

                        pixel_val = interpolated_segment->GetPixel(source_idx);
                        interpolation->SetPixel(idx, pixel_val);

                        ++It_source;
                        ++It_destination;
                    }
                    It_source.NextLine();
                    It_destination.NextLine();
                }
                It_source.NextSlice();
                It_destination.NextSlice();
            }

            // Copy probability map
//            std::cout << "Copying probability map to original image space" << std::endl;
            itk::ImageSliceConstIteratorWithIndex<T2DdoubleImage> It_sourceprob(pmap_segment,pmap_segment->GetLargestPossibleRegion() );
            It_sourceprob.SetFirstDirection(0);
            It_sourceprob.SetSecondDirection(1);

            itk::ImageSliceConstIteratorWithIndex<TInputImage> It_destinationprob( input, input->GetRequestedRegion() );
            It_destinationprob.SetFirstDirection(m_firstdirection);
            It_destinationprob.SetSecondDirection(m_seconddirection);

            It_destinationprob.GoToBegin();
            It_sourceprob.GoToBegin();
            typename TOutputImage::IndexType idx_prob;
            typename T2DdoubleImage::IndexType source_idx_prob;
            double prob_val;

            while( !It_sourceprob.IsAtEnd() )
            {
                while( !It_sourceprob.IsAtEndOfSlice() )
                {
                    while( !It_sourceprob.IsAtEndOfLine() )
                    {

                        idx_prob  = It_destinationprob.GetIndex();
                        source_idx_prob = It_sourceprob.GetIndex();
                        idx_prob[0] = idx_prob[0] +  m_SegmentationIndices[i] + intermediate_slice -1;
                        source_idx_prob[0] = sourceDim[0] - source_idx_prob[0] -1;
                        source_idx_prob[1] = sourceDim[1] - source_idx_prob[1] -1;

                        prob_val = pmap_segment->GetPixel(source_idx_prob);
                        probabilitymap->SetPixel(idx_prob, prob_val);

                        ++It_sourceprob;
                        ++It_destinationprob;
                    }

                    It_sourceprob.NextLine();
                    It_destinationprob.NextLine();

                }
                It_sourceprob.NextSlice();
                It_destinationprob.NextSlice();
            }
        }

    } // End of loop through slices

} // End of GenerateData()

} // End of itk namespace

#endif // itkBWAandRFinterpolation_txx

