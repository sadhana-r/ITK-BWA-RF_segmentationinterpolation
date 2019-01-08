#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBWAandRFinterpolation.h"
//#include <ctime>
#include <string>


int main(int argc, char *argv[])
{
    if(argc != 4){
        std::cout << "Incorrect number of imput arguments" << std::endl;

    }
    else{
        std::string intensity_image(argv[1]);
        std::string segmentation_image(argv[2]);
        bool m_intermediateslices = argv[3];
    
    std::cout << intensity_image << std::endl;
    std::cout << segmentation_image << std::endl;

    //int start_s=clock();

    // Setup types
    typedef itk::Image<int, 3> ThreeDImageType;
    typedef itk::Image<double, 3> ThreeDImageDoubleType;

    typedef itk::ImageFileWriter< ThreeDImageDoubleType > ThreeDDoubleWriterType;
    typedef itk::ImageFileWriter< ThreeDImageType > ThreeDWriterType;
    typedef itk::ImageFileReader< ThreeDImageType > ReaderType;
    typedef itk::CombineBWAandRFFilter < ThreeDImageDoubleType, ThreeDImageType > FilterType;

    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(segmentation_image);
    reader->Update();
    ThreeDImageType::Pointer SegmentationImage = reader->GetOutput();

    typedef itk::ImageFileReader< ThreeDImageDoubleType > ReaderDoubleType;
    ReaderDoubleType::Pointer reader_intensityimg = ReaderDoubleType::New();
    reader_intensityimg->SetFileName(intensity_image);
    reader_intensityimg->Update();
    ThreeDImageDoubleType::Pointer IntensityImage = reader_intensityimg->GetOutput();

    FilterType::Pointer filter = FilterType::New();
    filter->SetIntensityImage(IntensityImage);
    filter->SetSegmentationImage(SegmentationImage);
    filter->SetIntermediateSlicesOnly(m_intermediateslices);
    filter->Update();

    ThreeDWriterType::Pointer writer = ThreeDWriterType::New();
    writer->SetInput(filter->GetInterpolation());
    writer->SetFileName("OutputInterpolation.nii" );
    writer->Update();

    ThreeDDoubleWriterType::Pointer writer2 = ThreeDDoubleWriterType::New();
    writer2->SetInput(filter->GetProbabilityMap());
    writer2->SetFileName("OutputProbabilityMap.nii" );
    writer2->Update();

   // int stop_s=clock();
   // std::cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << std::endl;

    }
    return EXIT_SUCCESS;
}


