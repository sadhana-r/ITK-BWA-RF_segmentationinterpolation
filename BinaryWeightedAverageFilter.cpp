#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBWAandRFinterpolation.h"
#include <string>


int main(int argc, char *argv[])
{
    if(argc != 4){
        std::cout << "Incorrect number of input arguments. Please enter 1) Intensity Image, "
                     "2) Segmentation Image and 3) true/false for intermediate slices vs whole volume" << std::endl;

    }
    else{
        std::string intensity_image(argv[1]);
        std::string segmentation_image(argv[2]);
        bool m_intermediateslices = argv[3];
    
    std::cout << intensity_image << std::endl;
    std::cout << segmentation_image << std::endl;


    // Setup types
    typedef itk::Image<short, 3>            LabelType;
    typedef itk::VectorImage<short,3>       GreyVectorType;
    typedef itk::Image<short, 3>            GreyScalarType;
    typedef itk::Image<double,3>            ProbabilityType;

    typedef itk::ImageFileWriter< ProbabilityType > ProbabilityWriterType;
    typedef itk::ImageFileWriter< LabelType > LabelWriterType;
    typedef itk::ImageFileReader< LabelType > ReaderType;
    typedef itk::CombineBWAandRFFilter < GreyScalarType, GreyVectorType, LabelType > FilterType;

    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(segmentation_image);
    reader->Update();
    LabelType::Pointer SegmentationImage = reader->GetOutput();

    ReaderType::Pointer reader_intensityimg = ReaderType::New();
    reader_intensityimg->SetFileName(intensity_image);
    reader_intensityimg->Update();
    GreyScalarType::Pointer IntensityImage = reader_intensityimg->GetOutput();

    FilterType::Pointer filter = FilterType::New();
    filter->AddScalarImage(IntensityImage);
    filter->SetSegmentationImage(SegmentationImage);
    filter->SetIntermediateSlicesOnly(m_intermediateslices);
    filter->Update();

    LabelWriterType::Pointer writer = LabelWriterType::New();
    writer->SetInput(filter->GetInterpolation());
    writer->SetFileName("output_Interpolation.nii" );
    writer->Update();

    ProbabilityWriterType::Pointer writer2 = ProbabilityWriterType::New();
    writer2->SetInput(filter->GetProbabilityMap());
    writer2->SetFileName("output_ProbabilityMap.nii" );
    writer2->Update();

    }
    return EXIT_SUCCESS;
}


