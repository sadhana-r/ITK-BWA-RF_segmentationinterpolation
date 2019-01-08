#ifndef itkBWAandRFinterpolation_h
#define itkBWAandRFinterpolation_h

#include "itkBWAfilter.h"
#include "itkRFLabelMap.h"
#include "itkRandomForest.h"
#include "itkAddImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImageSliceConstIteratorWithIndex.h"
#include "itkImageToImageFilter.h"

namespace itk
{
template< class TInputImage, class TOutputImage>
class CombineBWAandRFFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
    /** Standard class typedefs. */
    typedef CombineBWAandRFFilter             Self;
    typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
    typedef SmartPointer< Self >                 Pointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(CombineBWAandRFFilter , ImageToImageFilter);

    void SetIntensityImage(const TInputImage* image);

    void SetSegmentationImage(const TOutputImage* mask);

    void SetIntermediateSlicesOnly(bool flag)
    {
        m_intermediateslices = flag;
    }

    TOutputImage* GetInterpolation();
    TInputImage* GetProbabilityMap();

    typedef itk::BinaryWeightedAveragingFilter< TOutputImage, TInputImage >  BWAFilterType;
    typedef itk::RFLabelMap<TOutputImage, TInputImage> RFLabelFilterType;
    typedef itk::RandomForest< TInputImage > RandomForestFilterType;
    typedef itk::BinaryThresholdImageFilter<TInputImage, TOutputImage> BinaryThresholdFilterType;
    typedef itk::MultiplyImageFilter< TInputImage > MultiplyImageFilterType;
    typedef itk::AddImageFilter<TInputImage, TInputImage> AddImageFilterType;

protected:
    CombineBWAandRFFilter ();
    ~CombineBWAandRFFilter (){}

    typename TOutputImage::ConstPointer GetSegmentationImage();
    typename TInputImage::ConstPointer GetIntensityImage();

    DataObject::Pointer MakeOutput(unsigned int idx);

    /** Does the real work. */
    virtual void GenerateData();

private:

    CombineBWAandRFFilter (const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
    std::vector<int> m_SegmentationIndices;
    int m_SlicingAxis;
    bool m_intermediateslices;
};
} //namespace ITK

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBWAandRFinterpolation.txx"
#endif

#endif // itkBWAandRFinterpolation_h

