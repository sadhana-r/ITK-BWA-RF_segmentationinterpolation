#ifndef itkBWAfilter_H
#define itkBWAfilter_H

#include "itkImageToImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkAndImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkOrImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkTileImageFilter.h"
#include "itkComputeInterpolation.h"
#include "itkMaximumImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"

namespace itk
{
template< class TInputImage, class TOutputImage>
class BinaryWeightedAveragingFilter:public ImageToImageFilter< TInputImage, TInputImage >
{
public:
    /** Standard class typedefs. */
    typedef BinaryWeightedAveragingFilter Self;
    typedef ImageToImageFilter< TInputImage, TInputImage > Superclass;
    typedef SmartPointer< Self > Pointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(BinaryWeightedAveragingFilter, ImageToImageFilter);

    TInputImage* GetInterpolation();
    TOutputImage* GetProbabilityMap();

    /** Image type information. */
    typedef itk::Image< int, 2 > T2DintImage;
    typedef itk::Image< double, 2 > T2DdoubleImage;

    /** Extract some information from the image types.  Dimensionality
   * of the two images is assumed to be the same. */
    typedef typename TInputImage::PixelType          InputPixelType;
    typedef typename TInputImage::InternalPixelType  InputInternalPixelType;

    /** Initialize ITK filters required */
    typedef itk::AndImageFilter<T2DintImage,T2DintImage,T2DintImage>AndImageFilterType;
    typedef itk::SubtractImageFilter<T2DintImage,T2DintImage,T2DintImage>SubtractImageFilterType;
    typedef itk::OrImageFilter<T2DintImage,T2DintImage,T2DintImage>OrImageFilterType;
    typedef itk::BinaryThresholdImageFilter<T2DintImage, T2DintImage>BinaryThresholdImageFilterType;
    typedef itk::ConnectedComponentImageFilter <T2DintImage,T2DintImage >ConnectedComponentImageFilterType;
    typedef itk::SignedMaurerDistanceMapImageFilter<T2DintImage, T2DdoubleImage> SignedDistanceMapFilterType;
    typedef itk::TileImageFilter< T2DintImage, TInputImage > TilerType;
    typedef itk::TileImageFilter< T2DdoubleImage, TOutputImage > DoubleTilerType;
    typedef itk::OrImageFilter<TInputImage,TInputImage,TInputImage>ThreeDOrImageFilterType;
    typedef itk::ComputeInterpolation< T2DdoubleImage, T2DintImage>  InterpolationFilterType;
    typedef itk::MaximumImageFilter< TOutputImage>  MaximumImageFilterType;
    typedef itk::MaximumImageFilter< T2DdoubleImage>  T2DMaximumImageFilterType;

    void SetIntermediateSlicesOnly(bool flag)
    {
        m_intermediateslices = flag;
    }

    void SetSegmentationIndices(std::vector<int> seg_index){
        m_SegmentationIndices = seg_index;
    }

    void SetBoundingBox(typename TInputImage::RegionType bbox){
        m_boundingbox = bbox;
    }

    void SetSlicingAxis(int SlicingAxis){
        m_slicingaxis = SlicingAxis;
    }

protected:
    BinaryWeightedAveragingFilter();
    ~BinaryWeightedAveragingFilter(){}

    /** Does the real work. */
    virtual void GenerateData();

    /** Create the Output */
    DataObject::Pointer MakeOutput(unsigned int idx);

private:
    BinaryWeightedAveragingFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
    bool m_intermediateslices;
    std::vector<int> m_SegmentationIndices;
    int m_slicingaxis;
    int m_firstdirection;
    int m_seconddirection;
    typename TInputImage::RegionType m_boundingbox;

    template <typename TImage>
    void CreateImage(TImage* const image, int numRows, int numCols)
    {
      // Create an image with 2 connected components
      typename TImage::IndexType corner = {{0,0}};

      typename TImage::SizeType size = {{numRows, numCols}};
      typename TImage::RegionType region(corner, size);

      image->SetRegions(region);
      image->Allocate();

      image->FillBuffer(0);
    }

};
} //namespace ITK

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBWAfilter.hxx"
#endif

#endif // BINARYWEIGHTEDAVERAGEFILTER_H
