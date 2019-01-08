#ifndef itkRandomForest_h
#define itkRandomForest_h

#include "itkImageToImageFilter.h"
#include "RFLibrary/RFTrain.h"
#include "RFLibrary/RandomForestClassifyImageFilter.h"
#include "RFLibrary/ImageCollectionToImageFilter.h"
#include "itkVectorImage.h"
#include "RFLibrary/Library/classifier.h"
#include "RFLibrary/Library/classification.h"
#include "RFLibrary/RandomForestClassifier.h"

namespace itk
{
template< class TInputImage>
class RandomForest : public ImageToImageFilter< TInputImage, TInputImage >
{
public:
    /** Standard class typedefs. */
    typedef RandomForest            Self;
    typedef ImageToImageFilter< TInputImage, TInputImage > Superclass;
    typedef SmartPointer< Self >                 Pointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(RandomForest, ImageToImageFilter);

    /** The image to be inpainted in regions where the mask is white.*/
    void SetInputImage(const TInputImage* image);

    /** The mask to be inpainted. White pixels will be inpainted, black pixels will be passed through to the output.*/
    void SetLabelMap(const TInputImage* mask);

    typedef typename TInputImage::PixelType TPixel;
    typedef itk::VectorImage<TPixel,3 > VectorImageType;
    typedef RandomForestClassifier<TPixel, TPixel,3> RFClassifierType;
    // Iterator for grouping images into a multi-component image
    typedef ImageCollectionConstRegionIteratorWithIndex< TInputImage, VectorImageType> CollectionIter;

    typedef itk::Image< double, 2 > T2DdoubleImage;
    typedef itk::Image<int,3>T3DintImage;

    void SetSegmentationIndices(std::vector<int> SegmentationIndices)
    {
        m_SegmentationIndices = SegmentationIndices;
    }

    void SetBoundingBox(typename T3DintImage::RegionType bbox){
        m_boundingbox = bbox;
    }

    void SetIntermediateSlicesOnly(bool flag)
    {
        m_intermediateslices = flag;
    }

    void SetSlicingAxis(int SlicingAxis){
        m_slicingaxis = SlicingAxis;
    }

protected:
    RandomForest();
    ~RandomForest(){}

    typename TInputImage::Pointer GetInputImage();
    typename TInputImage::Pointer GetLabelMap();

    DataObject::Pointer MakeOutput(unsigned int idx);

    /** Does the real work. */
    virtual void GenerateData();

private:

    RandomForest(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
    std::vector<int> m_SegmentationIndices;
    typename T3DintImage::RegionType m_boundingbox;
    bool m_intermediateslices;
    int m_slicingaxis;

};
} //namespace ITK


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRandomForest.txx"
#endif


#endif // itkRandomForest_h
