/****************************************************************************
**
** Copyright (C) 2014 Tim Zaman.
** Contact: http://www.timzaman.com/ , timbobel@gmail.com
**  Original version by G. Bradski 2013
**
****************************************************************************/
#ifndef ZXINGWRAPPER_H
#define ZXINGWRAPPER_H

#pragma once


#include <opencv2/opencv.hpp>

//////////////ZXING BARCODE READER////////////////////////////////////////////////////
#include <zxing/LuminanceSource.h>
#include <zxing/MultiFormatReader.h>
#include <zxing/oned/OneDReader.h>
#include <zxing/oned/EAN8Reader.h>
#include <zxing/oned/EAN13Reader.h>
#include <zxing/oned/Code39Reader.h>
#include <zxing/oned/Code128Reader.h>
#include <zxing/datamatrix/DataMatrixReader.h>
#include <zxing/qrcode/QRCodeReader.h>
#include <zxing/aztec/AztecReader.h>
#include <zxing/common/GlobalHistogramBinarizer.h>
#include <zxing/Exception.h>
#include <zxing/common/GreyscaleLuminanceSource.h>

using namespace zxing;
using namespace oned;
using namespace datamatrix;
using namespace qrcode;
using namespace aztec;

class zxingwrapper : public LuminanceSource
{
private:
    cv::Mat m_pImage;

public:
    zxingwrapper(const cv::Mat &image)
    : LuminanceSource(image.cols, image.rows)
    {
        m_pImage = image.clone();
    }

    ~zxingwrapper()
    {
    }

    int getWidth() const { return m_pImage.cols; }
    int getHeight() const { return m_pImage.rows; }

    ArrayRef<char> getRow(int y, ArrayRef<char> row) const
    {
        int width_ = getWidth();
        if (!row)
            row = ArrayRef<char>(width_);
        const char *p = m_pImage.ptr<char>(y);
        for(int x = 0; x<width_; ++x, ++p)
            row[x] = *p;
        return row;
    }

    ArrayRef<char> getMatrix() const
    {
        int width_ = getWidth();
        int height_ =  getHeight();
        ArrayRef<char> matrix = ArrayRef<char>(width_*height_);
        for (int y = 0; y < height_; ++y)
        {
            const char *p = m_pImage.ptr<char>(y);
            for(int x = 0; x < width_; ++x, ++p)
            {
                matrix[y*width_ + x] = *p;
            }
        }
        return matrix;
    }
    /*
    // The following methods are not supported by this demo (the DataMatrix Reader doesn't call these methods)
    bool isCropSupported() const { return false; }
    Ref<LuminanceSource> crop(int left, int top, int width, int height) {}
    bool isRotateSupported() const { return false; }
    Ref<LuminanceSource> rotateCounterClockwise() {}
    */
};
/*
void decode_image(Reader *reader, cv::Mat &image)
{
    try
    {
        Ref<zxingwrapper> source(new zxingwrapper(image));
        Ref<Binarizer> binarizer(new GlobalHistogramBinarizer(source));
        Ref<BinaryBitmap> bitmap(new BinaryBitmap(binarizer));
        Ref<Result> result(reader->decode(bitmap, DecodeHints(DecodeHints::TRYHARDER_HINT)));//+DecodeHints::DEFAULT_HINT)));
        cout << result->getText()->getText() << endl;
    }
    catch (zxing::Exception& e)
    {
        cerr << "Error: " << e.what() << endl;
    }
}
*/

//ZXINGWRAPPER_H
#endif