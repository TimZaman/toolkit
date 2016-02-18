/**
--------------------------------------------------------------------------------
-   Module      :   zxingwrapper.h
-   Description :   Provides ability for easy use of OpenCV with Zxing-C++
-   Author      :   G. Bradski, 2013
-                   Tim Zaman, 18-FEB-2016
--------------------------------------------------------------------------------
*/

/*

Copyright (c) 2013 G. Bradski
Copyright (c) 2016 Tim Zaman

Permission to use, copy, modify, distribute, and sell this software
for any purpose is hereby granted without fee, provided
that (i) the above copyright notices and this permission notice appear in
all copies of the software and related documentation, and (ii) the names of
Mike Johnson and BancTec may not be used in any advertising or
publicity relating to the software.

THE SOFTWARE IS PROVIDED "AS-IS" AND WITHOUT WARRANTY OF ANY KIND,
EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY
WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.
 
*/

#ifndef ZXINGWRAPPER_H
#define ZXINGWRAPPER_H

#pragma once

#include <opencv2/opencv.hpp>

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
	zxingwrapper(const cv::Mat &image) : LuminanceSource(image.cols, image.rows) {
		m_pImage = image.clone();
	}

	~zxingwrapper(){

	}

	int getWidth() const { 
		return m_pImage.cols; 
	}
	int getHeight() const {
		return m_pImage.rows;
	}

	ArrayRef<char> getRow(int y, ArrayRef<char> row) const {
		int width_ = getWidth();
		if (!row){
			row = ArrayRef<char>(width_);
		}
		const char *p = m_pImage.ptr<char>(y);
		for(int x = 0; x<width_; ++x, ++p){
			row[x] = *p;
		}
		return row;
	}

	ArrayRef<char> getMatrix() const {
		int width_ = getWidth();
		int height_ =  getHeight();
		ArrayRef<char> matrix = ArrayRef<char>(width_*height_);
		for (int y = 0; y < height_; ++y){
			const char *p = m_pImage.ptr<char>(y);
			for(int x = 0; x < width_; ++x, ++p){
				matrix[y*width_ + x] = *p;
			}
		}
		return matrix;
	}
};

//ZXINGWRAPPER_H
#endif