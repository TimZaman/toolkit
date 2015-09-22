/****************************************************************************
**
** Copyright (C) 2015 Tim Zaman.
** Contact: http://www.timzaman.com/ , timbobel@gmail.com
**
**
****************************************************************************/

#pragma once


#include <iostream>
#include <stdio.h> //file streams
#include <string>
#include <clocale>

#include <fcntl.h>
#include <pthread.h> //Multithreading (GTK)

#include <time.h>

#include <sys/statvfs.h> //Free disk space checker



#include <boost/program_options.hpp> //libboost-dev, libboost-1.53-all-dev
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp> //int to std:string conversion
#include <boost/format.hpp>
#include <boost/variant.hpp>
#include <boost/regex.hpp>	
#include <boost/foreach.hpp>



#include <sys/types.h>
#include <sys/stat.h>

#include <unistd.h>
#include <pwd.h>

#include <iomanip>
#include <cstdlib>
#include <sys/stat.h> //mkdir
#include <dirent.h> //dir listing

#include <opencv2/opencv.hpp>

#include <vector>

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
#include <zxing/common/IllegalArgumentException.h>
#include <zxing/common/Counted.h>
#include <zxing/Binarizer.h>
#include <zxing/MultiFormatReader.h>
#include <zxing/Result.h>
#include <zxing/ReaderException.h>
#include <zxing/common/GlobalHistogramBinarizer.h>
#include <zxing/common/HybridBinarizer.h>
#include <exception>
#include <zxing/Exception.h>
#include <zxing/common/IllegalArgumentException.h>
#include <zxing/BinaryBitmap.h>
#include <zxing/DecodeHints.h>
#include <zxing/qrcode/QRCodeReader.h>
#include <zxing/multi/qrcode/QRCodeMultiReader.h>
#include <zxing/multi/ByQuadrantReader.h>
#include <zxing/multi/MultipleBarcodeReader.h>
#include <zxing/multi/GenericMultipleBarcodeReader.h>

#include <zxing/datamatrix/decoder/Decoder.h>
#include <zxing/common/reedsolomon/ReedSolomonDecoder.h>
#include <zxing/common/Counted.h>
#include <zxing/common/Array.h>
#include <zxing/common/DecoderResult.h>
#include <zxing/common/BitMatrix.h>
#include <zxing/ResultPoint.h>
#include <zxing/ChecksumException.h>

#include "zxingwrapper.h" //ZXING

#ifdef __APPLE__
	#include <asl.h> //Apple System Logger API
#endif

#include "toolkit/sort.h"
#include "toolkit/utils_opencv.h"
#include "toolkit/utils_general.h"


using zxing::Ref;
using zxing::ArrayRef;
using zxing::LuminanceSource;
using namespace zxing;
using namespace datamatrix;
using namespace zxing::datamatrix;



namespace bc {


	std::string decode_pure_barcode(cv::Mat matImage);
	std::string decode_image_barcode(const cv::Mat &matImage, std::vector<int> vecBarcodeTypes, int numwiggles=1);

	std::string readQR(cv::Mat matImage, double dpi);
	std::string readDMTX(cv::Mat matImage, double dpi);

	// **************** //
};	// END NAMESPACE BC //
	// **************** //









