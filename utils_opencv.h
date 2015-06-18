/* Name :       utils*.h [part of tim zaman's utility functions]
 * Version:     v1.0
 * Date :       2012-05-07
 * Developer :  Tim Zaman, Pixelprisma BV (timbobel@gmail.com)
 * Distribution in any form is prohibited.
 */

#ifndef UTILS_OPENCV_TIM_H
#define UTILS_OPENCV_TIM_H

#include <iostream>
#include <stdio.h> //file streams
#include <string>
#include <clocale>

#include <fcntl.h>
#include <pthread.h> //Multithreading (GTK)

#include <time.h>

#include <sys/statvfs.h> //Free disk space checker


#ifndef Q_MOC_RUN
	#include <boost/program_options.hpp> //libboost-dev, libboost-1.53-all-dev
	#include <boost/filesystem.hpp>
	#include <boost/lexical_cast.hpp> //int to std:string conversion
	#include <boost/format.hpp>
	//#include <boost/thread/thread.hpp> //for this::thread::sleep
	#include <boost/variant.hpp>
	#include <boost/regex.hpp>	
	#include <boost/foreach.hpp>
#endif


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

#include <libexif/exif-data.h> //libexif-dev
#include <exiv2/exiv2.hpp> //libexiv2-dev

#ifdef __APPLE__
	#include <asl.h> //Apple System Logger API
#endif


using namespace std;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

namespace util{

	//OpenCV related
	cv::Mat crop(cv::Mat src, cv::RotatedRect rRect);
	void rot90(cv::Mat &matImage, int rotflag);
	void autoClipBrighten(cv::Mat &matImage, double percentile_lower, double percentile_upper);
	cv::Rect constrainRectInSize(cv::Rect rCrop, cv::Size sImage);
	double pointDist(cv::Point pt1, cv::Point pt2);
	double pointDist(cv::Point2f pt1, cv::Point2f pt2);
	void rotate(cv::Mat& src, double angle, cv::Mat& dst);

};




#endif





