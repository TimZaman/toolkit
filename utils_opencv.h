/* Name :       utils*.h [part of tim zaman's utility functions]
 * Version:     v1.0
 * Date :       2012-05-07
 * Developer :  Tim Zaman, Pixelprisma BV (timbobel@gmail.com)
 * Distribution in any form is prohibited.
 */

#ifndef UTILS_OPENCV_TIM_H
#define UTILS_OPENCV_TIM_H

#include <iomanip>
#include <iostream>
#include <stdio.h> //file streams
#include <string>
#include <vector>
#include <clocale>
#include <fcntl.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h> //mkdir
#include <sys/statvfs.h> //Free disk space checker
#include <unistd.h>
#include <pwd.h>
#include <cstdlib>
#include <dirent.h> //dir listing

#ifndef Q_MOC_RUN
	#include <boost/filesystem.hpp>
	#include <boost/lexical_cast.hpp> //int to std:string conversion
	#include <boost/format.hpp>
	#include <boost/variant.hpp>
	#include <boost/regex.hpp>	
	#include <boost/foreach.hpp>
#endif


#include <opencv2/opencv.hpp>

#include <libexif/exif-data.h> //libexif-dev
#include <exiv2/exiv2.hpp> //libexiv2-dev

#ifdef __APPLE__
	#include <asl.h> //Apple System Logger API
#endif

#include "toolkit/utils_general.h"

namespace util{
	//OpenCV related

	bool detectClipping(cv::Mat, int, double, int, double, std::string &);
	void expand(cv::Rect & , double );
	std::string matToJpgString(cv::Mat);
	cv::Rect constrainRectInSize(cv::Rect rCrop, cv::Size sImage);
	cv::RotatedRect scale(cv::RotatedRect, double);
	cv::Rect findBiggestBlob(cv::Mat & matImage);
	cv::RotatedRect findBiggestBlobRot(cv::Mat & matImage);
	cv::Rect retainCenterBlob(cv::Mat & matImage, int);
	cv::Point2f ptMove(cv::Point2f pt, double dist, double angle_deg);
	std::vector<cv::Point> vecrotrect2vecpt(std::vector<cv::RotatedRect> vecRotRect);
	std::vector<std::vector<int> > groupPoints(std::vector<cv::Point> vecPts, double mindist, int mingroupsize);
	cv::Point rect2cp(cv::Rect rRect);
	cv::RotatedRect fixRotatedRect(cv::RotatedRect );
	void rectangle(cv::Mat matImage, cv::RotatedRect rRect, cv::Scalar color, int thickness);
	cv::Point2f rotatePoint(const cv::Point2f& inPoint, const cv::Point2f& center, const double& angRad);
	cv::Point2f rotate2d(const cv::Point2f& inPoint, const double& angRad);
	double pts2angleRad(cv::Point, cv::Point);
	double pts2angleDeg(cv::Point, cv::Point);
	cv::Mat crop(cv::Mat src, cv::RotatedRect rRect);
	void rot90(cv::Mat &matImage, int rotflag);
	void autoClipBrighten(cv::Mat &matImage, double percentile_lower, double percentile_upper);
	cv::Rect constrainRectInSize(cv::Rect rCrop, cv::Size sImage);
	double pointDist(cv::Point pt1, cv::Point pt2);
	double pointDist(cv::Point2f pt1, cv::Point2f pt2);
	void rotate(cv::Mat& src, double angle, cv::Mat& dst);
	cv::Size getFitSize(cv::Size sizeIn, cv::Size sizeOut);
	static void rotatingCalipers( const cv::Point2f* points, int n, float* out );
	cv::RotatedRect minAreaSquare( cv::InputArray _points );
	cv::Mat correctGamma(cv::Mat &img, cv::Vec3d gamma , double contrast , bool forSaving);
	cv::Mat correctGamma(cv::Mat &img, double gammaR,  double gammaG ,  double gammaB , double contrast , bool forSaving);
};




#endif





