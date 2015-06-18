/* Name :       utils.h [part of pixelprisma utility functions]
 * Version:     v1.0
 * Date :       2012-05-07
 * Developer :  Tim Zaman, Pixelprisma BV (timbobel@gmail.com)
 * Distribution in any form is prohibited.
 */

#ifndef UTILS_TIM_H
#define UTILS_TIM_H

#ifdef QT_VERSION //Built with QT support
	#include <QObject>
	#include <QThread>
	#include <QMutex>
	#include <QString>
	#include <QGraphicsPixmapItem>
	#include <QtGui>
	#include <QList>
	#include <QMetaType>
#endif 


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

#ifdef CV_VERSION //If built with OpenCV support
	#include <opencv2/opencv.hpp>
#endif

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

	enum{ GPTIM_RAW=0, GPTIM_JPEG=1, GPTIM_TIFF=2};

	// === CONSTANTS === //
	//enum { T_DZ_QWALLI=0, T_GTD=1, T_GTO_FULL=2, T_UTT=3, T_NN=4, T_SG=5, T_CC_CRC=6,
	//			T_CC_MINI=7, T_CC_PASS=8, T_QA62_1=9, T_QA62_2=10, T_KGS_2007=11,
	//			T_KGS_2000=12, T_KGS_Q14=13, T_KCCP_13=14, T_KCCP_14=15, T_COUNT=16,T_GTO_HALF=17, T_DZ_DYN=18, T_DZ_IT8=19, T_DZ_SFR=20};

	enum errType { ERR_NONE = 0, ERR_GTO = 1 , ERR_BARCODE = 2, ERR_WRITE = 3, ERR_TARGET = 4, ERR_GTD = 5  };
	enum enumProfile { PROF_SRGB=0, PROF_ADOBE98=1, PROF_ECIRGBV2=2, PROF_GRAYGAMMA22=3, PROF_NONE=4 };
	enum enumSavetype{ SAVE_TIFF=0, SAVE_JPEG=1, SAVE_JP2=2, SAVE_RAW=3};

	//const float ITU709weight[3]={0.2126, 0.7152, 0.0722}; //In R G B

	//const float dpiboost_range = 1.07; //+-[p] was 1.06 for deltae and naturalis



	// === FUNCTIONS === //
	bool isValidURL(std::string);

	void logASL(std::string);
	std::string escapeRegex(std::string);
	std::vector<std::string> regexReplaceInVector(std::vector<std::string> , std::string, std::string);
	std::vector<std::string> correlateFileVectorFormat(std::vector<std::string> , std::string , int , int &, std::vector<std::string> &);
	std::vector<std::string> folderFilesToVector(string folder);
	std::map<std::string, std::string> relateFormatAndFile(std::string, std::string);
	vector<string> getRegexMatches(std::string strRegex, std::string strMatches);
	std::string ReplaceAll(std::string str, const std::string& from, const std::string& to);
	std::string fileformatToRegex(std::string);
	std::string regex_escape(const std::string&);

	int xfilelength(int );

	double interpolate(double x, vector< pair<double, double> > &table);
	void makeBezier(double gamma, double contrast, int N_SEG, int lutX[], int lutY[]);
	void makeBezier(double gamma, double contrast, int N_SEG, vector<int> &, vector<int> &);

	ExifEntry* init_tag(ExifData *, ExifIfd, ExifTag );
	ExifEntry* create_tag(ExifData *, ExifIfd, ExifTag , size_t );

	//OpenCV related
	#ifdef CV_VERSION //If built with OpenCV support
		cv::Mat crop(cv::Mat src, cv::RotatedRect rRect);
		void rot90(cv::Mat &matImage, int rotflag);
		void autoClipBrighten(cv::Mat &matImage, double percentile_lower, double percentile_upper);
		cv::Rect constrainRectInSize(cv::Rect rCrop, cv::Size sImage);
		double pointDist(cv::Point pt1, cv::Point pt2);
		double pointDist(cv::Point2f pt1, cv::Point2f pt2);
		void rotate(cv::Mat& src, double angle, cv::Mat& dst);
		#ifdef QT_VERSION //Built with QT support
			QImage Mat2QImage(const cv::Mat3b &);
		#endif
	#endif

};




#endif





