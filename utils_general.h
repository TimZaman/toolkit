/* Name :       utils*.h [part of tim zaman's utility functions]
 * Version:     v1.0
 * Date :       2012-05-07
 * Developer :  Tim Zaman, Pixelprisma BV (timbobel@gmail.com)
 * Distribution in any form is prohibited.
 */

#ifndef UTILS_GENERAL_TIM_H
#define UTILS_GENERAL_TIM_H

#include <iomanip>
#include <iostream>
#include <stdio.h> //file streams
#include <string>
#include <vector>
#include <clocale>
#include <fcntl.h>
//#include <pthread.h> //Multithreading (GTK)
#include <time.h>
#include <sys/statvfs.h> //Free disk space checker
#include <sys/types.h>
#include <sys/stat.h> //mkdir
#include <dirent.h> //dir listing
#include <unistd.h>
#include <pwd.h>
#include <cstdlib>

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

#include <opencv2/opencv.hpp>

#include <libexif/exif-data.h> 
#include <exiv2/exiv2.hpp>  //is this needed?

#ifdef __APPLE__
	#include <asl.h> //Apple System Logger API
#endif


template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}



namespace util{

	enum{ GPTIM_RAW=0, GPTIM_JPEG=1, GPTIM_TIFF=2};

	enum errType { ERR_NONE = 0, ERR_GTO = 1 , ERR_BARCODE = 2, ERR_WRITE = 3, ERR_TARGET = 4, ERR_GTD = 5  }; //Is this used?
	enum enumProfile { PROF_SRGB=0, PROF_ADOBE98=1, PROF_ECIRGBV2=2, PROF_GRAYGAMMA22=3, PROF_NONE=4 };
	enum enumSavetype{ SAVE_TIFF=0, SAVE_JPEG=1, SAVE_JP2=2, SAVE_RAW=3};

	// === FUNCTIONS === //
	std::string urlencode(const std::string &s);
	double calcMean(std::vector<double> scores);
	double calcMedian(std::vector<double> scores);
	bool isValidURL(std::string);

	void logASL(std::string);
	std::string escapeRegex(std::string);
	std::vector<std::string> regexReplaceInVector(std::vector<std::string>, std::string, std::string);
	std::vector<std::string> correlateFileVectorFormat(std::vector<std::string> , std::string , int , int &, std::vector<std::string> &);
	std::vector<std::string> folderFilesToVector(std::string folder);
	std::map<std::string, std::string> relateFormatAndFile(std::string, std::string);
	std::vector<std::string> getRegexMatches(std::string strRegex, std::string strMatches);
	std::string ReplaceAll(std::string str, const std::string& from, const std::string& to);
	std::string fileformatToRegex(std::string);
	std::string regex_escape(const std::string&);

	int xfilelength(int );

	double interpolate(double x, std::vector< std::pair<double, double> > &table);
	void makeBezier(double gamma, double contrast, int N_SEG, int lutX[], int lutY[]);
	void makeBezier(double gamma, double contrast, int N_SEG, std::vector<int> &, std::vector<int> &);

	ExifEntry* init_tag(ExifData *, ExifIfd, ExifTag );
	ExifEntry* create_tag(ExifData *, ExifIfd, ExifTag , size_t );

};




#endif





