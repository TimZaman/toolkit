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
#include <sstream>
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
#include <math.h>
#include <limits>

// @TODO(tzaman): replace all boost with C++11
//#include <regex>

#ifdef __APPLE__
    #include <asl.h> //Apple System Logger API
#endif

template <typename Tp>
std::string to_string_with_precision(const Tp a_value, const int n = 6) {
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
    std::vector<std::string> split(std::string, std::string);
    
    bool anySubstringInString(std::vector<std::string> &, std::string );

    std::string urlencode(const std::string &s);

    template <typename T>
    double calcMean(std::vector<T> scores) { // Calculates mean
        if (scores.size() == 0) {
            std::cerr << "Warning! Empty vector in calcMedian!" << std::endl;
            return 0;
        }
        double sum = 0;
        for (int i = 0; i < scores.size(); i++) {
            sum += scores[i];
        }
        return sum / double(scores.size());
    }

    double calcMedian(std::vector<double>);

    double calcMeanOfQuarterAndThreeQuarterPercentile(std::vector<double>);
    
    void logASL(std::string);
    std::string escapeRegex(std::string);

 #if defined BOOST_RE_REGEX_H || defined BOOST_RE_REGEX_HPP
    std::vector<std::string> regexReplaceInVector(std::vector<std::string>, std::string, std::string);
    bool isValidURL(std::string);
    std::vector<std::string> getRegexMatches(std::string, std::string );
    void splitDoubleRegex(std::vector<std::string> &, std::vector<std::string> &);
    std::map<std::string, std::string> relateFormatAndFile(std::string, std::string);
#ifdef BOOST_FORMAT_HPP
    std::string regex_escape(const std::string&);
    std::vector<std::string> correlateFileVectorFormat(std::vector<std::string>, std::string , int , int &, std::vector<std::string> &);
#endif // BOOST_FORMAT_HPP
#endif // BOOST_RE_REGEX_H || BOOST_RE_REGEX_HPP   


#ifdef BOOST_FILESYSTEM_FILESYSTEM_HPP
    std::vector<std::string> folderFilesToVector(std::string folder);
    std::string changeFileExtension(std::string, std::string);
#endif //BOOST_FILESYSTEM_FILESYSTEM_HPP


    
    std::string ReplaceAll(std::string, const std::string&, const std::string&);
    

    int xfilelength(int );

    double interpolate(double, std::vector< std::pair<double, double> > &);
    void makeBezier(double, double, int, int*, int*);
    void makeBezier(double, double, int, std::vector<int> &, std::vector<int> &);

#ifdef __EXIF_DATA_H__
    ExifEntry* init_tag(ExifData *, ExifIfd, ExifTag);
    ExifEntry* create_tag(ExifData *, ExifIfd, ExifTag, size_t);
#endif

};

#endif
