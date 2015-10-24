/****************************************************************************
**
** Copyright (C) 2015 Tim Zaman.
** Contact: http://www.timzaman.com/ , timbobel@gmail.com
**
**
****************************************************************************/

#ifndef CODE1D_H
#define CODE1D_H

#pragma once



#include <iostream>
#include <stdio.h> //file streams
#include <string>
#include <clocale>

#include <fcntl.h>
#include <pthread.h> //Multithreading (GTK)

#include <time.h>



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


#include "sort.h"
#include "utils_opencv.h"
#include "utils_general.h"

#define C39_SENTINEL '*'
#define C39_SENTINEL_STRING "nwnnwnwnn"
#define C39_CHARACTERS 44

static const char C39_Characters[C39_CHARACTERS] = {'0','1','2','3','4','5','6','7',
                                                    '8','9','A','B','C','D','E','F',
                                                    'G','H','I','J','K','L','M','N',
                                                    'O','P','Q','R','S','T','U','V',
                                                    'W','X','Y','Z','-','.',' ','$',
                                                    '/','+','%', C39_SENTINEL};

static const char* C39_Strings[C39_CHARACTERS] = {"nnnwwnwnn", "wnnwnnnnw", "nnwwnnnnw",
                                                  "wnwwnnnnn", "nnnwwnnnw", "wnnwwnnnn",
                                                  "nnwwwnnnn", "nnnwnnwnw", "wnnwnnwnn",
                                                  "nnwwnnwnn", "wnnnnwnnw", "nnwnnwnnw",
                                                  "wnwnnwnnn", "nnnnwwnnw", "wnnnwwnnn",
                                                  "nnwnwwnnn", "nnnnnwwnw", "wnnnnwwnn",
                                                  "nnwnnwwnn", "nnnnwwwnn", "wnnnnnnww",
                                                  "nnwnnnnww", "wnwnnnnwn", "nnnnwnnww",
                                                  "wnnnwnnwn", "nnwnwnnwn", "nnnnnnwww",
                                                  "wnnnnnwwn", "nnwnnnwwn", "nnnnwnwwn",
                                                  "wwnnnnnnw", "nwwnnnnnw", "wwwnnnnnn",
                                                  "nwnnwnnnw", "wwnnwnnnn", "nwwnwnnnn",
                                                  "nwnnnnwnw", "wwnnnnwnn", "nwwnnnwnn",
                                                  "nwnwnwnnn", "nwnwnnnwn", "nwnnnwnwn",
                                                  "nnnwnwnwn", C39_SENTINEL_STRING};



static const std::map<std::string, char> generateDecodingMap(){
    std::map<std::string, char> mapping;
    for(int i = 0; i < C39_CHARACTERS; i++){
        mapping[C39_Strings[i]] = C39_Characters[i];
    }
    return mapping;
}

static const std::map<std::string, char> decoding = generateDecodingMap();

struct stripeCode {
  cv::RotatedRect rotRect;
  int barcodeType;
  std::string str;
  std::string strRaw;
};
 

namespace bc1D {

	std::vector<stripeCode> readStripeCode(cv::Mat, double);
	void cpRansac_barcode(std::vector<cv::Point> vecPtsIn, int min_inliers, double max_px_dist, std::vector< std::vector<int> > & vecVecInlierIdx, std::vector<cv::Vec4f> & vecLines, cv::Mat image);
	
	
};	// END NAMESPACE bc1D //
	



#endif





