/* Name :       utils*.h [part of tim zaman's utility functions]
 * Version:     v1.0
 * Date :       2012-05-07
 * Developer :  Tim Zaman, Pixelprisma BV (timbobel@gmail.com)
 * Distribution in any form is prohibited.
 */

#ifndef UTILS_QT_OPENCV_TIM_H
#define UTILS_QT_OPENCV_TIM_H


#include <QObject>
#include <QString>
#include <QGraphicsPixmapItem>
//#include <QtGui>

#include <iostream>
#include <string>
#include <clocale>

#include <sys/types.h>
#include <sys/stat.h>

#include <unistd.h>
#include <pwd.h>

#include <iomanip>
#include <cstdlib>

#include <opencv2/opencv.hpp>




namespace util{
	QImage Mat2QImage(const cv::Mat3b &); //Opencv and Qt
};


#endif





