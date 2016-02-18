/**
--------------------------------------------------------------------------------
-   Module      :   keepalive.h
-   Description :   A wrapper to keep a QT app alive, written in C++ and C#,
-                   so made for the OSX platform.
-   Author      :   Tim Zaman, 18-FEB-2016
--------------------------------------------------------------------------------
*/

/*

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

#ifndef KEEPALIVE_H
#define KEEPALIVE_H

#include <QDebug>

class KeepAlive
{
public:

	KeepAlive(){
		qDebug() << "KeepAlive::KeepAlive()";
		KeepAliveMM();
	}

	~KeepAlive(){
		qDebug() << "KeepAlive::~KeepAlive()";
		KeepAliveDestructorMM();
	}

#ifdef __APPLE__ 
	void KeepAliveMM();
	void KeepAliveDestructorMM();
#else
	//?
#endif 

};

#endif