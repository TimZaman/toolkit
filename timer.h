/**
--------------------------------------------------------------------------------
-   Module      :   timer.h
-   Description :   A nice timer class.
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

#ifndef TIM_TIMER_h
#define TIM_TIMER_h

#pragma once

#include <iostream>
#include <chrono>
#include <ctime>

/*! 
 * \brief Timer
 *
 * Timer function that can be used to benchmark code speed.
 * \author Tim Zaman
 *
 */

class timer {

	public:
		timer(){
			this->start = std::chrono::system_clock::now();
		}

		void tic(){
			this->start = std::chrono::system_clock::now();
		}

		double toc(){
			end = std::chrono::system_clock::now();
			std::chrono::duration<double> elapsed_seconds = end-start;
 			return double (elapsed_seconds);
		}

		void toc(std::string strLabel){
			std::cout << strLabel << " took " << toc() << " ms." << std::endl;
		}

	private:
		 std::chrono::time_point<std::chrono::system_clock> start, end;


}