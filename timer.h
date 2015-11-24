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