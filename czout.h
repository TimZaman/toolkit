/**
--------------------------------------------------------------------------------
-   Module      :   czout.h
-   Description :   A stdout wrapper with a mutex to print sentences entirely
-                   so that in a multithreaded environment sentences do not 
-                   mingle with each character.
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

#ifndef CZOUT_h
#define CZOUT_h


#include <mutex>
#include <string>
#include <iostream> //cout, cerr, cin


#define RESET   "\033[0m"
#define BOLD    "\033[1m"
#define BLACK   "\033[30m"
#define RED     "\033[31m"
#define GREEN   "\033[32m"
#define YELLOW  "\033[33m"
#define BLUE    "\033[34m"
#define MAGENTA "\033[35m"
#define CYAN    "\033[36m"
#define WHITE   "\033[37m"

#define BGBLACK   "\033[40m"
#define BGRED     "\033[41m"
#define BGGREEN   "\033[42m"
#define BGBROWN   "\033[43m"
#define BGBLUE    "\033[44m"
#define BGMAGENTA "\033[45m"
#define BGCYAN    "\033[46m"
#define BGWHITE   "\033[47m"

#define COLRESET "\033[0m"
#define endlr "\033[0m\r\n"
#define WARNING "\033[1m\033[31mWARNING:\033[0m "

/**
* Interface structure that can be used for printing with mutex locks so that one can make sure entire strings are 
* printed without streams intertwining instead of just the characters as in 'cout'.
*/

struct czout_printer {
	std::mutex mutex;

    /**
    * Prints a complete uninterruptable string without line ending
    */
    void print(std::string str){
        std::lock_guard<std::mutex> guard(mutex);
        std::cout << str;
    }

	/**
	* Prints a complete uninterruptable string with line ending
	*/
    void println(std::string str, int verbosity){
        std::lock_guard<std::mutex> guard(mutex);
        if (verbosity<=verbose){
	        std::cout << str << std::endl;
	    }
    }

    int verbose = 1; //Default verbosity of a call (0=very important)
};


/*! 
 * \brief czout
 *
 * Interface for clean verbose printing, need a printer to make this work, see czout_printer
 * There can be many czout classes hooked into one printer.
 * \author Tim Zaman
 *
 */
class czout {

	public:

		czout(){
			//Nothing. Make sure to set the printer though!
		}

		//Constructor that sets the printer immediatelly
		czout(czout_printer * my_pzout_printer){
			setPrinter(my_pzout_printer);
		}

		void setPrinter(czout_printer * my_pzout_printer){
			this->t_pzout_printer = my_pzout_printer;
		}

		int verbose = 1; //Default verbosity of a call (0=very important)


	
		//template<typename T> czout& operator << (T&& x)  {
		//	//TODO
		//	cout << RED << " TODO: " << x << COLRESET << endl;
		//	return *this;
		//}

		czout& operator << (std::string data)  {
			allStr += data;
			return *this;
		};

		czout& operator << (void * data)  {
			//allStr += std::to_string(data);
			allStr += " ptr <0xTODO!>";
			return *this;
		};

		czout& operator << (int data)  {
			allStr += std::to_string(data);
			return *this;
		};

		czout& operator << (unsigned int data)  {
			allStr += std::to_string(data);
			return *this;
		};

		czout& operator << (unsigned long data)  {
			allStr += std::to_string(data);
			return *this;
		};

		czout& operator << (float data)  {
			allStr += std::to_string(data);
			return *this;
		};

		czout& operator << (double data)  {
			allStr += std::to_string(data);
			return *this;
		};

		//Call to 'std::endl' or 'endl'
		czout& operator << (std::ostream&(*pManip)(std::ostream&)){
			//Write and flush it out
			//Construct the full string
			std::string line = color + prefix + COLRESET + allStr;
			t_pzout_printer->println(line, verbose);
			//And reset the data
			allStr = "";
		}

		void setPrefix(std::string myprefix, int id = -1){
			this->prefix = myprefix;
			if (id>=0){
				setColor(id);
			}
		}


		void setColor(int id){
			std::vector<std::string> VECCOLSTR; //TODO: this can be more static, global and made generally less cumbersome
			VECCOLSTR.resize(5);
			VECCOLSTR[0] = CYAN;
			VECCOLSTR[1] = GREEN;
			VECCOLSTR[2] = YELLOW;
			VECCOLSTR[3] = BLUE;
			VECCOLSTR[4] = MAGENTA;
			this->color = VECCOLSTR[id % VECCOLSTR.size()]; //Pick sequentially changing color per ID
		}

	private:
		czout_printer * t_pzout_printer = NULL; //The interface that actually prints
		std::string allStr; //All data is concatenated here
		std::string prefix = "[UNSET]"; //The thing that's written before each string
		std::string color = "";
};




//END CZOUT_h
#endif 