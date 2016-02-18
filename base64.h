/**
--------------------------------------------------------------------------------
-   Module      :   base64.h
-   Description :   Just another Base64 encoding/decoding class
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

#ifndef BASE64_TIM_H
#define BASE64_TIM_H

#include <iostream>
#include <stdint.h>
#include <stdlib.h>
#include <string>

class b64 {
	public:
		b64();
		~b64();
		std::string base64_encode(unsigned char const* , unsigned int len);
		std::string base64_encode(std::string);
		std::string base64_decode(char const* , unsigned int len);
		std::string base64_decode(std::string);
	private:
		char *decoding_table = NULL;

		void base64_cleanup();
		void build_decoding_table();

		char * base64_encode(const unsigned char *data, unsigned int input_length, unsigned int *output_length);
		unsigned char * base64_decode(const char *data, unsigned int input_length, unsigned int *output_length);

};

#endif