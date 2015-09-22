/* Name :       utils*.h [part of tim zaman's utility functions]
 * Version:     v1.0
 * Date :       2012-05-07
 * Developer :  Tim Zaman, Pixelprisma BV (timbobel@gmail.com)
 * Distribution in any form is prohibited.
 */

#ifndef UTILS_JSONCPP_TIM_H
#define UTILS_JSONCPP_TIM_H


#include <iostream>
#include <string>
#include <clocale>

#include <sys/types.h>
#include <sys/stat.h>

#include <unistd.h>
#include <pwd.h>

#include <iomanip>
#include <cstdlib>

#include <json/json.h>
#include <json/writer.h>
#include <json/reader.h>

#include <algorithm>


namespace util{
	template <typename Iterable>
	Json::Value iterable2json(Iterable const& cont);
	void PrintJSONValue( Json::Value val );
	bool PrintJSONTree( Json::Value &root, unsigned short depth /* = 0 */);
	
};


#endif





