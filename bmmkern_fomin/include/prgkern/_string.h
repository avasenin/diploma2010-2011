/***************************************************************************
 *   Copyright (C) 2008 by Eduard S. Fomin                                 *
 *   fomin@bionet.nsc.ru                                                   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef __STRING__FC42EF6_3B55_5d9a_1B59_514158A70601__H
#define __STRING__FC42EF6_3B55_5d9a_1B59_514158A70601__H

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#define ba boost::algorithm

#include "prgkern/_prgconfig.h"

#include <string>
#include <sstream>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <complex>

namespace prgkern
{

	class prg_exception : public std::exception {};

	#undef _S
	#define _S   std::string

	#undef __STRING
	#define __STRING(expr) #expr

	#undef  __LOCATE__
	#define __LOCATE__  _S(__FILE__) +_S("(") + _S(itoa(__LINE__)) + _S(") ")

	#define PRINT_MESSAGE(msg) \
	{ \
		COUT_MUTEX_LOCK \
		std::cout << msg << std::endl; \
		COUT_MUTEX_UNLOCK \
	}

	#define PRINT_ERR(msg) \
	{ \
		COUT_MUTEX_LOCK \
		std::cout << "[ERROR] " << msg << std::endl; \
		COUT_MUTEX_UNLOCK \
		throw std::exception(); \
	}

	#define PRINT_BREAK(msg) \
	{ \
		COUT_MUTEX_LOCK \
		std::cout << msg << std::endl; \
		COUT_MUTEX_UNLOCK \
		throw prg_exception(); \
	}

	#ifdef USE_LINUX_SPECIFIC_CODE

		#define SYSTEM_ERR(operation, error) \
		{ \
			SYSTEM_MUTEX_LOCK; \
			operation; \
			int errno_ = error; \
			SYSTEM_MUTEX_UNLOCK; \
			if (errno_ != 0) \
			{ \
				std::string msg = _S("[SYSTEM ERROR] ") + _S(__FILE__) + _S("(") \
					+ _S(itoa(__LINE__)) + _S("): ") + _S(::strerror(errno_)); \
				PRINT_MESSAGE(msg); \
				throw std::exception(); \
			} \
		}

		#define SYSTEM_WARN(operation, error) \
		{ \
			SYSTEM_MUTEX_LOCK; \
			operation; \
			int errno_ = error; \
			SYSTEM_MUTEX_UNLOCK; \
			if (errno_ != 0) \
			{ \
				std::string msg = _S("[SYSTEM WARNING] ") + _S(__FILE__) \
					+ _S("(") + _S(itoa(__LINE__)) \
					+ _S("): ") + _S(::strerror(errno_)); \
				PRINT_MESSAGE(msg); \
			} \
		}

	#endif

	const int MAX_STRING_LEN  = 256; // length of text file strings (*.dat)
	const int MAX_LINE_LEN = 256; // length of text file strings (*.dat)

	template <typename T>
	INLINE std::string make_string(const char *format, T v)
	{
		char line[1200];
		::sprintf(line, format, v);
		return std::string(line);
	}

	template <typename T1, typename T2>
	INLINE std::string make_string(const char *format, T1 v1, T2 v2)
	{
		char line[1200];
		::sprintf(line, format, v1, v2);
		return std::string(line);
	}

	template <typename T1, typename T2, typename T3>
	INLINE std::string make_string(const char *format, T1 v1, T2 v2, T3 v3)
	{
		char line[1200];
		::sprintf(line, format, v1, v2, v3);
		return std::string(line);
	}

	template <typename T1, typename T2, typename T3, typename T4>
	INLINE std::string make_string(const char *format, T1 v1, T2 v2, T3 v3, T4 v4)
	{
		char line[1200];
		::sprintf(line, format, v1, v2, v3, v4);
		return std::string(line);
	}

	template <typename T1, typename T2, typename T3, typename T4, typename T5>
	INLINE std::string make_string(const char *format, T1 v1, T2 v2, T3 v3, T4 v4, T5 v5)
	{
		char line[1200];
		::sprintf(line, format, v1, v2, v3, v4, v5);
		return std::string(line);
	}

	template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
	INLINE std::string make_string(const char *format, T1 v1, T2 v2, T3 v3, T4 v4, T5 v5, T6 v6)
	{
		char line[1200];
		::sprintf(line, format, v1, v2, v3, v4, v5, v6);
		return std::string(line);
	}

	INLINE std::string make_string(const char *s) { return std::string(s); }
	INLINE std::string make_string(const std::string &s) { return s; }
	INLINE std::string make_string(bool s) { return s ? _S("true"): _S("false"); }

#define MAKE_STRING(type, format) \
	INLINE std::string make_string(type v) { \
		char line[120]; \
		::sprintf(line, format, v); \
		return std::string(line); \
	}

	MAKE_STRING(int,      "%d")
	MAKE_STRING(unsigned, "%u")
	MAKE_STRING(long, 		"%ld")
	MAKE_STRING(unsigned long, "%lu")
	MAKE_STRING(float,    "%e")
	MAKE_STRING(double,   "%le")
	MAKE_STRING(void*,    "%p")

#undef MAKE_STRING

#define MAKE_STRING(type, format) \
	INLINE std::string make_string(type v) { \
		char line[120]; \
		::sprintf(line, format, v.real(), v.imag()); \
		return std::string(line); \
	}
	MAKE_STRING(std::complex<double>, "[%11.3e %11.3e]")
	MAKE_STRING(std::complex<float>, "[%11.3e %11.3e]")
	MAKE_STRING(std::complex<int>, "[%d %d]")

#undef MAKE_STRING

	INLINE std::string make_string(const char *format, const std::string &v)
	{
		char line[1200];
		::sprintf(line, format, v.c_str());
		return std::string(line);
	}

	INLINE std::string itoa(int v)    { return make_string(v); }
	INLINE std::string ftoa(double v) { return make_string(v); }
	INLINE std::string ptoa(void *v)  { return make_string(v); }
	INLINE std::string ctoa(char v)   { return std::string(1, v); }

	//--------------------------------------------------------------------------
	// When parsing the input from a user, strings usually have unwanted leading
	// or trailing characters. To get rid of them, we need trim functions.
	//--------------------------------------------------------------------------
	INLINE void ltrim(std::string &s) { ba::trim_left(s); }
	INLINE void rtrim(std::string &s) { ba::trim_right(s); }
	INLINE void trim(std::string &s) { ba::trim(s); }
	INLINE std::string ltrim(const std::string &s) { return ba::trim_left_copy(s); }
	INLINE std::string rtrim(const std::string &s) { return ba::trim_right_copy(s); }
	INLINE std::string trim(const std::string &s) { return ba::trim_copy(s); }

	INLINE void trim(std::string &s, const char *subs)
	{ ba::trim_if(s, ba::is_any_of(subs)); }
	INLINE std::string trim(const std::string &s, const char *subs)
	{ return ba::trim_copy_if(s, ba::is_any_of(subs)); }


	//--------------------------------------------------------------------------
	// STL has a nice way of converting character case. Unfortunately, it works
	// only for a single character and we want to convert a string then ...
	//--------------------------------------------------------------------------
	INLINE void to_lower(std::string &s) { ba::to_lower(s); }
	INLINE void to_upper(std::string &s) { ba::to_upper(s); }
	INLINE std::string to_lower(const std::string &s) { return ba::to_lower_copy(s); }
	INLINE std::string to_upper(const std::string &s) { return ba::to_upper_copy(s); }

	/// return true in case of the string contains the space symbols only
	INLINE bool isspace(const std::string &s) { return ba::trim_copy(s).empty(); }

}
#undef ba
#endif
