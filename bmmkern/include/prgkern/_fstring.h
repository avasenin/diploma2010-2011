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
#ifndef _FSTRING__FC42EF6_3B55_5d9a_1B59_514158A70601__H
#define _FSTRING__FC42EF6_3B55_5d9a_1B59_514158A70601__H
#include "prgkern/_prgconfig.h"

#include <string>
#include <cctype>
#include <cstdlib>
#include "prgkern/_assert.h"
#include "prgkern/_string.h"

namespace prgkern
{

	/**
	* @brief short fast string optimized for assign & compare
	* Пробелы по краям строки незначимы, потому должны быть удалены при печати.
	*/
	template <typename A>
	class fstring_
	{
	#define VALUE(s)  (*(A*)&s)

	public:

		fstring_() { ::memset(&v_, ' ', sizeof(A)); }
		fstring_(const char *s) { init_(s); }
		fstring_(const std::string &s) { init_(s.c_str()); }
		fstring_(const fstring_ &s) { VALUE(v_) = VALUE(s.v_); }
		fstring_ &operator=(const char *s) { init_(s); return *this; }
		fstring_ &operator=(const std::string &s) { init_(s.c_str()); return *this; }
		fstring_ &operator=(const fstring_ &s) { VALUE(v_) = VALUE(s.v_); return *this; }
		char operator[](unsigned i) const { assert(_LT(i, (unsigned)sizeof(A))); return v_[i]; }
		char &operator[](unsigned i) { assert(_LT(i, (unsigned)sizeof(A))); return v_[i]; }
		unsigned size() const { return sizeof(A); }

	protected:

		void init_(const char *s)
		{
			assert(_LE(::strlen(s), sizeof(A)));
			::memset(&v_, ' ', sizeof(A));
			::memcpy(&v_, s, std::min(sizeof(A), ::strlen(s)));
		}

		char v_[sizeof(A)];

	#undef VALUE
	};
	typedef fstring_<int> fstring;

	INLINE bool operator<(fstring s1, fstring s2) { return *(int*)&s1 < *(int*)&s2; }
	INLINE bool operator>(fstring s1, fstring s2) { return *(int*)&s1 > *(int*)&s2; }
	INLINE bool operator==(fstring s1, fstring s2) { return *(int*)&s1 == *(int*)&s2; }
	INLINE bool operator!=(fstring s1, fstring s2) { return *(int*)&s1 != *(int*)&s2; }

	template <typename A>
	INLINE std::string make_string(const fstring_<A> &v)
	{
		int n = sizeof(A);
		std::ostringstream oss_convert;
		for (int i=0; i<n; i++) oss_convert << v[i];
		return trim(oss_convert.str());
	}

	template <typename A>
	INLINE std::ostream &operator<<(std::ostream &os, const fstring_<A> &v)
	{ return os << make_string(v); }

	template <typename A>
	INLINE void swap(fstring_<A> &s1, fstring_<A> &s2)
	{ fstring_<A> s = s1; s1 = s2; s2 = s; }

	INLINE void to_lower(fstring &s)
	{
		s[0] = (char)tolower(s[0]);
		s[1] = (char)tolower(s[1]);
		s[2] = (char)tolower(s[2]);
		s[3] = (char)tolower(s[3]);
	}

	INLINE void to_upper(fstring &s)
	{
		s[0] = (char)toupper(s[0]);
		s[1] = (char)toupper(s[1]);
		s[2] = (char)toupper(s[2]);
		s[3] = (char)toupper(s[3]);
	}

	INLINE fstring to_lower(const fstring &s)
	{
		fstring r;
		r[0] = (char)tolower(s[0]);
		r[1] = (char)tolower(s[1]);
		r[2] = (char)tolower(s[2]);
		r[3] = (char)tolower(s[3]);
		return r;
	}

	INLINE fstring to_upper(const fstring &s)
	{
		fstring r;
		r[0] = (char)toupper(s[0]);
		r[1] = (char)toupper(s[1]);
		r[2] = (char)toupper(s[2]);
		r[3] = (char)toupper(s[3]);
		return r;
	}

	INLINE std::string make_string(const char *format, const fstring &v)
	{
		char line[120];
		::sprintf(line, format, make_string(v).c_str());
		return std::string(line);
	}

}
#endif
