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
#ifndef _OS__1C505E08_33B0_4526_8B7F_E892B34B0A53__H
#define _OS__1C505E08_33B0_4526_8B7F_E892B34B0A53__H

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#define fs boost::filesystem

#include <fstream>

#include "prgkern/_prgconfig.h"
#include "prgkern/_string.h"

namespace prgkern
{

	const unsigned _1K   = 0x00000400;
	const unsigned _2K   = 0x00000800;
	const unsigned _4K   = 0x00001000;
	const unsigned _8K   = 0x00002000;
	const unsigned _16K  = 0x00004000;
	const unsigned _32K  = 0x00008000;
	const unsigned _64K  = 0x00010000;
	const unsigned _128K = 0x00020000;
	const unsigned _256K = 0x00040000;
	const unsigned _512K = 0x00080000;
	const unsigned _1M   = 0x00100000;
	const unsigned _2M   = 0x00200000;

	/// default max length of any string which is read from textual file
	const std::streamsize DEFAULT_GETLINE_LEN = 120;

	INLINE std::string make_dirname(const std::string &name)
	{

	#ifdef USE_LINUX_SPECIFIC_CODE
		if (name.empty()) return _S("./");
		if (name[name.size() - 1] == '/') return name;
		return name + _S("/");
	#endif

	#ifdef USE_WINDOWS_SPECIFIC_CODE
		if (name.empty()) return _S(".\\");
		if (name[name.size() - 1] == '\\') return name;
		return name + _S("\\");
	#endif
	}

	#define MAKE_DIRNAME(name) make_dirname(name)

	/** @brief safety reading of file line
	* @note allow to mix the WINDOWS and UNIX text file in same projects
	*/
	INLINE void getline(std::istream &file, char* buffer,
		std::streamsize num = DEFAULT_GETLINE_LEN)
	{

	#ifdef USE_LINUX_SPECIFIC_CODE
		file.getline(buffer, num);
	#endif
	#ifdef USE_WINDOWS_SPECIFIC_CODE
		file.getline(buffer, num, 0x0A);
	#endif

		int pos = ::strlen(buffer) - 1;
		if (pos > 0 && buffer[pos] == 0x0D) buffer[pos] = '\0';
	}

	INLINE std::string extension(const std::string &s)
	{ return fs::extension(fs::path(s)); }

	INLINE std::string basename(const std::string &s)
	{ return fs::basename(fs::path(s)); }

	class directory_iterator : public fs::directory_iterator
	{
		typedef fs::directory_iterator _base;
	public:

		/// iterator begin()
		directory_iterator(const std::string &dirname)
		: _base(fs::path(dirname)) {}

		/// iterator end()
		directory_iterator() : _base() {}

		/// iterator of extraction of file name
//		std::string operator*() { return (*this)->leaf(); }
		std::string operator*() { return ((*(_base*)this)->path()).filename();/*.generic_string();*/ }

		/// iterator of moving to next file name (directories will bw skipped)
		directory_iterator &operator++()
		{
			do { _base::operator++(); }
			while ( (*this) != _base() && fs::is_directory(**(_base*)this) );
			return *this;
		}

		bool operator==(const directory_iterator &it) const
		{ return (const fs::directory_iterator&)(*this) == (const fs::directory_iterator&)it; }

		bool operator!=(const directory_iterator &it) const
		{ return !(*this == it); }
	};
}

#undef fs
#endif
