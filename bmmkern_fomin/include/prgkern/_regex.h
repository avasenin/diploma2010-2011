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
#ifndef _REGEX__0077A726_56A4_5b9b_E3D3_D1451DA80100__H
#define _REGEX__0077A726_56A4_5b9b_E3D3_D1451DA80100__H
#include <boost/regex.hpp>

#include "prgkern/_prgconfig.h"
#include "prgkern/_string.h"

namespace prgkern
{

	class RegEx
	{
	public:
		RegEx(const std::string &mask) : compiled_regex_(mask) {}
		bool match(const std::string &s)
		{ return boost::regex_match(s, compiled_regex_); }
	private:
		boost::regex compiled_regex_;
	};

}
#endif

