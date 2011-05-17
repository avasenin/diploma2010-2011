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

	class no_match_exception {};

	class RegEx
	{
	public:

		typedef std::pair<std::string, std::string>  match_type;

		RegEx(const std::string &mask) : compiled_regex_(mask) {}

		/*
		 * Проверяет наличие заданного образца в строке s
		 */
		bool is_match(const std::string &s) const
		{ return boost::regex_match(s, compiled_regex_); }

		/**
		 *
		 * @param results совокупность подстрок, удовлетворящих шаблонам подвыражений
		 * @param s строка для анализа
		 * @return
		 */
		bool get_match(std::vector<std::string> &results, const std::string &s) const
		{
			boost::match_results<std::string::const_iterator> what;
			if (boost::regex_match(s, what, compiled_regex_))
			{
				unsigned cnt = what.size() - 1;
				results.resize(cnt);
				for (unsigned i=0; i<cnt; i++) results[i].assign(what[i + 1].first, what[i + 1].second);
				return true;
			}
			return false;
		}

		std::string get_match(unsigned n, const std::string &s) const
		{
			boost::match_results<std::string::const_iterator> what;
			if (boost::regex_match(s, what, compiled_regex_))
			{
				assert(_LE(n, what.size()));
				return _S(what[n + 1].first, what[n + 1].second);
			}
			return _S("");
		}

		std::string get_match(const std::string &s) const
		{
			boost::match_results<std::string::const_iterator> what;
			if (boost::regex_match(s, what, compiled_regex_))
			{
				return _S(what[1].first, what[1].second);
			}
			return _S("");
		}

		std::pair<std::string, std::string> get_match2(const std::string &s) const
		{
			typedef std::pair<std::string, std::string> _Pair;

			boost::match_results<std::string::const_iterator> what;
			if (boost::regex_match(s, what, compiled_regex_))
			{
				return _Pair(_S(what[1].first, what[1].second), _S(what[2].first, what[2].second));
			}
			return _Pair(_S(""), _S(""));
		}

	private:
		boost::regex compiled_regex_;
	};

}
#endif

