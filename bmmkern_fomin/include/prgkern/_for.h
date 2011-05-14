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
#ifndef _FOR__0077A726_AC88_55c0_4A44_CE4344760601__H
#define _FOR__0077A726_AC88_55c0_4A44_CE4344760601__H

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/if.hpp>

namespace prgkern
{

	template <int N1, int N2>
	struct for_
	{
		template <class F>
		static void expand(F f)
		{ f(); for_<N1+1, N2>::expand(f); }

		template <class F, typename T1>
		static void expand(F f, T1 *s)
		{ f(s[N1]); for_<N1+1, N2>::expand(f, s); }

		template <class F, typename T1, typename T2>
		static void expand(F f, T1 *s, T2 *p)
		{ f(s[N1], p[N1]); for_<N1+1, N2>::expand(f, s, p); }

		template <class F, typename T1, typename T2, typename T3>
		static void expand(F f, T1 *s, T2 *p, T3 *q)
		{ f(s[N1], p[N1], q[N1]); for_<N1+1, N2>::expand(f, s, p, q); }
	};

	template <int N2>
	struct for_<N2, N2>
	{
		template <class F>
		static void expand(F) {}

		template <class F, typename T1>
		static void expand(F, T1 *) {}

		template <class F, typename T1, typename T2>
		static void expand(F, T1*, T2*) {}

		template <class F, typename T1, typename T2, typename T3>
		static void expand(F, T1*, T2*, T3*) {}
	};

}
#endif
