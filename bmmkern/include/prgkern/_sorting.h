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
#ifndef _SORTING__F9ED1116_D518_53ef_5E7E_DE43FDAB0911__H
#define _SORTING__F9ED1116_D518_53ef_5E7E_DE43FDAB0911__H

#include <queue>

#include "prgkern/_prgconfig.h"
#include "prgkern/_assert.h"

namespace prgkern
{
	// __sync_fetch_and_add( &global_int, 1 );


	template <unsigned N, typename T>
	void radix_sorting(unsigned n, T *out, const T *in)
	{
		unsigned mask = 0x0001;
		std::queue<T> queue0, queue1;
		for (unsigned i=0; i<n; i++)
		{
			if (in[i] & mask) queue1.push(in[i]); else queue0.push(in[i]);
		}

		for (unsigned iter=1; iter<N; iter++)
		{
			mask <<= 1;
			unsigned n0 = queue0.size();
			unsigned n1 = queue1.size();
			for (unsigned i=0; i<n0; i++)
			{
				T s = queue0.front(); queue0.pop();
				if (s & mask) queue1.push(s); else queue0.push(s);
			}
			for (unsigned i=0; i<n1; i++)
			{
				T s = queue1.front(); queue1.pop();
				if (s & mask) queue1.push(s); else queue0.push(s);
			}
		}

		for (unsigned i=0; i<n; i++)
		{
			while (!queue0.empty()) { *out++ = queue0.front(); queue0.pop(); }
			while (!queue1.empty()) { *out++ = queue1.front(); queue1.pop(); }
		}
	}


}
#endif

