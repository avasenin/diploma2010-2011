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
#ifndef _ASSERT__09267106_F08A_58bd_BC7E_BF415E780A02__H
#define _ASSERT__09267106_F08A_58bd_BC7E_BF415E780A02__H
#include "prgkern/_prgconfig.h"

#include <cassert>
#include <iostream>
#include "prgkern/_string.h"

#define STATIC_ASSERT(expr) { char _[(expr)? 1 : 0]; _[0] = 0; }

namespace prgkern
{

	const std::string ASSERT_MSG("[ASSERTION FAILED] : ");

	#define NO_IMPLEMENTATION	\
	{ \
		std::string msg = _S("%[ERROR] No implemantation : ") + __LOCATE__; \
		PRINT_ERR(msg); \
	}

	#ifndef _DEBUG
		#define ASSERT_IMPL(op, t1, t2) { return true; }
	#else
		#define ASSERT_IMPL(op, t1, t2) \
		{ \
			if (!(t1 op t2)) \
			{	\
				std::string msg = ASSERT_MSG + _S(make_string(t1)) \
					+ _S(__STRING(op)) + _S(make_string(t2)); \
				PRINT_MESSAGE(msg); \
				return false; \
			} \
			return true; \
		}
	#endif

	#undef _EQ
	template <typename T1, typename T2> INLINE bool _EQ(T1 t1, T2 t2) ASSERT_IMPL(==, t1, t2)

	#undef _NE
	template <typename T1, typename T2> INLINE bool _NE(T1 t1, T2 t2) ASSERT_IMPL(!=, t1, t2)

	#undef _GT
	template <typename T1, typename T2> INLINE bool _GT(T1 t1, T2 t2) ASSERT_IMPL(>, t1, t2)

	#undef _GE
	template <typename T1, typename T2> INLINE bool _GE(T1 t1, T2 t2) ASSERT_IMPL(>=, t1, t2)

	#undef _LT
	template <typename T1, typename T2> INLINE bool _LT(T1 t1, T2 t2) ASSERT_IMPL(<, t1, t2)

	#undef _LE
	template <typename T1, typename T2>  INLINE bool _LE(T1 t1, T2 t2) ASSERT_IMPL(<=, t1, t2)

}
#endif
