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
#ifndef _IEEE754__F9ED1116_D518_53ef_5E7E_DE43FDAB0900__H
#define _IEEE754__F9ED1116_D518_53ef_5E7E_DE43FDAB0900__H

#include "prgkern/_prgconfig.h"
#include "prgkern/_assert.h"

namespace prgkern
{

	/**
	* @brief bit format(IEEE754) of T type
	*/
	template <typename T> struct IEEE754_traits {};

	/**
	* @brief bit format of double
	* @note
	* [63      62..52        51..0] bits of double
	* sign  exponent+1023 mantissa (without first 1)
	*/
	template <> struct IEEE754_traits<double>
	{
		typedef uint64_t allocator_type;

		static const int binary_code_digits       = 64;
		static const int exponent_shift           = 52;
		static const int exponent_bias            = 1023;
		static const allocator_type sign_mask     = 0x8000000000000000LL;
		static const allocator_type exponent_mask = 0x7ff0000000000000LL;
		static const allocator_type mantissa_mask = 0x000fffffffffffffLL;
		static const allocator_type mantissa_unit = 0x0010000000000000LL;
	};

	/**
	* @brief bit format of float
	* @note
	* [31      30..23        22..0] bits of float
	* sign  exponent+127  mantissa (without first 1)
	*/
	template <> struct IEEE754_traits<float>
	{
		typedef unsigned allocator_type;

		static const int binary_code_digits       = 32;
		static const int exponent_shift           = 23;
		static const int exponent_bias            = 127;
		static const allocator_type sign_mask     = 0x80000000L;
		static const allocator_type exponent_mask = 0x7f800000L;
		static const allocator_type mantissa_mask = 0x007fffffL;
		static const allocator_type mantissa_unit = 0x00800000L;
	};

	/**
	* @brief Получение двоичного кода вещественного числа, находящегося в области (0, 1). Код ноля
	* 	всегда равен нулю, но данная функция дает неверное число для нуля.
	* @note необходимо помнить о длине получаемого кода, поскольку нули в начале кода ЗНАЧИМЫ
	* @param s - real (float, double ..)
	* @return code as integer (unsigned int, int64t ..)
	*/
	template <typename I, typename S>
	INLINE I binary_code(S s, int len=CHAR_BIT*sizeof(I))
	{
		assert(_GT(s, (S)0) && _LT(s, (S)1));
		assert(_GE((unsigned)(CHAR_BIT*sizeof(I)), (unsigned)len));

		typedef typename IEEE754_traits<S>::allocator_type T;
		const int exponent_shift = IEEE754_traits<S>::exponent_shift;
		const T exponent_bias    = IEEE754_traits<S>::exponent_bias;
		const T exponent_mask    = IEEE754_traits<S>::exponent_mask;
		const T mantissa_mask    = IEEE754_traits<S>::mantissa_mask;
		const T mantissa_unit    = IEEE754_traits<S>::mantissa_unit;

		int exp_ = (int)(((*(T*)&s & exponent_mask) >> exponent_shift)
			- exponent_bias);
		T mantissa_ = (*(T*)&s & mantissa_mask) + mantissa_unit;
		int rsh = len - exponent_shift + exp_;

		return (I)(rsh >= 0 ? mantissa_ << rsh : mantissa_ >> -rsh);
	}

}
#endif

