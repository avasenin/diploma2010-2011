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
#ifndef _PREPROCESSOR__0077A726_7A2B_5be2_F165_CE43ECC00900__H
#define _PREPROCESSOR__0077A726_7A2B_5be2_F165_CE43ECC00900__H

#define STATIC_DIMENSION(array, type)  sizeof(array)/sizeof(type)

/// стандартный размер выровненной структуры в памяти
#define DEFAULT_PADDING_SIZE  (16)
	// для некоторых архитектур (CELL) такое выравниваение должно делаться вручную
	// программистом

#define PADDED_STRUCT_START(name) \
	struct name##80F12B46_2FC7_5d01_0DD6_C54A41090B01 {

#define PADDED_STRUCT_END(name, size) }; \
	struct name \
	: public name##80F12B46_2FC7_5d01_0DD6_C54A41090B01 { \
		static const unsigned base_size_80F12B46_2FC7_5d01_0DD6_C54A41090B00 = \
			sizeof(name##80F12B46_2FC7_5d01_0DD6_C54A41090B01); \
		static const unsigned pad_size_80F12B46_3C0F_563c_3ED6_C54A897A0800 = \
			(base_size_80F12B46_2FC7_5d01_0DD6_C54A41090B00 % size) ? \
			(size - (base_size_80F12B46_2FC7_5d01_0DD6_C54A41090B00 % size)) : 0; \
		char pad80F12B46_5474_538c_7DD6_C54A71650700[pad_size_80F12B46_3C0F_563c_3ED6_C54A897A0800]; \
	};

#define PADDED_STRUCT(name) \
	struct name##80F12B46_2FC7_5d01_0DD6_C54A41090B01

#define INSERT_PADDING(name, size) }; \
	struct name \
	: public name##80F12B46_2FC7_5d01_0DD6_C54A41090B01 { \
		static const unsigned base_size_80F12B46_2FC7_5d01_0DD6_C54A41090B00 = \
			sizeof(name##80F12B46_2FC7_5d01_0DD6_C54A41090B01); \
		static const unsigned pad_size_80F12B46_3C0F_563c_3ED6_C54A897A0800 = \
			(base_size_80F12B46_2FC7_5d01_0DD6_C54A41090B00 % size) ? \
			(size - (base_size_80F12B46_2FC7_5d01_0DD6_C54A41090B00 % size)) : 0; \
		char pad80F12B46_5474_538c_7DD6_C54A71650700[pad_size_80F12B46_3C0F_563c_3ED6_C54A897A0800];

#define OPERATOR_EQ(class_, m_) \
	bool operator==(const class_ &_) const { return m_ == _.m_; }

#define OPERATOR_NE(class_, m_) \
	bool operator!=(const class_ &_) const { return !(*this == _); }

#define OPERATOR_LT(class_, m_) \
	bool operator<(const class_ &_) const { return m_ < _.m_; }

#define OPERATOR_PMEQ(class_, m_) \
	bool operator==(const class_ &_) const { return *m_ == *_.m_; }

#define OPERATOR_PMNE(class_, m_) \
	bool operator!=(const class_ &_) const { return !(*this == _); }

#define OPERATOR_PMLT(class_, m_) \
	bool operator<(const class_ &_) const { return *m_ < *_.m_; }

#define OPERATOR_EQ2(class_, m_, m2_) \
	bool operator==(const class_ &_) const { return m_ == _.m_ && m2_ == _.m2_; }

#define OPERATOR_NE2(class_, m_, m2_) \
	bool operator!=(const class_ &_) const { return !(*this == _); }

#define OPERATOR_LT2(class_, m_, m2_) \
	bool operator<(const class_ &_) const \
	{ return m_ != _.m_ ? m_ < _.m_ : m2_ < _.m2_; }

#define PRINT_LINE(len, symbol) \
	{ \
		for (unsigned i=0; i<len; i++) std::cout << symbol; \
		std::cout << std::endl; \
	} \

namespace prgkern
{
	/** @brief conversion of expression to value
	* @param expr expression
	* @return value of expression
	* @note the function is used with preprocessor defines
	*/
	template <typename T> inline T expr2value(T expr) { return expr; }

	template <int N> class Guard_
	{
		std::string msg_;
		static int count_;

	protected:
		std::string indent() const
		{
			_S indent = _S("");
			for (unsigned i=0; i<count_; i++) indent += _S(" ");
			return indent;
		}

	public:
		Guard_(const char *message) : msg_(message)
		{
			++count_;
			_S msg = indent() + _S("[ START  ] ") + msg_;
			PRINT_MESSAGE(msg);
		}
		~Guard_()
		{
			_S msg = indent() + _S("[ FINISH ] ") + msg_;
			PRINT_MESSAGE(msg);
			--count_;
		}

	};
	template <int N> int Guard_<N>::count_ = -1;
}

#ifdef _DEBUG
#define _DEBUG_GUARD(msg)  Guard_<0> guard(msg);
#else
#define _DEBUG_GUARD(msg)
#endif
#define DEBUG_GUARD(msg)

#endif
