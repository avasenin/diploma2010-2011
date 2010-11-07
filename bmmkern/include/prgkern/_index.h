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
#ifndef _INDEX__0077A726_7DF9_5da9_5082_D74308840100__H
#define _INDEX__0077A726_7DF9_5da9_5082_D74308840100__H

#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include "prgkern/_for.h"
#define bll boost::lambda

#include "prgkern/_prgconfig.h"
#include "prgkern/_string.h"
#include "prgkern/_assert.h"

namespace prgkern
{

	/// invalid value (used as array index)
	const int nill = -1;

	/**
	* @brief complex index_ to data container
	* @param N - number of indexes
	*/
	template <unsigned N, typename T=int>
	class index_
	{
	#define _A(i)  assert(_LT((unsigned)i, N))

	public:

		index_() { clear(); }

		index_(T i0) { clear(); ndx_[0] = i0; }

		index_(T i0, T i1)
		{ clear(); ndx_[0] = i0; ndx_[1] = i1; }

		index_(T i0, T i1, T i2)
		{ clear(); ndx_[0] = i0; ndx_[1] = i1; ndx_[2] = i2; }

		index_(T i0, T i1, T i2, T i3)
		{ clear(); ndx_[0] = i0; ndx_[1] = i1; ndx_[2] = i2; ndx_[3] = i3; }

		index_(T i0, T i1, T i2, T i3, T i4)
		{ clear(); ndx_[0] = i0; ndx_[1] = i1; ndx_[2] = i2; ndx_[3] = i3; ndx_[4] = i4; }

		index_(T i0, T i1, T i2, T i3, T i4, T i5)
		{
			clear(); ndx_[0] = i0; ndx_[1] = i1; ndx_[2] = i2; ndx_[3] = i3;
			ndx_[4] = i4; ndx_[5] = i5;
		}

		index_(T i0, T i1, T i2, T i3, T i4, T i5, T i6)
		{
			clear(); ndx_[0] = i0; ndx_[1] = i1; ndx_[2] = i2; ndx_[3] = i3;
			ndx_[4] = i4; ndx_[5] = i5; ndx_[6] = i6;
		}

		index_(const index_ &ndx) { ::memcpy(this, &ndx, sizeof(ndx)); }

		index_ &operator=(const index_ &ndx)
		{ ::memcpy(this, &ndx, sizeof(ndx)); return *this; }

		T operator[](unsigned i) const { _A(i); return ndx_[i]; }
		T &operator[](unsigned i) { _A(i); return ndx_[i]; }

		bool operator==(const index_ &ndx) const
		{ return ::memcmp(this, &ndx, sizeof(ndx)); }
		bool operator!=(const index_ &ndx) const { return !(*this==ndx); }

		void clear() { ::memset(&ndx_[0], 0, sizeof(ndx_)); }

	protected:

		T ndx_[N];

	#undef _A
	};

	/*
	* @note this way produces bad ordering with negative indexes
	* template <unsigned N, typename T>
	* INLINE bool operator<(const index_<N, T> &ndx1, const index_<N, T> &ndx2)
	* { return ::memcmp(&ndx1, &ndx2, sizeof(index_<N, T>)) < 0; }
	*/

	template <typename T>
	INLINE bool operator<(const index_<1, T> &ndx1, const index_<1, T> &ndx2)
	{ return ndx1[0] < ndx2[0]; }

	template <typename T>
	INLINE bool operator<(const index_<2, T> &ndx1, const index_<2, T> &ndx2)
	{
		if (ndx1[0] < ndx2[0]) return true;
		if (ndx1[0] > ndx2[0]) return false;
		return (ndx1[1] < ndx2[1]);
	}

	template <typename T>
	INLINE bool operator<(const index_<3, T> &ndx1, const index_<3, T> &ndx2)
	{
		if (ndx1[0] < ndx2[0]) return true;
		if (ndx1[0] > ndx2[0]) return false;
		if (ndx1[1] < ndx2[1]) return true;
		if (ndx1[1] > ndx2[1]) return false;
		return (ndx1[2] < ndx2[2]);
	}

	template <typename T>
	INLINE bool operator<(const index_<4, T> &ndx1, const index_<4, T> &ndx2)
	{
		if (ndx1[0] < ndx2[0]) return true;
		if (ndx1[0] > ndx2[0]) return false;
		if (ndx1[1] < ndx2[1]) return true;
		if (ndx1[1] > ndx2[1]) return false;
		if (ndx1[2] < ndx2[2]) return true;
		if (ndx1[2] > ndx2[2]) return false;
		return (ndx1[3] < ndx2[3]);
	}

	template <typename T>
	INLINE void swap(index_<2,T> &ndx)
	{ std::swap(ndx[0], ndx[1]); }

	template <typename T>
	INLINE void swap(index_<3,T> &ndx)
	{ std::swap(ndx[0], ndx[2]); }

	template <typename T>
	INLINE void swap(index_<4,T> &ndx)
	{
		std::swap(ndx[0], ndx[3]);
		std::swap(ndx[1], ndx[2]);
	}

	template <typename T>
	INLINE void swap(index_<2,T> &n1, index_<2,T> &n2)
	{
		std::swap(n1[0], n2[0]);
		std::swap(n1[1], n2[1]);
	}

	template <typename T>
	INLINE void swap(index_<3,T> &n1, index_<3,T> &n2)
	{
		std::swap(n1[0], n2[0]);
		std::swap(n1[1], n2[1]);
		std::swap(n1[2], n2[2]);
	}

	template <typename T>
	INLINE void swap(index_<4,T> &n1, index_<4,T> &n2)
	{
		std::swap(n1[0], n2[0]);
		std::swap(n1[1], n2[1]);
		std::swap(n1[2], n2[2]);
		std::swap(n1[3], n2[3]);
	}

	template <typename T>
	INLINE index_<1,T> mod(const index_<1, T> &n, const index_<1, T> &k)
	{ return index_<1,T>(mod(n[0], k[0])); }

	template <typename T>
	INLINE index_<2,T> mod(const index_<2, T> &n, const index_<2, T> &k)
	{ return index_<2,T>(mod(n[0], k[0]), mod(n[1], k[1])); }

	template <typename T>
	INLINE index_<3,T> mod(const index_<3, T> &n, const index_<3, T> &k)
	{ return index_<3,T>(mod(n[0], k[0]), mod(n[1], k[1]), mod(n[2], k[2])); }

	template <typename T>
	INLINE index_<4,T> mod(const index_<4, T> &n, const index_<4, T> &k)
	{ return index_<4,T>(mod(n[0], k[0]), mod(n[1], k[1]), mod(n[2], k[2]), mod(n[3], k[3])); }

	template <typename T>
	INLINE index_<1,T> operator-(const index_<1, T> &ndx)
	{ return index_<1,T>(-ndx[0]); }

	template <typename T>
	INLINE index_<2,T> operator-(const index_<2, T> &ndx)
	{ return index_<2,T>(-ndx[0], -ndx[1]); }

	template <typename T>
	INLINE index_<3,T> operator-(const index_<3, T> &ndx)
	{ return index_<3,T>(-ndx[0], -ndx[1], -ndx[2]); }

	template <typename T>
	INLINE index_<4,T> operator-(const index_<4, T> &ndx)
	{ return index_<4,T>(-ndx[0], -ndx[1], -ndx[2], -ndx[3]); }

	template <typename T, typename S>
	INLINE index_<1,T> operator*(const index_<1, T> &ndx, S s)
	{ return index_<1,T>((T)(ndx[0]*s)); }

	template <typename T, typename S>
	INLINE index_<2,T> operator*(const index_<2, T> &ndx, S s)
	{ return index_<2,T>((T)(ndx[0]*s), (T)(ndx[1]*s)); }

	template <typename T, typename S>
	INLINE index_<3,T> operator*(const index_<3, T> &ndx, S s)
	{ return index_<3,T>((T)(ndx[0]*s), (T)(ndx[1]*s), (T)(ndx[2]*s)); }

	template <typename T, typename S>
	INLINE index_<4,T> operator*(const index_<4, T> &ndx, S s)
	{ return index_<4,T>((T)(ndx[0]*s), (T)(ndx[1]*s), (T)(ndx[2]*s), (T)(ndx[3]*s)); }

#define IMPLEMENT_OPERATOR(op) \
	template <typename T> \
	INLINE index_<1, T> operator op(const index_<1, T> &ndx1, const index_<1, T> &ndx2) \
	{ return index_<1, T>(ndx1[0] op ndx2[0]); }

	IMPLEMENT_OPERATOR(+)
	IMPLEMENT_OPERATOR(-)

#undef IMPLEMENT_OPERATOR

#define IMPLEMENT_OPERATOR(op) \
	template <typename T> \
	INLINE index_<2, T> operator op(const index_<2, T> &ndx1, const index_<2, T> &ndx2) \
	{ return index_<2, T>(ndx1[0] op ndx2[0], ndx1[1] op ndx2[1]); }

	IMPLEMENT_OPERATOR(+)
	IMPLEMENT_OPERATOR(-)

#undef IMPLEMENT_OPERATOR

#define IMPLEMENT_OPERATOR(op) \
	template <typename T> \
	INLINE index_<3, T> operator op(const index_<3, T> &ndx1, const index_<3, T> &ndx2) \
	{ return index_<3, T>(ndx1[0] op ndx2[0], ndx1[1] op ndx2[1], ndx1[2] op ndx2[2]); }

	IMPLEMENT_OPERATOR(+)
	IMPLEMENT_OPERATOR(-)

#undef IMPLEMENT_OPERATOR

#define IMPLEMENT_OPERATOR(op) \
	template <typename T> \
	INLINE index_<4, T> operator op(const index_<4, T> &ndx1, const index_<4, T> &ndx2) \
	{ return index_<4, T>(ndx1[0] op ndx2[0], ndx1[1] op ndx2[1], ndx1[2] op ndx2[2], ndx1[3] op ndx2[3]); }

	IMPLEMENT_OPERATOR(+)
	IMPLEMENT_OPERATOR(-)

#undef IMPLEMENT_OPERATOR

	/**
	* Возвращает смещение в памяти элемента задаваемого многомерным индексом.
	* Массив раполагается в памяти так, что последний индекс меняется быстрее
	* остальных. При расчете предполагается то, что элементы массива располагаются
	* друг за другом без промежутков. Если это не так, то для других вариантов
	* нужно умножить результат на расстояние между элементами.
	* @param ndx - многомерный целочисленный индекс
	* @param sz - размерности многомерного массива
	* @return смещение в памяти
	*/
	template <typename _Integer>
	INLINE _Integer memory_offset(const index_<1, _Integer> &ndx,
		const index_<1, _Integer> &sz) { return ndx[0]; }

	/**
	* см. комментарий для memory_offset(const index_<1, _Integer> &ndx, ..
	*/
	template <typename _Integer>
	INLINE _Integer memory_offset(const index_<2, _Integer> &ndx,
		const index_<2, _Integer> &sz)
	{ return ndx[0] * sz[1] + ndx[1]; }

	/**
	* см. комментарий для memory_offset(const index_<1, _Integer> &ndx, ..
	*/
	template <typename _Integer>
	INLINE _Integer memory_offset(const index_<3, _Integer> &ndx,
		const index_<3, _Integer> &sz)
	{ return ( ndx[0] * sz[1] + ndx[1] ) * sz[2] + ndx[2]; }

	/**
	* см. комментарий для memory_offset(const index_<1, _Integer> &ndx, ..
	*/
	template <typename _Integer>
	INLINE _Integer memory_offset(const index_<4, _Integer> &ndx,
		const index_<4, _Integer> &sz)
	{ return ( ( ndx[0] * sz[1] + ndx[1] ) * sz[2] + ndx[2] ) * sz[3] + ndx[3]; }

	/**
	* Возвращает многомерный индекс по смещению в памяти элемента.
	* Массив раполагается в памяти так, что последний индекс меняется быстрее
	* остальных. При расчете предполагается то, что элементы массива располагаются
	* друг за другом без промежутков. Если это не так, то для других вариантов
	* нужно умножить результат на расстояние между элементами.
	*
	* Учитывает цикличность, то есть смещение может выйти за границы диапозона.
	*
	* @param n - смещение в памяти
	* @param sz - размерности многомерного массива
	* @return многомерный индекс
	*/
	template <typename _Integer>
	INLINE index_<1, _Integer> memory_offset(unsigned n, const index_<1, _Integer> &sz)
	{ return index_<1, _Integer>(n % sz[0]); }

	/**
	* см. комментарий для memory_offset(unsigned n, const index_<1, _Integer> &ndx, ..
	*/
	template <typename _Integer>
	INLINE index_<2, _Integer> memory_offset(unsigned n, const index_<2, _Integer> &sz)
	{ return index_<2, _Integer>((n / sz[1]) % sz[0], n % sz[1]); }

	/**
	* см. комментарий для memory_offset(unsigned n, const index_<1, _Integer> &ndx, ..
	*/
	template <typename _Integer>
	INLINE index_<3, _Integer> memory_offset(unsigned n, const index_<3, _Integer> &sz)
	{
		unsigned n__ = n / sz[2];
		return index_<3, _Integer>((n__ / sz[1]) % sz[0], n__ % sz[1], n % sz[2]);
	}

	/**
	* см. комментарий для memory_offset(unsigned n, const index_<1, _Integer> &ndx, ..
	*/
	template <typename _Integer>
	INLINE index_<4, _Integer> memory_offset(unsigned n, const index_<4, _Integer> &sz)
	{
		unsigned n__ = n / sz[3];
		unsigned m__ = n__ / sz[2];
		return index_<4, _Integer>((m__ / sz[1]) % sz[0], m__ % sz[1], n__ % sz[2], n % sz[3]);
	}

	template <unsigned N, typename T>
	INLINE T scalar_product(const index_<N, T> &v1, const index_<N, T> &v2)
	{
		T s = 0;
		for (unsigned i=0; i<N; i++) s += v1[i] * v2[i];
		return s;
	}

	template <typename T>
	INLINE T scalar_product(const index_<1, T> &v1, const index_<1, T> &v2)
	{ return v1[0] * v2[0]; }

	template <typename T>
	INLINE T scalar_product(const index_<2, T> &v1, const index_<2, T> &v2)
	{ return v1[0] * v2[0] + v1[1] * v2[1]; }

	template <typename T>
	INLINE T scalar_product(const index_<3, T> &v1, const index_<3, T> &v2)
	{ return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]; }

	template <typename T>
	INLINE T scalar_product(const index_<4, T> &v1, const index_<4, T> &v2)
	{ return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3]; }

	template <typename T>
	INLINE T component_product(const index_<1, T> &v) { return v[0]; }

	template <typename T>
	INLINE T component_product(const index_<2, T> &v) { return v[0] * v[1]; }

	template <typename T>
	INLINE T component_product(const index_<3, T> &v) { return v[0] * v[1] * v[2]; }

	template <typename T>
	INLINE T component_product(const index_<4, T> &v)
	{ return v[0] * v[1] * v[2] * v[3]; }

	template <unsigned N, typename T>
	INLINE std::string make_string(const index_<N, T> &v)
	{
		std::ostringstream oss_convert;
		oss_convert << "\'" << v[0];
		for (unsigned i=1; i<N; i++) oss_convert << ", " << v[i];
		oss_convert << "\'";
		return oss_convert.str();
	}

	template <unsigned N, typename T>
	INLINE std::ostream &operator<<(std::ostream &os, const index_<N, T> &v)
	{ os << make_string(v); return os; }


	/**
	* @brief index_ with ordering components
	* @param N - dimension of index_
	* @param T - type of index_ component
	*/
	template <unsigned N, typename T=int>
	struct sym_index_ : public index_<N, T>
	{
		typedef index_<N, T>  _Base;

		sym_index_() {}
		sym_index_(T a0)
		{
			_Base::operator[](0) = a0;
		}
		sym_index_(T a0, T a1)
		{
			if (a1 < a0) swap(a0, a1);
			_Base::operator[](0) = a0;
			_Base::operator[](1) = a1;
		}
		sym_index_(T a0, T a1, T a2)
		{
			if (a2 < a0) swap(a0, a2);
			_Base::operator[](0) = a0;
			_Base::operator[](1) = a1;
			_Base::operator[](2) = a2;
		}
		sym_index_(T a0, T a1, T a2, T a3)
		{
			if (a3 < a0 || (a3 == a0 && a2 < a1)) { swap(a0, a3); swap(a1, a2); }
			_Base::operator[](0) = a0;
			_Base::operator[](1) = a1;
			_Base::operator[](2) = a2;
			_Base::operator[](3) = a3;
		}
	};

	/**
	* @brief cast multidimentional index to offset in linear memory
	* @param I1 - start value of index
	* @param I2 - afterlast value of index
	* @return index of one dimentional array
	*/
	template <unsigned N> class index_cast_
	{
	protected:

		typedef index_<1, int>  _Index1;
		typedef index_<2, int>  _Index2;
		typedef index_<3, int>  _Index3;
		typedef index_<N, int>  _Index;

		_Index b_; ///< нижняя граница области (по умалчиванию 0)
		_Index s_; ///< размеры области

	public:
		typedef _Index  multi_index_type;

		index_cast_() {}
		index_cast_(const _Index &up) : b_(0), s_(up) {}
		index_cast_(const _Index &bottom, const _Index &up) : b_(bottom), s_(up - bottom) {}
		index_cast_(const index_cast_ &c) : b_(c.b_), s_(c.s_) {}
		index_cast_ &operator=(const index_cast_ &c)
		{ b_ = c.b_; s_ = c.s_; return *this; }

		unsigned operator[](const index_<N, int> &ndx) const
		{ return make_single_index_(ndx - b_); }

		_Index operator[](unsigned n) const
		{ return make_multi_index_(s_, n) + b_; }
			// use dummy argument "s_" for choose overload function

		unsigned size() const { return component_product(s_); }
		_Index dimension() const { return s_; }
		_Index bottom() const { return b_; }
		_Index top() const { return b_ + s_; }
		_Index center() const { return b_ + s_ * 0.5; }

	private:

		unsigned make_single_index_(const _Index1 &ndx) const
		{
			assert(_LT((unsigned)ndx[0], (unsigned)s_[0]));
			return ndx[0];
		}

		unsigned make_single_index_(const _Index2 &ndx) const
		{
			assert(_LT((unsigned)ndx[0], (unsigned)s_[0]));
			assert(_LT((unsigned)ndx[1], (unsigned)s_[1]));
			return ndx[0] * s_[1] + ndx[1];
		}

		unsigned make_single_index_(const _Index3 &ndx) const
		{
			assert(_LT((unsigned)ndx[0], (unsigned)s_[0]));
			assert(_LT((unsigned)ndx[1], (unsigned)s_[1]));
			assert(_LT((unsigned)ndx[2], (unsigned)s_[2]));
			return ( ndx[0] * s_[1] + ndx[1] ) * s_[2] + ndx[2];
		}

		_Index1 make_multi_index_(const _Index1&, unsigned n) const
		{ return _Index1(n); }

		_Index2 make_multi_index_(const _Index2&, unsigned n) const
		{ return _Index2(n/s_[1], n%s_[1]); }

		_Index3 make_multi_index_(const _Index3&, unsigned n) const
		{
			unsigned n__ = n / s_[2];
			return _Index3(n__ / s_[1], n__ % s_[1], n % s_[2]);
		}
	};

	/**
	* @brief converts the multidimentional index to offset in linear memory by cycle method
	* @param I1 - start value of index
	* @param I2 - afterlast value of index
	* @return index of one dimentional array
	*/
	template <unsigned N> class cycle_index_cast : public index_cast_<N>
	{
		typedef index_cast_<N>  _Base;
		typedef typename _Base::_Index1  _Index1;
		typedef typename _Base::_Index2  _Index2;
		typedef typename _Base::_Index3  _Index3;
		typedef typename _Base::_Index   _Index;

	public:

		typedef _Index  multi_index_type;
		using _Base::operator[];
		using _Base::dimension;
		using _Base::size;
		using _Base::bottom;
		using _Base::top;
		using _Base::center;

		cycle_index_cast() {}
		cycle_index_cast(const _Index &up) : _Base(up) {}
		cycle_index_cast(const cycle_index_cast& c) : _Base(c) {}
		cycle_index_cast &operator=(const cycle_index_cast &c)
		{ _Base::operator=(c); return *this; }

		unsigned operator[](const _Index1 &ndx) const
		{
			_Index1 ndx__(mod(ndx[0], _Base::top()[0]));
			return _Base::operator[](ndx__);
		}

		unsigned operator[](const _Index2 &ndx) const
		{
			_Index2 ndx__(
				mod(ndx[0], _Base::top()[0]),
				mod(ndx[1], _Base::top()[1])
			);
			return _Base::operator[](ndx__);
		}

		unsigned operator[](const _Index3 &ndx) const
		{
			_Index3 ndx__(
				mod(ndx[0], _Base::top()[0]),
				mod(ndx[1], _Base::top()[1]),
				mod(ndx[2], _Base::top()[2])
			);
			return _Base::operator[](ndx__);
		}
	};

#undef bll
}
#endif
