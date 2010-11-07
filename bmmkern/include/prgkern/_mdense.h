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
#ifndef _MDENSE___0077A726_2963_5882_057A_CE43D46B0000__H
#define _MDENSE___0077A726_2963_5882_057A_CE43D46B0000__H
#include "prgkern/_prgconfig.h"

#include "prgkern/_type.h"
#include "prgkern/_math.h"
#include "prgkern/_dense.h"
#include "prgkern/_vdense.h"
#include "prgkern/_blas1.h"

#include <vector>

#define _T typename
namespace prgkern
{

	template <typename T>
	class mdense_<UNLIMITED_, UNLIMITED_, T> : protected basic_dense_<2, T>
	{
		INTERNAL_DENSE_TYPEDEFS(2, T, _Base);

	#define _ASSERT(i, j) \
		assert(_LT((unsigned)i, (unsigned)_Base::dim_[0])); \
		assert(_LT((unsigned)j, (unsigned)_Base::dim_[1]))

		mdense_(const mdense_ &m);
		mdense_ &operator=(const mdense_ &m);

	public:

		EXTERNAL_DENSE_TYPEDEFS;
		using _Base::size;
		using _Base::erase;
		using _Base::clear;

		mdense_(unsigned n1=0, unsigned n2=0) { resize(n1, n2); }
		unsigned resize(unsigned n1, unsigned n2) { return _Base::resize_(_Index(n1, n2)); }

		const T& operator()(unsigned i, unsigned j) const { return (*this)[_Index(i,j)]; }
		T& operator()(unsigned i, unsigned j) { return (*this)[_Index(i,j)]; }

		const T& operator[](const _Index &ndx) const
		{
			unsigned i = ndx[0] * _Base::dim_[1] + ndx[1];
			return _Base::operator[](i);
		}

		T& operator[](const _Index &ndx)
		{
			unsigned i = ndx[0] * _Base::dim_[1] + ndx[1];
			return _Base::operator[](i);
		}

		_Index dimension() const { return _Base::dim_; }

	#define ASSERT_DIMENSIONS \
		if ( _Base::dim_[0] != _Base::dim_[1] ) \
		{ \
			std::string msg = _S("[ERROR] bad matrix dimensions"); \
			PRINT_ERR(msg); \
		}

		template <_T S> mdense_ &operator=(S s)
		{
			ASSERT_DIMENSIONS;
			unsigned n = _Base::dim_[0] * _Base::dim_[1];
			T *x = &(*this)[0];
			VECTOR_EXPRESSION_0(T, n, x, =, (T)0.); // обнулили матрицу

			n = _Base::dim_[1];
			for (unsigned i=0; i<n; i++, x+=n+1) *x = (T)s;
			return *this;
		}

		template <_T S> mdense_ &operator+=(S s)
		{
			ASSERT_DIMENSIONS;
			int n = _Base::dim_[1]; T *x = &(*this)[0];
			for (unsigned i=0; i<n; i++, x+=n+1) *x += s; // забили диагональ
			return *this;
		}

		template <_T S> mdense_ &operator-=(S s)
		{
			ASSERT_DIMENSIONS;
			int n = _Base::dim_[1]; T *x = &(*this)[0];
			for (unsigned i=0; i<n; i++, x+=n+1) *x -= s;
			return *this;
		}

		template <_T S> mdense_ &operator*=(S s)
		{
			ASSERT_DIMENSIONS;
			int n = _Base::dim_[0] * _Base::dim_[1];
			T *x = &(*this)[0];
			VECTOR_EXPRESSION_1(T, n, x, =, s, *, x);
			return *this;
		}

		template <_T S> mdense_ &operator/=(S s)
		{
			ASSERT_DIMENSIONS;
			unsigned n = _Base::dim_[0] * _Base::dim_[1];
			T *x = &(*this)[0]; T q = 1./s;
			VECTOR_EXPRESSION_1(T, n, x, =, q, *, x);
			return *this;
		}
	#undef ASSERT_DIMENSIONS
	};

	/**
	* @brief get U(-1) from upper triangle matrix
	* @note result matrix replaces the origin matrix
	*   can be used for small matrix 'cos the scaling is O(N**3)
	*/
	template <typename T>
	INLINE void invert(mdense_<UNLIMITED_, UNLIMITED_, T> &U)
	{
		typedef typename mdense_<UNLIMITED_, UNLIMITED_, T>::index_type  _Index;

		_Index ndx = U.dimension();
		assert(_EQ(ndx[0], ndx[1]));

		unsigned n = ndx[0];
		for (unsigned i=0; i<n; ++i)
		{
			U(i, i) = (T) 1. / U(i, i);
			for (unsigned j=i+1; j<n; ++j)
			{
				typename extended<T>::type s = 0.;
				for (unsigned k=i; k<j; ++k)
					s += U(i, k) * U(k, j);
				U(i, j) = (T)(-s / U(j, j));
			}
		}
	}

}

#undef _T
#endif
