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
#ifndef _DENSE__0077A726_4743_5b7f_B566_CE43E5B50700__H
#define _DENSE__0077A726_4743_5b7f_B566_CE43E5B50700__H
#include "prgkern/_prgconfig.h"
#include "prgkern/_math.h"
#include "prgkern/_index.h"
#include "prgkern/_blas1.h"
#include <vector>

#define IMPLEMENT_OPERATOR_I1(R, A1, op, A2) \
	INLINE R operator op(const A1 &a1, const A2 &a2) \
	{ R tmp(a1); tmp op##= a2; return tmp; }

#define IMPLEMENT_OPERATOR_I2(R, A1, op, A2) \
	INLINE R operator op(const A1 &a1, const A2 &a2) \
	{ R tmp(a2); tmp op##= a1; return tmp; }

namespace prgkern
{

	enum
	{
		UNIT_ = 1,
		ORTOGONAL_UNIT_ = 2
	};

	/// count of elements in matrix row to output
	const size_t MATRIX_DEBUG_PRINT_COUNT = 10;
	const unsigned UNLIMITED_ = 0; // unlimited dimention of vector or matrix

	template <unsigned N, typename T>
	class basic_dense_ : public std::vector<T>
	{
	#define _ASSERT(i) assert(_LT((unsigned)i, (unsigned)_Base::size()))

		basic_dense_ &operator=(const basic_dense_ &s);

	protected:

		typedef std::vector<T>                  _Base;
		typedef index_<N, unsigned>             _Index;
		typedef T                               _Value;
		typedef typename _Base::iterator        iterator;
		typedef typename _Base::const_iterator  const_iterator;

		using _Base::size;
		using _Base::resize;
		using _Base::assign;
		using _Base::push_back;
		using _Base::erase;

		basic_dense_() {}
		void clear() {
			unsigned n = size();
			T *x = &(*this)[0];
			VECTOR_EXPRESSION_0(T, n, x, =, (T)0.);
		}

		const T& operator[](unsigned i) const { return (*(_Base*)this)[i]; }
		T& operator[](unsigned i) { return (*(_Base*)this)[i]; }

		unsigned resize_(const _Index &n)
		{
			unsigned sz = 1;
			dim_ = n;
			for (unsigned i=0; i<N; i++) sz *= n[i];
			_Base::resize(sz, T());
			return sz;
		}

		_Index dim_; // верхние границы массивов

	#undef _ASSERT
	};

#define INTERNAL_DENSE_TYPEDEFS(N, T, _Base) \
	typedef basic_dense_<N, T>      _Base; \
	typedef typename _Base::_Index  _Index; \
	typedef typename _Base::_Value  _Value; \

#define EXTERNAL_DENSE_TYPEDEFS \
	typedef _Index     index_type; \
	typedef _Value     value_type; \
	typedef typename _Base::iterator        iterator; \
	typedef typename _Base::const_iterator  const_iterator; \

	template <unsigned N, typename T> class vdense_;
		// if you want have large (>4) dimension than use N=UNLIMITED

	template <unsigned N, unsigned M, typename T>
	class mdense_;
		// if you want have float dimensions than use N,M=UNLIMITED

}
#endif
