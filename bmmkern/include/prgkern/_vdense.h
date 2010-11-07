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
#ifndef _VDENSE___0077A726_9715_5da8_D372_CE433D550D00__H
#define _VDENSE___0077A726_9715_5da8_D372_CE433D550D00__H
#include "prgkern/_prgconfig.h"

#include "prgkern/_type.h"
#include "prgkern/_math.h"
#include "prgkern/_dense.h"
#include "prgkern/_blas1.h"

namespace prgkern
{

	template <typename T>
	class vdense_<UNLIMITED_, T> : public basic_dense_<1, T>
	{
	#define _A(i) assert(_LT((unsigned)i, (unsigned)_Base::size()))
	#define _N (int)size()

		INTERNAL_DENSE_TYPEDEFS(1, T, _Base);

	public:
		EXTERNAL_DENSE_TYPEDEFS;
		using _Base::size;
		using _Base::resize;
		using _Base::assign;
		using _Base::push_back;
		using _Base::erase;
		using _Base::clear;

		vdense_() {}
		explicit vdense_(unsigned n) { resize(n, T()); }

		vdense_ &operator=(const vdense_ &s)
		{
			resize(s.size());
			T *x = &(*this)[0], *y = &s[0];
			VECTOR_EXPRESSION_1(T, _N, x, =, , , y);
			return *this;
		}

		vdense_(_I2T<UNIT_>, const vdense_ &dir) : _Base(dir.size())
		{
			T s = dir.length();
			assert(_GT(s, (T)0.));
			T *x = &(*this)[0]; T q = (T)1. / s;
			VECTOR_EXPRESSION_1(T, _N, x, =, q, *, x);
		}

		vdense_(_I2T<ORTOGONAL_UNIT_>, const vdense_ &f, const vdense_ &unit) : _Base(f)
		{
			// проверка на единичность, так как для эффективности алгоритм определен для таких векторов
			assert(_LT(fabs(unit.length2() - 1.), (T)(100. * std::numeric_limits<T>::epsilon())));
			T s = scalar_product(_N, &(*this)[0], &unit[0]);
			T *x = &(*this)[0], *y = &unit[0];
			VECTOR_EXPRESSION_1(T, _N, x, -=, s, *, y)
		}

		vdense_(unsigned n, const T *x) : _Base(&x[0], &x[n]) {}
		vdense_(const T *x1, const T *x2) : _Base(x1, x2) {}

		const T& operator[](unsigned i) const { return _Base::operator[](i); }
		T& operator[](unsigned i) { return _Base::operator[](i); }

		template <typename S> vdense_ &operator=(S s)
		{
			T *x = &(*this)[0];
			VECTOR_EXPRESSION_0(T, _N, x, =, (T)s);
			return *this;
		}

		template <typename S> vdense_ &operator+=(S s)
		{
			T *x = &(*this)[0];
			VECTOR_EXPRESSION_0(T, _N, x, +=, s)
			return *this;
		}

		template <typename S> vdense_ &operator-=(S s)
		{
			T *x = &(*this)[0];
			VECTOR_EXPRESSION_0(T, _N, x, -=, s)
			return *this;
		}

		template <typename S> vdense_ &operator*=(S s)
		{
			T *x = &(*this)[0];
			VECTOR_EXPRESSION_1(T, _N, x, =, (T)s, *, x);
			return *this;
		}

		template <typename S> vdense_ &operator/=(S s)
		{
			T *x = &(*this)[0]; S q = (S)1. / s;
			VECTOR_EXPRESSION_1(T, _N, x, =, q, *, x);
			return *this;
		}

		vdense_ &operator+=(const vdense_ &s)
		{
			T *x = &(*this)[0]; T *y = &s[0];
			VECTOR_EXPRESSION_1(T, _N, x, +=, , , y)
			return *this;
		}

		vdense_ &operator-=(const vdense_ &s)
		{
			T *x = &(*this)[0]; T *y = &s[0];
			VECTOR_EXPRESSION_1(T, _N, x, -=, , , y)
			return *this;
		}

		T length2() const { return prgkern::scalar_product(_N, &(*this)[0], &(*this)[0]); }
		T length() const { return sqrt(length2()); }

		/**
		* @brief normalize to new_length & return previous length of vector
		*/
		T normalize(T new_length=(T)1.)
		{
			T s = length();
			if (s <= 0)
			{
				std::string msg = _S("zero vector length");
				PRINT_ERR(msg);
			}
			T q = new_length/s;
			T *x = &(*this)[0];
			VECTOR_EXPRESSION_1(T, _N, x, =, q, *, x);
			return s;
		}

	#undef _N
	#undef _A
	};

	template <typename T>
	INLINE T scalar_product(const vdense_<UNLIMITED_, T> &v1,
		const vdense_<UNLIMITED_, T> &v2)
	{
		assert(_EQ(v1.size(), v2.size()));
		return prgkern::scalar_product((int)v1.size(), &v1[0], &v2[0]);
	}

	// Операторы излишне затратны (создание временного вектора), потому не должны использоваться.

	//	template <typename T>
	//	INLINE vdense_<UNLIMITED_, T> operator+(const vdense_<UNLIMITED_, T> &a1, T a2)
	//	{ vdense_<UNLIMITED_, T> tmp(a1); tmp += a2; return tmp; }
	//
	//	template <typename T>
	//	INLINE vdense_<UNLIMITED_, T> operator-(const vdense_<UNLIMITED_, T> &a1, T a2)
	//	{ vdense_<UNLIMITED_, T> tmp(a1); tmp -= a2; return tmp; }
	//
	//	template <typename T>
	//	INLINE vdense_<UNLIMITED_, T> operator*(const vdense_<UNLIMITED_, T> &a1, T a2)
	//	{ vdense_<UNLIMITED_, T> tmp(a1); tmp *= a2; return tmp; }
	//
	//	template <typename T>
	//	INLINE vdense_<UNLIMITED_, T> operator+(T a1, const vdense_<UNLIMITED_, T> &a2)
	//	{ vdense_<UNLIMITED_, T> tmp(a2); tmp += a1; return tmp; }
	//
	//	template <typename T>
	//	INLINE vdense_<UNLIMITED_, T> operator-(T a1, const vdense_<UNLIMITED_, T> &a2)
	//	{ vdense_<UNLIMITED_, T> tmp(a2); tmp -= a1; return tmp; }
	//
	//	template <typename T>
	//	INLINE vdense_<UNLIMITED_, T> operator*(T a1, const vdense_<UNLIMITED_, T> &a2)
	//	{ vdense_<UNLIMITED_, T> tmp(a2); tmp *= a1; return tmp; }
	//
	//	template <typename T>
	//	INLINE vdense_<UNLIMITED_, T> operator+(const vdense_<UNLIMITED_, T> &a1,
	//		const vdense_<UNLIMITED_, T> &a2)
	//	{ vdense_<UNLIMITED_, T> tmp(a1); tmp += a2; return tmp; }
	//
	//	template <typename T>
	//	INLINE vdense_<UNLIMITED_, T> operator-(const vdense_<UNLIMITED_, T> &a1,
	//		const vdense_<UNLIMITED_, T> &a2)
	//	{ vdense_<UNLIMITED_, T> tmp(a1); tmp -= a2; return tmp; }

	template <typename T>
	INLINE std::string make_string(const vdense_<UNLIMITED_, T> &v)
	{
		std::ostringstream oss_convert;
		oss_convert << "[ " ;
		for (unsigned i=0, sz=v.size(); i<sz; i++)
			oss_convert << v[i] << " ";
		oss_convert << "]";
		return oss_convert.str();
	}

	template <typename T>
	INLINE std::ostream &operator<<(std::ostream &os, const vdense_<UNLIMITED_, T> &v)
	{ os << make_string(v); return os; }

	template <typename T>
	INLINE bool equal(const vdense_<UNLIMITED_, T> &v1, const vdense_<UNLIMITED_, T> &v2,
		T accuacy=std::numeric_limits<T>::epsilon())
	{
		bool is_equal = true;
		unsigned n = v1.size();
		if (n!=v2.size())
		{
			std::string msg = _S("[ERROR] different dimensions of vectors");
			PRINT_ERR(msg);
		}
		for (unsigned i=0; i<n; ++i)
			if (!equal(v1[i], v2[i], accuacy))
			{
				std::string msg = _S("[WARNING] i=") + make_string(i) + _S(" ")
					+ make_string(v1[i]) + _S(" vs. ") + make_string(v2[i]);
				PRINT_MESSAGE(msg);
				is_equal = false;
			}

		return is_equal;
	}

	template <unsigned N> struct extended<vdense_<N, float> > { typedef vdense_<N, double> type; };
	template <unsigned N> struct extended<vdense_<N, double> > { typedef vdense_<N, double> type; };

}
#endif
