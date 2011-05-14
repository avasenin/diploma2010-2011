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
#ifndef _V3DENSE___0077A726_E6E3_58a1_C16D_CE436AB90500__H
#define _V3DENSE___0077A726_E6E3_58a1_C16D_CE436AB90500__H

#include <complex>
#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>

#include "prgkern/_prgconfig.h"
#include "prgkern/_type.h"
#include "prgkern/_math.h"
#include "prgkern/_random.h"
#include "prgkern/_dense.h"

#include "prgkern/_for.h"
#define bll boost::lambda

namespace prgkern
{

	template <unsigned N, typename T>
	class vdense_
	{
	#define _A(i) assert(_LT((unsigned)i, (unsigned)N))
	#define _assign bll::_1 = bll::_2

	public:

		typedef T value_type;
		typedef unsigned dimension_type;

		vdense_() {} // undefined vector, used for efficiency in implementation
		            // of the expressions such as v = v1 + v2
		vdense_(T s) { for_<0, N>::expand(bll::_1 = s, v_); }

	#define MACRO(_z, _n, _text)  _text[_n]=s##_n;
		vdense_(T s0, T s1)
		{ assert(_EQ(N, (unsigned)2)); BOOST_PP_REPEAT(2, MACRO, v_); }
		vdense_(T s0, T s1, T s2)
		{ assert(_EQ(N, (unsigned)3)); BOOST_PP_REPEAT(3, MACRO, v_); }
		vdense_(T s0, T s1, T s2, T s3)
		{ assert(_EQ(N, (unsigned)4)); BOOST_PP_REPEAT(4, MACRO, v_); }
	#undef MACRO

		template <typename S> vdense_ &operator=(const vdense_<N, S> &v)
		{ for_<0, N>::expand(_assign, v_, &v[0]); return *this; }

		vdense_ &operator=(T s)
		{ for_<0, N>::expand(bll::_1 = s, v_); return *this; }

		T &operator[](size_t i) { _A(i); return v_[i]; }
		const T &operator[](size_t i) const { _A(i); return v_[i]; }
			// return const T& to get direct access to memory position of &v[0]

		T &operator()(size_t i) { _A(i); return v_[i]; }
		T operator()(size_t i) const { _A(i); return v_[i]; }

		void clear() { for_<0, N>::expand(bll::_1 = T(), v_);  }

	#define IMPLEMENT_SCALAR_OPERATOR(op) \
		template <typename __T> \
		vdense_ &operator op(__T s) { \
			for_<0, N>::expand(bll::_1 op s, v_); \
			return *this; \
		}
		IMPLEMENT_SCALAR_OPERATOR(+=)
		IMPLEMENT_SCALAR_OPERATOR(-=)
		IMPLEMENT_SCALAR_OPERATOR(*=)
	#undef IMPLEMENT_SCALAR_OPERATOR

	#define IMPLEMENT_VECTOR_OPERATOR(op) \
		template <typename __T> \
		vdense_ &operator op(const vdense_<N, __T> &s) { \
			for_<0, N>::expand(bll::_1 op bll::_2, v_, &s[0]); \
			return *this; \
		}
		IMPLEMENT_VECTOR_OPERATOR(+=)
		IMPLEMENT_VECTOR_OPERATOR(-=)
	#undef IMPLEMENT_VECTOR_OPERATOR

		T length2() const
		{ T s = 0; for_<0, N>::expand(bll::var(s) += bll::_1 * bll::_1, v_); return s; }

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

			T q = new_length/s; for_<0, N>::expand(bll::_1 *= q, v_);
			return s; // return previous length of vector
		}

		unsigned size() const { return N; }

		dimension_type dimension() const { return N; }
		const T *data_() const { return &v_[0]; }
		T *data_() { return &v_[0]; }

	protected:
		T v_[N];

	#undef _A
	};

	/**
	* @note negate operator
	*/
	template <unsigned N, typename T>
	INLINE vdense_<N, T> operator-(const vdense_<N, T> &v)
	{
		vdense_<N, T> tmp;
		for_<0, N>::expand(bll::_1 = -bll::_2, &tmp.v_[0], &v[0]);
		return tmp;
	}

	template <unsigned N, typename T>
	INLINE T distance2(const vdense_<N, T> &v1, const vdense_<N, T> &v2)
	{
		T s = 0;
		for_<0, N>::expand(bll::var(s) += (bll::_1 - bll::_2) * (bll::_1 - bll::_2),
			&v1[0], &v2[0]);
		return s;
	}
	#define _choose(op) \
		bll::if_then_else(bll::_2 op bll::_3, bll::_1 = bll::_2, bll::_1 = bll::_3)

	/**
	* @note minimize vector V on place such way that
	*       all coordinates <= coordinates of two vectors V1 & V2
	* @param v - new vector
	* @param v1,v2 - arguments
	*/
	template <unsigned N, typename T>
	INLINE void minimize(vdense_<N, T> &r, const vdense_<N, T> &a1, const vdense_<N, T> &a2)
	{ for_<0, N>::expand(_choose(<), &r[0], &a1[0], &a2[0]); }

	/**
	* @note maximize vector V on place such way that
	*       all coordinates >= coordinates of two vectors V1 & V2
	* @param v - new vector
	* @param v1,v2 - arguments
	*/
	template <unsigned N, typename T>
	INLINE void maximize(vdense_<N, T> &r, const vdense_<N, T> &a1, const vdense_<N, T> &a2)
	{ for_<0, N>::expand(_choose(>), &r[0], &a1[0], &a2[0]); }

	#undef choose

	#define IMPLEMENT_SCALAR_OPERATOR(op) \
		template <unsigned N, typename T, typename S> \
		INLINE vdense_<N, T> operator op(const vdense_<N, T> &v, S s) { \
			vdense_<N, T> tmp; \
			for_<0, N>::expand(bll::_1 = bll::_2 op s, &tmp[0], &v[0]); \
			return tmp; \
		} \
		template <unsigned N, typename T, typename S> \
		INLINE vdense_<N, T> operator op(S s, const vdense_<N, T> &v) \
		{ return v op s; }
	IMPLEMENT_SCALAR_OPERATOR(+) // V = V1 + s and V = s + V1
	IMPLEMENT_SCALAR_OPERATOR(-) // V = V1 - s ..
	IMPLEMENT_SCALAR_OPERATOR(*) // V = V1 * s ..
	#undef IMPLEMENT_SCALAR_OPERATOR

	#define IMPLEMENT_VECTOR_OPERATOR(op) \
	template <unsigned N, typename T, typename S> \
	INLINE vdense_<N, typename maxtype<T,S>::type> operator op(const vdense_<N, T> &v1, const vdense_<N, S> &v2) \
	{ \
		vdense_<N, typename maxtype<T,S>::type> tmp; \
		for_<0, N>::expand(bll::_1 = bll::_2 op bll::_3, &tmp[0], &v1[0], &v2[0]); \
		return tmp; \
	}
	IMPLEMENT_VECTOR_OPERATOR(+)
	IMPLEMENT_VECTOR_OPERATOR(-)
	IMPLEMENT_VECTOR_OPERATOR(*)
	#undef IMPLEMENT_VECTOR_OPERATOR

	/// max abs value of vector proection
	template <unsigned N, typename T>
	INLINE T max_proection(const vdense_<N, T> &v)
	{
		T s = fabs(v[0]);
		for_<1, N>::expand(bll::if_then(bll::var(s) < bll::_1, bll::var(s) = bll::_1), &v[1]);
		return s;
	}

	template <unsigned N, typename T>
	INLINE vdense_<N, T> min(const vdense_<N, T> &v1, const vdense_<N, T> &v2)
	{
		vdense_<N, T> tmp(v1);
		for_<0, N>::expand(bll::if_then(bll::_1 > bll::_2, bll::_1 = bll::_2), &tmp[0], &v2[0]);
		return tmp;
	}

	template <unsigned N, typename T>
	INLINE vdense_<N, T> max(const vdense_<N, T> &v1, const vdense_<N, T> &v2)
	{
		vdense_<N, T> tmp(v1);
		for_<0, N>::expand(bll::if_then(bll::_1 < bll::_2, bll::_1 = bll::_2), &tmp[0], &v2[0]);
		return tmp;
	}

	/// copy operation to avoid problems with alignment in memory
	template <unsigned N, typename T>
	INLINE void memcpy(T *to, const vdense_<N, T> &from)
	{
		for_<0, N>::expand(bll::_1 = bll::_2, &to[0], &from[0]);
	}

	/// copy operation to avoid problems with alignment in memory
	template <unsigned N, typename T>
	INLINE void memcpy(vdense_<N, T> &to, const T *from)
	{
		for_<0, N>::expand(bll::_1 = bll::_2, &to[0], &from[0]);
	}

	/// vector algebra operations s = (v1, v2)
	template <unsigned N, typename T1, typename T2>
	INLINE typename maxtype<T1, T2>::type
	scalar_product(const vdense_<N, T1> &v1, const vdense_<N, T2> &v2)
	{
		typename maxtype<T1, T2>::type s = 0.;
		for_<0, N>::expand(bll::var(s) += bll::_1 * bll::_2, &v1[0], &v2[0]);
		return s;
	}

	template <typename T>
	INLINE bool operator<(const vdense_<2, T> &v1, const vdense_<2, T> &v2)
	{
		if (v1[0] >= v2[0]) return false;
		if (v1[1] >= v2[1]) return false;
		return true;
	}

	template <typename T>
	INLINE bool operator<(const vdense_<3, T> &v1, const vdense_<3, T> &v2)
	{
		if (v1[0] >= v2[0]) return false;
		if (v1[1] >= v2[1]) return false;
		if (v1[2] >= v2[2]) return false;
		return true;
	}

	template <typename T>
	INLINE bool operator<(const vdense_<4, T> &v1, const vdense_<4, T> &v2)
	{
		if (v1[0] >= v2[0]) return false;
		if (v1[1] >= v2[1]) return false;
		if (v1[2] >= v2[2]) return false;
		if (v1[3] >= v2[3]) return false;
		return true;
	}

	template <typename T>
	INLINE bool operator<=(const vdense_<2, T> &v1, const vdense_<2, T> &v2)
	{
		if (v1[0] > v2[0]) return false;
		if (v1[1] > v2[1]) return false;
		return true;
	}

	template <typename T>
	INLINE bool operator<=(const vdense_<3, T> &v1, const vdense_<3, T> &v2)
	{
		if (v1[0] > v2[0]) return false;
		if (v1[1] > v2[1]) return false;
		if (v1[2] > v2[2]) return false;
		return true;
	}

	template <typename T>
	INLINE bool operator<=(const vdense_<4, T> &v1, const vdense_<4, T> &v2)
	{
		if (v1[0] > v2[0]) return false;
		if (v1[1] > v2[1]) return false;
		if (v1[2] > v2[2]) return false;
		if (v1[3] > v2[3]) return false;
		return true;
	}

	template <typename T>
	INLINE void swap(vdense_<2, T> &v1, vdense_<2, T> &v2)
	{
		std::swap(v1[0], v2[0]);
		std::swap(v1[1], v2[1]);
	}

	template <typename T>
	INLINE void swap(vdense_<3, T> &v1, vdense_<3, T> &v2)
	{
		std::swap(v1[0], v2[0]);
		std::swap(v1[1], v2[1]);
		std::swap(v1[2], v2[2]);
	}

	template <typename T>
	INLINE void swap(vdense_<4, T> &v1, vdense_<4, T> &v2)
	{
		std::swap(v1[0], v2[0]);
		std::swap(v1[1], v2[1]);
		std::swap(v1[2], v2[2]);
		std::swap(v1[3], v2[3]);
	}

	template <unsigned N, typename T>
	INLINE T distance1(const vdense_<N, T> &v1, const vdense_<N, T> &v2)
	{ return ::sqrt(distance2(v1, v2)); }

	/**
	* @note construct new vector V with all coordinates <= than
	*       coordinates of two vectors V1 & V2
	* @param v - new vector
	* @param v1,v2 - arguments
	*/
	template <unsigned N, typename T>
	INLINE vdense_<N, T> minimize(const vdense_<N, T> &a1, const vdense_<N, T> &a2)
	{ vdense_<N, T> tmp; minimize(tmp, a1, a2); return tmp; }

	/**
	* @note construct new vector V with all coordinates >= than
	*       coordinates of two vectors V1 & V2
	* @param v - new vector
	* @param v1,v2 - arguments
	*/
	template <unsigned N, typename T>
	INLINE vdense_<N, T> maximize(const vdense_<N, T> &a1, const vdense_<N, T> &a2)
	{ vdense_<N, T> tmp; maximize(tmp, a1, a2); return tmp; }

	template <typename T>
	INLINE bool equal(const vdense_<3, T> &v1, const vdense_<3, T> &v2,
		T accuacy=std::numeric_limits<T>::epsilon())
	{
		bool s0 = equal(v1[0], v2[0], accuacy);
		bool s1 = equal(v1[1], v2[1], accuacy);
		bool s2 = equal(v1[2], v2[2], accuacy);
		if (s0 && s1 && s2) return true;
		return false;
	}

	template <typename T>
	INLINE std::string make_string(const vdense_<2, T> &v)
	{
		char line[128];
		if (fabs(v[0]) + fabs(v[1]) > 1e3)
			sprintf(line, "[%10.3e%10.3e]", v[0], v[1]);
		else sprintf(line, "[%6.3f %6.3f]", (float)v[0], (float)v[1]);
		return std::string(line);
	}

	template <typename T>
	INLINE std::string make_string(const vdense_<3, T> &v)
	{
		char line[128];
		if (fabs(v[0]) + fabs(v[1]) + fabs(v[2]) > 1e3)
			sprintf(line, "[%10.3e%10.3e%10.3e]", v[0], v[1], v[2]);
		else sprintf(line, "[%6.3f %6.3f %6.3f]", (float)v[0], (float)v[1], (float)v[2]);
		return std::string(line);
	}

	template <typename T>
	INLINE std::string make_string(const vdense_<4, T> &v)
	{
		char line[128];
		if (fabs(v[0]) + fabs(v[1]) + fabs(v[2]) + fabs(v[3]) > 1e3)
			sprintf(line, "[%10.3e%10.3e%10.3e%10.3e]", v[0], v[1], v[2], v[3]);
		else sprintf(line, "[%6.3f %6.3f %6.3f %6.3f]", (float)v[0], (float)v[1],
			(float)v[2], (float)v[3]);
		return std::string(line);
	}

	INLINE std::string make_string(const vdense_<2, int> &v)
	{
		char line[128];
		sprintf(line, "[%d %d]", v[0], v[1]);
		return std::string(line);
	}

	INLINE std::string make_string(const vdense_<3, int> &v)
	{
		char line[128];
		sprintf(line, "[%d %d %d]", v[0], v[1], v[2]);
		return std::string(line);
	}

	INLINE std::string make_string(const vdense_<4, int> &v)
	{
		char line[128];
		sprintf(line, "[%d %d %d %d]", v[0], v[1], v[2], v[3]);
		return std::string(line);
	}

	template <unsigned N, typename T>
	INLINE std::ostream &operator<<(std::ostream &os, const vdense_<N, T> &v)
	{ os << make_string(v); return os; }

	/// vector algebra operations s = [v1 x v2]
	template <typename T1, typename T2, typename T3>
	INLINE void vector_product(vdense_<3, T1> &r,
		const vdense_<3, T2> &v1, const vdense_<3, T3> &v2)
	{
		r[0] = v1[1]*v2[2] - v1[2]*v2[1];
		r[1] = v1[2]*v2[0] - v1[0]*v2[2];
		r[2] = v1[0]*v2[1] - v1[1]*v2[0];
	}

	/// vector algebra operations s = [v1 x v2]
	template <typename T1, typename T2>
	INLINE vdense_<3, typename maxtype<T1, T2>::type>
	vector_product(const vdense_<3, T1> &v1, const vdense_<3, T2> &v2)
	{
		vdense_<3, typename maxtype<T1, T2>::type> s;
		vector_product(s, v1, v2);
		return s;
	}

	template <unsigned N, typename T>
	INLINE T cos_angle(const vdense_<N, T> &v1, const vdense_<N, T> &v2)
	{
		T s = v1.length2() * v2.length2();
		if (s <= 0.)
		{
			std::string msg = _S("[ERROR] zero vector length");
			PRINT_ERR(msg);
				// it's an invalid argument for next operation
		}
		return scalar_product(v1, v2) / sqrt(s);
	}

	template <typename T>
	INLINE T sin_angle(const vdense_<3, T> &v1, const vdense_<3, T> &v2)
	{
		T s = v1.length2() * v2.length2();
		if (s <= 0.)
		{
			std::string msg = _S("[ERROR] zero vector length");
			PRINT_ERR(msg);
				// it's an invalid argument for next operation
		}

		vdense_<3, T> tmp = vector_product(v1, v2);
		return sqrt(tmp.length2() / s);
	}

	template <unsigned N, typename T>
	INLINE T get_angle(const vdense_<N, T> &a, const vdense_<N, T> &b)
	{
		T phi = cos_angle(a, b);
		return safe_acos(phi);
	}

	/**
	* @brief makes random ortogonal vector to given vector
	*/
	template <typename T>
	INLINE void ortogonal(vdense_<3, T> &ort, const vdense_<3, T> &v,
		T accuracy=std::numeric_limits<T>::epsilon())
	{
		vdense_<3, T> ranx(ran0<T>(), ran0<T>(), ran0<T>());

		//<linear dependence test>
		int max_regenerate_count = 2; // max regenerate count
		while (max_regenerate_count--
			&& fabs(ranx[0] * v[1] - v[0] * ranx[1]) < accuracy
			&& fabs(ranx[1] * v[2] - v[1] * ranx[2]) < accuracy
		)
		{
			// the using of x/x0 == y/y0 == z/z0
			// produces x*y0*z0 == y*x0*z0 == z*x0*y0
			// and then x*y0 == y*x0 && y*z0 == z*y0
			ranx = vdense_<3, T>(ran0<T>(), ran0<T>(), ran0<T>()); // regenerate
		}
		if (max_regenerate_count < 0)
		{
			std::string msg = _S("[ERROR] bad generation of random vector");
			PRINT_ERR(msg);
		}
		//</linear dependence test>

		vdense_<3, T> ran__ = v * (scalar_product(ranx, v) / v.length2());
		ort = ranx - ran__;
	}

	template <typename T>
	INLINE vdense_<3, T> ortogonal(const vdense_<3, T> &v)
	{
		vdense_<3, T> ort;
		ortogonal(ort, v);
		return ort;
	}

	/**
	*  Генерит случайный вектор в заданном диапазоне.
	*/
	template <unsigned N, typename T>
	INLINE vdense_<N,T> ran0(const vdense_<N,T> &max,
		const vdense_<N,T> &min=vdense_<N,T>(0.))
	{
		vdense_<N,T> v;
		for (unsigned i=0; i<N; i++)
			v[i] = ran0(max[i], min[i]);
		return v;
	}

	template <unsigned N, typename T> struct extended<vdense_<N, T> >
	{ typedef vdense_<N, typename extended<T>::type > type; };

}
#undef bll
#endif
