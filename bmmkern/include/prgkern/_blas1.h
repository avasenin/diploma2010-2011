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
#ifndef _BLAS1__0077A726_4743_5b7f_B566_CE43E5B50700__H
#define _BLAS1__0077A726_4743_5b7f_B566_CE43E5B50700__H
#include <boost/preprocessor/iteration/local.hpp>

#include "prgkern/_prgconfig.h"
#include "prgkern/_math.h"
#include "prgkern/_type.h"
#include "prgkern/_debug.h"
#include "prgkern/_sse.h"

namespace prgkern
{

	#define SCALAR_EXPRESSION_1(T, n, s, op, a1, o1, x1)                                                   \
		{                                                                                                \
			typedef vecreal_<4, T> vreal_t;                                                                \
			const unsigned N__ = CYCLE_EXPANSION_SIZE;                                                     \
			unsigned n__ = n % N__;                                                                        \
			for (unsigned i=0; i<n__; i++) s op a1 o1 x1[i];                                               \
			vreal_t result = 0.f;                                                                          \
			const T *x1__ = x1 + n__;                                                                      \
			for (unsigned i=n__; i<n; i+=N__, x__+=N__, x1__+=N__)                                         \
				result op a1 o1 vreal_t(x1__);                                                               \
			s op summarize(result);                                                                        \
		}

	#define VECTOR_EXPRESSION_0(T, n, x, op, a1)                                                           \
		{                                                                                                \
			typedef vecreal_<4, T> vreal_t;                                                                \
			const unsigned N__ = CYCLE_EXPANSION_SIZE;                                                     \
			unsigned n__ = n % N__;                                                                        \
			for (unsigned i=0; i<n__; i++) x[i] op a1;                                                     \
			vreal_t result; T *x__ = x + n__;                                                              \
			for (unsigned i=n__; i<n; i+=N__, x__+=N__)                                                    \
			{ result op a1; result.store(x__); }                                                           \
		}

	#define VECTOR_EXPRESSION_1(T, n, x, op, a1, o1, x1)                                               \
		{                                                                                                \
			typedef vecreal_<4, T> vreal_t;                                                                \
			const unsigned N__ = CYCLE_EXPANSION_SIZE;                                                     \
			unsigned n__ = n % N__;                                                                        \
			for (unsigned i=0; i<n__; i++) x[i] op a1 o1 x1[i];                                            \
			vreal_t result; T *x__ = x + n__;                                                              \
			const T *x1__ = x1 + n__;                                                                      \
			for (unsigned i=n__; i<n; i+=N__, x__+=N__, x1__+=N__)                                         \
			{ result op a1 o1 vreal_t(x1__); result.store(x__); }                                          \
		}

	#define VECTOR_EXPRESSION_2(T, n, x, op, a1, o1, x1, op1, a2, o2, x2)                              \
		{                                                                                                \
			typedef vecreal_<4, T> vreal_t;                                                                \
			const unsigned N__ = CYCLE_EXPANSION_SIZE;                                                     \
			unsigned n__ = n % N__;                                                                        \
			for (unsigned i=0; i<n__; i++) x[i] op (a1 o1 x1[i]) op1 (a2 o2 x2[i]);                        \
			vreal_t result; T *x__ = x + n__;                                                              \
			const T *x1__ = x1 + n__;                                                                      \
			const T *x2__ = x2 + n__;                                                                      \
			for (unsigned i=n__; i<n; i+=N__, x__+=N__, x1__+=N__, x2__+=N__)                              \
			{                                                                                              \
				result op (a1 o1 vreal_t(x1__)) op1 (a2 o2 vreal_t(x2__));                                   \
				result.store(x__);                                                                           \
			}                                                                                              \
		}

	#define VECTOR_EXPRESSION_3(T, n, x, op, a1, o1, x1, op1, a2, o2, x2, op2, a3, o3, x3)             \
		{                                                                                                \
			typedef vecreal_<4, T> vreal_t;                                                                \
			const unsigned N__ = CYCLE_EXPANSION_SIZE;                                                     \
			unsigned n__ = n % N__;                                                                        \
			for (unsigned i=0; i<n__; i++) x[i] op (a1 o1 x1[i]) op1 (a2 o2 x2[i]) op2 (a3 o3 x3[i]);      \
			vreal_t result; T *x__ = x + n__;                                                              \
			const T *x1__ = x1 + n__;                                                                      \
			const T *x2__ = x2 + n__;                                                                      \
			const T *x3__ = x3 + n__;                                                                      \
			for (unsigned i=n__; i<n; i+=N__, x__+=N__, x1__+=N__, x2__+=N__, x3__+=N__)                   \
			{                                                                                              \
				result op (a1 o1 vreal_t(x1__)) op1 (a2 o2 vreal_t(x2__)) op2 (a3 o3 vreal_t(x3__));         \
				result.store(x__);                                                                           \
			}                                                                                              \
		}

	template <typename T>
	INLINE T scalar_product(unsigned n, const T *x, const T *y)
	{
		typedef vecreal_<4, T> vreal_t;
		const unsigned N = CYCLE_EXPANSION_SIZE;
		T s = (T)0.f;
		unsigned n__ = n % N;
		for (unsigned i=0; i<n__; i++) s += x[i] * y[i];

		vreal_t s__ = (T)0.f;
		const T *x__ = x + n__, *y__ = y + n__;
		for (unsigned i=n__; i<n; i+=N, x__+=N, y__+=N)
			s__ += vreal_t(x__) * vreal_t(y__);

		return s + summarize(s__);
	}

	template <typename T>
	INLINE T scalar_product(unsigned n, const T *x, const T *y, const T *z)
	{
		typedef vecreal_<4, T> vreal_t;
		const unsigned N = CYCLE_EXPANSION_SIZE;
		T s = (T)0.f;
		unsigned n__ = n % N;
		for (unsigned i=0; i<n__; i++) s += x[i] * y[i] * z[i];

		vreal_t s__ = (T)0.f;
		const T *x__ = x + n__, *y__ = y + n__, *z__ = z + n__;
		for (unsigned i=n__; i<n; i+=N, x__+=N, y__+=N, z__+=N)
			s__ += vreal_t(x__) * vreal_t(y__) * vreal_t(z__);

		return s + summarize(s__);
	}

	INLINE float max(unsigned n, const float *x)
	{
		typedef vecreal_<4, float> vreal_t;
		const unsigned N = CYCLE_EXPANSION_SIZE;

		float s = x[0];
		unsigned n__ = n % N;
		for (unsigned i=0; i<n__; i++) s = std::max(s, x[i]);

		vreal_t s__ = vreal_t(s);
		const float *x__ = x + n__;
		for (unsigned i=n__; i<n; i+=N, x__+=N)
			s__ = max(s, vreal_t(x__));

		return max(s__[0], s__[1], s__[2], s__[3]);
	}

	INLINE float min(unsigned n, const float *x)
	{
		typedef vecreal_<4, float> vreal_t;
		const unsigned N = CYCLE_EXPANSION_SIZE;

		float s = x[0];
		unsigned n__ = n % N;
		for (unsigned i=0; i<n__; i++) s = std::min(s, x[i]);

		vreal_t s__ = vreal_t(s);
		const float *x__ = x + n__;
		for (unsigned i=n__; i<n; i+=N, x__+=N)
			s__ = min(s, vreal_t(x__));

		return min(s__[0], s__[1], s__[2], s__[3]);
	}

	template <typename T>
	INLINE T maxabs(unsigned n, const T *x)
	{
		typedef vecreal_<4, T> vreal_t;
		const unsigned N = CYCLE_EXPANSION_SIZE;

		T s = x[0], p = x[0];
		unsigned n__ = n % N;
		for (unsigned i=0; i<n__; i++) { s = std::max(s, x[i]); p = std::min(p, x[i]); }

		vreal_t max__ = vreal_t(s), min__ = vreal_t(p);
		const T *x__ = x + n__;
		for (unsigned i=n__; i<n; i+=N, x__+=N)
		{
			vreal_t v(x__);
			max__ = max(max__, v);
			min__ = min(min__, v);
		}
		s = max(s, max__[0], max__[1], max__[2], max__[3]);
		p = min(p, min__[0], min__[1], min__[2], min__[3]);
		return std::max(fabs(s), fabs(p));
	}

	template <typename T>
	INLINE T normalize(unsigned n, T *x)
	{
		T s = sqrt(scalar_product(n, x, x));
		T s__ = 1.f / s;
		VECTOR_EXPRESSION_1(T, n, x, =, s__, *, x);
		return s;
	}

	/**
	* @brief normalize vector & returns previous norma
	*/
	template <typename T>
	INLINE T normalize(unsigned n, T *x, T new_length=1.)
	{
		T s = sqrt(scalar_product(n, x, x));
		if (s <= 0)
		{
			std::string msg = _S("zero vector length");
			PRINT_ERR(msg);
		}
		T q = new_length / s;
		VECTOR_EXPRESSION_1(T, n, x, =, q, *, x);

		return s; // return previous length of vector
	}

	/**
	* @brief fast operation s = max(V', V'')
	*/
	template <typename T>
	INLINE T max_difference(unsigned n, const T *v1, const T *v2)
	{
		T s = 0, s__;
		unsigned n__ = n % CYCLE_EXPANSION_SIZE;
		for (unsigned i=0; i<n__; i++) { s__ = fabs(v1[i] - v2[i]); if (s < s__) s = s__; }

		#define BOOST_PP_LOCAL_MACRO(__n)  T s##__n, p##__n = (T)0;
		#define BOOST_PP_LOCAL_LIMITS      (0, CYCLE_EXPANSION_SIZE - 1)
		#include BOOST_PP_LOCAL_ITERATE()
		for (unsigned i=n__; i<n; i+=CYCLE_EXPANSION_SIZE)
		{
		#define BOOST_PP_LOCAL_MACRO(__n)  \
			s##__n = fabs(v1[i+ __n] - v2[i+ __n]); if (p##__n < s##__n) p##__n = s##__n;
	  #define BOOST_PP_LOCAL_LIMITS      (0, CYCLE_EXPANSION_SIZE - 1)
		#include BOOST_PP_LOCAL_ITERATE()
		}
		#define BOOST_PP_LOCAL_MACRO(__n)  if (s < p##__n) s =  p##__n;
		#define BOOST_PP_LOCAL_LIMITS      (0, CYCLE_EXPANSION_SIZE - 1)
		#include BOOST_PP_LOCAL_ITERATE()
		return s;
	}

	/**
	* @brief fast operation s = sqrt( (1/n) * sum(i) (sqr(V'[i]- V''[i]))
	*/
	template <typename T>
	INLINE T rmsd(unsigned n, const T *v1, const T *v2)
	{
		typename extended<T>::type s = (typename extended<T>::type)0;
		unsigned n__ = n % CYCLE_EXPANSION_SIZE;
		for (unsigned i=0; i<n__; i++) s += sqr(v1[i] - v2[i]);

		#define BOOST_PP_LOCAL_MACRO(__n)  T s##__n = (T)0;
		#define BOOST_PP_LOCAL_LIMITS      (0, CYCLE_EXPANSION_SIZE - 1)
		#include BOOST_PP_LOCAL_ITERATE()
		for (unsigned i=n__; i<n; i+=CYCLE_EXPANSION_SIZE)
		{
		#define BOOST_PP_LOCAL_MACRO(__n)  s##__n += sqr(v1[i+ __n] - v2[i+ __n]);
	  #define BOOST_PP_LOCAL_LIMITS      (0, CYCLE_EXPANSION_SIZE - 1)
		#include BOOST_PP_LOCAL_ITERATE()
		}

		#define BOOST_PP_LOCAL_MACRO(__n)  s +=  s##__n;
		#define BOOST_PP_LOCAL_LIMITS      (0, CYCLE_EXPANSION_SIZE - 1)
		#include BOOST_PP_LOCAL_ITERATE()
		return (T)::sqrt(s / n);
	}

	/**
	* @brief fast operation s = (A, B) / (|A| * |B|)
	*/
	template <typename T>
	INLINE T cos_(unsigned n, const T *v1, const T *v2)
	{
		typename extended<T>::type s  = (typename extended<T>::type)0;
		typename extended<T>::type p = s, q = s;

		unsigned n__ = n % CYCLE_EXPANSION_SIZE;
		for (unsigned i=0; i<n__; i++)
		{
			s += sqr(v1[i]);
			p += sqr(v2[i]);
			q += v1[i] * v2[i];
		}
		#define BOOST_PP_LOCAL_MACRO(__n)  T s##__n = (T)0, p##__n = (T)0, q##__n = (T)0;
		#define BOOST_PP_LOCAL_LIMITS      (0, CYCLE_EXPANSION_SIZE - 1)
		#include BOOST_PP_LOCAL_ITERATE()
		for (unsigned i=n__; i<n; i+=CYCLE_EXPANSION_SIZE)
		{
		#define BOOST_PP_LOCAL_MACRO(__n)  \
			s##__n += sqr(v1[i + __n]); \
			p##__n += sqr(v2[i + __n]); \
			q##__n += v1[i + __n] * v2[i + __n];
	  #define BOOST_PP_LOCAL_LIMITS      (0, CYCLE_EXPANSION_SIZE - 1)
		#include BOOST_PP_LOCAL_ITERATE()
		}
		#define BOOST_PP_LOCAL_MACRO(__n)  s +=  s##__n; p +=  p##__n; q +=  q##__n;
		#define BOOST_PP_LOCAL_LIMITS      (0, CYCLE_EXPANSION_SIZE - 1)
		#include BOOST_PP_LOCAL_ITERATE()
		return (T)(q / ::sqrt(s * p));
	}

	/**
	* @brief fast operation s = |A|
	*/
	template <typename T> INLINE T norma(unsigned n, const T *v)
	{ return ::sqrt(scalar_product(n, v, v)); }

	/**
	* @brief ortonormalize the vector system by Gram-Schmidt method
	* @param m - number of vectors to ortogonalize
	* @param n - dimension of vector
	* @param v - array of vectors
	*/
	template <typename T>
	INLINE void ortonormalize(unsigned m, unsigned n, T **v)
	{
		for (unsigned i=0; i<m; i++)
		{
			for (unsigned j=0; j<i; j++)
			{
				T sp = scalar_product(n, &v[i][0], &v[j][0]);
				T *x = &v[i][0], *y = &v[j][0];
				VECTOR_EXPRESSION_1(T, n, x, -=, sp, *, y)
			}
			normalize(n, &v[i][0]);
		}
	}

	template <typename T>
	INLINE T cos_angle(unsigned n, const T *a, const T *b)
	{
		T s = scalar_product(n, a, a) * scalar_product(n, b, b);
		if (s <= 0)
		{
			std::string msg = _S("[ERROR] zero vector length");
			PRINT_ERR(msg);
		}
		return scalar_product(n, a, b) / sqrt(s);
	}

	template <typename T>
	INLINE T get_angle(unsigned n, const T *a, const T *b)
	{
		T cosphi = cos_angle(a, b);
		return safe_acos(cosphi);
	}


}
#endif
