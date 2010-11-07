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
#ifndef _M3X3DENSE___0077A726_ABC2_5180_A573_CE435B530E00__H
#define _M3X3DENSE___0077A726_ABC2_5180_A573_CE435B530E00__H
#include "prgkern/_prgconfig.h"

#include "prgkern/_math.h"
#include "prgkern/_dense.h"
#include "prgkern/_v3dense.h"

namespace prgkern
{

	template <unsigned N, typename T>
	class mdense_<N, N, T>
	{
	#define _A(i) assert(_LT((unsigned)i,(unsigned)N))

	public:
		typedef T value_type;
		typedef vdense_<N, T> row_type;
		enum { dim = N, dim1 = N, dim2 = N };

		mdense_() { scalar_assign_equal((_I2T<N> *)0, (T)0, (T)0); }
		template <typename S> mdense_(S s) { scalar_assign_equal((_I2T<N> *)0, (T)s, (T)0); }
		template <typename S> mdense_(S s00, S s01, S s10, S s11)
		{
			assert(_EQ(N, 2));
			v_[0][0] = (T)s00; v_[0][1] = (T)s01;
			v_[1][0] = (T)s10; v_[1][1] = (T)s11;
		}
		template <typename S> mdense_(S s00, S s01, S s02, S s10, S s11, S s12, S s20, S s21, S s22)
		{
			assert(_EQ(N, 3));
			v_[0][0] = (T)s00; v_[0][1] = (T)s01; v_[0][2] = (T)s02;
			v_[1][0] = (T)s10; v_[1][1] = (T)s11; v_[1][2] = (T)s12;
			v_[2][0] = (T)s20; v_[2][1] = (T)s21; v_[2][2] = (T)s22;
		}
		template <typename S> mdense_(S s00, S s01, S s02, S s03, S s10, S s11, S s12, S s13,
			S s20, S s21, S s22, S s23, S s30, S s31, S s32, S s33)
		{
			assert(_EQ(N, 4));
			v_[0][0] = (T)s00; v_[0][1] = (T)s01; v_[0][2] = (T)s02; v_[0][3] = (T)s03;
			v_[1][0] = (T)s10; v_[1][1] = (T)s11; v_[1][2] = (T)s12; v_[1][3] = (T)s13;
			v_[2][0] = (T)s20; v_[2][1] = (T)s21; v_[2][2] = (T)s22; v_[2][3] = (T)s23;
			v_[3][0] = (T)s30; v_[3][1] = (T)s31; v_[3][2] = (T)s32; v_[3][3] = (T)s33;
		}

		mdense_(const mdense_& v) { matrix_assign_equal((_I2T<N> *)0, v); }

		void clear() { scalar_assign_equal((_I2T<N> *)0, (T)0, (T)0); }

		row_type &operator[](unsigned n) { _A(n); return v_[n]; }
		const row_type &operator[](unsigned n) const { _A(n); return v_[n]; }
		T operator()(unsigned i, unsigned j) const { return (*this)[i][j]; }
		T &operator()(unsigned i, unsigned j) { return (*this)[i][j]; }

	// implementation with non-diagonal elemenrs = 0
	#define IMPLEMENT_SCALAR_OPERATOR(op, name) \
		template <typename S> mdense_ &operator op(S s) \
		{ scalar_assign_##name((_I2T<N> *)0, (T)s, (T)0); return *this; }
		IMPLEMENT_SCALAR_OPERATOR(=, equal)
		IMPLEMENT_SCALAR_OPERATOR(+=, plus_equal)
		IMPLEMENT_SCALAR_OPERATOR(-=, minus_equal)
	#undef IMPLEMENT_SCALAR_OPERATOR

	// implementation with non-diagonal elemenrs != 0
	#define IMPLEMENT_SCALAR_OPERATOR(op, name) \
		template <typename S> mdense_ &operator op(S s) \
		{ scalar_assign_##name((_I2T<N> *)0, (T)s, (T)s); return *this; }
		IMPLEMENT_SCALAR_OPERATOR(*=, multiplies_equal)
		// IMPLEMENT_SCALAR_OPERATOR(/=) unefficient
		//   you must use *=(1/s) instead
	#undef IMPLEMENT_SCALAR_OPERATOR

	#define IMPLEMENT_MATRIX_OPERATOR(op, name) \
		mdense_ &operator op(const mdense_ &v) \
		{ matrix_assign_##name((_I2T<N> *)0, v); return *this; }
		IMPLEMENT_MATRIX_OPERATOR(=, equal)
			// neglect self assignment A=A control
		IMPLEMENT_MATRIX_OPERATOR(+=, plus_equal)
		IMPLEMENT_MATRIX_OPERATOR(-=, minus_equal)
	#undef IMPLEMENT_MATRIX_OPERATOR

	protected:

	#define ASSIGN(op, d, s, name) \
		void scalar_assign_##name(_I2T<2> *, T d, T s) { \
			v_[0][0] op d; v_[0][1] op s; v_[1][0] op s; v_[1][1] op d; \
		}
		ASSIGN(=, d, s, equal)
		ASSIGN(+=, d, s, plus_equal)
		ASSIGN(-=, d, s, minus_equal)
		ASSIGN(*=, d, s, multiplies_equal)
	#undef ASSIGN

	#define ASSIGN(op, d, s, name) \
		void scalar_assign_##name(_I2T<3> *, T d, T s) { \
			v_[0][0] op d; v_[0][1] op s; v_[0][2] op s; \
			v_[1][0] op s; v_[1][1] op d; v_[1][2] op s; \
			v_[2][0] op s; v_[2][1] op s; v_[2][2] op d; \
		}
		ASSIGN(=, d, s, equal)
		ASSIGN(+=, d, s, plus_equal)
		ASSIGN(-=, d, s, minus_equal)
		ASSIGN(*=, d, s, multiplies_equal)
	#undef ASSIGN

	#define ASSIGN(op, d, s, name) \
		void scalar_assign_##name(_I2T<4> *, T d, T s) { \
			v_[0][0] op d; v_[0][1] op s; v_[0][2] op s; v_[0][3] op s; \
			v_[1][0] op s; v_[1][1] op d; v_[1][2] op s; v_[1][3] op s; \
			v_[2][0] op s; v_[2][1] op s; v_[2][2] op d; v_[2][3] op s; \
			v_[3][0] op s; v_[3][1] op s; v_[3][2] op s; v_[3][3] op d; \
		}
		ASSIGN(=, d, s, equal)
		ASSIGN(+=, d, s, plus_equal)
		ASSIGN(-=, d, s, minus_equal)
		ASSIGN(*=, d, s, multiplies_equal)
	#undef ASSIGN

	#define ASSIGN(op, m, name) \
		void matrix_assign_##name(_I2T<2> *, const mdense_ &m) { \
			v_[0] op m.v_[0]; v_[1] op m.v_[1]; \
		}
		ASSIGN(=, m, equal)
		ASSIGN(+=, m, plus_equal)
		ASSIGN(-=, m, minus_equal)
	#undef ASSIGN

	#define ASSIGN(op, m, name) \
		void matrix_assign_##name(_I2T<3> *, const mdense_ &m) { \
			v_[0] op m.v_[0]; v_[1] op m.v_[1]; v_[2] op m.v_[2]; \
		}
		ASSIGN(=, m, equal)
		ASSIGN(+=, m, plus_equal)
		ASSIGN(-=, m, minus_equal)
	#undef ASSIGN

	#define ASSIGN(op, m, name) \
		void matrix_assign_##name(_I2T<4> *, const mdense_ &m) { \
			v_[0] op m.v_[0]; v_[1] op m.v_[1]; v_[2] op m.v_[2]; v_[3] op m.v_[3]; \
		}
		ASSIGN(=, m, equal)
		ASSIGN(+=, m, plus_equal)
		ASSIGN(-=, m, minus_equal)
	#undef ASSIGN

	protected:
		row_type v_[N];
	#undef _A
	};

	template <unsigned N, typename T>
	INLINE mdense_<N, N, T> operator*(const mdense_<N, N, T> &a1, T a2)
	{ mdense_<N, N, T> tmp(a1); tmp *= a2; return tmp; }

	template <unsigned N, typename T>
	INLINE mdense_<N, N, T> operator*(T a1, const mdense_<N, N, T> &a2)
	{ mdense_<N, N, T> tmp(a2); tmp *= a1; return tmp; }

	template <unsigned N, typename T>
	INLINE mdense_<N, N, T> operator+(const mdense_<N, N, T> &a1, const mdense_<N, N, T> &a2)
	{ mdense_<N, N, T> tmp(a1); tmp += a2; return tmp; }

	template <unsigned N, typename T>
	INLINE mdense_<N, N, T> operator-(const mdense_<N, N, T> &a1, const mdense_<N, N, T> &a2)
	{ mdense_<N, N, T> tmp(a1); tmp -= a2; return tmp; }


	template <typename T>
	INLINE void multiplies(mdense_<2, 2, T> &R, const mdense_<2, 2, T> &A,
		const mdense_<2, 2, T> &B)
	{
	#define M(i,j) A[i][0]*B[0][j] + A[i][1]*B[1][j]
		R[0][0] = M(0,0); R[0][1] = M(0,1);
		R[1][0] = M(1,0); R[1][1] = M(1,1);
	#undef M
	}

	template <typename T>
	INLINE void multiplies(mdense_<3, 3, T> &R, const mdense_<3, 3, T> &A,
		const mdense_<3, 3, T> &B)
	{
	#define M(i,j) A[i][0]*B[0][j] + A[i][1]*B[1][j] + A[i][2]*B[2][j]
		R[0][0] = M(0,0); R[0][1] = M(0,1); R[0][2] = M(0,2);
		R[1][0] = M(1,0); R[1][1] = M(1,1); R[1][2] = M(1,2);
		R[2][0] = M(2,0); R[2][1] = M(2,1); R[2][2] = M(2,2);
	#undef M
	}

	template <typename T>
	INLINE void multiplies(mdense_<4, 4, T> &R, const mdense_<4, 4, T> &A, const mdense_<4, 4, T> &B)
	{
	#define M(i,j) A[i][0]*B[0][j] + A[i][1]*B[1][j] + A[i][2]*B[2][j] + A[i][3]*B[3][j]
		R[0][0] = M(0,0); R[0][1] = M(0,1); R[0][2] = M(0,2); R[0][3] = M(0,3);
		R[1][0] = M(1,0); R[1][1] = M(1,1); R[1][2] = M(1,2); R[1][3] = M(1,3);
		R[2][0] = M(2,0); R[2][1] = M(2,1); R[2][2] = M(2,2); R[2][3] = M(2,3);
		R[3][0] = M(3,0); R[3][1] = M(3,1); R[3][2] = M(3,2); R[3][3] = M(3,3);
	#undef M
	}

	template <unsigned N, typename T>
	INLINE mdense_<N, N, T> operator*(const mdense_<N, N, T> &A, const mdense_<N, N, T> &B)
	{ mdense_<N, N, T> R; multiplies(R, A, B); return R; }

	template <unsigned N, typename T>
	INLINE std::string make_string(const mdense_<N, N, T> &v)
	{
		std::ostringstream oss_convert;
		oss_convert << "[" ;
		for (unsigned i=0; i<N-1; i++) oss_convert << make_string(v[i]) << ", ";
		oss_convert << make_string(v[N-1]) << "]";
		return oss_convert.str();
	}

	template <unsigned N, typename T>
	INLINE std::ostream &operator<<(std::ostream &os, const mdense_<N, N, T> &v)
	{ os << make_string(v); return os; }

	template <typename T>
	INLINE void multiplies(vdense_<2, T> &R, const mdense_<2, 2, T> &M, const vdense_<2, T> &V)
	{
		R[0] = M[0][0]*V[0] + M[0][1]*V[1];
		R[1] = M[1][0]*V[0] + M[1][1]*V[1];
	}

	template <typename T>
	INLINE void multiplies(vdense_<3, T> &R, const mdense_<3, 3, T> &M, const vdense_<3, T> &V)
	{
		R[0] = M[0][0]*V[0] + M[0][1]*V[1] + M[0][2]*V[2];
		R[1] = M[1][0]*V[0] + M[1][1]*V[1] + M[1][2]*V[2];
		R[2] = M[2][0]*V[0] + M[2][1]*V[1] + M[2][2]*V[2];
	}

	template <typename T>
	INLINE void multiplies(vdense_<4, T> &R, const mdense_<4, 4, T> &M, const vdense_<4, T> &V)
	{
		R[0] = M[0][0]*V[0] + M[0][1]*V[1] + M[0][2]*V[2] + M[0][3]*V[3];
		R[1] = M[1][0]*V[0] + M[1][1]*V[1] + M[1][2]*V[2] + M[1][3]*V[3];
		R[2] = M[2][0]*V[0] + M[2][1]*V[1] + M[2][2]*V[2] + M[2][3]*V[3];
		R[3] = M[3][0]*V[0] + M[3][1]*V[1] + M[3][2]*V[2] + M[3][3]*V[3];
	}

	template <typename T>
	INLINE void multiplies(vdense_<2, T> &R, const vdense_<2, T> &V, const mdense_<2, 2, T> &M)
	{
		R[0] = M[0][0]*V[0] + M[1][0]*V[1];
		R[1] = M[0][1]*V[0] + M[1][1]*V[1];
	}

	template <typename T>
	INLINE void multiplies(vdense_<3, T> &R, const vdense_<3, T> &V, const mdense_<3, 3, T> &M)
	{
		R[0] = M[0][0]*V[0] + M[1][0]*V[1] + M[2][0]*V[2];
		R[1] = M[0][1]*V[0] + M[1][1]*V[1] + M[2][1]*V[2];
		R[2] = M[0][2]*V[0] + M[1][2]*V[1] + M[2][2]*V[2];
	}

	template <typename T>
	INLINE void multiplies(vdense_<4, T> &R, const vdense_<4, T> &V, const mdense_<4, 4, T> &M)
	{
		R[0] = M[0][0]*V[0] + M[1][0]*V[1] + M[2][0]*V[2] + M[3][0]*V[3];
		R[1] = M[0][1]*V[0] + M[1][1]*V[1] + M[2][1]*V[2] + M[3][1]*V[3];
		R[2] = M[0][2]*V[0] + M[1][2]*V[1] + M[2][2]*V[2] + M[3][2]*V[3];
		R[3] = M[0][3]*V[0] + M[1][3]*V[1] + M[2][3]*V[2] + M[3][3]*V[3];
	}

	template <unsigned N, typename T>
	INLINE vdense_<N, T> operator*(const mdense_<N, N, T> &M, const vdense_<N, T> &V)
	{ vdense_<N, T> R; multiplies(R, M, V); return R; }

	template <unsigned N, typename T>
	INLINE vdense_<N, T> operator*(const vdense_<N, T> &V, const mdense_<N, N, T> &M)
	{ vdense_<N, T> R; multiplies(R, V, M); return R; }

	template <typename T>
	INLINE void direct_product(mdense_<2, 2, T> &m, const vdense_<2, T> &X1,
		const vdense_<2, T> &X2)
	{
		m[0][0] = X1[0] * X2[0]; m[0][1] = X1[0] * X2[1];
		m[1][0] = X1[1] * X2[0]; m[1][1] = X1[1] * X2[1];
	}

	template <typename T>
	INLINE void direct_product(mdense_<3, 3, T> &m, const vdense_<3, T> &X1,
		const vdense_<3, T> &X2)
	{
		m[0][0] = X1[0] * X2[0]; m[0][1] = X1[0] * X2[1]; m[0][2] = X1[0] * X2[2];
		m[1][0] = X1[1] * X2[0]; m[1][1] = X1[1] * X2[1]; m[1][2] = X1[1] * X2[2];
		m[2][0] = X1[2] * X2[0]; m[2][1] = X1[2] * X2[1]; m[2][2] = X1[2] * X2[2];
	}

	template <typename T>
	INLINE void direct_product(mdense_<4, 4, T> &m, const vdense_<4, T> &X1,
		const vdense_<4, T> &X2)
	{
	#define _M(i,j) m[i][j] = X1[i] * X2[j];
	_M(0,0) _M(0,1) _M(0,2) _M(0,3)
	_M(1,0) _M(1,1) _M(1,2) _M(1,3)
	_M(2,0) _M(2,1) _M(2,2) _M(2,3)
	_M(3,0) _M(3,1) _M(3,2) _M(3,3)
	#undef _M
	}

	template <unsigned N, typename T>
	INLINE mdense_<N, N, T> direct_product(const vdense_<N, T> &X1, const vdense_<N, T> &X2)
	{ mdense_<N, N, T> m; direct_product(m, X1, X2); return m; }

	template <typename T>
	INLINE T determinant(const mdense_<2, 2, T> &m)
	{ return m[0][0] * m[1][1]  - m[0][1] * m[1][0]; }

	template <typename T>
	INLINE T determinant(const mdense_<3, 3, T> &m)
	{
		return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
			+ m[0][1] * (m[1][2] * m[2][0] - m[1][0] * m[2][2])
			+ m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
	}

	/**
	* @brief get U(-1) from 3x3 matrix
	*/
	template <typename T>
	INLINE mdense_<3, 3, T> invert(const mdense_<3, 3, T> &u)
	{
		T u00 = u[1][1] * u[2][2] - u[1][2] * u[2][1];
		T u01 = u[1][2] * u[2][0] - u[1][0] * u[2][2];
		T u02 = u[1][0] * u[2][1] - u[1][1] * u[2][0];
		T det__ = 1. / (u[0][0] * u00 + u[0][1] * u01 + u[0][2] * u02);
		mdense_<3, 3, T> m;
		m[0][0] = det__ * u00;
		m[1][0] = det__ * u01;
		m[2][0] = det__ * u02;
		m[0][1] = det__ * (u[0][2]*u[2][1] - u[0][1]*u[2][2]);
		m[0][2] = det__ * (u[0][1]*u[1][2] - u[0][2]*u[1][1]);
		m[1][1] = det__ * (u[0][0]*u[2][2] - u[0][2]*u[2][0]);
		m[1][2] = det__ * (u[0][2]*u[1][0] - u[0][0]*u[1][2]);
		m[2][1] = det__ * (u[0][1]*u[2][0] - u[0][0]*u[2][1]);
		m[2][2] = det__ * (u[0][0]*u[1][1] - u[0][1]*u[1][0]);
		return m;
	}

}
#endif
