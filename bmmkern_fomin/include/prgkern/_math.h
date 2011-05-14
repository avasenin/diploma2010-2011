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
#ifndef _MATH__0077A726_728A_56f7_9465_CE437EC80000__H
#define _MATH__0077A726_728A_56f7_9465_CE437EC80000__H

#include "prgkern/_prgconfig.h"
#include "prgkern/_assert.h"
#include "prgkern/_type.h"

#include "prgkern/_for.h"
#define bll boost::lambda

#include <cstddef>
#include <limits>
#include <memory>
#include <cmath>

#ifdef USE_LINUX_SPECIFIC_CODE
	#include <sys/time.h>
	#include <sys/utsname.h>
#endif

#undef INFINITY
#undef max
#undef min

namespace prgkern
{

	/// constants
	const double M_SQRT_1  = 1.;
	const double M_SQRT_2  = ::sqrt(2.);
	const double M_SQRT_3  = ::sqrt(3.);
	const double M_SQRT_4  = 2.;
	const double M_SQRT_5  = ::sqrt(5.);
	const double M_SQRT_6  = ::sqrt(6.);
	const double M_SQRT_7  = ::sqrt(7.);
	const double M_SQRT_8  = ::sqrt(8.);
	const double M_SQRT_9  = 3.;
	const double M_SQRT_10 = ::sqrt(10.);
	const double M_SQRT_11 = ::sqrt(11.);
	const double M_SQRT_12 = ::sqrt(12.);

	const double M_1_SQRT_2 = 1. / M_SQRT_2;
	const double M_1_SQRT_3 = 1. / M_SQRT_3;

	const double M_2PI = M_PI + M_PI;
	const double M_4PI = M_2PI + M_2PI;
	const double M_SQRT_PI = ::sqrt(M_PI);
	const double M_SQRT_2PI = ::sqrt(M_2PI);
	const double M_REVERSE_SQRT_PI = 1. / M_SQRT_PI;
	const double M_REVERSE_SQRT_2PI = 1. / M_SQRT_2PI;
	const double M_REVERSE_LN2 = 1. / ::log(2.);

	const double M_INV_2 = 1./ 2.;
	const double M_INV_3 = 1./ 3.;
	const double M_INV_4 = 1./ 4.;
	const double M_INV_5 = 1./ 5.;
	const double M_INV_6 = 1./ 6.;
	const double M_INV_7 = 1./ 7.;
	const double M_INV_8 = 1./ 8.;
	const double M_INV_9 = 1./ 9.;

	/// const to translate degree -> radians
	const double DEG2RAD = M_PI / 180.;
	const double RAD2DEG = 180. / M_PI;

	const double M_COS_0  = 1.;
	const double M_COS_30 = 0.5 * M_SQRT_3;
	const double M_COS_45 = 0.5 * M_SQRT_2;
	const double M_COS_60 = 0.5;
	const double M_COS_90 = 0.;

	const double M_SIN_0  = 0.;
	const double M_SIN_30 = 0.5;
	const double M_SIN_45 = 0.5 * M_SQRT_2;
	const double M_SIN_60 = 0.5 * M_SQRT_3;
	const double M_SIN_90 = 1.;

	/// infinities
	template <typename T> struct infinity
	{ operator T() { return std::numeric_limits<T>::max(); } };

	/// maxima
	template <typename T> INLINE T max(T t1, T t2, T t3)
	{ return std::max(t1, std::max(t2, t3)); }

	template <typename T> INLINE T max(T t1, T t2, T t3, T t4)
	{ return std::max(t1, max(t2, t3, t4)); }

	template <typename T> INLINE T max(T t1, T t2, T t3, T t4, T t5)
	{ return std::max(t1, max(t2, t3, t4, t5)); }

	/// minima
	template <typename T> INLINE T min(T t1, T t2, T t3)
	{ return std::min(t1, std::min(t2, t3)); }

	template <typename T> INLINE T min(T t1, T t2, T t3, T t4)
	{ return std::min(t1, min(t2, t3, t4)); }

	template <typename T> INLINE T min(T t1, T t2, T t3, T t4, T t5)
	{ return std::min(t1, min(t2, t3, t4, t5)); }

	/// minimaxima
	template <typename T> INLINE T minmax(T t1, T t2, T t3)
	{ return min(std::max(t1, t2), std::max(t2, t3), std::max(t3, t1)); }

	template <typename T> INLINE T minmax(T t1, T t2, T t3, T t4)
	{ return min(std::max(t1, t2), std::max(t2, t3), std::max(t3, t4), std::max(t4, t1)); }

	/// maximinima
	template <typename T> INLINE T maxmin(T t1, T t2, T t3)
	{ return max(std::min(t1, t2), std::min(t2, t3), std::min(t3, t1)); }

	template <typename T> INLINE T maxmin(T t1, T t2, T t3, T t4)
	{ return max(std::min(t1, t2), std::min(t2, t3), std::min(t3, t4), std::min(t4, t1)); }

	/// polinomes
	template <typename T> INLINE
	T polynome(T c0, T x) { return c0; }

	template <typename T> INLINE
	T polynome(T c1, T c0, T x) { return c1 * x + polynome(c0, x); }

	template <typename T> INLINE
	T polynome(T c2, T c1, T c0, T x) { return polynome(c2, c1, x) * x + polynome(c0, x); }

	template <typename T> INLINE
	T polynome(T c3, T c2, T c1, T c0, T x)
	{ return polynome(c3, c2, c1, x) * x + polynome(c0, x); }

	template <typename T> INLINE
	T polynome(T c4, T c3, T c2, T c1, T c0, T x)
	{ return polynome(c4, c3, c2, c1, x) * x + polynome(c0, x); }

	template <typename T> INLINE
	T polynome(T c5, T c4, T c3, T c2, T c1, T c0, T x)
	{ return polynome(c5, c4, c3, c2, c1, x) * x + polynome(c0, x); }

	template <typename T> INLINE
	T polynome(T c6, T c5, T c4, T c3, T c2, T c1, T c0, T x)
	{ return polynome(c6, c5, c4, c3, c2, c1, x) * x + polynome(c0, x); }

	template <typename T> INLINE
	T polynome(T c7, T c6, T c5, T c4, T c3, T c2, T c1, T c0, T x)
	{ return polynome(c7, c6, c5, c4, c3, c2, c1, x) * x + polynome(c0, x); }

	template <typename T> INLINE
	T polynome(T c8, T c7, T c6, T c5, T c4, T c3, T c2, T c1, T c0, T x)
	{ return polynome(c8, c7, c6, c5, c4, c3, c2, c1, x) * x + polynome(c0, x); }

	template <typename T> INLINE
	T polynome(T c9, T c8, T c7, T c6, T c5, T c4, T c3, T c2, T c1, T c0, T x)
	{ return polynome(c9, c8, c7, c6, c5, c4, c3, c2, c1, x) * x + polynome(c0, x); }

	template <typename T> INLINE
	T polynome(T c10, T c9, T c8, T c7, T c6, T c5, T c4, T c3, T c2, T c1, T c0, T x)
	{ return polynome(c10, c9, c8, c7, c6, c5, c4, c3, c2, c1, x) * x + polynome(c0, x); }

	/// integrate of polinomes
	template <typename T> INLINE
	T integral(T c0, T a, T b) { return c0 * (b - a); }

	template <typename T> INLINE
	T integral(T c1, T c0, T a, T b)
	{
		return polynome(c1 * M_INV_2, c0, b) * b - polynome(c1 * M_INV_2, c0, a) * a;
	}

	template <typename T> INLINE
	T integral(T c2, T c1, T c0, T a, T b)
	{
		return polynome(c2 * M_INV_3, c1 * M_INV_2, c0, b) * b
			- polynome(c2 * M_INV_3, c1 * M_INV_2, c0, a) * a;
	}

	template <typename T> INLINE
	T integral(T c3, T c2, T c1, T c0, T a, T b)
	{
		return polynome(c3 * M_INV_4, c2 * M_INV_3, c1 * M_INV_2, c0, b) * b
			- polynome(c3 * M_INV_4, c2 * M_INV_3, c1 * M_INV_2, c0, a) * a;
	}

	/// sign functions
	template <typename T> INLINE T abs(T a) { return a >= 0 ? a : -a; }
	template <typename T, typename T2> INLINE T sign(T a, T2 b) { return b >= 0 ? abs(a) : -abs(a); }
	template <typename T> INLINE int sign(T a) { return a >= 0 ? (a > 0 ? 1 : 0) : -1; }

	template <typename T> INLINE void sort(T &n0, T &n1)
	{ if (n0 > n1) std::swap(n0, n1); }

	template <typename T> INLINE void sort(T &n0, T &n1, T &n2)
	{
		if (n0 > n2) std::swap(n0, n2);
		if (n0 > n1) std::swap(n0, n1);
		if (n1 > n2) std::swap(n1, n2);
	}

	template <typename T> INLINE void sort(T &n0, T &n1, T &n2, T &n3)
	{
		if (n0 > n1) std::swap(n0, n1);
		if (n2 > n3) std::swap(n2, n3);
		if (n0 > n2) std::swap(n0, n2);
		if (n1 > n3) std::swap(n1, n3);
		if (n1 > n2) std::swap(n1, n2);
	}

	template <typename T> INLINE T multiple_not(bool b, T a) { return !b ? a : (T)0; }
	template <typename T> INLINE T multiple    (bool b, T a) { return  b ? a : (T)0; }

	/**
	* @brief modulo of integer value
	* @note safety for using with negative values vs c++ operator%
	* @param m - base of modulo
	* @param x - value for calculate
	* @return modulo(x)
	*/
	template <typename T> INLINE T mod(T x, T m)
	{ return x<0 ? m - 1 - ((-x) - 1)%m : x%m; }

	INLINE long double mod(long double x, long double m)
	{ return x<0 ? x + m - int(x/m) * m : x - int(x/m) * m; }

	INLINE double mod(double x, double m)
	{ return x<0 ? x + m - int(x/m) * m : x - int(x/m) * m; }

	INLINE float mod(float x, float m)
	{ return x<0 ? x + m - int(x/m) * m : x - int(x/m) * m; }

	INLINE int celling(int x, int radix=CYCLE_EXPANSION_SIZE)
	{ return (x - 1) + radix - (x - 1)%radix; }

	/// degrees functions
	template <typename T> INLINE T sqr(T s) { return s * s; }
	template <typename T> INLINE T cube(T s) { return s * s * s; }
	template <typename T> INLINE T rsqrt(T a) { return T(1.f) / sqrt(a); }
	template <typename T> INLINE T rcp(T a) { return T(1.f) / a; }

	template <unsigned N, typename T> INLINE T degree(T s)
	{
		T result = 1;
		for_<0, N>::expand(bll::var(result) *= s);
		return result;
	}

	template <int S, unsigned DEGREE> struct Degree
	{ enum { value = S * Degree<S, DEGREE-1>::value }; };
	template <int S> struct Degree<S, 0> { enum { value = 1 }; };

	/// round function
	template <unsigned digits, unsigned radix>
	INLINE double round(double x)
	{
		const int k = degree<digits + 1>(radix);
			// support given digits after points ^^^

		int exponent;
		double mantissa = ::frexp(x, &exponent);
		unsigned round_mantissa = size_t((mantissa + mantissa) * k);
			// [0.5, 1) -> [1, 2)                   ^^^
			//   'cos integer division must be applied to unit interval

		double factor = 0.5 / k;
		mantissa = round_mantissa * factor;
		return ldexp(mantissa, exponent);
	}

	/// comparison functions
	template <typename T> INLINE bool equal(T a, T b,
		T accuracy=std::numeric_limits<T>::epsilon())
	{
		T s = abs(a) + abs(b);
		if (s < accuracy) return true; // absolute comparison with 0
		if (abs(a - b) < accuracy * s) return true; // relative comparison
		return false;
	}

	/// extansion	acos upto [-inf, +inf] to avoid round errors
	template <typename T> INLINE T safe_acos(T cosphi)
	{
		if (cosphi <= -1.) return T(M_PI);
		if (cosphi >=  1.) return T(0.);
		return (T) acos(cosphi);
	}


	template <typename T> INLINE T mult2(T s) { return s + s; }
	INLINE unsigned mult2(unsigned s) { return s << 1; }

	template <int N> struct _FLOG2 { enum { result = 1 + _FLOG2<(N >> 1)>::result }; };
	template <> struct _FLOG2<1> { enum { result = 0 }; };

	/** @brief zero of function
	* @note build linear approximation f() ~ B * x + C
	* @param xa - left bound of interval
	* @param xb - right bound of interval
	* @param fa - f(xa)
	* @param fb - f(xb)
	* @return position of zero
	*/
	template <typename T>
	INLINE T zero_linear(T xa, T xb, T fa, T fb)
	{
		return (xa * fb - xb * fa) /  (fb - fa);
	}

	/** @brief zero of function
	* @note build linear approximation f() = A * x + B
	* @param xa - left bound of interval
	* @param xb - right bound of interval
	* @param fa - f(xa)
	* @param fb - f(xb)
	* @return position of zero
	*/
	template <typename T>
	INLINE T zero_linear(T &A, T &B, T xa, T xb, T fa, T fb)
	{
		T _1dx = (T)1. / (xb - xa);
		A = (fb - fa) * _1dx;
		B = (fa *xb - fb * xa) * _1dx;
		return -B / A;
	}

	/**
	 * Интерполяция фунции по трем точкам на интервале [xa, xc] с целью получения пересечения
	 * функции с 0. Функция специально выведена в двойную точность (не используется шаблон),
	 * так как для float часто происходит nan при счете, так как аргуметны имеют высокие степени.
	 */
	INLINE double zero_quadratic(double xa, double xb, double xc, double fa, double fb, double fc,
		double accuracy=10*std::numeric_limits<double>::epsilon())
	{
		assert(_LT(xa, xb)); assert(_LT(xb, xc));
		assert(_LT(fa, 0.)); assert(_GT(fc, 0.));

		double fa__ = fa - fb; double fc__ = fc - fb;
		double xa__ = xa - xb;	double xc__ = xc - xb;
		double fcxa = fc__ * xa__;
		double faxc = fa__ * xc__;

		if (fabs(fcxa - faxc) < accuracy) return zero_linear(xa, xc, fa, fc);
			// при вырождении используем линейную аппроксимацию

		double coef = (double)1. / (xa__ * xc__ * (xc - xa));
		double A = (fcxa - faxc) * coef;
		double B = (faxc * xc__ - fcxa * xa__) * coef;
		double D = sqr(B) - 4 * A * fb;
		assert(_GT(D, (double)0.));

		D = sqrt(D);
		A = (double)(-1.) / (A + A);
		double x1 = xb + (B + D) * A;
		double x2 = xb + (B - D) * A;

		double minx = std::min(x1, x2);
		double maxx = std::max(x1, x2);

		// отбрасываем значения вне интервала
		if (minx < xa) return maxx;
		if (maxx > xc) return minx;

		// возвращаем всегда решение, которое ближе к левому краю
		return minx;
	}

	/** @brief zero of function
	* @note build trigonometric approximation f() ~ B * cos(x) + C * sin(x)
	*  for perioic functions and finding the zero of it
	* @param xa - left bound of interval
	* @param xb - right bound of interval
	* @param fa - f(xa)
	* @param fb - f(xb)
	* @param n  - порядок цикличности (число максимумов на интервале [0..2pi])
	* @return position of zero
	*/
	template <typename T>
	INLINE T zero_trigonometric(T xa, T xb, T fa, T fb, unsigned n=2)
	{
		xa *= n; xb *= n;
		return (T)atan( (fa * sin(xb) - fb * sin(xa) / (fa * cos(xb) - fb * cos(xa)))) / n;
	}

	/** @brief extreme of function
	* @note build quadratic approximation f() ~ A * x ** 2 + B * x + C
	* @param x[out] position of extreme
	* @param f[out] estimate of value in extreme
	* @param xa - left bound of interval
	* @param xb - right bound of interval
	* @param fa - f(xa)
	* @param fb - f(xb)
	* @param ga - f'(xa)
	* @param gb - f'(xb)
	* @return number of roots found (0 - degeneracy)
	*/
	template <typename T>
	INLINE int extreme_sqr(T *x, T *f, T xa, T xb, T fa, T fb, T ga, T gb,
		T accuracy=std::numeric_limits<T>::epsilon())
	{
		T h = xb - xa;
		T C = h * ga;
		T B = fb - fa + C; // B -> B / h
		if (equal(B, (T)0., accuracy))
		{
		#ifdef DEGENERACY_PRINT
			{
				std::string msg = _S("[WARNING] extreme degeneracy");
				PRINT_MESSAGE(msg);
			}
		#endif
			return 0;
		}
		*x = -ga / mult2(B);
		*f = polynome(B, C, fa, *x);
		*x = h * (*x) + xa; // restore x -> x * h + xa
		return 1;
	}

	/** @brief extreme(s) of function
	* @note build cubic approximation f() ~ A * x**3 + B * x**2 + C * x + D
	* @param x[out] position of extremes
	* @param f[out] estimate of value in extremes
	* @param xa - left bound of interval
	* @param xb - right bound of interval
	* @param fa - f(xa)
	* @param fb - f(xb)
	* @param ga - f'(xa)
	* @param gb - f'(xb)
	* @return number of roots found (0 - degeneracy)
	*  return roots are ordered by increasing of f()
	*  so in case of one root it may be minimum (or maximum)
	*/
	template <typename T>
	INLINE int extreme_cube(T *xmin, T *xmax, T *fmin, T *fmax,
		T xa, T xb, T fa, T fb, T ga, T gb,
		T accuracy=std::numeric_limits<T>::epsilon())
	{
		T h = xb - xa;
		T A = mult2(fa - fb) + h * (ga + gb); // #A = A * h**3

		if (equal(A, 0., accuracy)) // // probe the quadratic approximation
			return extreme_sqr(xmin, fmin, xa, xb, fa, fb, ga, gb, accuracy);

		T _3A = 3 * (-A); // A -> -3*A
		T B = 3 * (fb - fa) - h * (gb + mult2(ga)); // #B = B * h**2
		T C = h * ga; // #C = C * h
		T det = (sqr(B) + _3A * C);

		if (det <= 0) return 0;
		det = ::sqrt(det);

		*xmin = (B - det) / _3A;
		*xmax = (B + det) / _3A;
		*fmin = polynome(A, B, C, fa, *xmin);
		*fmax = polynome(A, B, C, fa, *xmax);

		*xmin = h * (*xmin) + xa;
		*xmax = h * (*xmax) + xa;
		if (*fmin > *fmax)
		{
			std::swap(*xmin, *xmax);
			std::swap(*fmin, *fmax);
		}
		return 2;
	}

	/**
	*   Линейная интерполяция функции заданной на 2D-сетке.
	* @note Предполагается, что координаты узлов равны 0 и 1.
	*   Нумерация индексов такова, что быстрее всего меняется последний индекс.
	*   Для координат (x,y,z) также быстрее меняется последний индекс z.
	* @param f00,f01,f10,f11 значения функции в узлах сетки
	* @param x0,y0 относительное расстояние от точки до левого нижнего узла
	* @return аппроксимированное значение функции
	*/
	template <typename _Real>
	INLINE _Real linear_interpolation(_Real x0, _Real y0,
		_Real f00, _Real f01, _Real f10, _Real f11)
	{
		_Real x1 = 1. - x0;
		_Real y1 = 1. - y0;
		return ( x1 * f00 + x0 * f10 ) * y1
		     + ( x1 * f01 + x0 * f11 ) * y0;
	}

	/**
	*   линейная интерполяция функции заданной на 3D-сетке
	* @note Предполагается, что координаты узлов равны 0 и 1.
	*   Нумерация индексов такова, что быстрее всего меняется последний индекс.
	*   Для координат (x,y,z) также быстрее меняется последний индекс z.
	* @param f000,f001,f010,f011... значения функции в узлах сетки
	* @param x0,y0,z0 относительное расстояние от точки до левого нижнего узла
	* @return аппроксимированное значение функции
	*/
	template <typename _Real>
	INLINE _Real linear_interpolation(_Real x0, _Real y0, _Real z0,
		_Real f000, _Real f001, _Real f010, _Real f011,
		_Real f100, _Real f101, _Real f110, _Real f111)
	{
		_Real x1 = 1. - x0;
		_Real y1 = 1. - y0;
		_Real z1 = 1. - z0;
		return ( ( x1 * f000 + x0 * f100 ) * y1 ) * z1
		     + ( ( x1 * f010 + x0 * f110 ) * y0 ) * z1
		     + ( ( x1 * f001 + x0 * f101 ) * y1 ) * z0
		     + ( ( x1 * f011 + x0 * f111 ) * y0 ) * z0;
	}

	/**
	*   Линейная интерполяция набора из n-функции заданных на 2D-сетке
	*   Для n>2 данная функция быстрее, чем n-кратный вызов функции интерполяции
	*   для одиночной функции.
	* @note Предполагается, что координаты узлов равны 0 и 1.
	*   Нумерация индексов такова, что быстрее всего меняется последний индекс.
	*   Для координат (x,y) также быстрее меняется последний индекс y.
	* @param n число функций
	* @param f набор из n аппроксимированных значений
	* @param f00,f01,f10,f11 значения функции в узлах сетки
	* @param x0,y0 относительное расстояние от точки до левого нижнего узла
	* @return аппроксимированное значение функции
	*/
	template <typename _Real>
	INLINE void linear_interpolation(unsigned n, _Real *f, _Real x0, _Real y0,
		const _Real *f00, const _Real *f01, const _Real *f10, const _Real *f11)
	{
		_Real x1 = 1. - x0;
		_Real y1 = 1. - y0;
		_Real xy00 = x0 * y0;
		_Real xy01 = x0 * y1;
		_Real xy10 = x1 * y0;
		_Real xy11 = x1 * y1;
		for (unsigned i=0; i<n; i++) f[i]
			= f00[i] * xy11 + f01[i] * xy10 + f10[i] * xy01 + f11[i] * xy00;
	}

	/**
	*   Линейная интерполяция набора из n-функции заданных на 3D-сетке
	*   Для n>2 данная функция быстрее, чем n-кратный вызов функции интерполяции
	*   для одиночной функции.
	* @note Предполагается, что координаты узлов равны 0 и 1.
	*   Нумерация индексов такова, что быстрее всего меняется последний индекс.
	*   Для координат (x,y,z) также быстрее меняется последний индекс z.
	* @param n число функций
	* @param f набор из n аппроксимированных значений
	* @param f000,f001,f010,f011... значения функции в узлах сетки
	* @param x0,y0,z0 относительное расстояние от точки до левого нижнего узла
	* @return аппроксимированное значение функции
	*/
	template <typename _Real>
	INLINE void linear_interpolation(unsigned n, _Real *f, _Real x0, _Real y0, _Real z0,
		const _Real *f000, const _Real *f001, const _Real *f010, const _Real *f011,
		const _Real *f100, const _Real *f101, const _Real *f110, const _Real *f111)
	{
		_Real x1 = 1. - x0;
		_Real y1 = 1. - y0;
		_Real z1 = 1. - z0;
		_Real xy00 = x0 * y0;
		_Real xy01 = x0 * y1;
		_Real xy10 = x1 * y0;
		_Real xy11 = x1 * y1;
		_Real xyz000 = xy00 * z0;
		_Real xyz001 = xy00 * z1;
		_Real xyz010 = xy01 * z0;
		_Real xyz011 = xy01 * z1;
		_Real xyz100 = xy10 * z0;
		_Real xyz101 = xy10 * z1;
		_Real xyz110 = xy11 * z0;
		_Real xyz111 = xy11 * z1;
		for (unsigned i=0; i<n; i++) f[i]
			= f000[i] * xyz111 + f001[i] * xyz110 + f010[i] * xyz101 + f011[i] * xyz100
		  + f100[i] * xyz011 + f101[i] * xyz010 + f110[i] * xyz001 + f111[i] * xyz000;
	}

	/**
	*  Рассчитывает интерполяцию функции вперед на равномерной сетке на один шаг,
	*  используя полином первого порядка, созданный по значениям функции и ее
	*  первой производной в точке x=0. Точность интерполяции o(h**2).
	* @param f0 значение функции в точке x=0
	* @param p0 значение первой производной функции в точке x=0
	* @param h шаг сетки
	* @return интерполированное значение функции в точке x=h
	*/
	template <typename S, typename _Real>
	INLINE S forward_interpolation(S f0, S p0, _Real h)
	{ return f0 + p0 * h; }

	/**
	*  Рассчитывает интерполяцию набора функции (в том числе и векторных)
	*  вперед на равномерной сетке на один шаг, используя полином первого порядка,
	*  созданный по значениям функции и ее первой производной в точке x=0.
	*  Точность интерполяции o(h**2).
	* @param n число функций
	* @param f[out] массив интерполированных значений функции в точке x=h
	* @param f0 массив значений функций в точке x=0
	* @param p0 массив значений первой производной функций в точке x=0
	* @param h шаг сетки
	*/
	template <typename S, typename _Real>
	INLINE void forward_interpolation(unsigned n, S *f,
		const S *f0, const S *p0, _Real h)
	{
		for (unsigned i=0; i<n; i++) f[i] = f0[i] + p0[i] * h;
	}

	/**
	*  Рассчитывает интерполяцию функции вперед на равномерной сетке на один шаг,
	*  используя полином третьего порядка, созданный по значениям функции и ее
	*  первой производной в точках x=-h и x=0. Точность интерполяции o(h**4).
	* @param f_1 значение функции в точке x=-h
	* @param f0 значение функции в точке x=0
	* @param p_1 значение первой производной функции в точке x=-h
	* @param p0 значение первой производной функции в точке x=0
	* @param h шаг сетки
	* @return интерполированное значение функции в точке x=h
	*/
	template <typename S, typename _Real>
	INLINE S forward_interpolation(S f_1, S f0, S p_1, S p0, _Real h)
	{
		return f_1 + 4 * (f_1 - f0) + (p_1 + p0 + p0) * (h + h);
	}

	/**
	*  Рассчитывает интерполяцию набора функции вперед на равномерной сетке на один шаг,
	*  используя полином третьего порядка, созданный по значениям функции и ее
	*  первой производной в точках x=-h и x=0. Точность интерполяции o(h**4).
	* @param n число функций
	* @param f[out] массив интерполированных значений функции в точке x=h
	* @param f_1 массив значений функций в точке x=-h
	* @param f0 массив значений функций в точке x=0
	* @param p_1 массив значений первой производной функций в точке x=-h
	* @param p0 массив значений первой производной функций в точке x=0
	* @param h шаг сетки
	*/
	template <typename S, typename _Real>
	INLINE void forward_interpolation(unsigned n, S *f,
		const S *f_1, const S *f0, const S *p_1, const S *p0, _Real h)
	{
		_Real hh = h + h;
		for (unsigned i=0; i<n; i++) f[i] =
			f_1[i] + (f_1[i] + -f0[i]) * 4 + (p_1[i] + p0[i] + p0[i]) * hh;
	}

	/**
	*  Рассчитывает интерполяцию функции вперед на равномерной сетке на один шаг,
	*  используя полином пятого порядка, созданный по значениям функции и ее
	*  первой производной в точках x={-2h,-h, 0}. Точность интерполяции o(h**6).
	* @param f_2 значение функции в точке x=-2h
	* @param f_1 значение функции в точке x=-h
	* @param f0 значение функции в точке x=0
	* @param p_2 значение первой производной функции в точке x=-2h
	* @param p_1 значение первой производной функции в точке x=-h
	* @param p0 значение первой производной функции в точке x=0
	* @param h шаг сетки
	* @return интерполированное значение функции в точке x=h
	*/
	template <typename S, typename _Real>
	INLINE S forward_interpolation(S f_2, S f_1, S f0, S p_2, S p_1, S p0, _Real h)
	{
		return f_2 + 9 * (f_2 + f_1 - mult2(f0))
			+ (p_2 + 3 * (p_1 + p_1 + p0)) * 3 * h;
	}

	/**
	*  Рассчитывает интерполяцию функции вперед на равномерной сетке на один шаг,
	*  используя полином пятого порядка, созданный по значениям функции и ее
	*  первой производной в точках x={-2h,-h, 0}. Точность интерполяции o(h**6).
	* @param n число функций
	* @param f[out] массив интерполированных значений функции в точке x=h
	* @param f_2 массив значений функций в точке x=-2h
	* @param f_1 массив значений функций в точке x=-h
	* @param f0 массив значений функций в точке x=0
	* @param p_2 массив значений первой производной функций в точке x=-2h
	* @param p_1 массив значений первой производной функций в точке x=-h
	* @param p0 массив значений первой производной функций в точке x=0
	* @param h шаг сетки
	*/
	template <typename S, typename _Real>
	INLINE void forward_interpolation(unsigned n, S *f,
		const S *f_2, const S *f_1, const S *f0,
		const S *p_2, const S *p_1, const S *p0, _Real h)
	{
		_Real h3 = h + h + h;
		for (unsigned i=0; i<n; i++) f[i] =
			f_2[i] + 9 * (f_2[i] + f_1[i] - mult2(f0[i]))
			+ (p_2[i] + 3 * (p_1[i] + p_1[i] + p0[i])) * h3;
	}

	template <typename T>
	struct Range_
	{
		Range_() : min(T()), max(infinity<T>()) {}
		Range_(T a, T b) : min(a), max(b) {}
		Range_(const Range_ &range) : min(range.min), max(range.max) {}
		T min, max;
	};

}
#undef bll
#endif
