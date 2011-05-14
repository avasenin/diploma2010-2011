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
#ifndef _MINIMIZE__0077A726_478A_5b23_91D7_D9448B790B00__H
#define _MINIMIZE__0077A726_478A_5b23_91D7_D9448B790B00__H

#include "prgkern/_prgconfig.h"
#include "prgkern/_type.h"

//#define OPTIMIZATION_DEBUG
//#define DEGENERACY_PRINT

#ifdef OPTIMIZATION_DEBUG
#define OPTIMIZATION_PRINT(msg) \
	{ \
		_S msg__ = make_string(" x[%12.5e %12.5e]  f[%12.5e %12.5e]  g[%12.5e %12.5e]", \
				xa, xb, fa / n, fb / n, ga / n, gb / n); \
		PRINT_MESSAGE(_S(msg) + _S(msg__)); \
	}
#else
	#define OPTIMIZATION_PRINT(msg)
#endif


namespace prgkern
{

	/// after OPTIMIZE_PRINT_FREQUENCY will be printed all information
	const int OPTIMIZE_PRINT_FREQUENCY = 1;
	const int MAX_ITERATION = 2000;

	/// число итераций выполняемых самым STEEP методом вначале работы любого минимизатора
	const int DEFAULT_STEEP_ITERATIONS = 10;

	/// optimizer step for first line search
	const float MINIMIZER_STEP_START = 1e-6;

	template <int TYPE> class Minimizer_;
	enum { SCITBX_LBFGS_, GSL_SIMPLEX_, STEEP_, LMBFGS_, ROUND_ };

	/// минимально разрешенное смещение в x-пространстве
	const float DEFAULT_STPMIN = 0.001;

	/// максимально разрешенное смещение в x-пространстве
	const float DEFAULT_STPMAX = 0.100;

	/// Принимает как допустимый минимум, если "падение" энергии оказалось в пределах DEFAULT_WOLFE1
	const float DEFAULT_WOLFE1 = 0.1;

	/// Принимает как допустимый минимум, если "падение" градиента оказалось в пределах DEFAULT_WOLFE2
	const float DEFAULT_WOLFE2 = 0.1;

	/// Максимальное число вызовов внешней функции, которая считает энергию и градиент
	const unsigned DEFAULT_MAXFEV = 200; // need large value ~200 for first iteration

	/// число последовательных попыток прервать итерационный процесс
	const unsigned DEFAULT_MAX_HALTS = 3;

	/**
	* @param m The number of corrections used in the BFGS update.
	*   Values of m less than 3 are not recommended;
	*   large values of m will result in excessive computing time.
	*   3 <= m <= 7 is recommended by authors, but you can use upto 80.
	*/
	const unsigned DEFAULT_M = 7;

	/// предельная точность локализации X
	template <typename T>
	INLINE T DEFAULT_XTOL() { return (T) 10. * std::numeric_limits<T>::epsilon(); }

	/// предельная точность локализации энергии
	template <typename T>
	INLINE T DEFAULT_FTOL() { return (T) 10. * std::numeric_limits<T>::epsilon(); }

	/// предельная точность локализации градиента
	template <typename T>
	INLINE T DEFAULT_GTOL() { return 1e-3; }

	/**
	 * @param fn - внешняя функция для расчета энергии и градиента
	 * @param param - параметры внешней функции
	 * @param n  - число степеней свободы
	 * @param x_ - конечные координаты системы
	 * @param g_ - градиент в конечных координатах
	 * @param x0 - стартовые координаты системы (центр окружности)
	 * @param g0 - градиент в стартовых координатах системы
	 * @param d  - направление движения
	 * @param nfev - число оценок энергии
	 * @param wolfe1 - условие Вольфа завершения по значению функции
	 * @param wolfe2 - условие Вольфа завершения по значению градиента
	 * @param xtol   - предел согласования по координатам
	 * @param ftol   - предел согласования по значению функции
	 * @param gtol   - предел согласования по значению градиента
	 */
	#define SHARE_PARAM_DEF(T)  T (*fn)(unsigned, const T *, T *, void *param), void *param, \
		unsigned n, T *x_, T *g_, const T *x, const T *g, T *d, \
		int &nfev, T wolfe1=DEFAULT_WOLFE1, T wolfe2=DEFAULT_WOLFE2, \
		T xtol=DEFAULT_XTOL<T>(), T ftol=DEFAULT_FTOL<T>(), T gtol=DEFAULT_GTOL<T>()

	#define SHARE_PARAM_VAL   fn, param, n, x_, g_, x, g, d, \
		nfev, wolfe1, wolfe2, xtol, ftol, gtol

	#define MAKE_EQ(a, b) { x##a = x##b; f##a = f##b; g##a = g##b; }

	/// line search optimizer
	class Line_minimizer_
	{

	protected:

		enum { FAIL, INTERVAL_OK, FUNCTION_DIFFERENCE_OK, WOLFE_OK, TOLERANCE_OK };

		/**
		 * Сосчитать пробную точку на заданной прямой
		 * @param xu - координаты пробной точки
		 * @param fu - значение энергии в точке
		 * @param gu - проекция градиента в точке на заданную прямую
		 */
		template <typename T> void calculate_fgt(bool find_minimum, T &xu, T &fu, T &gu, SHARE_PARAM_DEF(T))
		{
			VECTOR_EXPRESSION_2(T, n, x_, =, , , x, +, xu, *, d); // x_ = x + b * d
			fu = (*fn)(n, x_, g_, param);
			if (--nfev < 0) { PRINT_BREAK(_S("[ERROR] exceed the limit of f() evaluations")); }
			gu = scalar_product(n, g_, d);
			if (!find_minimum) { fu = -fu; gu = -gu; }
		}

		/**
		*  Позволяет избежать точного поиска минимума на заданном направлении,
		*  что значительно (в десяток раз) ускоряет расчеты за счет пренебрежения
		*  вычислениями в областях далеких от минимума.
		*  Минимум определяется с точностью, которая задается двумя параметрами:
		*   wolfe1 - точность нахождения значения функции
		*   wolfe2 - точность нахождения производной
		*  Это не настоящий WOLFE. Я изменил принцип работы с wolfe1 параметром.
		*  Мой вариант оценивает положение минимума и значение функции в нем.
		*  Условие удовлетворяется, если в точке xb значение функции fb ближе к
		*  минимуму в 1/wolfe1 раз, чем значение fa
		*
		* @param n degree of freedom
		* @param nfev number of f() evaluations
		* @param ax,bx,cx interval of minimum localization
		* @param fa,fb,fc function values on interval
		* @param ga,gb,gc derivation values on interval
		* @param wolfe1,wolfe2 Wolfe conditions of exit of line search
		* @param tol max accuracy of minima localization
		* @param f function value at start
		* @param gd scalar_product(grad, direction) at start
		*/
		template <typename T>
		bool test_Wolfe_(unsigned n, T xa, T fa, T ga, T xb, T fb, T gb, T wolfe1, T wolfe2) const
		{
			if (fabs(gb) < wolfe2 * fabs(ga))
			{
			#ifdef OPTIMIZATION_DEBUG
				PRINT_MESSAGE("wolf2 ok");
			#endif
				return true;
			}

			if ( (fa > 0.f && fb > 0.f && fb < fa * wolfe1)
				|| (fa < 0.f && fb < 0.f && fb * wolfe1 < fa)
				|| (fa > 0.f && fb < 0.f)
			)
			{
			#ifdef OPTIMIZATION_DEBUG
				PRINT_MESSAGE("wolf0 ok");
			#endif
				return true;
			}

			if (xa < xb && fabs(ga) > fabs(gb))
			{
				T A, B;
				T xu = zero_linear(A, B, xa, xb, ga, gb);
				T dfa = (T)0.5 * A * (sqr(xu) - sqr(xa)) + B * (xu - xa);

				if (fb <= fa + dfa * ((T)1. - wolfe1))
				{
				#ifdef OPTIMIZATION_DEBUG
					PRINT_MESSAGE("wolf1 ok");
				#endif
					return true;
				}
			}

			return false;
		}

		template <typename T>
		bool test_tolerance_(unsigned n, T xa, T fa, T ga, T xb, T fb, T gb, T xtol, T ftol, T gtol) const
		{
			if ( fabs(xb - xa) <= xtol
				|| fabs(fb - fa) <= std::max(1., fabs(fa)) * ftol
				|| fabs(gb - ga) <= std::max(1., fabs(ga)) * gtol
			)
			{
			#ifdef OPTIMIZATION_DEBUG
				PRINT_MESSAGE("tolerance ok");
			#endif
				return true;
			}
			return false;
		}

		template <typename T>
		int dec_interval_(bool find_minimum, T &beta, T &f0, T &g0, T &xa, T &fa, T &ga, T &xb, T &fb, T &gb,
			SHARE_PARAM_DEF(T))
		{
			bool good = (xa < xb) && (ga < 0.f) && (fa < fb);
			if (!good)
			{
				PRINT_MESSAGE(" test of condition: (xa < xb) && (ga < 0) && (fa < fb)");
				OPTIMIZATION_PRINT("range-:");
				PRINT_BREAK(_S("[ERROR] bad start values"));
			}

			T xu, fu, gu; const T SHIFT = 0.1;
			while (fa < fb)
			{
				xu = xa + SHIFT * (xb - xa);
				calculate_fgt(find_minimum, xu, fu, gu, SHARE_PARAM_VAL);

				MAKE_EQ(b, u);
				OPTIMIZATION_PRINT("range-:");

				if (test_Wolfe_(n, T(0.), f0, g0, xb, fb, gb, wolfe1, wolfe2))
				{
					beta = xb; f0 = fb; g0 = gb;
					return WOLFE_OK;
				}
				if (sign(ga) != sign(gb)) return INTERVAL_OK;
			}
			return INTERVAL_OK;
		}

		template <typename T>
		int inc_interval_(bool find_minimum, T &beta, T &f0, T &g0, T &xa, T &fa, T &ga, T &xb, T &fb, T &gb,
			SHARE_PARAM_DEF(T))
		{
			bool good = (xa < xb) && (ga < 0) && (gb < 0);
			if (!good)
			{
				PRINT_MESSAGE(" test of condition: (xa < xb) && (ga < 0) && (gb < 0)");
				OPTIMIZATION_PRINT("range+:");
				PRINT_BREAK(_S("[ERROR] bad start values"));
			}

			const T SHIFT_LARGE = 10.f, SHIFT_SMALL = 2.f;
			T xu, fu, gu;
			while (gb < 0.f)
			{
				xu = zero_linear(xa, xb, ga, gb); // линейная аппроксимация нуля

				T ca = xb; // - (xa = 0);
				T sca = ca * SHIFT_SMALL;
				T lca = ca * SHIFT_LARGE;
				if (xu - xb < sca) xu = xb + sca;
				if (xu - xb > lca) xu = xb + lca;

				calculate_fgt(find_minimum, xu, fu, gu, SHARE_PARAM_VAL);

				MAKE_EQ(a, b);
				MAKE_EQ(b, u);
				OPTIMIZATION_PRINT("range+:");

				while (fu > fa && gu < 0.f)
				{
					xu = 0.5f * (xa + xu);
					calculate_fgt(find_minimum, xu, fu, gu, SHARE_PARAM_VAL);
					MAKE_EQ(b, u);
					OPTIMIZATION_PRINT("rang++:");

					if (test_Wolfe_(n, T(0.), f0, g0, xb, fb, gb, wolfe1, wolfe2))
					{
						beta = xb; f0 = fb; g0 = gb;
						return WOLFE_OK;
					}
				}

				if (test_Wolfe_(n, T(0.), f0, g0, xb, fb, gb, wolfe1, wolfe2))
				{
					beta = xb; f0 = fb; g0 = gb;
					return WOLFE_OK;
				}
			}
			return INTERVAL_OK;
		}

		template <typename T>
		void locate_zero_derivation_(bool find_minimum, T &beta, T &f0, T &g0, T &xa, T &fa, T &ga, T &xb, T &fb, T &gb,
			SHARE_PARAM_DEF(T))
		{
			bool good = (xa < xb) && (ga <= 0) && (gb >= 0);
			if (!good)
			{
				PRINT_MESSAGE(" test of condition: (xa < xb) && (ga <= 0) && (gb >= 0)");
				OPTIMIZATION_PRINT("gzero :");
				PRINT_BREAK(_S("[ERROR] bad start values"));
			}

			const T SHIFT = 0.1;
			T xu, fu, gu, xl, xr, fl, gl, fr, gr, ca, ca__;
			MAKE_EQ(l, a);
			MAKE_EQ(r, b);

			xu = 0.5 * (xa + xb);
			calculate_fgt(find_minimum, xu, fu, gu, SHARE_PARAM_VAL);

			if (gu < 0.f) { MAKE_EQ(l, a); MAKE_EQ(a, u); }
			else { MAKE_EQ(r, b); MAKE_EQ(b, u); }
			OPTIMIZATION_PRINT("gzero :");

			while (true)
			{
				if (gu < 0.f) xu = zero_quadratic(xl, xa, xb, gl, ga, gb);
				else					xu = zero_quadratic(xa, xb, xr, ga, gb, gr);

				ca = xb - xa;
				ca__ = ca * SHIFT;
				if (xu - xa < ca__) xu = xa + ca__;
				if (xb - xu < ca__) xu = xb - ca__;

				calculate_fgt(find_minimum, xu, fu, gu, SHARE_PARAM_VAL);
				if (gu < 0.f) { MAKE_EQ(l, a); MAKE_EQ(a, u); }
				else { MAKE_EQ(r, b); MAKE_EQ(b, u); }
				OPTIMIZATION_PRINT("gzero :");

				if ( test_Wolfe_(n, T(0.), f0, g0, xu, fu, gu, wolfe1, wolfe2)
					|| test_tolerance_(n, xa, fa, ga, xb, fb, gb, xtol, ftol, gtol)
				)
				{
					beta = xu;
					f0 = fu;
					g0 = gu;
					break;
				}
			}
		}

	public:

		/**
		* @brief runs optimization
		* @param x - general coords at start[in] & at finish[out]
		* @param d - direction
		* @param g - gradient at start[in] & at finish[out]
		* @return f(x)
		*/
		template <typename T>
		T operator()(bool find_minimum, T &beta, T &f0, T &g0, SHARE_PARAM_DEF(T))
		{
			if (!find_minimum) { f0 = -f0; g0 = -g0; }

			// смена направления в сторону уменьшения функции
			if (g0 > 0.f)
			{
				VECTOR_EXPRESSION_1(T, n, d, =, -1.f, *, d);
				g0 = -g0; // при смене направления также меняется знак градиента
			}

			T xa = 0.f, fa = f0, ga=g0;
			T xb = beta, fb, gb;

			xb = beta;
			calculate_fgt(find_minimum, xb, fb, gb, SHARE_PARAM_VAL);
			OPTIMIZATION_PRINT("start :");

			int result = FAIL;
			if (fb > fa)
			{
				result = dec_interval_(find_minimum, beta, f0, g0, xa, fa, ga, xb, fb, gb, SHARE_PARAM_VAL);
				if (result == WOLFE_OK)
				{
					if (!find_minimum) { f0 = -f0; g0 = -g0; }
					return f0;
				}
			}

			if (fb < fa && gb < 0.)
			{
				result = inc_interval_(find_minimum, beta, f0, g0, xa, fa, ga, xb, fb, gb, SHARE_PARAM_VAL);
				if (result == WOLFE_OK)
				{
					if (!find_minimum) { f0 = -f0; g0 = -g0; }
					return f0;
				}
			}

			locate_zero_derivation_(find_minimum, beta, f0, g0, xa, fa, ga, xb, fb, gb, SHARE_PARAM_VAL);

			if (!find_minimum) { f0 = -f0; g0 = -g0; }
			return f0;
		}

	};

	#undef SHARE_PARAM_DEF
	#undef SHARE_PARAM_VAL
	#undef MAKE_EQ
}
#endif
