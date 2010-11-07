#ifndef _OPTIMIZER_ROUND_0077A726_0C5C_5a4e_DCD8_D144196A0C01__H
#define _OPTIMIZER_ROUND_0077A726_0C5C_5a4e_DCD8_D144196A0C01__H

#include "prgkern/_prgconfig.h"
#include "prgkern/_math.h"
#include "prgkern/_minimize.h"
#include "prgkern/_blas1.h"
#include <vector>
#include <limits>

namespace prgkern
{

	/// round search optimizer
	class Round_minimizer_
	{

		enum { FAIL, INTERVAL_OK, WOLFE_OK, TOLERANCE_OK };

		#define MAKE_EQ(a, b) { x##a = x##b; f##a = f##b; g##a = g##b; }

		/**
		 * @param fn - внешняя функция для расчета энергии и градиента
		 * @param param - параметры внешней функции
		 * @param n  - число степеней свободы
		 * @param x_ - конечные координаты системы
		 * @param g_ - градиент в конечных координатах
		 * @param r_ - всмогательный массив для хранения радиус векторов
		 * @param x0 - стартовые координаты системы (центр окружности)
		 * @param g0 - градиент в стартовых координатах системы
		 * @param r  - длина вектора, определяющего окружность
		 * @param ru - единичный вектор, указывающий в точку 0 окружности
		 * @param rn - ортогональный к ru единичный вектор
		 * @param nfev - число оценок энергии
		 * @param wolfe1 - условие Вольфа завершения по значению функции
		 * @param wolfe2 - условие Вольфа завершения по значению градиента
		 * @param xtol   - предел согласования по координатам
		 * @param ftol   - предел согласования по значению функции
		 * @param gtol   - предел согласования по значению градиента
		 */
		#define SHARE_PARAM_DEF(T)  T (*fn)(unsigned, const T *, T *, void *param), void *param, \
			unsigned n, T *X_, T *G_, T *R_, \
			const T *X0, const T *G0, T rlen, const T *Ru, const T *Rn, \
			int &nfev, T wolfe1=DEFAULT_WOLFE1, T wolfe2=DEFAULT_WOLFE2, \
			T xtol=DEFAULT_XTOL<T>(), T ftol=DEFAULT_FTOL<T>(), T gtol=DEFAULT_GTOL<T>()

		#define SHARE_PARAM_VAL   fn, param, n, X_, G_, R_, X0, G0, rlen, Ru, Rn, \
			nfev, wolfe1, wolfe2, xtol, ftol, gtol

			/**
			*  Позволяет избежать точного поиска минимума на заданном направлении, что значительно
			*  (в десяток раз) ускоряет расчеты за счет пренебрежения вычислениями в областях далеких
			*  от минимума. Минимум определяется с точностью, которая задается двумя параметрами:
			*   wolfe1 - точность нахождения значения функции
			*   wolfe2 - точность нахождения производной
			*  Это не настоящий WOLFE. Я изменил принцип работы с wolfe1 параметром.
			*  Мой вариант оценивает положение минимума и значение функции в нем.
			*  Условие удовлетворяется, если в точке xb значение функции fb ближе к
			*  минимуму в 1/wolfe1 раз, чем значение fa
			*
			* @param x0,f0,g0 параметры с которыми сравниваются найденные значения
			* @param x,f,g найденные значения
			* @param wolfe1,wolfe2 условия Вольфа
			*/
			template <typename T>
			bool test_Wolfe_(unsigned n, T xa, T fa, T ga, T xb, T fb, T gb, T wolfe1, T wolfe2) const
			{
				if (fa < fb || fabs(ga) < fabs(gb)) return false;

				// стратегия основывается только на проверке градиента
				if (fabs(gb)  < wolfe2 * fabs(ga))
				{
				#ifdef OPTIMIZATION_DEBUG
					PRINT_MESSAGE("wolf2 ok");
				#endif
					return true;
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
					OPTIMIZATION_PRINT("toler :");
					return true;
				}
				return false;
			}

		/**
		 * Сосчитать компоненту градиента вдоль окружности
		 * @param xu - угол поворота вдоль окружности
		 * @param fu - значение энергии в точке xu
		 * @param gu - касательная компонента градиента в точке равной углу поворота
		 */
		template <typename T> void calculate_fgt(T &xu, T &fu, T &gu, SHARE_PARAM_DEF(T))
		{
			T sb = sin(xu), cb = cos(xu);
			T rsb = rlen * sb, rcb = rlen * cb;

			VECTOR_EXPRESSION_2(T, n, R_, =, rcb, *, Ru, +, rsb, *, Rn);  // r_ = r * (cb * ru + sb * rn)
				// поворот единичного вектора ru

			VECTOR_EXPRESSION_2(T, n, X_, =, , , X0, +, , , R_);        // x_ = x0 + r_
				// конечная координата системы

			fu = (*fn)(n, X_, G_, param);
			if (--nfev < 0) { PRINT_BREAK(_S("[ERROR] exceed the limit of f() evaluations")); }

			VECTOR_EXPRESSION_2(T, n, R_, =, cb, *, Rn, -, sb, *, Ru);  // r_ = cb * rn - sb * ru
				// поворот ортогонального вектора rn

			gu = scalar_product(n, G_, R_);
				// касательная компонента градиента в конечной точке
		}

		/**
		 * Находит знакопеременный интервал
		 * @param x0,f0,g0 - значения координаты, функции и градиента для сравнения по Вольфу
		 * @param xa,fa,ga - значения в левой точке интервала
		 * @param xa,fa,ga - значения в правой точке интервала
		 * @return критерий обнаружения
		 */
		template <typename T>
		int inc_interval_(T &beta, T &f0, T &g0, T &xa, T &fa, T &ga, T &xb, T &fb, T &gb,
			SHARE_PARAM_DEF(T))
		{
			T xu, fu, gu;
			const T SHIFT_LARGE = 2.f, SHIFT_SMALL = 0.1f;
			int sign_beta = sign(beta);

			while (true)
			{
				xu = zero_linear(xa, xb, ga, gb); // линейная аппроксимация нуля

				T ca = fabs(xa - xb);
				T sca = ca * SHIFT_SMALL;
				T lca = ca * SHIFT_LARGE;

				if (fabs(xb - xu) < sca ) xu = xb + sign_beta * sca;
				if (fabs(xb - xu) > lca ) xu = xb + sign_beta * lca;
				if (fabs(xb - xa) > M_PI) xu = xa + sign_beta * M_PI;

				calculate_fgt(xu, fu, gu, SHARE_PARAM_VAL);
				MAKE_EQ(a, b); MAKE_EQ(b, u);
				OPTIMIZATION_PRINT("range+:");

				if (test_Wolfe_(n, 0.f, f0, g0, xu, fu, gu, wolfe1, wolfe2))
				{
					beta = xu; f0 = fu; g0 = gu;
					return WOLFE_OK;
				}

				if (sign(ga) != sign(gb)) return INTERVAL_OK;
			}

			return FAIL;
		}

		template <typename T>
		void locate_zero_derivation_(T &beta, T &f0, T &g0, T &xa, T &fa, T &ga, T &xb, T &fb, T &gb,
			SHARE_PARAM_DEF(T))
		{
			if (xa > xb) { std::swap(xa, xb); std::swap(fa, fb); std::swap(ga, gb); }

			bool good = sign(ga) != sign(gb);
			if (!good)
			{
				PRINT_MESSAGE(" test of condition: sign(ga) != sign(gb)");
				OPTIMIZATION_PRINT("gzero :");
				PRINT_BREAK(_S("[ERROR] bad start values"));
			}

			const T SHIFT = 0.1;
			T xu, fu, gu, xl, xr, fl, gl, fr, gr, ca, ca__;
			MAKE_EQ(l, a);
			MAKE_EQ(r, b);

			xu = 0.5 * (xa + xb);
			calculate_fgt(xu, fu, gu, SHARE_PARAM_VAL);

			if (sign(gu) == sign(ga)) { MAKE_EQ(l, a); MAKE_EQ(a, u); }
			else { MAKE_EQ(r, b); MAKE_EQ(b, u); }
			OPTIMIZATION_PRINT("gzero :");

			while (true)
			{
				if (sign(gu) == sign(ga)) xu = zero_quadratic(xl, xa, xb, gl, ga, gb);
				else					            xu = zero_quadratic(xa, xb, xr, ga, gb, gr);

				ca = xb - xa;
				ca__ = ca * SHIFT;
				if (xu - xa < ca__) xu = xa + ca__;
				if (xb - xu < ca__) xu = xb - ca__;

				calculate_fgt(xu, fu, gu, SHARE_PARAM_VAL);
				if (sign(gu) == sign(ga)) { MAKE_EQ(l, a); MAKE_EQ(a, u); }
				else                      { MAKE_EQ(r, b); MAKE_EQ(b, u); }
				OPTIMIZATION_PRINT("gzero :");

				if ( test_Wolfe_(n, 0.f, f0, g0, xu, fu, gu, wolfe1, wolfe2)
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
		 * Выполнить минимизацию на окружности заданного радиуса
		 * @param xu,fu,gu - значения координаты, функции и касательной градиента в найденной точке
		 */
		template <typename T>
		void operator()(T &beta, T &f0, T &g0, SHARE_PARAM_DEF(T))
		{
			T xa = 0.f, fa = f0, ga = g0;
			if (sign(g0 * beta) > 0) beta = -beta;
				// смещаемся в сторону убывания функции, меняя знак beta, иначе возможно локализовать
				// максимум на окружности вместо минимума

			T xb = beta, fb, gb;
			calculate_fgt(xb, fb, gb, SHARE_PARAM_VAL);
			OPTIMIZATION_PRINT("start :");

			int found = FAIL;
			if (sign(ga) == sign(gb))
			{
				found = inc_interval_(beta, f0, g0, xa, fa, ga, xb, fb, gb, SHARE_PARAM_VAL);
				if (found == WOLFE_OK) return;
			}
			locate_zero_derivation_(beta, f0, g0, xa, fa, ga, xb, fb, gb, SHARE_PARAM_VAL);
		}
	};


	/**
	* @brief Minimizer
	*/
	template <> class Minimizer_<ROUND_>
	{
		template <typename T>
		T round_minimize_(T &bx, T &f0, T &g0, T (*fn)(unsigned, const T *, T *, void *param), void *param,
			unsigned n, T *X_, T *G_, T *D_, const T *X0, const T *G0,
			unsigned max_iter=MAX_ITERATION, T stpmin=DEFAULT_STPMIN, T stpmax=DEFAULT_STPMAX,
			unsigned maxfev=DEFAULT_MAXFEV, unsigned maxhalt=DEFAULT_MAX_HALTS,
			T wolfe1=DEFAULT_WOLFE1, T wolfe2=DEFAULT_WOLFE2,
			T xtol=DEFAULT_XTOL<T>(), T ftol=DEFAULT_FTOL<T>(), T gtol=DEFAULT_GTOL<T>()
		)
		{
			std::vector<T> Ru_(n), Rn_(n);
			T *Ru = &Ru_[0], *Rn = &Rn_[0];

			//--------------------------------------------------------------------------------
			// строим едининый вектор Ru, определяющий начальную точку окружности минимизации
			//--------------------------------------------------------------------------------
			T s = 1.f / sqrt(scalar_product(n, G0, G0));
			VECTOR_EXPRESSION_1(T, n, Ru, =, s, *, G0);

			// установим длину вектора такую, чтобы максимальное смещение в x-пространстве
			// было равно MAX_X
			T rlen = 0.f;
			for (unsigned i=0; i<n; i++)
				if (rlen < fabs(Ru[i])) rlen = fabs(Ru[i]);
			rlen = stpmax / rlen;

			//--------------------------------------------------------------------------
			// строим касательный вектор Rn к начальной точке окружности минимизации
			//--------------------------------------------------------------------------
			VECTOR_EXPRESSION_2(T, n, X_, =, , , X0, +, rlen, *, Ru);
			f0 = (*fn)(n, X_, Rn, param);

			s = scalar_product(n, Rn, Ru);
			VECTOR_EXPRESSION_2(T, n, Rn, =, , , Rn, -, s, *, Ru);
			g0 = normalize(n, Rn);

			{
				_S msg = _S("\n        ***** _____ start of round minimization ___ *****\n");
				msg += make_string("    iteration    0  f() = %12.5e    |f'()| = %12.5e  nfev =  1",
						f0 / n, g0 / n);
				PRINT_MESSAGE(msg);
			}

			int nfev, full_nfev = 0; // number of evaluation of f()
			T fprev = T(), gprev = T();

			unsigned finish_flag = 0; // определяет завершение работы, если условия
				// завершения случаются дважды (трижды) на двух последующих итерациях

			unsigned iter__ = 0;
			T beta = M_PI * 0.1;
			for(;;)
			{
				// save previous value
				fprev = f0; gprev = g0;
				nfev = maxfev;

				if (beta >  M_PI * 0.1) beta =  M_PI * 0.1;
				if (beta < -M_PI * 0.1) beta = -M_PI * 0.1;

				Round_minimizer_()(beta, f0, g0, fn, param, n, X_, G_, D_, X0, G0, rlen, Ru, Rn,
					nfev, wolfe1, wolfe2, xtol, ftol, gtol);
				full_nfev += maxfev - nfev;

				if (iter__ % OPTIMIZE_PRINT_FREQUENCY == 0)
				{
					_S msg = make_string("    iteration %4u  f() = %12.5e    |f'()| = %12.5e  nfev = %2d",
							iter__, f0 / n, g0 / n, (int)(maxfev - nfev));
					PRINT_MESSAGE(msg);
				}

				// рассчитаем условия завершения и если они повторятся два раза подряд
				// завершим оптимизацию
				bool halt_flag = (fabs(beta) < xtol)
					|| (fabs(f0 - fprev) < std::max(1., fabs(fprev)) * ftol)
					|| (fabs(g0 - gprev) < std::max(1., fabs(gprev)) * gtol);

				if (halt_flag) ++finish_flag; else finish_flag = 0;
				if (finish_flag > DEFAULT_MAX_HALTS) break;
				if (++iter__ >= max_iter) break;

				//--------------------------------------------------------------------------
				// строим вектор Ru, определяющий начальную точку окружности минимизации
				//--------------------------------------------------------------------------
				VECTOR_EXPRESSION_2(T, n, Ru, =, , , X_, -, , , X0); // Ru = X_ - X0
				normalize(n, Ru);

				//--------------------------------------------------------------------------
				// строим касательный вектор Rn к начальной точке окружности минимизации
				//--------------------------------------------------------------------------
				VECTOR_EXPRESSION_2(T, n, X_, =, , , X0, +, rlen, *, Ru);
				f0 = (*fn)(n, X_, Rn, param);

				s = scalar_product(n, Rn, Ru);
				VECTOR_EXPRESSION_2(T, n, Rn, =, , , Rn, -, s, *, Ru);

				g0 = normalize(n, Rn);
					// усредненное на число степеней свободы
			}

			{
				_S msg = _S("\n        ***** _______ finish of minimization ______ *****\n");
				msg += make_string("    iteration %4u  f() = %12.5e    |f'()| = %12.5e  full nfev = %4d\n",
						iter__, f0 /n, g0 /n, full_nfev);
				msg += make_string("    |dx| = %12.5e  |df| = %12.5e  |dg| = %12.5e\n",
						fabs(beta), fabs(f0 - fprev) / n, fabs(g0 - gprev) / n);
				PRINT_MESSAGE(msg);
			}

			bx = rlen;
			return f0;
		}

		template <typename T>
		T n_1_minimize_(T &bx, T &f0, T &g0, T (*fn)(unsigned, const T *, T *, void *param), void *param,
			unsigned n, T *X_, T *G_, T *D0, T *X0, T *G0,
			unsigned max_iter=MAX_ITERATION, T stpmin=DEFAULT_STPMIN, T stpmax=DEFAULT_STPMAX,
			unsigned maxfev=DEFAULT_MAXFEV, unsigned maxhalt=DEFAULT_MAX_HALTS,
			T wolfe1=DEFAULT_WOLFE1, T wolfe2=DEFAULT_WOLFE2,
			T xtol=DEFAULT_XTOL<T>(), T ftol=DEFAULT_FTOL<T>(), T gtol=DEFAULT_GTOL<T>()
		)
		{
			normalize(n, D0);

			std::vector<T> D__(n); T *D = &D__[0];

			unsigned iter__ = 0;
			{
				_S msg = _S("\n        ***** _______ start of n - 1 minimization _______ *****\n");
				msg += make_string("    iteration    0  f() = %12.5e    |f'()| = %12.5e  nfev =  1",
						f0 / n, g0 / n);
				PRINT_MESSAGE(msg);
			}

			int nfev, full_nfev = 0; // number of evaluation of f()
			T fprev = f0, gprev = g0;
			T maxd = 0.f;

			unsigned finish_flag = 0; // определяет завершение работы, если условия
				// завершения случаются дважды (трижды) на двух последующих итерациях

			for(;;)
			{
				fprev = f0; gprev = g0;
				nfev = maxfev;

				T s = scalar_product(n, G0, D0);
				VECTOR_EXPRESSION_2(T, n, D, =, , , G0, -, s, *, G0);
				g0 = scalar_product(n, G0, D);

				// масштабируем bx таким образом, чтобы произведение bx * D[i]
				// не превысило stpmax (максимальное смещение атома по любому направлению)
				maxd = sqrt(scalar_product(n, D, D) / n);
				if (maxd <= (T) 10. * std::numeric_limits<T>::min()) // нет никакого направления движения
				{
					_S msg = _S("Can't find direction for minimization");
					PRINT_MESSAGE(msg);
					break;
				}
				if (stpmax < bx * maxd) bx = stpmax / maxd;
				if (stpmin > bx * maxd) bx = stpmin / maxd;

				f0 = Line_minimizer_()(true, bx, f0, g0, fn, param, n, X_, G_, X0, G0, D,
					nfev, wolfe1, wolfe2, xtol, ftol, gtol);

				if (++iter__ % OPTIMIZE_PRINT_FREQUENCY == 0)
				{
					_S msg = make_string("    iteration %4u  f() = %12.5e    |f'()| = %12.5e  nfev = %2d",
						iter__, f0 / n, g0 / n, (int)(maxfev - nfev));
					PRINT_MESSAGE(msg);
				}

				if (fabs(bx) <= (T) 10. * std::numeric_limits<T>::min())
				{
					_S msg = _S("Can't move from given point");
					PRINT_MESSAGE(msg);
					break;
				}

				VECTOR_EXPRESSION_1(T, n, X0, =, , , X_);
				VECTOR_EXPRESSION_1(T, n, G0, =, , , G_);

				full_nfev += maxfev - nfev;

				// рассчитаем условия завершения и если они повторятся два раза подряд
				// завершим оптимизацию
				bool halt_flag = (fabs(maxd * bx) < stpmin)
					|| (fabs(f0 - fprev) / n < std::max(1., fabs(f0/n)) * ftol)
					|| (fabs(g0 - gprev) / n < std::max(1., fabs(g0/n)) * gtol);
				if (halt_flag) ++finish_flag; else finish_flag = 0;
				if (finish_flag > maxhalt) break;
				if (iter__ >= max_iter) break;
			}
			{
				_S msg = _S("\n        ***** _______ finish of minimization ______ *****\n");
				msg += make_string("    iteration %4u  f() = %12.5e    |f'()| = %12.5e  full nfev = %4d\n",
						++iter__, f0 /n, g0 / n, full_nfev);
				msg += make_string("    |dx| = %12.5e  |df| = %12.5e  |dg| = %12.5e",
						maxd * bx, fabs(f0 - fprev) / n, fabs(g0 - gprev) / n);
				PRINT_MESSAGE(msg);
			}
			return f0;
		}

	public:

		/**
		* @brief runs optimization
		* @param fn - object function f(n, x, g, param) for calculation dE__dX at point X
		* @param x - general coords at start pos[input] & finish pos[output]
		* @param g - general gradient at start pos[input] & finish pos[output]
		* @param param - help variable (to pass object, for example)
		* @param iteration - criteria of finish using iterations
		* @return f(x)
		*/
		template <typename T>
		T operator()(T (*fn)(unsigned, const T *, T *, void *param), void *param,
			unsigned n, T *X, T *G, unsigned max_iter=MAX_ITERATION,
			T stpmin=DEFAULT_STPMIN, T stpmax=DEFAULT_STPMAX,
			unsigned maxfev=DEFAULT_MAXFEV, unsigned maxhalt=DEFAULT_MAX_HALTS,
			T wolfe1=DEFAULT_WOLFE1, T wolfe2=DEFAULT_WOLFE2,
			T xtol=DEFAULT_XTOL<T>(), T ftol=DEFAULT_FTOL<T>(), T gtol=DEFAULT_GTOL<T>(),
			unsigned m=DEFAULT_M, int steepit=DEFAULT_STEEP_ITERATIONS
		)
		{
			if (n == 0) return 0;
			unsigned full_nfev = 0;
			int nfev = maxfev;

			T f0 = (*fn)(n, X, G, param);
			T g0 = scalar_product(n, G, G);
			{
				_S msg = _S("\n        ******************************************************\n");
				msg +=     _S("        *****   start of transition state minimization   *****\n");
				msg +=     _S("        ******************************************************\n");
				msg += make_string("    iteration    0  f() = %12.5e    |f'()| = %12.5e  nfev =  1",
						f0 / n, g0 / n);
				PRINT_MESSAGE(msg);
			}

			std::vector<T> D__(n), X__(n), G__(n);
			T *D_ = &D__[0], *X_ = &X__[0], *G_ = &G__[0];

			T bx = 0.f;
			f0 = round_minimize_(bx, f0, g0, fn, param, n, X_, G_, D_, X, G, max_iter, stpmin, stpmax,
				nfev, maxhalt, wolfe1, wolfe2, xtol, ftol, gtol);
			VECTOR_EXPRESSION_2(T, n, D_, =, , , X_, -, , , X);
				// в D лежит направление, по которому выстроился димер

			full_nfev += maxfev - nfev;

			unsigned iter__ = 0;
			for(;;)
			{
				nfev = maxfev;
				f0 = (*fn)(n, X, G, param);
				g0 = scalar_product(n, G, D_);

				// масштабируем bx таким образом, чтобы произведение bx * D[i]
				// не превысило stpmax (максимальное смещение атома по любому направлению)
				T maxd = sqrt(scalar_product(n, D_, D_) / n);
				if (maxd <= (T) 10. * std::numeric_limits<T>::min()) // нет никакого направления движения
				{
					_S msg = _S("Can't find direction for minimization");
					PRINT_MESSAGE(msg);
				}
				if (stpmax < bx * maxd) bx = stpmax / maxd;
				if (stpmin > bx * maxd) bx = stpmin / maxd;

				f0 = Line_minimizer_()(false, bx, f0, g0, fn, param, n, X_, G_, X, G, D_,
					nfev, wolfe1, wolfe2, xtol, ftol, gtol);

				full_nfev += maxfev - nfev;

				nfev = maxfev;
				f0 = n_1_minimize_(bx, f0, g0, fn, param, n, X_, G_, D_, X, G,
					max_iter, stpmin, stpmax, nfev, maxhalt, wolfe1, wolfe2, xtol, ftol, gtol);

				full_nfev += maxfev - nfev;

				iter__++;
				if (iter__ >= max_iter) break;
			}

			{
				_S msg = _S("\n        ******************************************************\n");
				msg     += _S("        *****   finish of transition state minimization  *****\n");
				msg     += _S("        ******************************************************\n");
				msg += make_string("    iteration %4u  f() = %12.5e    |f'()| = %12.5e  full nfev = %4d\n",
						++iter__, f0 /n, g0 / n, full_nfev);
				PRINT_MESSAGE(msg);
			}
			return f0;
		}
	};

	#undef MAKE_EQ
	#undef SHARE_PARAM_DEF
	#undef SHARE_PARAM_VAL

}
#endif

