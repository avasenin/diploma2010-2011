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

#include <string>
#include <stdexcept>

//#define OPTIMIZATION_DEBUG
//#define DEGENERACY_PRINT

#define _S std::string
#define PRINT_ERR(msg)   { std::cout << msg << std::endl; throw std::exception(); }
#define PRINT_MSG(msg)   { std::cout << msg << std::endl; }
#define PRINT_BREAK(msg) { std::cout << msg << std::endl; throw std::exception(); }

namespace prgkern
{
	const double M_HIGH_ACCURACY = 1e-9;
	const double M_SAFETY_HIGH   = 1e-9;

	template <typename T> inline void set(unsigned n, T s, T *x)
	{ for (unsigned i_=0; i_<n; i_++) x[i_] = s; }

	template <typename T> inline void copy(unsigned n, const T *x, T *y)
	{ for (unsigned i_=0; i_<n; i_++) y[i_] = x[i_]; }

	template <typename T> inline void scal(unsigned n, T alpha, T *x)
	{ for (unsigned i_=0; i_<n; i_++) x[i_] *= alpha; }

	template <typename T>
	static void _xpby(unsigned n, const T *x, T beta, const T *y, T *z)
	{ for (unsigned i_=0; i_<n; i_++) z[i_] = x[i_] + beta * y[i_]; }
	
	template <typename T> inline T scalar_product(unsigned n, const T *x, const T *y)
	{ T s = (T)0; for (unsigned i_=0; i_<n; i_++) s += x[i_] * y[i_]; return s; }
	
	template <typename T> inline T scalar_product(unsigned n, const T *x, const T *y, const T *z)
	{ T s = (T)0; for (unsigned i_=0; i_<n; i_++) s += x[i_] * y[i_] * z[i_]; return s; }
	
	template <typename T> inline T sqr(T s) { return s * s; }
	template <typename T> inline T cube(T s) { return s * s * s; }
	
	template <typename T> bool equal(T a, T b, T accuracy=M_HIGH_ACCURACY)
	{
		T s = std::abs(a) + std::abs(b);
		if (s < accuracy) return true; // absolute comparison with 0
		if (std::abs(a - b) < accuracy * s) return true; // relative comparison
		return false;
	}
	
	template <typename T> inline int sign(T a) { return a >= 0 ? (a > 0 ? 1 : 0) : -1; }
	
	/** @brief zero of function
	* @note build linear approximation f() ~ B * x + C
	* @param xa - left bound of interval
	* @param xb - right bound of interval
	* @param fa - f(xa)
	* @param fb - f(xb)
	* @return position of zero
	*/
	template <typename T>
	inline T zero_linear(T xa, T xb, T fa, T fb)
	{
		return (xa * fb - xb * fa) /  (fb - fa);
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
	inline int extreme_sqr(T *x, T *f, T xa, T xb, T fa, T fb, T ga, T gb,
		T accuracy=M_HIGH_ACCURACY)
	{
		T h = xb - xa;
		T C = h * ga;
		T B = fb - fa + C; // B -> B / h
		if (equal(B, 0., accuracy))
		{
		#ifdef DEGENERACY_PRINT
			{
				std::string msg = _S("[WARNING] extreme degeneracy");
				PRINT_MSG(msg);
			}
		#endif
			return 0;
		}
		*x = -ga / (B + B);
		*f = polynome(B, C, fa, *x);
		*x = h * (*x) + xa; // restore x -> x * h + xa
		return 1;
	}
	
	template <typename T> inline
	T polynome(T c0, T x) { return c0; }

	template <typename T> inline
	T polynome(T c1, T c0, T x) { return c1 * x + polynome(c0, x); }

	template <typename T> inline
	T polynome(T c2, T c1, T c0, T x) { return polynome(c2, c1, x) * x + polynome(c0, x); }
}

namespace prgkern
{
	
	/// after OPTIMIZE_PRINT_FREQUENCY will be printed all information
	const int OPTIMIZE_PRINT_FREQUENCY = 1;
	const int MAX_ITERATION = 2000;
	
	/// optimizer step for first line search
	const double MINIMIZER_STEP_START = 1e-6;
	
	template <int TYPE> class Minimizer_;
	enum { SCITBX_LBFGS_, GSL_SIMPLEX_, STEEP_, LMBFGS_ };
	
	/// line search optimizer
	class Line_minimizer_
	{
		
		enum { FAIL, INTERVAL_OK, FUNCTION_DIFFERENCE_OK, WOLFE_OK, TOLERANCE_OK };
		
	public:
		
		/**
		* @brief runs optimization
		* @param x - general coords at start[in] & at finish[out]
		* @param d - direction
		* @param g - gradient at start[in] & at finish[out]
		* @return f(x)
		*/
		template <typename T>
		T operator()(T (*fn)(unsigned, const T *, T *, void *param), void *param,
			unsigned n, T *x_, T *g_, T &beta, T f, const T *x, const T *g, const T *d,
			int &maxfev, T wolfe1, T wolfe2,
			T xtol=M_HIGH_ACCURACY, T ftol=M_HIGH_ACCURACY, T gtol=M_HIGH_ACCURACY);
		
	protected:
		/**
		* @brief parameters description
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
		bool test_Wolfe_(T ax, T bx, T fa, T fb, T ga, T gb, T wolfe1, T wolfe2) const
		{
			if (fb - fa < wolfe1 * (bx - ax) * ga && fabs(gb) < wolfe2 * fabs(ga))
			{
			#ifdef OPTIMIZATION_DEBUG
				PRINT_MSG(_S("    ************** test Wolfe is OK ************** "));
			#endif
				return true;
			}
			return false;
		}
		
		template <typename T>
		bool test_tolerance_(T fa, T fb, T ftol=M_HIGH_ACCURACY) const
		{
			if (fabs(fb - fa) <= std::max(1., fabs(fa)) * ftol)
			{
	#ifdef OPTIMIZATION_DEBUG
		PRINT_MSG(_S("    ************** test tolerance is OK ************** "));
	#endif
				return true;
			}
			return false;
		}
		
		template <typename T>
		unsigned locate_zero_derivation_(T (*fn)(unsigned, const T *, T *, void *param), void *param,
			unsigned n, T *x_, T *g_, T &xa, T &xb, T &fa, T &fb, T &ga, T &gb, const T *x, const T *d,
			int &nfev, T wolfe1, T wolfe2,
			T xtol=M_HIGH_ACCURACY, T ftol=M_HIGH_ACCURACY, T gtol=M_HIGH_ACCURACY);
			
		template <typename T>
		unsigned dec_interval_(T (*fn)(unsigned, const T *, T *, void *param), void *param,
			unsigned n, T *x_, T *g_, T &xa, T &xb, T &fa, T &fb, T &ga, T &gb, const T *x, const T *d,
			int &nfev, T wolfe1, T wolfe2,
			T xtol=M_HIGH_ACCURACY, T ftol=M_HIGH_ACCURACY, T gtol=M_HIGH_ACCURACY);
		
		template <typename T>
		unsigned inc_interval_(T (*fn)(unsigned, const T *, T *, void *param), void *param,
			unsigned n, T *x_, T *g_, T &xa, T &xb, T &fa, T &fb, T &ga, T &gb, const T *x, const T *d,
			int &nfev, T wolfe1, T wolfe2,
			T xtol=M_HIGH_ACCURACY, T ftol=M_HIGH_ACCURACY, T gtol=M_HIGH_ACCURACY);
		
		template <typename T>
		unsigned undefinite_interval_(T (*fn)(unsigned, const T *, T *, void *param), void *param,
			unsigned n, T *x_, T *g_, T &xa, T &xb, T &fa, T &fb, T &ga, T &gb, const T *x, const T *d,
			int &nfev, T wolfe1, T wolfe2,
			T xtol=M_HIGH_ACCURACY, T ftol=M_HIGH_ACCURACY, T gtol=M_HIGH_ACCURACY);
	};
	
	template <typename T>
	inline T Line_minimizer_
	::operator()(T (*fn)(unsigned, const T *, T *, void *param), void *param,
		unsigned n, T *x_, T *g_, T &beta, T fa, const T *x, const T *g, const T *d,
		int &nfev, T wolfe1, T wolfe2, T xtol, T ftol, T gtol)
	{
		if (n == 0) return fa;
		T ga = scalar_product(n, g,  d);
		if (ga >= 0) { PRINT_ERR(_S("[ERROR] Bad direction in line minimizer")); }
		
		T xa = 0; T xb = beta;
		_xpby(n, x, xb, d, x_);
		T fb = (*fn)(n, x_, g_, param); nfev--;
		T gb = scalar_product(n, g_, d);
		if (test_Wolfe_(xa, xb, fa, fb, ga, gb, wolfe1, wolfe2)) { return fb; }
		
		unsigned found = FAIL;
		if (gb <= 0)
		{
			if (fabs(fb - fa) <= M_SAFETY_HIGH)
			{
				found = undefinite_interval_(fn, param, n, x_, g_, xa, xb, fa, fb, ga, gb, x, d,
					nfev, wolfe1, wolfe2, xtol, ftol, gtol);
				if (found == WOLFE_OK || found == TOLERANCE_OK) goto LABEL_FINISH;
				if (found == INTERVAL_OK) goto LABEL_LOCATE;
			}
			if (fb > fa)
			{
				found = dec_interval_(fn, param, n, x_, g_, xa, xb, fa, fb, ga, gb, x, d,
					nfev, wolfe1, wolfe2, xtol, ftol, gtol);
				if (found == WOLFE_OK || found == TOLERANCE_OK) goto LABEL_FINISH;
			}
			if (fb < fa)
			{
				found = inc_interval_(fn, param, n, x_, g_, xa, xb, fa, fb, ga, gb, x, d,
					nfev, wolfe1, wolfe2, xtol, ftol, gtol);
				if (found == WOLFE_OK || found == TOLERANCE_OK) goto LABEL_FINISH;
			}
		}
		
	LABEL_LOCATE:
		if (gb > 0)
			locate_zero_derivation_(fn, param, n, x_, g_, xa, xb, fa, fb, ga, gb, x, d,
				nfev, wolfe1, wolfe2, xtol, ftol, gtol);
		
	LABEL_FINISH:
		_xpby(n, x, xb, d, x_);
		fb = (*fn)(n, x_, g_, param); nfev--;
		
		beta = xb;
		return fb;
	}
	
	
	template <typename T>
	inline unsigned Line_minimizer_
	::locate_zero_derivation_(T (*fn)(unsigned, const T *, T *, void *param), void *param,
		unsigned n, T *x_, T *g_, T &xa, T &xb, T &fa, T &fb, T &ga, T &gb, const T *x, const T *d,
		int &nfev, T wolfe1, T wolfe2, T xtol, T ftol, T gtol)
	{
	#ifdef OPTIMIZATION_DEBUG
		PRINT_MSG("locate_zero_derivation_:");
		{
			char line[140];
			::sprintf(line, "[%17.10e %17.10e]  f[%17.10e %17.10e]  g[%17.10e %17.10e]",
				xa, xb, fa, fb, ga, gb);
			PRINT_MSG(_S(line));
		}
	#endif
		bool good = (xa < xb) && (ga <= 0) && (gb >= 0);
		if (!good)
		{
			PRINT_MSG(" test of condition: (xa < xb) && (ga <= 0) && (gb >= 0)");
			char line[140];
			::sprintf(line, "[%17.10e %17.10e]  f[%17.10e %17.10e]  g[%17.10e %17.10e]",
				xa, xb, fa, fb, ga, gb);
			PRINT_BREAK(_S("[ERROR] locate_zero_derivation_ bad start values"));
		}
		const double SHIFT = 0.1;
		
		T f0 = fa, x0 = xa, g0 = ga;
		T xu, fu, gu, ca = xb - xa; T sign_edge; T shift;
		while (ca > xtol)
		{
			xu = zero_linear(xa, xb, ga, gb);
			// make shift of xu to catch the large edge
			if (xb - xu > xu - xa) { shift = SHIFT * (xb - xu); sign_edge = sign(gb); }
			else { shift = SHIFT * (xa - xu); sign_edge = sign(ga); }
			xu += shift;
			
			_xpby(n, x, xu, d, x_);
			fu = (*fn)(n, x_, g_, param); nfev--;
			gu = scalar_product(n, g_, d);
			// try to cut the large edge
			if (gu < 0.) { xa = xu; fa = fu; ga = gu;}
			else { xb = xu; fb = fu; gb = gu; }
			
		#ifdef OPTIMIZATION_DEBUG
			{
				char line[140];
				::sprintf(line, "[%17.10e %17.10e]  f[%17.10e %17.10e]  g[%17.10e %17.10e]",
					xa, xb, fa, fb, ga, gb);
				PRINT_MSG(_S(line));
			}
		#endif
			
			// fail to cut large edge so cut the another edge
			if (sign(gu) != sign_edge)
			{
				xu = 0.5 * (xa + xb);
				_xpby(n, x, xu, d, x_);
				fu = (*fn)(n, x_, g_, param); nfev--;
				gu = scalar_product(n, g_, d);
				if (gu < 0.) { xa = xu; fa = fu; ga = gu;}
				else { xb = xu; fb = fu; gb = gu; }
			#ifdef OPTIMIZATION_DEBUG
				{
					char line[140];
					::sprintf(line, "[%17.10e %17.10e]  f[%17.10e %17.10e]  g[%17.10e %17.10e]",
						xa, xb, fa, fb, ga, gb);
					PRINT_MSG(_S(line));
				}
			#endif
			}
			
			ca = xb - xa;
			if (test_Wolfe_(x0, xb, f0, fb, g0, gb, wolfe1, wolfe2)) return WOLFE_OK;
			if (test_tolerance_(ga, gb, gtol)) return TOLERANCE_OK;
			if (nfev < 0) { PRINT_BREAK(_S("[ERROR] exceed the limit of f() evaluations")); }
		}
		xb = xa + 0.5 * ca;
		return TOLERANCE_OK;
	}
	
	template <typename T>
	inline unsigned Line_minimizer_
	::inc_interval_(T (*fn)(unsigned, const T *, T *, void *param), void *param,
		unsigned n, T *x_, T *g_, T &xa, T &xb, T &fa, T &fb, T &ga, T &gb, const T *x, const T *d,
		int &nfev, T wolfe1, T wolfe2, T xtol, T ftol, T gtol)
	{
	#ifdef OPTIMIZATION_DEBUG
		PRINT_MSG("inc_interval_:");
		{
			char line[140];
			::sprintf(line, "[%17.10e %17.10e]  f[%17.10e %17.10e]  g[%17.10e %17.10e]",
				xa, xb, fa, fb, ga, gb);
			PRINT_MSG(_S(line));
		}
	#endif
		bool good = (xa < xb) && (ga < 0) && (gb < 0);
		if (!good)
		{
			PRINT_MSG(" test of condition: (xa < xb) && (ga < 0) && (gb < 0)");
			char line[140];
			::sprintf(line, "[%17.10e %17.10e]  f[%17.10e %17.10e]  g[%17.10e %17.10e]",
				xa, xb, fa, fb, ga, gb);
			PRINT_BREAK(_S("[ERROR] inc_interval_ bad start values"));
		}
		
		T xu, fu, gu, hmax, xc = xb - xa;
		xtol = 1. / xtol; // reverse for increasing of intreval
		while (xc < xtol)
		{
			hmax = 8 * (xb - xa); // upper limit of interval
			unsigned num_roots = extreme_sqr(&xu, &fu, xa, xb, fa, fb, ga, gb);
			
			if (num_roots == 0 || (num_roots == 1 && (xu < xa || xu > xa + hmax)))
			{
			#ifdef DEGENERACY_PRINT
				if (num_roots == 0) PRINT_MSG(_S("[WARNING] inc_interval_ degeneracy"));
			#endif
				xu = xa + hmax;
			}
			
			_xpby(n, x, xu, d, x_);
			fu = (*fn)(n, x_, g_, param); nfev--;
			gu = scalar_product(n, g_, d);
		#ifdef OPTIMIZATION_DEBUG
			{
				char line[140];
				::sprintf(line, "[%17.10e %17.10e]  f[%17.10e %17.10e]  g[%17.10e %17.10e]",
					xa, xu, fa, fu, ga, gu);
				PRINT_MSG(_S(line));
			}
		#endif
			
			xa = xb; fa = fb; ga = gb;
			xb = xu; fb = fu; gb = gu;
			
			if (test_Wolfe_(xa, xb, fa, fb, ga, gb, wolfe1, wolfe2)) return WOLFE_OK;
			if (fu > fb || gu > 0.) return INTERVAL_OK;
			if (nfev < 0) { PRINT_BREAK(_S("[ERROR] exceed the limit of f() evaluations")); }
			xc = xb - xa;
		}
		return FAIL;
	}
	
	template <typename T>
	inline unsigned Line_minimizer_
	::dec_interval_(T (*fn)(unsigned, const T *, T *, void *param), void *param,
		unsigned n, T *x_, T *g_, T &xa, T &xb, T &fa, T &fb, T &ga, T &gb, const T *x, const T *d,
		int &nfev, T wolfe1, T wolfe2, T xtol, T ftol, T gtol)
	{
	#ifdef OPTIMIZATION_DEBUG
		PRINT_MSG("dec_interval_:");
		{
			char line[140];
			::sprintf(line, "[%17.10e %17.10e]  f[%17.10e %17.10e]  g[%17.10e %17.10e]",
				xa, xb, fa, fb, ga, gb);
			PRINT_MSG(_S(line));
		}
	#endif
		bool good = (xa < xb) && (ga < 0) && (gb < 0) && (fa < fb);
		if (!good)
		{
			PRINT_MSG(" test of condition: (xa < xb) && (ga < 0) && (gb < 0) && (fa < fb)");
			char line[140];
			::sprintf(line, "[%17.10e %17.10e]  f[%17.10e %17.10e]  g[%17.10e %17.10e]",
				xa, xb, fa, fb, ga, gb);
			PRINT_BREAK(_S("[ERROR] inc_interval_ bad start values"));
		}
		
		T xu, fu, gu; T ca = xb - xa;
		
		// save previous value for fail search
		T fa0 = fa, xa0 = xa, ga0 = ga;
		T fb0 = fb, xb0 = xb, gb0 = gb;
		
		while (ca > xtol)
		{
			xu = (xb + xa) * 0.5;
			_xpby(n, x, xu, d, x_);
			fu = (*fn)(n, x_, g_, param); nfev--;
			gu = scalar_product(n, g_, d);
		#ifdef OPTIMIZATION_DEBUG
			{
				char line[140];
				::sprintf(line, "[%17.10e %17.10e]  f[%17.10e %17.10e]  g[%17.10e %17.10e]",
					xa, xb, fa, fb, ga, gb);
				PRINT_MSG(_S(line));
			}
		#endif
			
			if (gu < 0.)
			{
				if (fu < fa) { xa = xu; fa = fu; ga = gu; }
				else { xb = xu; fb = fu; gb = gu; }
			}
			else
			{
				xb = xu; fb = fu; gb = gu;
				return INTERVAL_OK;
			}
			if (test_Wolfe_(xa0, xb, fa0, fb, ga0, gb, wolfe1, wolfe2)) return WOLFE_OK;
			if (test_tolerance_(ga0, gb, gtol)) return TOLERANCE_OK;
			if (nfev < 0) { PRINT_BREAK(_S("[ERROR] exceed the limit of f() evaluations")); }
			ca = xb - xa;
		}
		// restore previous values in fail case
		fa = fa0; fb = fb0;
		xa = xa0; xb = xb0;
		ga = ga0; gb = gb0;
		return FAIL;
	}
	
	template <typename T>
	inline unsigned Line_minimizer_
	::undefinite_interval_(T (*fn)(unsigned, const T *, T *, void *param), void *param,
		unsigned n, T *x_, T *g_, T &xa, T &xb, T &fa, T &fb, T &ga, T &gb, const T *x, const T *d,
		int &nfev, T wolfe1, T wolfe2, T xtol, T ftol, T gtol)
	{
	#ifdef OPTIMIZATION_DEBUG
		PRINT_MSG("undefinite_interval_:");
		{
			char line[140];
			::sprintf(line, "[%17.10e %17.10e]  f[%17.10e %17.10e]  g[%17.10e %17.10e]",
				xa, xb, fa, fb, ga, gb);
			PRINT_MSG(_S(line));
		}
	#endif
		bool good = (xa < xb) && (ga < 0) && (gb < 0);
		if (!good)
		{
			PRINT_MSG(" test of condition: (xa < xb) && (ga < 0) && (gb < 0)");
			char line[140];
			::sprintf(line, "[%17.10e %17.10e]  f[%17.10e %17.10e]  g[%17.10e %17.10e]",
				xa, xb, fa, fb, ga, gb);
			PRINT_BREAK(_S("[ERROR] undefinite_interval_ bad start values"));
		}
		
		T xu = 0.5 * (xb + xa);
		T xu__ = xa + 2 * (xb - xa);
		T xtol__ = 1. / xtol;
		
		unsigned count = 0;
		//while (xu >= xtol && xu__ <= xtol__)
		while (fabs(xu__ - xa) <= xtol__)
		{
			if (count++ % 2) { xb = xu; xu = 0.5 * (xb + xa); }
			else { xb = xu__; xu__ = xa + 2 * (xb - xa); }
			
			_xpby(n, x, xb, d, x_);
			fb = (*fn)(n, x_, g_, param); nfev--;
			gb = scalar_product(n, g_, d);
		#ifdef OPTIMIZATION_DEBUG
			{
				char line[140];
				::sprintf(line, "[%17.10e %17.10e]  f[%17.10e %17.10e]  g[%17.10e %17.10e]",
					xa, xb, fa, fb, ga, gb);
				PRINT_MSG(_S(line));
			}
		#endif
			if (test_Wolfe_(xa, xb, fa, fb, ga, gb, wolfe1, wolfe2)) return WOLFE_OK;
			if (sign(ga) != sign(gb)) return INTERVAL_OK;
			if (fabs(fb - fa) > M_SAFETY_HIGH) return FUNCTION_DIFFERENCE_OK;
			if (fabs(gb - ga) < std::max(1., fabs(ga)) * gtol) return TOLERANCE_OK;
			if (fabs(xb - xa) < xtol) return TOLERANCE_OK;
			if (nfev < 0) { PRINT_BREAK(_S("[ERROR] exceed the limit of f() evaluations")); }
		}
		PRINT_BREAK(_S("[ERROR] can't find interval"));
		return FAIL;
	}
	
}
#endif
