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
#ifndef _OPTIMIZER_LBFGS_0077A726_0C5C_5a4e_DCD8_D144196A0C00__H
#define _OPTIMIZER_LBFGS_0077A726_0C5C_5a4e_DCD8_D144196A0C00__H

#include "_minimize_line.h"
#include <vector>

namespace prgkern
{
	template <typename T>
	class mdense : public std::vector<T>
	{
		typedef std::vector<T> _Base;
	public:
		mdense(unsigned nn1 = 0, unsigned nn2 = 0) : _Base(nn1 * nn2), n2(nn2) {}
		mdense(_Base const& src) : _Base(((_Base&)(src)).size()), n2(src.n2) { std::copy(src.begin(), src.end(), this->begin()); }
		T &operator()(unsigned i1, unsigned i2) { return ((_Base&)(*this))[i1 * n2 + i2]; }
		void resize(unsigned nn1, unsigned nn2) {  std::vector<T>::resize(nn1 * nn2); n2 = nn2;}
		unsigned n2;
		int get_max_width() { return this-> size() / n2; }
		int get_max_height() { return n2;}
		_Base invert_diagonal() {
			_Base copy = *this;
			invert(copy);
			return copy;
		}
		_Base mul(_Base const& a, _Base const& b)
		{
			assert(a.get_max_height() == b.get_max_width());
			mdense<T> result = mdense<T>(a.get_max_width(), b.get_max_height());
			for (int i =0; i < a.get_max_width(); i++ )
			{
				for (int j=0; j < a.get_max_height(); j++)
				{
					for (int k=0; k < b.get_max_height; k++)
					{
						result(i,j) += a(i,j) * b(j,k);
					}
				}
			}
			return result;
		}
		_Base transpose()
		{
			mdense<T> result = mdense<T>(get_max_height(), get_max_width());
			for (int i =0; i < get_max_width(); i++ )
			{
				for (int j=0; j < get_max_height(); j++)
				{
					result(j,i) = (*this)(i,j);
				}
			}
		}
	};
	
	/**
	* @brief get U(-1) from upper triangle matrix
	* @note result matrix replaces the origin matrix
	*   can be used for small matrix 'cos the scaling is O(N**3)
	*/
	template <typename T>
	inline void invert(mdense<T> &U)
	{
		int n = U.n2;
		for (int i=0; i<n; ++i)
		{
			U(i, i) = 1. / U(i, i);
			for (int j=i+1; j<n; ++j)
			{
				double s = 0;
				for (int k=i; k<j; ++k)
					s += U(i, k) * U(k, j);
				U(i, j) = -s / U(j, j);
			}
		}
	}
	
}

namespace prgkern
{
	
	/**
	* @param DEFAULT_WOLFE1 Controls the accuracy of the line search for f().
	*/
	const double DEFAULT_WOLFE1 = 0.1;
	
	/**
	* @param DEFAULT_WOLFE1 Controls the accuracy of the line search for f'().
	*/
	const double DEFAULT_WOLFE2 = 0.1;
	
	/**
	* @param m The number of corrections used in the BFGS update.
	*   Values of m less than 3 are not recommended;
	*   large values of m will result in excessive computing time.
	*   3 <= m <= 7 is recommended by authors, but you can use upto 80.
	*/
	const int LBFGS_DEFAULT_M = 7;
	
	/**
	* @param maxfev Termination occurs when the number of evaluations
	*   of the objective function is at least maxfev by the end of an iteration.
	*/
	const int LBFGS_DEFAULT_MAXFEV = 40; // need large value for first iteration
	
	/**
	* @param [x,f,g]tol Controls the accuracy of the line search.
	*/
	const double LBFGS_DEFAULT_XTOL = 1.e-16;
	const double LBFGS_DEFAULT_FTOL = 1.e-16;
	const double LBFGS_DEFAULT_GTOL = 1e-6;

	/**
	* @param stp[min,max] Specifies the lower & upper bounds for the step in the line search.
	*/
	const double LBFGS_DEFAULT_STPMIN = 1.e-20;
	const double LBFGS_DEFAULT_STPMAX = 1.e20;
	
	/**
	* @brief direction updater
	* @param TYPE - updater type (STEEP_, LMBFGS_)
	*/
	template <int TYPE> class Updater;
	
	template <> class Updater<STEEP_>
	{
	public:
		
		/**
		* @brief ctor
		* @param nmax - max degrees of freedom
		* @param m - number of substeps
		*/
		Updater(unsigned nmax, unsigned m=0) : n_(nmax) {}
		
		/**
		* @brief update function
		* @param d - [out] direction of line search
		* @param f,x,g - new f(), coords & gradient
		*/
		void update(double *d, double f, const double *x, const double *g)
		{
			copy(n_, g, d); // d = -g
			scal(n_, -1., d);
		}
		
		/**
		* @brief name of minimizer
		*/
		static std::string name() { return _S("STEEP_"); }
		
	protected:
		unsigned n_;
	};
	
	/**
	* @brief compact form of modified L-BFGS
	*   Y.Yueting, X.Chengxian "A compact limited memory method
	*   for large scale unconstrained optimization", 2006,
	*   European Journal of operation reseach, in press
	* @note  save elements in order k%m to avoid reordering S, Y, D
	*   example [0, 1, 2, ..., m-1] or [m, m+1, 2, 3, ..., m-1]
	*   rewritten       --->          ^^^  ^^^
	*/
	template <> class Updater<LMBFGS_>
	{
		typedef std::vector<double>  _Vector;
		typedef mdense<double>       _Mdense;
		
	public:
		
		/**
		* @param nmax - max degrees of freedom
		* @param mmax - max sub iterations
		*/
		Updater(unsigned nmax, unsigned mmax=LBFGS_DEFAULT_M);
		
		/**
		* @brief update function
		* @param d - [out] direction of line search
		* @param f,x,g - f(), position & gradient
		*/
		void update(double *d, double f, const double *x, const double *g);
			
		/**
		* @brief name of minimizer
		*/
		static std::string name() { return _S("LMBFGS_"); }
		
	protected:
		
		unsigned n_; // the number of freedom degrees
		unsigned m_; // the subiterations limit
		int k_; // iterations number (update counter)
		
		std::vector<_Vector> Sk_; // X(k+1) - X(k)               - effective saving
		std::vector<_Vector> Yk_; // grad(k+1) - grad(k)         - effective saving
		_Mdense Rk_; // S(t)[i-1]*Y[j-1], if i<=j & 0 otherwise  - effective saving
		_Mdense Wk_; // Yk(T)*H0*Yk                              - effective saving
		_Mdense Tk_; // Rk**(-1)                - normal saving
		_Mdense Zk_; // [Sk(T)*qk, Yk(T)*H0*qk] - normal saving
		_Mdense Pk_; // MATR * Zk               - normal saving
		_Vector Dk_; // diag[S0(T)Y0, ..., Sk-1(T)Yk-1]
		_Vector H0_; // start diagonal Hessian**(-1)
		_Vector xk__; // previous x(k-1)
		_Vector gk__; // previous grad(k-1)
		double fk__; // previous f()
	};
	
	#define TEMPLATE_HEADER
	#define TEMPLATE_ARG      LMBFGS_
	
	TEMPLATE_HEADER
	Updater<TEMPLATE_ARG>
	::Updater(unsigned nmax, unsigned mmax) : n_(nmax), m_(mmax), k_(-1),
		Sk_(mmax), Yk_(mmax), Rk_(mmax, mmax), Wk_(mmax, mmax), Tk_(mmax, mmax),
		Zk_(2, mmax), Pk_(3, mmax), Dk_(mmax), H0_(nmax),
		xk__(nmax), gk__(nmax), fk__(0.)
	{
		for (int i=0; i<m_; ++i)
		{
			Sk_[i].resize(n_);
			Yk_[i].resize(n_);
		}
		set(n_, 1., &H0_[0]);
	}
	
	TEMPLATE_HEADER
	void Updater<TEMPLATE_ARG>
	::update(double *d, double f, const double *x, const double *g)
	{
	#define ADDR(m)  &(m)[0]
		if (k_ == -1)
		{
		 copy(n_, g, d);
		 scal(n_, -1., d);
		}
		else
		{
			int k__ = k_ % m_;
			
			_xpby(n_, x, -1., ADDR(xk__), ADDR(Sk_[k__]));
			_xpby(n_, g, -1., ADDR(gk__), ADDR(Yk_[k__]));
			
			// calc (g, Sk) & (g__, Sk)
			double gs   = scalar_product(n_, ADDR(Sk_[k__]), ADDR(g   )); // Sk * g
			double gs__ = scalar_product(n_, ADDR(Sk_[k__]), ADDR(gk__));
			double ys   = gs - gs__;
			double yy   = scalar_product(n_, ADDR(Yk_[k__]), ADDR(Yk_[k__])); // Yk * Yk
			
			double teta = 6 * (fk__ - f) + 3 * (gs__ + gs);
			if (teta < (M_SAFETY_HIGH - 1) * ys)
			{
			#ifdef DEGENERACY_PRINT
				{
					std::string msg = _S("[WARNING] LMBFGS H**(-1) correction teta: ") + ftoa(teta)
						+ _S(" -ys : ") + ftoa(-ys);
					PRINT_MSG(msg);
				}
			#endif
				// use safeguard strategy to garantee the positive definite of H matrix
				teta = (M_SAFETY_HIGH - 1) * ys;
			}
			
			// update Dk, H0_
			Dk_[k__] = ys / (1. + teta/ys);
			
			// set(n_, gs__/yy, ADDR(H0_));
				// this correction will produce unstable behauvour
				// so for LMBFGS sufficiently use only teta correction
			
			// update Wk & Rk
			int lmax = std::max(-1, (int)(k_ - m_));
			for (int i=k_; i>lmax; --i)
			{
				int i__ = i % m_;
				Wk_(k__, i__) = scalar_product(n_, ADDR(Yk_[k__]), ADDR(H0_), ADDR(Yk_[i__]));
				Wk_(i__, k__) = Wk_(k__, i__);
				Rk_(i__, k__) = scalar_product(n_, ADDR(Sk_[i__]), ADDR(Yk_[k__]));
			}
			Wk_(k__, k__) += Dk_[k__];
			
			// calc Tk = Rk**(-1)
			int lmin = std::min((int)(k_ + 1), (int)m_);
			for (int i=0; i<lmin; ++i)
			{
				int i__ = (k_ < m_) ? i : ((k_ + 1 + i) % m_);
				for (int j=i; j<lmin; ++j)
				{
					int j__ = (k_ < m_) ? j : ((k_ + 1 + j) % m_);
					Tk_(i, j) = Rk_(i__, j__);
				}
			}
			for (int i=lmin; i<m_; ++i) Tk_(i, i) = 1.;
			invert(Tk_);
			
			// update Zk
			for (int i=0; i<lmin; ++i)
			{
				int i__ = (k_ < m_) ? i : ((k_ + 1 + i) % m_);
				Zk_(0, i) = scalar_product(n_, ADDR(Sk_[i__]), ADDR(g));
				Zk_(1, i) = scalar_product(n_, ADDR(Yk_[i__]), ADDR(H0_), ADDR(g));
			}
			
			// update Pk
			for (int i=0; i<lmin; ++i)
			{
				double s = 0;
				for (int j=i; j<lmin; ++j) s -= Tk_(i, j) * Zk_(0, j);
				Pk_(1, i) = s;
			}
			for (int i=0; i<lmin; ++i)
			{
				int i__ = k_ < m_ ? i : (k_ + 1 + i) % m_;
				double s = 0;
				for (int j=0; j<lmin; ++j)
				{
					int j__ = (k_ < m_) ? j : ((k_ + 1 + j) % m_);
					s += Wk_(i__, j__) * Pk_(1, j);
				}
				Pk_(2, i) = s + Zk_(1, i);
			}
			for (int i=0; i<lmin; ++i)
			{
				double s = 0;
				for (int j=0; j<=i; ++j)
				{
					s -= Tk_(j, i) * Pk_(2, j);
				}
				Pk_(0, i) = s;
			}
			
			// update *d
			for (int i=0; i<n_; ++i)
			{
				double s = 0;
				for (int j=0; j<lmin; ++j)
				{
					int j__ = (k_ < m_) ? j : ((k_ + 1 + j) % m_);
					s += Sk_[j__][i] * Pk_(0, j);
					s += H0_[i] * Yk_[j__][i] * Pk_(1, j);
				}
				d[i] = -(H0_[i] * g[i] + s);
			}
		}
		
		fk__ = f;
		copy(n_, x, ADDR(xk__));
		copy(n_, g, ADDR(gk__));
		k_++;
	#undef ADDR
	}
	
	#undef TEMPLATE_HEADER
	#undef TEMPLATE_ARG
	
	/**
	* @brief Minimizer
	* @note All optimizers use the same line searcher but different updaters.
	*   Minimizer<STEEP_> uses Updater<STEEP_>.
	*   Minimizer<LMBFGS_> uses Updater<LMBFGS_> etc.
	*/
	template <int TYPE> class Minimizer_
	{
		typedef Updater<TYPE>  _Updater;
		
	public:
		
		/**
		* @brief name of minimizer
		*/
		std::string name() const { return Updater<TYPE>::name(); }
		
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
		T operator()(T (*fn)(unsigned, const T *, T *, void *param),
			unsigned n, T *x, T *g, void *param, unsigned max_iter=MAX_ITERATION,
			T wolfe1=DEFAULT_WOLFE1, T wolfe2=DEFAULT_WOLFE2,
			T xtol=LBFGS_DEFAULT_XTOL, T ftol=LBFGS_DEFAULT_FTOL, T gtol=LBFGS_DEFAULT_GTOL,
			unsigned m=LBFGS_DEFAULT_M, unsigned maxfev=LBFGS_DEFAULT_MAXFEV,
			T stpmin=LBFGS_DEFAULT_STPMIN, T stpmax=LBFGS_DEFAULT_STPMAX)
		{
			if (n == 0) return 0;
			std::vector<double> x__(n), g__(n), d(n);
			Updater<TYPE> dir_search(n, m);
			
			char line[120];
			T bx = std::max(MINIMIZER_STEP_START, 4 * xtol);
				// 4 * xtol it's heuristics to avoid problems with decreasing of interval
				
			T f = (*fn)(n, x, g, param);
			T df__dx = sqrt(scalar_product(n, g, g));
			int iter__ = 0;
			{
				::sprintf(line, "\n   ***** ____ start of minimization ____ *****");
				PRINT_MSG(line);
				::sprintf(line, "    iteration    0  f() = %22.12lf    |f'()| = %14.8lf  nfev =  1", f, df__dx);
				PRINT_MSG(line);
			}
			
			int nfev, full_nfev = 0; // number of evaluation of f()
			T f__, df__dx__; // previous values
			while (iter__ < max_iter)
			{
				if (bx > LBFGS_DEFAULT_STPMAX) bx = LBFGS_DEFAULT_STPMAX;
				if (bx < LBFGS_DEFAULT_STPMIN) bx = LBFGS_DEFAULT_STPMIN;
				
				 // save previous value
				f__ = f;
				df__dx__ = df__dx;
				nfev = maxfev;
				dir_search.update(&d[0], f, x, g);
				f = Line_minimizer_()(fn, param, n, &x__[0], &g__[0], bx, f, x, g, &d[0],
					nfev, wolfe1, wolfe2, xtol, ftol, gtol);
				copy(n, &x__[0], x);
				copy(n, &g__[0], g);
				df__dx = sqrt(scalar_product(n, g, g));
				
				full_nfev += maxfev - nfev;
				
				if (fabs(bx) < xtol) break;
				if (fabs(f - f__) < std::max(1., fabs(f)) * ftol) break;
				if (fabs(df__dx - df__dx__) < std::max(1., fabs(df__dx)) * gtol) break;

				if (++iter__ % OPTIMIZE_PRINT_FREQUENCY == 0)
				{
					::sprintf(line, "    iteration %4d  f() = %22.12lf    |f'()| = %14.8lf  nfev = %2d",
						iter__, f, df__dx, maxfev - nfev);
					PRINT_MSG(_S(line));
				}
			}
			{
				::sprintf(line, "   ***** ____ finish of minimization ____ *****  \n ");
				PRINT_MSG(line);
				::sprintf(line, "    iteration %4d  f() = %22.12lf    |f'()| = %14.8lf  full nfev = %4d",
						++iter__, f, df__dx, full_nfev);
				PRINT_MSG(line);
				::sprintf(line, "    |dx| = %12.5e with xtol = %12.5e", bx, xtol);
				PRINT_MSG(line);
				::sprintf(line, "    |df| = %12.5e with std::max(1., fabs(f)) * ftol = %12.5e",
					fabs(f - f__), std::max(1., fabs(f)) * ftol);
				PRINT_MSG(line);
				::sprintf(line, "    |dg| = %12.5e with std::max(1., fabs(g)) * gtol = %12.5e \n",
					fabs(df__dx - df__dx__), std::max(1., fabs(df__dx)) * gtol);
				PRINT_MSG(line);
			}
			return f;
		}
	};
	
}
#endif
