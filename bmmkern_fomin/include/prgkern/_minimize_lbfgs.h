#ifndef _OPTIMIZER_LBFGS_0077A726_0C5C_5a4e_DCD8_D144196A0C00__H
#define _OPTIMIZER_LBFGS_0077A726_0C5C_5a4e_DCD8_D144196A0C00__H

#include "prgkern/_prgconfig.h"
#include "prgkern/_math.h"
#include "prgkern/_minimize.h"
#include <vector>
#include <limits>

namespace prgkern
{
	/**
	* @brief direction updater
	* @param TYPE - updater type (STEEP_, LMBFGS_)
	*/
	template <int TYPE, typename _Real> class Updater;

	template <typename _Real> class Updater<STEEP_, _Real>
	{
	public:

		/**
		* @brief ctor
		* @param nmax - max degrees of freedom
		* @param m - number of substeps
		*/
		Updater(unsigned nmax, unsigned m=0, int=0) : n_(nmax) {}

		/**
		* @brief update function
		* @param d - [out] direction of line search
		* @param f,x,g - new f(), coords & gradient
		*/
		void update(_Real *d, _Real f, const _Real *x, const _Real *g)
		{
			VECTOR_EXPRESSION_1(_Real, n_, d, =, 0., -, g);
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
	template <typename _Real> class Updater<LMBFGS_, _Real>
	{
		typedef std::vector<_Real>                      _Vector;
		typedef mdense_<UNLIMITED_, UNLIMITED_, _Real>  mdense;

	public:

		/**
		* @param nmax - max degrees of freedom
		* @param mmax - max sub iterations
		*/
		Updater(unsigned nmax, unsigned mmax, int steep_iterations=DEFAULT_STEEP_ITERATIONS)
		: n_(nmax), m_(mmax), k_(-steep_iterations),
			Sk_(mmax), Yk_(mmax), Rk_(mmax, mmax), Wk_(mmax, mmax), Tk_(mmax, mmax),
			Zk_(2, mmax), Pk_(3, mmax), Dk_(mmax), H0_(nmax),
			xk__(nmax), gk__(nmax), fk__(0.)
		{
			for (unsigned i=0; i<m_; ++i)
			{
				Sk_[i].resize(n_);
				Yk_[i].resize(n_);
			}
			_Real *x = &H0_[0];
			VECTOR_EXPRESSION_0(_Real, n_, x, =, (_Real)1.);
		}

		/**
		* @brief update function
		* @param d - [out] direction of line search
		* @param f,x,g - f(), position & gradient
		*/
		void update(_Real *d, _Real f, const _Real *x, const _Real *g);

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
		mdense Rk_; // S(t)[i-1]*Y[j-1], if i<=j & 0 otherwise  - effective saving
		mdense Wk_; // Yk(T)*H0*Yk                              - effective saving
		mdense Tk_; // Rk**(-1)                - normal saving
		mdense Zk_; // [Sk(T)*qk, Yk(T)*H0*qk] - normal saving
		mdense Pk_; // MATR * Zk               - normal saving
		_Vector Dk_; // diag[S0(T)Y0, ..., Sk-1(T)Yk-1]
		_Vector H0_; // start diagonal Hessian**(-1)
		_Vector xk__; // previous x(k-1)
		_Vector gk__; // previous grad(k-1)
		_Real fk__; // previous f()
	};

	#define TEMPLATE_HEADER   template <typename _Real>
	#define TEMPLATE_ARG      LMBFGS_, _Real

	TEMPLATE_HEADER
	inline void Updater<TEMPLATE_ARG>
	::update(_Real *d, _Real f, const _Real *x, const _Real *g)
	{
	#define ADDR(m)  &(m)[0]
		if (k_ < 0)
		{
			VECTOR_EXPRESSION_1(_Real, n_, d, =, 0., -, g);
		}
		else
		{
			unsigned k__ = k_ % m_;

			_Real *xkk = &xk__[0], *gkk = &gk__[0], *skk = &Sk_[k__][0], *ykk = &Yk_[k__][0];
			VECTOR_EXPRESSION_2(_Real, n_, skk, =, , , x, -, , , xkk);
			VECTOR_EXPRESSION_2(_Real, n_, ykk, =, , , g, -, , , gkk);

			// calc (g, Sk) & (g__, Sk)
			_Real gs   = scalar_product(n_, ADDR(Sk_[k__]), ADDR(g   )); // Sk * g
			_Real gs__ = scalar_product(n_, ADDR(Sk_[k__]), ADDR(gk__));
			_Real ys   = gs - gs__;

			_Real teta = 6 * (fk__ - f) + 3 * (gs__ + gs);

			const _Real high_accuracy = (_Real) 10. * std::numeric_limits<_Real>::epsilon();
			if (teta < (high_accuracy - 1) * ys)
			{
			#ifdef DEGENERACY_PRINT
				{
					std::string msg = _S("[WARNING] LMBFGS H**(-1) correction teta: ") + ftoa(teta)
						+ _S(" -ys : ") + ftoa(-ys);
					PRINT_MESSAGE(msg);
				}
			#endif
				// use safeguard strategy to garantee the positive definite of H matrix
				teta = (high_accuracy - 1) * ys;
					// низкая точность отклонения использована, чтобы избежать ситуации,
					// M_SAFETY_LOW - 1 = -1 при float расчетах
			}

			// update Dk, H0_
			Dk_[k__] = sqr(ys) / (ys + teta);

			 //_Real yy   = scalar_product(n_, ADDR(Yk_[k__]), ADDR(Yk_[k__])); // Yk * Yk
			 //set(n_, gs__/yy, ADDR(H0_));
				// this correction will produce unstable behauvour
				// so for LMBFGS sufficiently use only teta correction

			// update Wk & Rk
			int lmax = -1 > (int)(k_ - m_) ? -1 : (int)(k_ - m_);
			for (int i=k_; i>lmax; --i)
			{
				int i__ = i % m_;
				Wk_(k__, i__) = scalar_product(n_, ADDR(Yk_[k__]), ADDR(H0_), ADDR(Yk_[i__]));
				Wk_(i__, k__) = Wk_(k__, i__);
				Rk_(i__, k__) = scalar_product(n_, ADDR(Sk_[i__]), ADDR(Yk_[k__]));
			}
			Wk_(k__, k__) += Dk_[k__];

			// calc Tk = Rk**(-1)
			unsigned lmin = k_ + 1;
			lmin = std::min(lmin, m_);

			for (unsigned i=0; i<lmin; ++i)
			{
				unsigned i__ = ((unsigned)k_ < m_) ? i : ((unsigned)(k_ + 1 + i) % m_);
				for (unsigned j=i; j<lmin; ++j)
				{
					unsigned j__ = ((unsigned)k_ < m_) ? j : ((unsigned)(k_ + 1 + j) % m_);
					Tk_(i, j) = Rk_(i__, j__);
				}
			}
			for (unsigned i=lmin; i<m_; ++i) Tk_(i, i) = 1.;
			invert(Tk_);

			// update Zk
			for (unsigned i=0; i<lmin; ++i)
			{
				unsigned i__ = ((unsigned)k_ < m_) ? i : ((unsigned)(k_ + 1 + i) % m_);
				Zk_(0, i) = scalar_product(n_, ADDR(Sk_[i__]), ADDR(g));
				Zk_(1, i) = scalar_product(n_, ADDR(Yk_[i__]), ADDR(H0_), ADDR(g));
			}

			// update Pk
			for (unsigned i=0; i<lmin; ++i)
			{
				_Real s = 0;
				for (unsigned j=i; j<lmin; ++j) s -= Tk_(i, j) * Zk_(0, j);
				Pk_(1, i) = s;
			}
			for (unsigned i=0; i<lmin; ++i)
			{
				unsigned i__ = (unsigned)k_ < m_ ? i : (unsigned)(k_ + 1 + i) % m_;
				_Real s = 0;
				for (unsigned j=0; j<lmin; ++j)
				{
					unsigned j__ = ((unsigned)k_ < m_) ? j : ((unsigned)(k_ + 1 + j) % m_);
					s += Wk_(i__, j__) * Pk_(1, j);
				}
				Pk_(2, i) = s + Zk_(1, i);
			}
			for (unsigned i=0; i<lmin; ++i)
			{
				_Real s = 0;
				for (unsigned j=0; j<=i; ++j)
				{
					s -= Tk_(j, i) * Pk_(2, j);
				}
				Pk_(0, i) = s;
			}

			// update *d
			for (unsigned i=0; i<n_; ++i)
			{
				_Real s = 0;
				for (unsigned j=0; j<lmin; ++j)
				{
					unsigned j__ = ((unsigned)k_ < m_) ? j : ((unsigned)(k_ + 1 + j) % m_);
					s += Sk_[j__][i] * Pk_(0, j);
					s += H0_[i] * Yk_[j__][i] * Pk_(1, j);
				}
				d[i] = -(H0_[i] * g[i] + s);
			}

			// [fomin] контроль на правильность направления
			if (scalar_product(n_, d, g) > 0.)
			{
				VECTOR_EXPRESSION_1(_Real, n_, d, =, 0., -, g);

				_Real *x = &H0_[0];
				VECTOR_EXPRESSION_0(_Real, n_, x, =, (_Real)1.);
				k_ = -1;
				_S msg = _S("[WARNING] minimizer has renewed update");
				PRINT_MESSAGE(msg);
			}
			normalize(n_, d, (_Real)1.);
		}

		fk__ = f;
		_Real *xkk = &xk__[0], *gkk = &gk__[0];
		VECTOR_EXPRESSION_1(_Real, n_, xkk, =, , , x);
		VECTOR_EXPRESSION_1(_Real, n_, gkk, =, , , g);
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
			T xtol=DEFAULT_XTOL<T>(), T ftol=DEFAULT_FTOL<T>(),
			T gtol=DEFAULT_GTOL<T>(), unsigned m=DEFAULT_M,
			int steepit=DEFAULT_STEEP_ITERATIONS
		)
		{
			if (n == 0) return 0;

			T f0 = (*fn)(n, X, G, param);
			T g0 = scalar_product(n, G, G);

			std::vector<T> X__(n), G__(n), D__(n);
			T *X_ = &X__[0], *G_ = &G__[0], *D = &D__[0];

			Updater<TYPE, T> dir_search(n, m, steepit);

			T maxd = 0, bx = std::max((T)MINIMIZER_STEP_START, (T)(4 * (xtol)));
				// 4 * xtol it's heuristics to avoid problems with decreasing of interval

			// Обычно при старте есть клеши, которые продуцируют резкий сдвиг для некоторых атомов.
			// Чтобы избавиться от него, самую большую компоненту выровняем на 1. При этом теряем
			// значение градиента в стартовой точке.
			T gmax = 1.f / maxabs(n, G);
			VECTOR_EXPRESSION_1(T, n, G, =, gmax, *, G);

			unsigned iter__ = 0;
			{
				_S msg = _S("\n        ***** _______ start of minimization _______ *****\n");
				msg += make_string("    iteration    0  f() = %12.5e    |f'()| = %12.5e  nfev =  1",
						f0 / n, g0 / n);
				PRINT_MESSAGE(msg);
			}

			int nfev, full_nfev = 0; // number of evaluation of f()
			T fprev = f0, gprev = g0;

			unsigned finish_flag = 0; // определяет завершение работы, если условия
				// завершения случаются дважды (трижды) на двух последующих итерациях

			for(;;)
			{
				if (iter__ >= max_iter) break;

				fprev = f0; gprev = g0;
				nfev = maxfev;

				dir_search.update(D, f0, X, G);
				g0 = scalar_product(n, G, D);

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

				f0 = Line_minimizer_()(true, bx, f0, g0, fn, param, n, X_, G_, X, G, D,
					nfev, wolfe1, wolfe2, xtol, ftol, gtol);

				if (++iter__ % OPTIMIZE_PRINT_FREQUENCY == 0)
				{
					_S msg = make_string("    iteration %4u  f() = %12.5e    |f'()| = %12.5e  nfev = %2d",
						iter__, f0 / n, g0 / n, (int)(maxfev - nfev));
					PRINT_MESSAGE(msg);
				}
				full_nfev += maxfev - nfev;

				VECTOR_EXPRESSION_1(T, n, X, =, , , X_);
				VECTOR_EXPRESSION_1(T, n, G, =, , , G_);

				if (fabs(bx) <= (T) 10. * std::numeric_limits<T>::min())
				{
					_S msg = _S("Can't move from given point");
					PRINT_MESSAGE(msg);
					break;
				}

				// рассчитаем условия завершения и если они повторятся два раза подряд
				// завершим оптимизацию
				bool halt_flag = (fabs(maxd * bx) < stpmin)
					|| (fabs(f0 - fprev) / n < std::max(1., fabs(f0/n)) * ftol)
					|| (fabs(g0 - gprev) / n < std::max(1., fabs(g0/n)) * gtol);
				if (halt_flag) ++finish_flag; else finish_flag = 0;
				if (finish_flag > maxhalt) break;
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
	};

}
#endif
