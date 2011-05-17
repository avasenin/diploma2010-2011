#ifndef _OPTIMIZE__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H
#define _OPTIMIZE__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H

#include "boost/tuple/tuple.hpp"
#include "boost/tuple/tuple_io.hpp"

#include "molkern/__moldefs.h"

namespace molkern
{
	using namespace prgkern;

	/**
	* @brief geometry optimizer for molecule (or complex)
	*/
	template <int OPTIMIZER_TYPE>
	class Optimizer_ : protected Minimizer_<OPTIMIZER_TYPE>
	{
	public:

		template <typename LPComplex, typename Param>
		real_t run(LPComplex *complex, const Param &param)
		{
			real_t energy = 0.;

			unsigned n = complex->count(FREEDOM);
			real_t *x_ = complex->get(XPOSITION);
			real_t *g_ = complex->get(GPOSITION);

		TIME_TESTING_START("Full cycle of optimization finished...", 1)

			energy = Minimizer_<OPTIMIZER_TYPE>::operator()(
				&dU__dX<LPComplex>, (void*)complex,
				n, &x_[0], &g_[0], param.maxiter,
				param.stpmin, param.stpmax,
				param.maxfev, param.maxhalt,
				param.wolfe1, param.wolfe2,
				param.xtol, param.ftol, param.gtol,
				param.m, param.steep
			);

		TIME_TESTING_FINISH;

		return energy;
		}

		template <typename LPComplex>
		static real_t dU__dX(unsigned n, const real_t *x, real_t *g, void *param)
		{
			PRINT_MESSAGE("BING");
			typedef typename LPComplex::atom_type _Atom;
			LPComplex* complex = (LPComplex*)param;
			_Atom *atoms = complex->get(ATOM);

			complex->read(POSITION, x, atoms);
			real_t energy = complex->dU__dX();
			complex->write(GRADIENT, g, atoms);

			return energy;
		}
	};

}
#endif
