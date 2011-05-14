#ifndef _BOND__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H
#define _BOND__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H

#include "molkern/__moldefs.h"
#include "molkern/forcefield/_atom.h"
#include "molkern/forcefield/_fparams.h"
#include "molkern/forcefield/_1interactions.h"

namespace molkern
{
	using namespace prgkern;

	template <int FORCEFIELD_TYPE> struct Bond_
	{
		typedef Params<BOND_, FORCEFIELD_TYPE>  param_type;
		typedef index_<2>  index_type;

		index_<2> ndx; // индексы атомов
		param_type param; // параметры связи

		/**
		* @brief энергия связи
		* @param R вектор XB-XA
		*/
		real_t U(const vector_t &R) const
		{
			real_t q = sqrt(scalar_product(R, R));
			return param.ke * sqr(q - param.q0);
		}

		/**
		* @brief расчет силы, действующих на атомы связи
		*   Вектор, определяющий расстояние, заменяется силой.
		* @param R[in,out] вектор XB-XA; сила, действующая на атомы
		* @return энергия
		*/
		real_t dU__dX(vector_t &R) const
		{
			real_t q = sqrt(scalar_product(R, R));
			real_t q__ = q - param.q0;
			R *= mult2(param.ke) * q__ / q;
			return param.ke * sqr(q__);
		}
	};

}
#endif
