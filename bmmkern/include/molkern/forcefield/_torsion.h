#ifndef _TORSION__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H
#define _TORSION__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H

#include "molkern/__moldefs.h"
#include "molkern/forcefield/_fparams.h"
#include "molkern/forcefield/_atom.h"
#include "molkern/forcefield/_1interactions.h"

namespace molkern
{
	using namespace prgkern;

	template <int FORCEFIELD_TYPE> struct Torsion_
	{
		typedef Params<TORSION_, FORCEFIELD_TYPE> param_type;

		index_<4> ndx; // индексы атомов
		param_type param;
	};

}
#endif
