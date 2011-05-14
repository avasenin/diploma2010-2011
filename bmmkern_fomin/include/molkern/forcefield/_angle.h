#ifndef _ANGLE__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H
#define _ANGLE__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H

#include "molkern/__moldefs.h"
#include "molkern/forcefield/_fparams.h"
#include "molkern/forcefield/_atom.h"
#include "molkern/forcefield/_1interactions.h"

namespace molkern
{
	using namespace prgkern;

	template <int FORCEFIELD_TYPE> struct Angle_
	{
		typedef Params<ANGLE_, FORCEFIELD_TYPE> param_type;

		index_<3> ndx; // индексы атомов
		param_type param;
	};

}
#endif
