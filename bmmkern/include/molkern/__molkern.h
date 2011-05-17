#ifndef MOLKERN__0077A726_C744_52f9_4B92_B244A9C70800__H
#define MOLKERN__0077A726_C744_52f9_4B92_B244A9C70800__H

#include "molkern/__moldefs.h"
	// должна быть включена первой, так как определяет директивы препроцессора

#include <algorithm>
#include <cassert>
#include <exception>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <set>
#include <stack>
#include <string>
#include <deque>
#include <vector>
#include <iterator>
#include <ctime>

#include "molkern/__config.h"
#include "molkern/forcefield/_rapple_goddard_params.h"
#include "molkern/forcefield/_fparams.h"
#include "molkern/forcefield/_forcefield.h"
#include "molkern/forcefield/_forcefield_amber.h"
#include "molkern/forcefield/_residue.h"
#include "molkern/forcefield/_residue_amber.h"
#include "molkern/forcefield/_residome.h"
#include "molkern/forcefield/_residome_amber.h"
#include "molkern/forcefield/_nuclear.h"
#include "molkern/forcefield/_atomdata.h"
#include "molkern/forcefield/_atom.h"
#include "molkern/forcefield/_bond.h"
#include "molkern/forcefield/_angle.h"
#include "molkern/forcefield/_torsion.h"
#include "molkern/forcefield/_1interactions.h"
#include "molkern/forcefield/_interactions.h"

#include "molkern/complex/_parallel.h"
#include "molkern/complex/_region.h"
#include "molkern/complex/_verlet.h"
#include "molkern/complex/_linkcell.h"
#include "molkern/complex/_protonization.h"
#include "molkern/complex/_water.h"
#include "molkern/complex/_molecule.h"
#include "molkern/complex/_complex.h"
#include "molkern/complex/_optimize.h"
#include "molkern/complex/_charge_optimize.h"
#include "molkern/complex/_geom_tool.h"
#include "molkern/complex/_phys_tool.h"
#include "molkern/complex/_ensemble.h"
#include "molkern/complex/_thermostat.h"
#include "molkern/complex/_mdintegrator.h"
#include "molkern/complex/_coulomb_params.h"
//---------------------------------------------------------------------------------------------------
//             снятие всех определений после использования в библиотеке MOLKERN
//---------------------------------------------------------------------------------------------------
#undef DEFINE_VECTOR_ACCESS_FUNCTION
#undef DEFINE_TAG

#endif

