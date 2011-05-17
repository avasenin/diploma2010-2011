#ifndef _NUCLEAR__0077A726_CBE6_5c54_73A2_F9434C9C0000__H
#define _NUCLEAR__0077A726_CBE6_5c54_73A2_F9434C9C0000__H

#include "molkern/__moldefs.h"

namespace molkern
{
	using namespace prgkern;

	const struct Nuclear
	{
		const char *full_name;
		int charge;
		fstring name;
		real_t mass;
	}
	nuclears[] = {
		{"unknown",    0, "  ", IntMass(0.0000)   },
		{"hydrogen",   1, "H ", IntMass(1.0079)   },
		{"helium",     2, "He", IntMass(4.00260)  },
		{"lithium",    3, "Li", IntMass(6.94)     },
		{"beryllium",  4, "Be", IntMass(9.01218)  },
		{"boron",      5, "B ", IntMass(10.81)    },
		{"carbon",     6, "C ", IntMass(12.011)   },
		{"nitrogen",   7, "N ", IntMass(14.0067)  },
		{"oxygen",     8, "O ", IntMass(15.999)   },
		{"fluorine",   9, "F ", IntMass(18.998403)},
		{"neon",      10, "Ne", IntMass(20.180)   },
		{"sodium",    11, "Na", IntMass(22.98977) },
		{"magnesium", 12, "Mg", IntMass(24.305)   },
		{"aluminium", 13, "Al", IntMass(26.98154) },
		{"silicon",   14, "Si", IntMass(28.086)   },
		{"phosphorus",15, "P ", IntMass(30.97376) },
		{"sulfur",    16, "S ", IntMass(32.07)    },
		{"chlorine",  17, "Cl", IntMass(35.453)   },
		{"argon",     18, "Ar", IntMass(39.948)   },
		{"potassium", 19, "K ", IntMass(39.0983)  },
		{"calcium",   20, "Ca", IntMass(40.08)    },
		{"scandium",  21, "Sc", IntMass(44.95591) },
		{"titanium",  22, "Ti", IntMass(47.867)   },
		{"vanadium",  23, "V ", IntMass(50.9415)  },
		{"chromium",  24, "Cr", IntMass(51.996)   },
		{"manganese", 25, "Mn", IntMass(54.93805) },
		{"iron",      26, "Fe", IntMass(55.84)    },
		{"cobalt",    27, "Co", IntMass(58.93320) },
		{"nickel",    28, "Ni", IntMass(58.693)   },
		{"copper",    29, "Cu", IntMass(63.55)    },
		{"zinc",      30, "Zn", IntMass(65.4)     },
		{"gallium",   31, "Ga", IntMass(69.723)   },
		{"germanium", 32, "Ge", IntMass(72.6)     },
		{"arsenic",   33, "As", IntMass(74.9216)  },
		{"selenium",  34, "Se", IntMass(79.0)     },
		{"bromine",   35, "Br", IntMass(79.904)   },
		{"krypton",   36, "Kr", IntMass(83.80)    },
		{"rubidium",  37, "Rb", IntMass(85.468)   },
		{"strontium", 38, "Sr", IntMass(87.62)    },
		{"yttrium",   39, "Y ", IntMass(88.9058)  },
		{"zirconium", 40, "Zr", IntMass(91.22)    },
		{"niobium",   41, "Nb", IntMass(92.9064)  },
		{"molybdenum",42, "Mo", IntMass(95.94)    },
		{"technetium",43, "Tc", IntMass(97.9072)  },
		{"ruthenium", 44, "Ru", IntMass(101.1)    },
		{"rhodium",   45, "Rh", IntMass(102.9055) },
		{"palladium", 46, "Pd", IntMass(106.42)   },
		{"silver",    47, "Ag", IntMass(107.868)  },
		{"cadmium",   48, "Cd", IntMass(112.41)   },
		{"indium",    49, "In", IntMass(114.82)   },
		{"tin",       50, "Sn", IntMass(118.71)   },
		{"antimony",  51, "Sb", IntMass(121.760)  },
		{"tellurium", 52, "Te", IntMass(127.6)    },
		{"iodine",    53, "I ", IntMass(126.9045) },
		{"xenon",     54, "Xe", IntMass(131.3)    },
		{"caesium",   55, "Cs", IntMass(132.9054) },
		{"barium",    56, "Ba", IntMass(137.33)   }
	};
	const unsigned MAX_NUCLEAR_INDEX = 56;

	/**
	* @brief finds nuclear with given mass
	* @note mass must be close to real nuclear mass otherwise will be exception
	*/
	inline unsigned find_nuclear(_I2T<MASS_>, real_t mass, real_t accuracy=IntMass(0.5))
	{
		assert(_GT(mass, IntMass(0.)));
		assert(_LT(mass, IntMass(137.331)));

		int ndx = (int)floor(0.5 * mass); // start approximation
		if (ndx > (int)MAX_NUCLEAR_INDEX) ndx = (int)MAX_NUCLEAR_INDEX;

		int h = nuclears[ndx].mass > mass ? -1 : 1;

		while (ndx >= 0 && ndx <= (int)MAX_NUCLEAR_INDEX
			&& fabs(nuclears[ndx].mass - mass) > accuracy)
			ndx += h;

		if (ndx > (int)MAX_NUCLEAR_INDEX || ndx < 0) ndx = 0; // unknown
			// it's possible due to such type as LP
		return (unsigned)ndx;
	}

	/**
	* @brief finds nuclear with given name
	*/
	INLINE unsigned find_nuclear(const std::string &name)
	{
		for (unsigned i=0; i<=MAX_NUCLEAR_INDEX; i++)
			if (trim(make_string(nuclears[i].name)) == trim(name)) return i;
		return 0;
	}

}
#endif
