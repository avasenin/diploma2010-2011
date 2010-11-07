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
		{"unknown",    0, "  ", 0.0000   },
		{"hydrogen",   1, "H ", 1.0079   },
		{"helium",     2, "He", 4.00260  },
		{"lithium",    3, "Li", 6.94     },
		{"beryllium",  4, "Be", 9.01218  },
		{"boron",      5, "B ", 10.81    },
		{"carbon",     6, "C ", 12.011   },
		{"nitrogen",   7, "N ", 14.0067  },
		{"oxygen",     8, "O ", 15.999   },
		{"fluorine",   9, "F ", 18.998403},
		{"neon",      10, "Ne", 20.180   },
		{"sodium",    11, "Na", 22.98977 },
		{"magnesium", 12, "Mg", 24.305   },
		{"aluminium", 13, "Al", 26.98154 },
		{"silicon",   14, "Si", 28.086   },
		{"phosphorus",15, "P ", 30.97376 },
		{"sulfur",    16, "S ", 32.07    },
		{"chlorine",  17, "Cl", 35.453   },
		{"argon",     18, "Ar", 39.948   },
		{"potassium", 19, "K ", 39.0983  },
		{"calcium",   20, "Ca", 40.08    },
		{"scandium",  21, "Sc", 44.95591 },
		{"titanium",  22, "Ti", 47.867   },
		{"vanadium",  23, "V ", 50.9415  },
		{"chromium",  24, "Cr", 51.996   },
		{"manganese", 25, "Mn", 54.93805 },
		{"iron",      26, "Fe", 55.84    },
		{"cobalt",    27, "Co", 58.93320 },
		{"nickel",    28, "Ni", 58.693   },
		{"copper",    29, "Cu", 63.55    },
		{"zinc",      30, "Zn", 65.4     },
		{"gallium",   31, "Ga", 69.723   },
		{"germanium", 32, "Ge", 72.6     },
		{"arsenic",   33, "As", 74.9216  },
		{"selenium",  34, "Se", 79.0     },
		{"bromine",   35, "Br", 79.904   },
		{"krypton",   36, "Kr", 83.80    },
		{"rubidium",  37, "Rb", 85.468   },
		{"strontium", 38, "Sr", 87.62    },
		{"yttrium",   39, "Y ", 88.9058  },
		{"zirconium", 40, "Zr", 91.22    },
		{"niobium",   41, "Nb", 92.9064  },
		{"molybdenum",42, "Mo", 95.94    },
		{"technetium",43, "Tc", 97.9072  },
		{"ruthenium", 44, "Ru", 101.1    },
		{"rhodium",   45, "Rh", 102.9055 },
		{"palladium", 46, "Pd", 106.42   },
		{"silver",    47, "Ag", 107.868  },
		{"cadmium",   48, "Cd", 112.41   },
		{"indium",    49, "In", 114.82   },
		{"tin",       50, "Sn", 118.71   },
		{"antimony",  51, "Sb", 121.760  },
		{"tellurium", 52, "Te", 127.6    },
		{"iodine",    53, "I ", 126.9045 },
		{"xenon",     54, "Xe", 131.3    },
		{"caesium",   55, "Cs", 132.9054 },
		{"barium",    56, "Ba", 137.33   }
	};
	const unsigned MAX_NUCLEAR_INDEX = 56;

	/**
	* @brief finds nuclear with given mass
	* @note mass must be close to real nuclear mass otherwise will be exception
	*/
	inline unsigned find_nuclear(_I2T<MASS_>, real_t mass, real_t accuracy=0.5)
	{
		assert(_GT(mass, 0.));
		assert(_LT(mass, 137.331));

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
