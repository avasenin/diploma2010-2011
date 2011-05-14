#ifndef _RESIDUE__F9ED1116_412E_59a7_5CEB_EA459CD20E00__H
#define _RESIDUE__F9ED1116_412E_59a7_5CEB_EA459CD20E00__H

#include "molkern/__moldefs.h"

namespace molkern
{

	using namespace prgkern;

	struct AminoAcidName
	{
		char _1name;
		fstring _3name;
		const char *name;
	};
	const AminoAcidName amino_acid_name[] =
	{
		{ 'A', "ALA", "ALANINE"    },
		{ 'C', "CYS", "CYSTEINE"   },
		{ 'D', "ASP", "ASPARAGINE ACID"},
		{ 'E', "GLU", "GLUTAMINE ACID" },
		{ 'F', "PHE", "PHENYLALANINE"  },
		{ 'G', "GLY", "GLYCINE"    },
		{ 'H', "HIS", "HISTIDINE"  },
		{ 'I', "ILE", "ISOLEUCINE" },
		{ 'K', "LYS", "LYSINE"     },
		{ 'L', "LEU", "LEUCINE"    },
		{ 'M', "MET", "METHIONINE" },
		{ 'N', "ASN", "ASPARAGINE" },
		{ 'P', "PRO", "PROLINE"    },
		{ 'Q', "GLN", "GLUTAMINE"  },
		{ 'R', "ARG", "ARGININE"   },
		{ 'S', "SER", "SERINE"     },
		{ 'T', "THR", "THREONINE"  },
		{ 'Y', "TYR", "TYROSINE"   },
		{ 'V', "VAL", "VALINE"     },
		{ 'W', "TRP", "TRIPTOPHANE"}
	};

	struct __ResidueAtomdata
	{
		typedef fstring index_type;

		fstring name; // топологическое имя
		fstring fftype; // тип силового поля
		int resx; // используется как идентификатор отдельной молекулы (в растворе)
		real_t charge; // заряд атома

		static const unsigned max_A = 3; // предел числа связей с НЕ-водородами
		static const unsigned max_H = 3; // предел числа связей с водородами
		unsigned na; // число связей с НЕ-водородами
		unsigned nh; // число связей с водородами

		int nid[max_A]; // идентификаторы атомов, связных с заданным
		int nhid[max_A]; // идентификаторы водородов, связных с заданным
		vector_t X0; // позиция атома в остатке

		/**
		* @param contacts[out] - all contacts of atom
		*/
		unsigned get_contacts(int *contacts) const
		{
			for (unsigned i=0; i<na; i++) contacts[i] = nid[i];
			for (unsigned i=0; i<nh; i++) contacts[i + na] = nhid[i];
			return na + nh;
		}
	};

	template <int RESIDOME_TYPE>
	class basic_residue_
	{
		typedef struct Rotamer_
		{
			std::pair<fstring, fstring> axis;
			std::vector<fstring> atom;
		} _Rotamer;

		typedef __ResidueAtomdata             _ParamT;
		typedef typename _ParamT::index_type  _Index;
		typedef std::pair<_Index, _ParamT>    _Pair;
		typedef std::vector<_Pair>            _VectorT;

	public:

		typedef __ResidueAtomdata atom_type;
		typedef _Index index_type;
		typedef _Pair value_type;
		typedef _Rotamer rotamer_type;

		basic_residue_() : name("UNK"), has_connects_(false),
			has_external_connects_(false)
		{
			ext_contacts_[0] = nill;
			ext_contacts_[1] = nill;
		}

		int get_data_ndx(const _Index &ndx) const
		{
			for (unsigned i=0, sz=residue_atomdata_.size(); i<sz; i++)
				if (residue_atomdata_[i].first == ndx) return i;
			return nill;
		}

		const __ResidueAtomdata *get_data(int n) const
		{
			assert(_LT((unsigned)n, (unsigned)residue_atomdata_.size()));
			return &residue_atomdata_[n].second;
		}

		const __ResidueAtomdata *get_data(const char *name) const
		{ return get_data(get_data_ndx(_Index(name))); }

		const __ResidueAtomdata *get(_I2T<ATOM_>, const char *name) const
		{
			int n = get_data_ndx(_Index(name));
			if (n == nill)
			{
				_S msg = _S("residue doesn't have the atom %s", name);
				PRINT_ERR(msg);
			}
			return get_data(n);
		}

		unsigned count(_I2T<ATOM_>) const { return residue_atomdata_.size(); }
		int get(_I2T<EXT_CONTACT_PREV_>) const { return ext_contacts_[0]; }
		int get(_I2T<EXT_CONTACT_NEXT_>) const { return ext_contacts_[1]; }
		unsigned count(_I2T<ROTAMER_>) const { return rotamers_.size(); }

		const _Rotamer *get(_I2T<ROTAMER_>, unsigned i) const { return &rotamers_[i]; }

		const box_t &box() const { return box_; }

		bool has_connects() const { return has_connects_; }
		bool has_external_connects() const { return has_external_connects_; }
		unsigned count(_I2T<CHAIN_>) const { return chain_count_; }

		/// name of residue
		std::string name;

	protected:

		_VectorT residue_atomdata_;
		box_t box_;
		bool has_connects_;
		bool has_external_connects_;
		int ext_contacts_[2];
		unsigned chain_count_;

		std::vector<_Rotamer> rotamers_;
	};

	template <int RESIDOME_TYPE>
	INLINE bool operator<(const basic_residue_<RESIDOME_TYPE> &a1,
		const basic_residue_<RESIDOME_TYPE> &a2)
	{ return a1.name < a2.name; }

	template <int RESIDOME_TYPE>
	INLINE bool operator==(const basic_residue_<RESIDOME_TYPE> &a1,
		const basic_residue_<RESIDOME_TYPE> &a2)
	{ return a1.name == a2.name; }

	template <int> class Residue_;

}
#endif
