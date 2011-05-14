#ifndef _RESIDOME__F9ED1116_2CEA_59cc_BFDF_EA45A58E0100__H
#define _RESIDOME__F9ED1116_2CEA_59cc_BFDF_EA45A58E0100__H

#include "molkern/__moldefs.h"
#include "molkern/forcefield/_residue.h"

namespace molkern
{
	using namespace prgkern;

	template <int RESIDOME_TYPE> class Residome_;

	/**
	* @brief Residue database template
	* @param RESIDOME_TYPE - residome type (AMBER, GROMACS etc)
  */
	template <int RESIDOME_TYPE>
	class basic_residome_
	{
	protected:

		typedef Residue_<RESIDOME_TYPE> _Residue;
		typedef std::vector<_Residue>   _VectorR;

	public:
		typedef _Residue residue_type;

		/**
		* @brief get access to residue data
		* @param name - residue name
		* @return NULL in case of the residue isn't found
		*/
		const _Residue *get_data(const std::string &name) const
		{ return get_data_ndx__(name); }

		/**
		* @brief defines whether residome includes residue name
		* @param name - residue name
		* @return false in case of the residome don't includes residue
		*/
		bool has_residue(const std::string &name) const
		{ return get_data_ndx__(name) != NULL; }

		unsigned count(_I2T<RESIDUE_>) const { return residues_.size(); }

		/**
		* @brief debug print
		*/
		void print() const;

	protected:

		basic_residome_() {}

		const _Residue *get_data_ndx__(const std::string &name) const
		{
			typedef typename _VectorR::const_iterator cit;
			_Residue tmp; tmp.name = name;
			typename std::pair<cit, cit> pit
				= std::equal_range(residues_.begin(), residues_.end(), tmp);

			if (pit.first == pit.second) { return NULL; }
			return &(*pit.first);
		}

		_VectorR residues_;
		std::string dirname_;
	};

	#define TEMPLATE_HEADER  template <int RESIDOME_TYPE>
	#define TEMPLATE_ARG     RESIDOME_TYPE

	TEMPLATE_HEADER
	inline void basic_residome_<TEMPLATE_ARG>
	::print() const
	{
		char line[MAX_LINE_LEN];
		std::string msg =
				_S("\n----------------------------------------")
			+ _S("\n  residues readed ...")
			+ _S("\n----------------------------------------");
		PRINT_MESSAGE(msg);

		std::set<fstring> atom_names;
		for (int i=0, sz=(int)residues_.size(); i<sz; i++)
		{
			::sprintf(line, "%4d %s", i, make_string(residues_[i].name).c_str());
			PRINT_MESSAGE(line);

			unsigned n = residues_[i].count(_I2T<ATOM_>());
			for (unsigned i__=0; i__<n; i__++) atom_names.insert(residues_[i].get_data(i__)->name);
		}

		msg = _S("\n----------------------------------------")
			+ _S("\n  topologied atoms readed ...")
			+ _S("\n----------------------------------------");
		PRINT_MESSAGE(msg);

		typename std::set<fstring>::const_iterator it = atom_names.begin(), ite = atom_names.end();
		for(; it!=ite; ++it) PRINT_MESSAGE(make_string(*it));
	}

	#undef TEMPLATE_HEADER
	#undef TEMPLATE_ARG

}
#endif
