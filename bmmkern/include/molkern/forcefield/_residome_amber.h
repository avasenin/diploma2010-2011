#ifndef _RESIDOME_AMBER__F9ED1116_EFB9_53a5_11E1_EA4512720100__H
#define _RESIDOME_AMBER__F9ED1116_EFB9_53a5_11E1_EA4512720100__H

#include "molkern/__moldefs.h"
#include "molkern/forcefield/_residue.h"
#include "molkern/forcefield/_residue_amber.h"
#include "molkern/forcefield/_residome.h"

namespace molkern
{
	using namespace prgkern;

	namespace amber
	{
		struct _Equi
		{
			static const int pH_MIN = 3;
			static const int pH_MAX = 9;
			static const int pH_RANGE = pH_MAX - pH_MIN + 1;

			fstring resi;
			fstring equi[pH_RANGE];
		};
		const _Equi equis_[] =
		{
			// pH       3      4      5      6      7      8      9
			{ "CYS", {"CYS", "CYS", "CYS", "CYS", "CYS", "CYS", "CYX"} },
			{ "ASP", {"ASH", "ASP", "ASP", "ASP", "ASP", "ASP", "ASP"} },
			{ "GLU", {"GLH", "GLH", "GLU", "GLU", "GLU", "GLU", "GLU"} },
			{ "HIS", {"HIP", "HIP", "HIP", "HID", "HID", "HID", "HID"} },
			{ "LYS", {"LYS", "LYS", "LYS", "LYS", "LYS", "LYS", "LYN"} }
		};

		/**
		* @brief subset of all residues set
		* @note these residues change their names on terminal positions
		*       for example, VAL -> NVAL
		*/
		const fstring basic_residue_name[] =
		{
			"ALA", "ARG", "ASH", "ASN", "ASP",
			"CYM", "CYS", "CYX",
			"GLH", "GLN", "GLU", "GLY",
			"HID", "HIE", "HIP",
			"ILE",
			"LEU", "LYN", "LYS",
			"MET",
			"PHE", "PRO",
			"SER",
			"THR", "TRP", "TYR",
			"VAL"
		};
		const unsigned residue_count = STATIC_DIMENSION(basic_residue_name, fstring);

		struct _Short
		{
			const char *full_name;
			const char *short_name;
		};
		const _Short short_residue_name_[] =
		{
			{"CHCL3BOX",   "CL3"},
			{"MEOHBOX",    "MOH"},
			{"NMABOX",     "NMA"},
			{"POL3BOX",    "PL3"},
			{"SPCBOX",     "SPC"},
			{"TIP3PBOX",   "TP3"},
			{"TIP4PBOX",   "TP4"},
			{"TIP4PEWBOX", "T4E"}
		};
	}

	/**
	* @brief residue AMBER database specialization
	*/
	template <> class Residome_<RESIDOME_AMBER_>
	: public basic_residome_<RESIDOME_AMBER_>
	{
		typedef basic_residome_<RESIDOME_AMBER_>   _Base;
		typedef Params<RESIDUE_, RESIDOME_AMBER_>  _Residata;
		typedef std::vector<_Residata>             _VectorR;

		/// number of residues in all *.lib files for reserving memory
		static const int DEFAULT_RESIDUE_COUNT = 120;

	public:
		typedef _Base::residue_type  residue_type;

		Residome_() {}
		Residome_(const std::string &dirname, const std::string &fileregex=_S(".*?\\.lib$"))
		{ load_dir(dirname, fileregex); }

		/**
		* @brief gets modified residue name at given pH
		* @param residue standard PDB residue name
		* @param pH pH
		* @return modified AMBER residue name
		*/
		std::string get_equivalent(_I2T<RESIDUE_>, const std::string &residue,
			int pH=NORMAL_PH_WATER) const
		{
			assert(_EQ(residue.length(), trim(residue).length()));
				// допускаем, что всегда приходит триммированное имя

			if (residue.length() > 4) return residue;
				// нельзя сделать преобразование для псевдоостатков (из *.lib файла)

			if (pH < amber::_Equi::pH_MIN) pH = amber::_Equi::pH_MIN;
			if (pH > amber::_Equi::pH_MAX) pH = amber::_Equi::pH_MAX;
			static int sz = sizeof(amber::equis_) / sizeof(amber::_Equi);

			fstring residue__(residue);
			for (int i=0; i<sz; i++)
			{
				if (amber::equis_[i].resi == residue__)
				{
					return make_string(amber::equis_[i].equi[pH - amber::_Equi::pH_MIN]);
				}
			}
			return residue;
		}

		/**
		* @brief gets standard residue name from modified residue name
		* @param residue modified AMBER residue name
		* @return standard PDB residue name
		*/
		std::string get_equivalent(_I2T<PDB_RESIDUE_NAME_>,
			const std::string &residue) const
		{
			assert(_EQ(residue.length(), trim(residue).length()));
				// допускаем, что всегда приходит триммированное имя

			if (residue.length() > 4) return residue;
				// нельзя сделать преобразование для псевдоостатков (из *.lib файла)

			static int sz = sizeof(amber::equis_) / sizeof(amber::_Equi);
			fstring residue__(residue);
			for (int i=0; i<sz; i++)
			{
				for (int i__=0; i__<amber::_Equi::pH_RANGE; i__++)
				{
					if (amber::equis_[i].equi[i__] == residue__)
						return make_string(amber::equis_[i].resi);
				}
			}
			return residue;
		}

		std::string get_short_name(_I2T<RESIDUE_>, const std::string &residue) const
		{
			assert(_EQ(residue.length(), trim(residue).length()));
				// допускаем, что всегда приходит триммированное имя

			if (residue.length() < 4) return residue;
				// незачем делать преобразование для коротких имен

			static int sz = sizeof(amber::short_residue_name_) / sizeof(amber::_Short);
			for (int i=0; i<sz; i++)
			{
				if (::strcmp(amber::short_residue_name_[i].full_name, residue.c_str()) == 0)
					return _S(amber::short_residue_name_[i].short_name);
			}
			return residue;
		}

		bool load_dir(const std::string &dirname, const std::string &fileregex=_S(".*?\\.lib$"));
		bool load(const std::string &filename);

		std::string name(_I2T<NTERM_>, const std::string &name) const
		{
			assert(_EQ(name.length(), trim(name).length()));
				// допускаем, что всегда приходит триммированное имя

			if (name.length() > 4) return name;
				// нельзя сделать преобразование для псевдоостатков (из *.lib файла)

			return has_basic_residue_(name)? make_term('N', name) : name; }

		std::string name(_I2T<CTERM_>, const std::string &name) const
		{
			assert(_EQ(name.length(), trim(name).length()));
				// допускаем, что всегда приходит триммированное имя

			if (name.length() > 4) return name;
				// нельзя сделать преобразование для псевдоостатков (из *.lib файла)
			return has_basic_residue_(name)? make_term('C', name) : name; }

	protected:

		bool has_basic_residue_(const std::string &name) const
		{
			assert(_EQ(name.length(), trim(name).length()));
				// допускаем, что всегда приходит триммированное имя

			if (name.length() > 4) return false;
				// нельзя сделать преобразование для псевдоостатков (из *.lib файла)

			fstring name__(name);
			for (unsigned i=0; i<amber::residue_count; i++)
				if (name__ == amber::basic_residue_name[i]) return true;
			return false;
		}

		std::string make_term(char symbol, const std::string &name) const
		{
			assert(_EQ(name.length(), trim(name).length()));
				// допускаем, что всегда приходит триммированное имя

			if (name.length() >= (unsigned)MAX_ATOM_NAME)
			{
				std::string msg = _S("incorrect name of residue ") + _S(name);
				PRINT_ERR(msg);
			}
			return _S(1, symbol) + name;
		}

	private:
		_VectorR residata_;
	};

	typedef Residome_<RESIDOME_AMBER_>  Residome;

	#define TEMPLATE_HEADER
	#define TEMPLATE_ARG     RESIDOME_AMBER_

	TEMPLATE_HEADER
	inline bool Residome_<TEMPLATE_ARG>
	::load(const std::string &filename)
	{
		std::ifstream src(filename.c_str());
		if (!src)
		{
			std::string msg = _S("can't open file ") + filename;
			PRINT_ERR(msg);
		}
		char line[MAX_STRING_LEN];
		src.getline(line, MAX_STRING_LEN);
		while (::strncmp(&line[0], "!entry", 6) != 0)
			src.getline(line, MAX_STRING_LEN);
		{ // INPUT FOR RESIDUES
			while (!src.eof())
			{
				residata_.push_back(_Residata());
				_Residata &data__ = residata_[residata_.size() - 1];
				data__.load(src);
			}
		}
		LOADED_OK_MESSAGE(filename);
		return true;
	}

	TEMPLATE_HEADER
	inline bool Residome_<TEMPLATE_ARG>
	::load_dir(const std::string &dirname, const std::string &fileregex)
	{
		residata_.reserve(DEFAULT_RESIDUE_COUNT);
		_Base::dirname_ = dirname; // save for something in future
		_Base::dirname_ += _S("AMBER/");
		std::string filename;

		RegEx regex(fileregex);
		directory_iterator it(_Base::dirname_), ite;
		for (; it!=ite; ++it)
		{
			filename = *it;
			if (regex.match(filename))
			{
				filename = _Base::dirname_ + filename;
				load(filename.c_str());
			}
		}
		std::sort(residata_.begin(), residata_.end());

		int sz = (int)residata_.size();
		_Base::residues_.resize(sz);
		for (int i=0; i<sz; i++)
			_Base::residues_[i].init(residata_[i]);

	#ifdef FULL_TOPOLOGY_DEBUG
		_Base::print();
	#endif
		return true;
	}

	#undef TEMPLATE_HEADER
	#undef TEMPLATE_ARG

}
#endif
