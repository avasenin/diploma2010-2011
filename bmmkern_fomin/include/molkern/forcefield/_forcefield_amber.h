#ifndef _FORCRFIELD_AMBER__0077A726_F6D3_5e58_EF3F_E84517630500__H
#define _FORCRFIELD_AMBER__0077A726_F6D3_5e58_EF3F_E84517630500__H

#include "molkern/__moldefs.h"
#include "molkern/forcefield/_forcefield.h"

namespace molkern
{

	#define SKIP_UNTIL_EMPTY_VALUE(src, NN, NA, NX) { \
		char buffer___[MAX_STRING_LEN]; \
		char buff___[NA + 1]; \
		bool stop = false; \
		do { \
			src.getline(buffer___, MAX_STRING_LEN); \
			if (::strlen(buffer___) < (NA + NX) * NN) break; \
			for (int i=0, k=0; i<NN; i++) { \
				::strncpy(buff___, &buffer___[k], NA); buff___[NA] = 0; \
				if (isspace(_S(buff___))) { stop = true; break; } \
				k += NA + NX; \
			} \
		} while (!stop); \
	}

	#ifdef FULL_TOPOLOGY_DEBUG
		#ifndef NDEBUG
		#define USES_EQUIVALENTS(ndx, index__) { \
			if (ndx != index__) { \
				std::string msg = _S("[WARNING] Database uses equivalent : ") \
					+ _S(make_string(ndx)) + _S(" -> ") \
					+ _S(make_string(index__)); \
				PRINT_MESSAGE(msg); \
			} \
		}
		#else
			#define USES_EQUIVALENTS(ndx, index__) {}
		#endif
	#else
		#define USES_EQUIVALENTS(ndx, index__)
	#endif

	template <> class Forcefield_<FORCEFIELD_AMBER_> :
		public __Forcefield<FORCEFIELD_AMBER_>
	{
		typedef __Forcefield <FORCEFIELD_AMBER_> _Base;

	public:
		using _Base::get_equivalent_data;

		Forcefield_() {}
		Forcefield_(const std::string &dirname) { load_dir(dirname); }

		bool load_dir(const std::string &dirname);
		bool load(const std::string &filename);

		/**
		* @brief extract scale data
		*/
		real_t get_scale(_I2T<COUL14_>) const { return AMBER_SCALE_COUL14; }
		real_t get_scale(_I2T<VDW14_>) const { return AMBER_SCALE_VDW14; }

	protected:

		/// @note we used that all atom symbols in *.dat file are ordered
		void make_equivalent_data_(_Base::_ParamE*, unsigned start_pos=0)
		{
			_Base::_PairE pair_node;
			_Base::_IndexE &node_index = pair_node.first;
			_Base::_ParamE &node_param = pair_node.second;
			node_param.n = 0;

			unsigned sz = _Base::symb_params_.size();
			if (start_pos >= sz) return;
				// no new equivalents

			fstring prev_name = _Base::symb_params_[start_pos].second.name;
			real_t prev_mass = _Base::symb_params_[start_pos].second.mass;
				// to avoid the pushing of empty first element
			real_t mass = 0.;

			for (unsigned i=start_pos; i<sz; i++)
			{
				_Base::_IndexS &symb_index = _Base::symb_params_[i].first;
				_Base::_ParamS &symb_param = _Base::symb_params_[i].second;
				mass = symb_param.mass;

				if (fabs(mass - prev_mass) > 0.5)
				{
					if (node_param.n > 1) // skip atoms without equivalents
					{
						// save atoms with several equivalents only
						int nuclear = find_nuclear(_I2T<MASS_>(), prev_mass);
						node_index = nuclears[nuclear].name;
						_Base::eqty_params_.push_back(pair_node);
					}
					node_param.n = 0;
				}

				if (node_param.n >= MAX_EQUIVALENT)
				{
					std::string msg = make_string("MAX_EQUIVALENT limit [%d] is exceed", MAX_EQUIVALENT);
					PRINT_ERR(msg);
				}
				node_param.eqv[node_param.n++] = symb_index[0];
				prev_mass = mass;
				prev_name = symb_param.name;
			}
		}
	};
	typedef Forcefield_<FORCEFIELD_AMBER_> ForceField;

	#define TEMPLATE_HEADER
	#define TEMPLATE_ARG     FORCEFIELD_AMBER_

	TEMPLATE_HEADER
	inline bool Forcefield_<TEMPLATE_ARG>::load(const std::string &filename)
	{
		char buffer[MAX_STRING_LEN];
		char buff__[MAX_STRING_LEN];
		std::ifstream src(filename.c_str());
		if (!src)
		{
			std::string msg = _S("[ERROR] can't open file ") + filename;
			PRINT_BREAK(msg);
		}
		SKIP_LINE(src);
			// A title for identification of the parameter set.
		{ // INPUT FOR ATOM SYMBOLS AND MASSES
			_Base::_PairS obj;
			while (extract_object(obj, &src) == CODE_SUCCESS)
				_Base::symb_params_.push_back(obj);
		}
		SKIP_UNTIL_EMPTY_VALUE(src, 20, 2, 2);
			// INPUT FOR ATOM SYMBOLS THAT ARE HYDROPHILIC
		{ // INPUT FOR BOND LENGTH PARAMETERS
			_Base::_PairB obj;
			while (extract_object(obj, &src) == CODE_SUCCESS)
				_Base::bond_params_.push_back(obj);
		}
		{ // INPUT FOR BOND ANGLE PARAMETERS
			_Base::_PairA obj;
			while (extract_object(obj, &src) == CODE_SUCCESS)
				_Base::angl_params_.push_back(obj);
		}
		{ // INPUT FOR DIHEDRAL PARAMETERS
			_Base::_PairD obj;
			while (extract_object(obj, &src) == CODE_SUCCESS)
				_Base::dihe_params_.push_back(obj);
		}
		SKIP_UNTIL_EMPTY_LINE(src);
			// INPUT FOR IMPROPER DIHEDRAL PARAMETERS
		SKIP_UNTIL_EMPTY_LINE(src);
			// INPUT FOR H-BOND 10-12 POTENTIAL PARAMETERS
		{ // INPUT FOR EQUIVALENCING ATOM SYMBOLS
			_Base::_PairE obj;
			while (extract_object(obj, &src) == CODE_SUCCESS)
				_Base::equi_params_.push_back(obj);
		}
		{ // INPUT FOR THE 6-12 POTENTIAL PARAMETERS
			src.getline(buffer, MAX_STRING_LEN);
			::strncpy(buff__, &buffer[10], 2); buff__[2] = 0;
			if (::strncmp(buff__, "RE", 2) != 0)
			{
				std::string msg = _S(" [NO_IMPLEMENTATION] buff_ = ") + _S(buff__);
				PRINT_ERR(msg);
			}
			_Base::_PairT obj;
			while (extract_object(obj, &src) == CODE_SUCCESS)
				_Base::atom_params_.push_back(obj);
		}
		LOADED_OK_MESSAGE(filename);
		return true;
	}

	TEMPLATE_HEADER
	inline bool Forcefield_<TEMPLATE_ARG>
	::load_dir(const std::string &dirname)
	{
		_Base::dirname_ = dirname; // save for something in future
		std::string filename = dirname + _S("AMBER/__parm99.dat");
		load(filename);
		make_equivalent_data_((_Base::_ParamE *)0);
			// make equi types before std::sort

		unsigned start_pos = _Base::symb_params_.size();
			// pos to treat new records in symb_params
		filename = dirname + _S("AMBER/__gaff.dat");
		load(filename);
		make_equivalent_data_((_Base::_ParamE *)0, start_pos);
			// make equi types before std::sort()
			// start pos for the making is equal to symb_params_.size()
			// the skipping of start pos may produce the code with the low efficiency

		std::sort(_Base::symb_params_.begin(), _Base::symb_params_.end());
		std::sort(_Base::atom_params_.begin(), _Base::atom_params_.end());
		std::sort(_Base::bond_params_.begin(), _Base::bond_params_.end());
		std::sort(_Base::angl_params_.begin(), _Base::angl_params_.end());
		std::sort(_Base::dihe_params_.begin(), _Base::dihe_params_.end());
		std::sort(_Base::equi_params_.begin(), _Base::equi_params_.end());

	#ifdef FULL_TOPOLOGY_DEBUG
		_Base::print();
	#endif
		return true;
	}

	#undef TEMPLATE_HEADER
	#undef TEMPLATE_ARG

	#undef SKIP_UNTIL_EMPTY_VALUE
	#undef USES_EQUIVALENTS

}
#endif

