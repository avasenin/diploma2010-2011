#ifndef FORCEFIELD__F9ED1116_E190_5bc1_48A9_DB43E7FD0700__H
#define FORCEFIELD__F9ED1116_E190_5bc1_48A9_DB43E7FD0700__H

#include "molkern/__moldefs.h"
#include "molkern/forcefield/_fparams.h"

namespace molkern
{

	using namespace prgkern;

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

	template <int FORCEFIELD_TYPE> class Forcefield_ {};

	/**
	* @brief ForceField template
	* @param FORCEFIELD_TYPE - identifier of forcefield
	* @note
	*  The given class is base class for different forcefield specializations,
	*  such as Forcefield<AMBER>, Forcefield<GROMACS> etc. The class is used
	*  as database. The properties of forcefields are describes in
	*  specializations.
	*/
	template <int FORCEFIELD_TYPE>
	class __Forcefield
	{

	public:

		typedef Params<NUCLEAR_, FORCEFIELD_TYPE> _ParamS; // params of nuclear Symbols
		typedef Params<ATOM_,    FORCEFIELD_TYPE> _ParamT; // params of aToms
		typedef Params<BOND_,    FORCEFIELD_TYPE> _ParamB; // params of Bonds
		typedef Params<ANGLE_,   FORCEFIELD_TYPE> _ParamA; // params of Angles
		typedef Params<TORSION_, FORCEFIELD_TYPE> _ParamD; // params of torsions (Dihedrals)
		typedef Params<EQUI_,    FORCEFIELD_TYPE> _ParamE; // params of Equivalents

		typedef typename _ParamS::index_type _IndexS;
		typedef typename _ParamT::index_type _IndexT;
		typedef typename _ParamA::index_type _IndexA;
		typedef typename _ParamB::index_type _IndexB;
		typedef typename _ParamD::index_type _IndexD;
		typedef typename _ParamE::index_type _IndexE;

		typedef std::pair<_IndexS, _ParamS> _PairS;
		typedef std::pair<_IndexT, _ParamT> _PairT;
		typedef std::pair<_IndexB, _ParamB> _PairB;
		typedef std::pair<_IndexA, _ParamA> _PairA;
		typedef std::pair<_IndexD, _ParamD> _PairD;
		typedef std::pair<_IndexE, _ParamE> _PairE;

		typedef std::vector<_PairS> _VectorS;
		typedef std::vector<_PairT> _VectorT;
		typedef std::vector<_PairB> _VectorB;
		typedef std::vector<_PairA> _VectorA;
		typedef std::vector<_PairD> _VectorD;
		typedef std::vector<_PairE> _VectorE;

		__Forcefield() {}

		/**
		* @brief get param data index
		* @note used for atoms only
		* @param ndx - index of param data (forcefield + atom types)
		* @return nill in case of index don't found
		*/
		int get_data_ndx(_I2T<ATOM_>, const _IndexT &ndx) const
		{ return get_data_ndx__(ndx, atom_params_); }

		/**
		* @brief get param data
		* @param param - adress to copy param data (or NULL)
		* @param ndx - index of param data (forcefield + atom types)
		* @return false in case of index don't found
		*/
		bool get_data(_ParamS *param, const _IndexS &ndx) const
		{ return get_data__(param, ndx, symb_params_); }

		bool get_data(_ParamT *param, const _IndexT &ndx) const
		{ return get_data__(param, ndx, atom_params_); }

		bool get_data(_ParamB *param, const _IndexB &ndx) const
		{ return get_data__(param, ndx, bond_params_); }

		bool get_data(_ParamA *param, const _IndexA &ndx) const
		{ return get_data__(param, ndx, angl_params_); }

		bool get_data(_ParamD *param, const _IndexD &ndx) const
		{ return get_data__(param, ndx, dihe_params_); }

		const _ParamT *get(_I2T<ATOM_>, int n) const
		{ return &atom_params_[n].second; }

		real_t get(_I2T<MAX_VDW_RADIUS_>) const;

		/**
		* @brief get equivalent param data
		* @param param - adress to copy param data (or NULL)
		* @param ndx - index of param data (forcefield + atom types)
		* @return nill in case of index don't found
		*/
		fstring get(_I2T<EQUI_>, const _IndexE &ndx) const;
		bool get_equivalent_data(_ParamB *param, const _IndexB &ndx) const;
		bool get_equivalent_data(_ParamA *param, const _IndexA &ndx) const;
		bool get_equivalent_data(_ParamD *param, const _IndexD &ndx) const;

		typename _VectorT::const_iterator begin(_PairT *) const { return atom_params_.begin(); }
		typename _VectorT::const_iterator end  (_PairT *) const { return atom_params_.end(); }

		/**
		* @brief debug print
		*/
		void print() const;

	protected:
		/**
		* @brief gets distance between equivalent type (for bond, angle, tors)
		*/
		int get(_I2T<EQTY_DISTANCE_>, const fstring &name,
			const fstring &ndx1, const fstring &ndx2) const;

		/**
		* @brief gets index of data in ordered array
		*/
		template <typename _Index, typename _Param>
		int get_data_ndx__(const _Index &ndx,
			const std::vector<std::pair<_Index, _Param> > &params__) const
		{
			typedef std::pair<_Index, _Param> _Pair;
			typedef typename std::vector<_Pair>::const_iterator _It;

			_Pair tmp; tmp.first = ndx;
			typename std::pair<_It, _It> pit
				= std::equal_range(params__.begin(), params__.end(), tmp);
			if (pit.first == pit.second) return nill;

			return (int)std::distance(params__.begin(), pit.first);
		}

		/**
		* @brief gets data from ordered array
		*/
		template <typename _Index, typename _Param>
		bool get_data__(_Param *param, const _Index &ndx,
			const std::vector<std::pair<_Index, _Param> > &params__) const
		{
			typedef std::pair<_Index, _Param> _Pair;
			typedef typename std::vector<_Pair>::const_iterator _It;

			_Pair tmp; tmp.first = ndx;
			typename std::pair<_It, _It> pit
				= std::equal_range(params__.begin(), params__.end(), tmp);
			if (pit.first == pit.second) return false;

			if (param != NULL) *param = (*pit.first).second;
			return true;
		}

		/**
		* @brief gets data from non-ordered array
		*/
		template <typename _Index, typename _Param>
		bool find_data__(_Param *param, const _Index &ndx,
			const std::vector<std::pair<_Index, _Param> > &params__) const
		{
			typedef std::pair<_Index, _Param> _Pair;
			typedef typename std::vector<_Pair>::const_iterator _It;

			_Pair tmp; tmp.first = ndx;
			_It it = std::find(params__.begin(), params__.end(), tmp);
			if (it == params__.end()) return false;

			if (param != NULL) *param = (*it).second;
			return true;
		}

		_VectorS symb_params_;
		_VectorT atom_params_;
		_VectorB bond_params_;
		_VectorA angl_params_;
		_VectorD dihe_params_;
		_VectorE equi_params_; // EQuivalent types (for vdw)
		_VectorE eqty_params_; // EQuivalent TYpes (for bonds, angls, tors)

		mutable _VectorB equi_bonds_;
		mutable _VectorA equi_angls_;
		mutable _VectorD equi_dihes_;

		std::string dirname_;
	};

	#define TEMPLATE_HEADER  template <int FORCEFIELD_TYPE>
	#define TEMPLATE_ARG     FORCEFIELD_TYPE

	TEMPLATE_HEADER
	INLINE real_t __Forcefield<TEMPLATE_ARG>
	::get(_I2T<MAX_VDW_RADIUS_>) const
	{
		int sz = (int)atom_params_.size();
		real_t max_radius = 0;
		for (int i=0; i<sz; i++)
		{
			const _ParamT &param = atom_params_[i].second;
			if (max_radius < param.radius) max_radius = param.radius;
		}
		return max_radius;
	}

	TEMPLATE_HEADER
	INLINE fstring __Forcefield<TEMPLATE_ARG>
	::get(_I2T<EQUI_>, const _IndexE &ndx) const
	{
		// very unefficient for large arrays BUT(!)
		// all arrays are very small
		for (int i=0, sz=(int)equi_params_.size(); i<sz; i++)
		{
			const _ParamE &param = equi_params_[i].second;
			for (int k=0; k<param.n; k++)
				if (param.eqv[k] == ndx[0]) return param.eqv[0];
		}
		return ndx[0];
	}

	TEMPLATE_HEADER
	INLINE int __Forcefield<TEMPLATE_ARG>
	::get(_I2T<EQTY_DISTANCE_>, const fstring &name,
		const fstring &ndx0, const fstring &ndx1) const
	{
		int pos0 = nill, pos1 = nill;
		for (unsigned  i=0, sz=(unsigned)eqty_params_.size(); i<sz; i++)
		{
			if (eqty_params_[i].first != name) continue;
			const _ParamE &param = eqty_params_[i].second;

			pos0 = nill, pos1 = nill;
			for (unsigned k=0; k<(unsigned)param.n; k++)
			{
				if (param.eqv[k] == ndx0) pos0 = k;
				if (param.eqv[k] == ndx1) pos1 = k;
			}
			if (pos0 != nill && pos1 != nill) return /*std::*/abs(pos1 - pos0);
		}
		return nill; // infinite difference
	}

	TEMPLATE_HEADER
	inline bool __Forcefield<TEMPLATE_ARG>
	::get_equivalent_data(_ParamB *param, const _IndexB &ndx) const
	{
		if (get_data__(param, ndx, equi_bonds_))
		{
		#ifdef FULL_TOPOLOGY_DEBUG
			std::string msg = _S("          The equivalent has been used");
			PRINT_MESSAGE(msg);
		#endif
			return true;
		}
		_IndexE ndxa0(ndx[0]), ndxa1(ndx[1]);
		unsigned pos0 = get_data_ndx__(ndxa0, symb_params_);
		unsigned pos1 = get_data_ndx__(ndxa1, symb_params_);

		fstring name0 = symb_params_[pos0].second.name;
		fstring name1 = symb_params_[pos1].second.name;

		int min_distance = mult2(MAX_EQUIVALENT + 1);
		int min_index = nill;

		for (int i=0,sz=bond_params_.size(); i<sz; i++)
		{
			const _IndexB &index = bond_params_[i].first;
			_IndexE ndxa0__(index[0]), ndxa1__(index[1]);

			unsigned pos__ = get_data_ndx__(ndxa0__, symb_params_);
			fstring name__ = symb_params_[pos__].second.name;
			if (name__ != name0) continue;

			pos__ = get_data_ndx__(ndxa1__, symb_params_);
			name__ = symb_params_[pos__].second.name;
			if (name__ != name1) continue;

			int n0 = get(_I2T<EQTY_DISTANCE_>(), name0, ndx[0], index[0]);
			int n1 = get(_I2T<EQTY_DISTANCE_>(), name1, ndx[1], index[1]);
			int distance = abs(n0) + abs(n1);
			if (min_distance > distance)
			{
				min_distance = distance;
				min_index = i;
			}
		}

		if (min_index != nill)
		{
			*param = bond_params_[min_index].second;
			equi_bonds_.push_back(_PairB(ndx, *param));
		#ifdef FULL_TOPOLOGY_DEBUG
			const _IndexB &index = bond_params_[min_index].first;
			{
				std::string msg = _S("          The equivalent ") + make_string(ndx)
					+ _S(" -> ") + make_string(index) + _S(" has been created.");
				PRINT_MESSAGE(msg);
			}
		#endif
			std::sort(equi_bonds_.begin(), equi_bonds_.end());
				// sort for fast search
			return true;
		}
		return false;
	}

	TEMPLATE_HEADER
	inline bool __Forcefield<TEMPLATE_ARG>
	::get_equivalent_data(_ParamA *param, const _IndexA &ndx) const
	{
		if (get_data__(param, ndx, equi_angls_))
		{
		#ifdef FULL_TOPOLOGY_DEBUG
			std::string msg = _S("          The equivalent has been used");
			PRINT_MESSAGE(msg);
		#endif
			return true;
		}
		_IndexE ndxa0(ndx[0]), ndxa1(ndx[1]), ndxa2(ndx[2]);

		unsigned pos0 = get_data_ndx__(ndxa0, symb_params_);
		unsigned pos1 = get_data_ndx__(ndxa1, symb_params_);
		unsigned pos2 = get_data_ndx__(ndxa2, symb_params_);
		fstring name0 = symb_params_[pos0].second.name;
		fstring name1 = symb_params_[pos1].second.name;
		fstring name2 = symb_params_[pos2].second.name;
		int min_distance = 3 * (MAX_EQUIVALENT + 1);
		int min_index = nill;

		for (int i=0,sz=angl_params_.size(); i<sz; i++)
		{
			const _IndexA &index = angl_params_[i].first;
			_IndexE ndxa0__(index[0]), ndxa1__(index[1]), ndxa2__(index[2]);
			unsigned pos__ = get_data_ndx__(ndxa0__, symb_params_);
			fstring name__ = symb_params_[pos__].second.name;
			if (name__ != name0) continue;

			pos__ = get_data_ndx__(ndxa1__, symb_params_);
			name__ = symb_params_[pos__].second.name;
			if (name__ != name1) continue;

			pos__ = get_data_ndx__(ndxa2__, symb_params_);
			name__ = symb_params_[pos__].second.name;
			if (name__ != name2) continue;

			int n0 = get(_I2T<EQTY_DISTANCE_>(), name0, ndx[0], index[0]);
			int n1 = get(_I2T<EQTY_DISTANCE_>(), name1, ndx[1], index[1]);
			int n2 = get(_I2T<EQTY_DISTANCE_>(), name2, ndx[2], index[2]);

			int distance = abs(n0) + abs(n1) + abs(n2);
			if (min_distance > distance)
			{
				min_distance = distance;
				min_index = i;
			}
		}

		if (min_index != nill)
		{
			if (param != NULL) *param = angl_params_[min_index].second;
			equi_angls_.push_back(_PairA(ndx, *param));
		#ifdef FULL_TOPOLOGY_DEBUG
			const _IndexA &index = angl_params_[min_index].first;
			{
				std::string msg = _S("          The equivalent ") + make_string(ndx)
					+ _S(" -> ") + make_string(index) + _S(" has been created.");
				PRINT_MESSAGE(msg);
			}
		#endif
			std::sort(equi_angls_.begin(), equi_angls_.end());
				// sort for fast search
			return true;
		}
		return false;
	}

	TEMPLATE_HEADER
	inline bool __Forcefield<TEMPLATE_ARG>
	::get_equivalent_data(_ParamD *param, const _IndexD &ndx) const
	{
		if (get_data__(param, ndx, equi_dihes_))
		{
		#ifdef FULL_TOPOLOGY_DEBUG
			std::string msg = _S("          The equivalent has been used");
			PRINT_MESSAGE(msg);
		#endif
			return true;
		}

		_IndexD index__;
		//<uses X-*-*-*> equivalents
		index__ = _IndexD(fstring("X"), ndx[1], ndx[2], ndx[3]);
		if (get_data__(param, index__, dihe_params_))
		{
		#ifdef FULL_TOPOLOGY_DEBUG
			std::string msg = _S("          The equivalent ")
				+ _S(make_string(index__)) + _S(" will be used.");
			PRINT_MESSAGE(msg);
		#endif
			return true;
		}
		index__ = _IndexD(ndx[0], ndx[1], ndx[2], fstring("X"));
		if (get_data__(param, index__, dihe_params_))
		{
		#ifdef FULL_TOPOLOGY_DEBUG
			std::string msg = _S("          The equivalent ")
				+ _S(make_string(index__)) + _S(" will be used.");
			PRINT_MESSAGE(msg);
		#endif
			return true;
		}

		//<uses X-*-*-X> equivalents
		index__ = _IndexD(fstring("X"), ndx[1], ndx[2], fstring("X"));
		if (get_data__(param, index__, dihe_params_))
		{
		#ifdef FULL_TOPOLOGY_DEBUG
			std::string msg = _S("          The equivalent ")
				+ _S(make_string(index__)) + _S(" will be used.");
			PRINT_MESSAGE(msg);
		#endif
			return true;
		}

		//<uses X-X-*-*> equivalents
		index__ = _IndexD(fstring("X"), fstring("X"), ndx[2], ndx[3]);
		if (get_data__(param, index__, dihe_params_))
		{
		#ifdef FULL_TOPOLOGY_DEBUG
			std::string msg = _S("          The equivalent ")
				+ _S(make_string(index__)) + _S(" will be used.");
			PRINT_MESSAGE(msg);
		#endif
			return true;
		}
		index__ = _IndexD(ndx[0], ndx[1], fstring("X"), fstring("X"));
		if (get_data__(param, index__, dihe_params_))
		{
		#ifdef FULL_TOPOLOGY_DEBUG
			std::string msg = _S("          The equivalent ")
				+ _S(make_string(index__)) + _S(" will be used.");
			PRINT_MESSAGE(msg);
		#endif
			return true;
		}

		_IndexB ndxb(ndx[1], ndx[2]);
			// auto ordering of indexes ('cos X-*-*-X might another order)
		_IndexE ndxa1(ndxb[0]), ndxa2(ndxb[1]);
			// use correct order

		unsigned pos1 = get_data_ndx__(ndxa1, symb_params_);
		unsigned pos2 = get_data_ndx__(ndxa2, symb_params_);
		fstring name1 = symb_params_[pos1].second.name;
		fstring name2 = symb_params_[pos2].second.name;
		int min_distance = mult2(MAX_EQUIVALENT + 1);
		int min_index = nill;

		for (int i=0,sz=dihe_params_.size(); i<sz; i++)
		{
			const _IndexD &index = dihe_params_[i].first;
			_IndexE ndxa1__(index[1]), ndxa2__(index[2]);

			unsigned pos__ = get_data_ndx__(ndxa1__, symb_params_);
			fstring name__ = symb_params_[pos__].second.name;
			if (name__ != name1) continue;

			pos__ = get_data_ndx__(ndxa2__, symb_params_);
			name__ = symb_params_[pos__].second.name;
			if (name__ != name2) continue;

			int n1 = get(_I2T<EQTY_DISTANCE_>(), name1, ndx[1], index[1]);
			int n2 = get(_I2T<EQTY_DISTANCE_>(), name2, ndx[2], index[2]);

			int distance = abs(n1) + abs(n2);
			if (min_distance > distance)
			{
				min_distance = distance;
				min_index = i;
			}
		}

		if (min_index != nill)
		{
			if (param != NULL) *param = dihe_params_[min_index].second;
			equi_dihes_.push_back(_PairD(ndx, *param));
		#ifdef FULL_TOPOLOGY_DEBUG
			const _IndexD &index = dihe_params_[min_index].first;
			{
				std::string msg = _S("          The equivalent ")
					+ make_string(index) + _S(" has been created.");
				PRINT_MESSAGE(msg);
			}
		#endif
			std::sort(equi_dihes_.begin(), equi_dihes_.end());
				// sort for fast search
			return true;
		}
		return false;
	}

	TEMPLATE_HEADER
	inline void __Forcefield<TEMPLATE_ARG>
	::print() const
	{
		std::string msg =
				_S("\n----------------------------------------")
			+ _S("\n  symbols readed ...")
			+ _S("\n----------------------------------------");
		PRINT_MESSAGE(msg);
		for (int i=0, sz=(int)symb_params_.size(); i<sz; i++)
		{
			msg = make_string(symb_params_[i].first);
			PRINT_MESSAGE(msg);
		}
		msg =
				_S("\n----------------------------------------")
			+ _S("\n  atoms readed ...")
			+ _S("\n----------------------------------------");
		PRINT_MESSAGE(msg);
		for (int i=0, sz=(int)atom_params_.size(); i<sz; i++)
		{
			msg = make_string(atom_params_[i].first);
			PRINT_MESSAGE(msg);
		}
		msg =
				_S("\n----------------------------------------")
			+ _S("\n  bonds readed ...")
			+ _S("\n----------------------------------------");
		PRINT_MESSAGE(msg);
		for (int i=0, sz=(int)bond_params_.size(); i<sz; i++)
		{
			msg = make_string(bond_params_[i].first);
			PRINT_MESSAGE(msg);
		}
		msg =
				_S("\n----------------------------------------")
			+ _S("\n  angles readed ...")
			+ _S("\n----------------------------------------");
		PRINT_MESSAGE(msg);
		for (int i=0, sz=(int)angl_params_.size(); i<sz; i++)
		{
			msg = make_string(angl_params_[i].first);
			PRINT_MESSAGE(msg);
		}
		msg =
				_S("\n----------------------------------------")
			+ _S("\n  dihedral readed ...")
			+ _S("\n----------------------------------------");
		PRINT_MESSAGE(msg);
		for (int i=0, sz=(int)dihe_params_.size(); i<sz; i++)
		{
			msg = make_string(dihe_params_[i].first);
			PRINT_MESSAGE(msg);
		}
		msg =
				_S("\n----------------------------------------")
			+ _S("\n  equivalent params generated ...")
			+ _S("\n----------------------------------------");
		PRINT_MESSAGE(msg);
		for (int i=0, sz=(int)equi_params_.size(); i<sz; i++)
		{
			msg = make_string(equi_params_[i].first);
			PRINT_MESSAGE(msg);
		}
		msg =
				_S("\n----------------------------------------")
			+ _S("\n  equi types params generated ...")
			+ _S("\n----------------------------------------");
		PRINT_MESSAGE(msg);

		for (int i=0, sz=(int)eqty_params_.size(); i<sz; i++)
		{
			std::string msg = make_string(eqty_params_[i].first) + _S(" : ");
			for (int i__=0; i__<eqty_params_[i].second.n; i__++)
			{
				msg += make_string(eqty_params_[i].second.eqv[i__]) + _S(" ");
			}
			PRINT_MESSAGE(msg);
		}
	}

	#undef TEMPLATE_HEADER
	#undef TEMPLATE_ARG

	#undef USES_EQUIVALENTS
}
#endif
