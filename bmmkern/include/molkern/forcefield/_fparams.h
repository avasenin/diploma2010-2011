#ifndef _FPARAMS__F9ED1116_B3CD_5e56_EF64_DB4362A40D00__H
#define _FPARAMS__F9ED1116_B3CD_5e56_EF64_DB4362A40D00__H

#include "molkern/__moldefs.h"
#include "molkern/forcefield/_nuclear.h"

namespace molkern
{

	/**
	* @brief access key to database of potential params
	*/
	typedef index_<1, fstring> ffindex1;
	typedef sym_index_<2, fstring> ffindex2;
	typedef sym_index_<3, fstring> ffindex3;
	typedef sym_index_<4, fstring> ffindex4;

	/**
	* @brief template of any params structure
	* @param PARAMS_TYPE - params type
	*/
	template <int PARAMS_TYPE, int FORCEFIELD_TYPE> struct Params;
		// по умалчиванию все вещественные поля real_t, так как внешние параметры
		// не могут быть точно определены. Кроме того, это уменьшает память
		// необходимую для хранения параметров, уменьшая промахи кэша.

	template <int X, int Y> inline bool operator<(
		const std::pair<typename Params<X, Y>::index_type, Params<X, Y> > &a1,
		const std::pair<typename Params<X, Y>::index_type, Params<X, Y> > &a2)
	{ return a1.first < a2.first; }

	template <int X, int Y> inline bool operator==(
		const std::pair<typename Params<X, Y>::index_type, Params<X, Y> > &a1,
		const std::pair<typename Params<X, Y>::index_type, Params<X, Y> > &a2)
	{ return a1.first == a2.first; }

	/**
	* @brief params of AMBER forcefield
	*/
	template <> struct Params<NUCLEAR_, FORCEFIELD_AMBER_>
	{
		typedef ffindex1 index_type;
		fstring name; // atom name (to efficient search)
		real_t mass; // effective mass of atom
		real_t polar; // atomic polarizability (какая размерность ????)
	};

	/**
	* @brief params of AMBER forcefield
	*/
	template <> struct Params<ATOM_, FORCEFIELD_AMBER_>
	{
		typedef ffindex1 index_type;
		union { real_t radius; real_t sigma; }; // effective radius of atom
		real_t eps; // used to calculate c6, c12 VdW coefficients
	};

	/**
	* @brief params of AMBER forcefield
	*/
	template <> struct Params<EQUI_, FORCEFIELD_AMBER_>
	{
		typedef ffindex1 index_type;
		fstring eqv[MAX_EQUIVALENT];
		int n; // full count of eqv
			// no copy operator
			// I think that using default copy operator will be faster due to
			//   cycle escaping
	};

	/**
	* @brief params of AMBER forcefield
	*/
	template <> struct Params<BOND_, FORCEFIELD_AMBER_>
	{
		typedef ffindex2 index_type;
		real_t ke;  // E''(q0_)
		real_t q0;  // equlibrium distance
	};

	/**
	* @brief params of AMBER forcefield
	*/
	template <> struct Params<ANGLE_, FORCEFIELD_AMBER_>
	{
		typedef ffindex3 index_type;
		real_t ke;  // E''(q0_)
		real_t q0;  // equlibrium distance
	};

	/**
	* @brief params of AMBER forcefield
	*/
	template <> struct Params<TORSION_, FORCEFIELD_AMBER_>
	{
		typedef ffindex4 index_type;
		static const int maxf = 4; // max count of subfuctions
		real_t v[maxf]; // potential
		real_t phi[maxf]; // angle phase

		int proper; // proper (1) & unproper(0) dihedral
		int nb; // number of paths in dihedral
		int curf; // current subfunctions
		int n[maxf]; // index subfunction
	};

	#define NDX_PARAM_PAIR(name) std::pair<Params<name, FORCEFIELD_AMBER_>::index_type, \
		Params<name, FORCEFIELD_AMBER_> >

	/**
	* @brief extract params from file stream
	*/
	inline int extract_object(NDX_PARAM_PAIR(NUCLEAR_) &p, std::ifstream *src)
	{
		typedef Params<NUCLEAR_, FORCEFIELD_AMBER_>  _Param;
		typedef _Param::index_type          _Index;

		char buffer[MAX_STRING_LEN]; // field to read the string
		char buff__[MAX_STRING_LEN]; // field to convert the string
		src->getline(buffer, MAX_STRING_LEN);

		if (isspace(_S(buffer))) return CODE_ERROR;
			// empty line is delimiter for AMBER format file

		if (buffer[0] == ';') return CODE_SUCCESS;
			// line which started with ';' is comment
		::strncpy(buff__, &buffer[0], 2); buff__[2] = 0; // atom name

		_Index &index = p.first;
		_Param &param = p.second;

		_S tmp = _S(buff__, 2);
		rtrim(tmp);
		const char *pname = tmp.c_str();
		bool gaff_forcefield = ::islower(pname[0]) ? true : false;
		index = _Index(pname);

		::strncpy(buff__, &buffer[3], 10); buff__[10] = 0; // mass
			// inconsistence here ^^^  with format description
		param.mass = IntMass((real_t)atof(buff__));

		int nuclear = find_nuclear(_I2T<MASS_>(), param.mass);
		if (gaff_forcefield) { param.name = to_lower(nuclears[nuclear].name); }
		else { param.name = to_upper(nuclears[nuclear].name); }
			// convert to according forcefield type

		::strncpy(buff__, &buffer[13], 10); buff__[10] = 0;
			// The atomic polarizability for each atom (in A**3)
		param.polar = (real_t)atof(buff__);

		return CODE_SUCCESS;
	}

	/**
	* @brief extract params from file stream
	*/
	inline int extract_object(NDX_PARAM_PAIR(ATOM_) &p, std::ifstream *src)
	{
		typedef Params<ATOM_, FORCEFIELD_AMBER_>  _Param;
		typedef _Param::index_type       _Index;

		char buffer[MAX_STRING_LEN]; // field to read the string
		char buff__[MAX_STRING_LEN]; // field to convert the string

		src->getline(buffer, MAX_STRING_LEN);
		if (isspace(_S(buffer))) return CODE_ERROR;
			// empty line is delimiter for AMBER format file

		if (buffer[0] == ';') return CODE_SUCCESS;
			// line which started with ';' is comment

		::strncpy(buff__, &buffer[2], 2); buff__[2] = 0; // atom name

		_Index &index = p.first;
		_Param &param = p.second;
		index = _Index(rtrim(_S(buff__, 2)));

		::strncpy(buff__, &buffer[14], 6); buff__[6] = 0; // VdW radius (angstroms)
		param.sigma = (real_t)atof(buff__);

	#ifdef USE_HARM_AVERAGE_OF_SIGMA
		if (fabs(param.sigma) < 1.) param.sigma = 1.;
			// приравниваем 1, чтобы при геометрическом усреднении происходил выбор
			// максимального sigma
	#endif

		::strncpy(buff__, &buffer[22], 9); buff__[9] = 0; // well depth (kcal/mol)
		param.eps = (real_t)atof(buff__) * CAL2J;

		return CODE_SUCCESS;
	}

	/**
	* @brief extract params from file stream
	*/
	inline int extract_object(NDX_PARAM_PAIR(EQUI_) &p, std::ifstream *src)
	{
		typedef Params<EQUI_, FORCEFIELD_AMBER_>  _Param;
		typedef _Param::index_type       _Index;

		char buffer[MAX_STRING_LEN]; // field to read the string
		char buff__[MAX_ATOM_NAME + 1]; // field to convert the string

		src->getline(buffer, MAX_STRING_LEN);

		if (isspace(_S(buffer))) return CODE_ERROR;
			// empty line is delimiter for AMBER format file

		if (buffer[0] == ';') return CODE_SUCCESS;
			// line which started with ';' is comment

		_Index &index = p.first;
		_Param &param = p.second;
		param.n = 0; // number of equivalents

		unsigned symbs = std::min((unsigned)MAX_EQUIVALENT, (unsigned)(::strlen(buffer)+3)/4);
		for (unsigned i=0; i<symbs; i++)
		{
			::strncpy(buff__, &buffer[4 * i], 2); buff__[2] = 0;
			fstring a = fstring(trim(_S(buff__, 2)));
			if (a == fstring("")) break;
			param.eqv[i] = a;
			param.n++;
		}

		fstring ndx = param.eqv[0];
		index = _Index(ndx); // make index

		return CODE_SUCCESS;
	}

	/**
	* @brief extract params from file stream
	*/
	inline int extract_object(NDX_PARAM_PAIR(BOND_) &p, std::ifstream *src)
	{
		typedef Params<BOND_, FORCEFIELD_AMBER_>  _Param;
		typedef _Param::index_type       _Index;

		char buffer[MAX_STRING_LEN]; // field to read the string
		char buff__[MAX_STRING_LEN]; // field to convert the string

		src->getline(buffer, MAX_STRING_LEN);
		if (isspace(_S(buffer))) return CODE_ERROR;
			// empty line is delimiter for AMBER format file

		if (buffer[0] == ';') return CODE_SUCCESS;
			// line which started with ';' is comment

		::strncpy(buff__, &buffer[0], 2); buff__[2] = 0; // atom name
		fstring a1(rtrim(_S(buff__, 2)));
		::strncpy(buff__, &buffer[3], 2); buff__[2] = 0; // atom name

		fstring a2(rtrim(_S(buff__, 2)));
		_Index &index = p.first;
		_Param &param = p.second;

		index = _Index(a1, a2);

		::strncpy(buff__, &buffer[5], 10); buff__[10] = 0;
			// harmonic force constant (kcal/(mol*A^2))
		param.ke = (real_t)atof(buff__) * CAL2J;

		::strncpy(buff__, &buffer[15], 10); buff__[10] = 0;
			// equilibrium bond length for the above bond in angstroms
		param.q0 = (real_t)atof(buff__);

		return CODE_SUCCESS;
	}

	/**
	* @brief extract params from file stream
	*/
	inline int extract_object(NDX_PARAM_PAIR(ANGLE_) &p, std::ifstream *src)
	{
		typedef Params<ANGLE_, FORCEFIELD_AMBER_>  _Param;
		typedef _Param::index_type        _Index;

		char buffer[MAX_STRING_LEN]; // field to read the string
		char buff__[MAX_STRING_LEN]; // field to convert the string

		src->getline(buffer, MAX_STRING_LEN);

		if (isspace(_S(buffer))) return CODE_ERROR;
			// empty line is delimiter for AMBER format file

		if (buffer[0] == ';') return CODE_SUCCESS;
			// line which started with ';' is comment

		::strncpy(buff__, &buffer[0], 2); buff__[2] = 0; // atom name
		fstring a1(rtrim(_S(buff__, 2)));

		::strncpy(buff__, &buffer[3], 2); buff__[2] = 0; // atom name
		fstring a2(rtrim(_S(buff__, 2)));

		::strncpy(buff__, &buffer[6], 2); buff__[2] = 0; // atom name
		fstring a3(rtrim(_S(buff__, 2)));

		_Index &index = p.first;
		_Param &param = p.second;

		index = _Index(a1, a2, a3);

		::strncpy(buff__, &buffer[8], 10); buff__[10] = 0;
			// harmonic force constant (kcal/(mol*rad^2))
		param.ke = (real_t)atof(buff__) * CAL2J;

		::strncpy(buff__, &buffer[18], 10); buff__[10] = 0;
			// equilibrium bond length for the above bond in rad
		param.q0 = (real_t) (atof(buff__) * DEG2RAD);

		if (param.q0 < 0.) param.q0 = (real_t) 0.;
		if (param.q0 > 180. * DEG2RAD) param.q0 = (real_t) (180. * DEG2RAD);
			// Нулевых и отрицательных углов быть не может,
			// также не могут быть углы > 180.

		return CODE_SUCCESS;
	}

	/**
	* @brief extract params from file stream
	*/
	inline int extract_object(NDX_PARAM_PAIR(TORSION_) &p, std::ifstream *src)
	{
		typedef Params<TORSION_, FORCEFIELD_AMBER_>  _Param;
		typedef _Param::index_type          _Index;

		char buffer[MAX_STRING_LEN]; // field to read the string
		char buff__[MAX_STRING_LEN]; // field to convert the string
		char ctrl__[9]; // control line
		src->getline(buffer, MAX_STRING_LEN);
		if (isspace(_S(buffer))) return CODE_ERROR;
			// empty line is delimiter for AMBER format file
		if (buffer[0] == ';') return CODE_SUCCESS;
			// line which started with ';' is comment

		::strncpy(ctrl__, &buffer[0], 9); buff__[9] = 0; // control line
		::strncpy(buff__, &buffer[0], 2); buff__[2] = 0; // atom name
		fstring a1(rtrim(_S(buff__, 2)));
		::strncpy(buff__, &buffer[3], 2); buff__[2] = 0; // atom name
		fstring a2(rtrim(_S(buff__, 2)));
		::strncpy(buff__, &buffer[6], 2); buff__[2] = 0; // atom name
		fstring a3(rtrim(_S(buff__, 2)));
		::strncpy(buff__, &buffer[9], 2); buff__[2] = 0; // atom name
		fstring a4(rtrim(_S(buff__, 2)));
		_Index &index = p.first;
		_Param &param = p.second;
		index = _Index(a1, a2, a3, a4);

		param.curf = 0; // init the counter
		::strncpy(buff__, &buffer[11], 4); buff__[4] = 0;
			// The factor by which the torsional barrier is divided.
		param.nb = (int)floor(atof(buff__));
		::strncpy(buff__, &buffer[15], 15); buff__[15] = 0;
			// The barrier height divided by a factor of 2.
		param.v[0] = (real_t)atof(buff__) * CAL2J / abs(param.nb);
		::strncpy(buff__, &buffer[30], 15); buff__[15] = 0;
			// The phase shift angle in the torsional function (degrees).
		param.phi[0] = (real_t)(atof(buff__) * DEG2RAD);
		::strncpy(buff__, &buffer[45], 15); buff__[15] = 0;
			// The periodicity of the torsional barrier.

		if (param.phi[0] < 0.) param.phi[0] = 0.;
		if (param.phi[0] > 180. * DEG2RAD) param.phi[0] = 180. * DEG2RAD;
			// Нулевых и отрицательных углов быть не может,
			// также не могут быть углы > 180.

		param.n[0] = atoi(buff__);

		while (param.n[param.curf] < 0) // read the next lines
		{
			param.n[param.curf] *= -1; // the only positive value is used
			param.curf++; // count of potential parts
			src->getline(buffer, MAX_STRING_LEN);
			if (isspace(_S(buffer))) return CODE_ERROR;
				// empty line is delimiter for AMBER format file
			if (buffer[0] == ';') continue;
				// line which started with ';' is comment
			if (::strncmp(buffer, ctrl__, 9) != 0)
			{
				std::string msg = _S("bad file format");
				PRINT_ERR(msg);
			}
			::strncpy(buff__, &buffer[11], 4); buff__[4] = 0;
				// The factor by which the torsional barrier is divided.
			param.nb = (int)floor(atof(buff__));
			::strncpy(buff__, &buffer[15], 15); buff__[15] = 0;
				// The barrier height divided by a factor of 2.
			param.v[param.curf] = (real_t)atof(buff__) * CAL2J / abs(param.nb);
			::strncpy(buff__, &buffer[30], 15); buff__[15] = 0;
				// The phase shift angle in the torsional function (degrees).
			param.phi[param.curf] = (real_t)(atof(buff__) * DEG2RAD);

			if (param.phi[0] < 0.) param.phi[param.curf] = 0.;
			if (param.phi[0] > 180. * DEG2RAD) param.phi[param.curf] = 180. * DEG2RAD;
				// Нулевых и отрицательных углов быть не может,
				// также не могут быть углы > 180.

			::strncpy(buff__, &buffer[45], 15); buff__[15] = 0;
				// The periodicity of the torsional barrier.
			param.n[param.curf] = (int)floor(atof(buff__));
		}
		param.curf++; // must be full number of cur functions, not index

		return CODE_SUCCESS;
	}
	#undef NDX_PARAM_PAIR
}
#endif
