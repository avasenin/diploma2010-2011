#ifndef _RESIDUE_AMBER__F9ED1116_23AF_57d7_F2EA_EA45CD1D0B00__H
#define _RESIDUE_AMBER__F9ED1116_23AF_57d7_F2EA_EA45CD1D0B00__H

#include "molkern/__moldefs.h"
#include "molkern/forcefield/_residue.h"

namespace molkern
{

	using namespace prgkern;

	template <> struct Params<RESIDUE_, RESIDOME_AMBER_>
	{
		typedef std::string index_type; // residue name

		// !entry.ALA.unit.atoms table  str name  str type  int typex  int resx  int flags
		//   int seq  int elmnt  dbl chg
		// "N" "N" 0 1 131072 1 7 -0.415700
		typedef struct // used partially
		{
			std::string name;
			fstring type;
			int typex;  // don't used
			int resx;  // используется, чтобы найти связные цепи (молекулы) в растворе
			int flags;  // don't used
			int seq;    // don't used (number of atom in seq)
			int elmnt;  // don't used
			real_t chg; // amber(and in gromacs the same) charge
		} _Atom;

		// !entry.ALA.unit.rotamers table  int count {str axisname str axisname int N str atomname_1 .. str atomname_N}(count)
		// 1
		// "CA" "CB" 3 "HB1" "HB2" "HB3"
		typedef struct // used
		{
			std::pair<fstring, fstring> axis; // атомы оси
			int atomnum; // число вращаемых атомов
			std::vector<fstring> atom;  // вращаемые атомы
		} _Rotamer;

		// !entry.ALA.unit.atomspertinfo table  str pname  str ptype  int ptypex
		//   int pelmnt  dbl pchg
		// "N" "N" 0 -1 0.0
		typedef struct // don't used
		{
			fstring pname;
			fstring ptype;
			int ptypex;
			int pelmnt;
			real_t pchg;
		} _Atompertinfo;

		// !entry.ALA.unit.boundbox array dbl
		// -1.000000
		typedef struct
		{
			real_t flag;
			real_t angle;
			real_t translation[3];

			void clear()
			{
				flag = 0.;
				angle = 0.;
				translation[0] = 0.;
				translation[1] = 0.;
				translation[2] = 0.;
			}
		} _BoundBox;

		// !entry.ALA.unit.childsequence single int
		// 2
		typedef int childsequence_type; // don't used

		// !entry.ALA.unit.connect array int
		// 1
		typedef int _Connect; // used as array of 2 elements ([0] - prev, [1] - next)

		// !entry.ALA.unit.connectivity table  int atom1x  int atom2x  int flags
		//  1 2 1
		typedef struct { int atom1x, atom2x, flags; } _Bond;

		// !entry.ALA.unit.hierarchy table  str abovetype  int abovex  str belowtype
		//   int belowx
		// "U" 0 "R" 1
		typedef struct // don't used
		{
			std::string abovetype;
			int abovex;
			std::string belowtype;
			int belowx;
		} _Hierarchy;

		// !entry.ALA.unit.name single str
		//  "ALA"
		typedef std::string _Name; // used

		// !entry.ALA.unit.positions table  dbl x  dbl y  dbl z
		// 3.325770 1.547909 -1.607257E-06
		typedef struct { real_t x, y, z; } _Position;

		// !entry.ALA.unit.residueconnect table  int c1x  int c2x  int c3x  int c4x
		//   int c5x  int c6x
		// 1 9 0 0 0 0
		typedef struct { int c1x, c2x, c3x, c4x, c5x, c6x; } _Residueconnect; // don't used

		// !entry.ALA.unit.residues table  str name  int seq  int childseq
		//   int startatomx  str restype  int imagingx
		// "ALA" 1 11 1 "p" 0
		typedef struct // don't used
		{
			std::string name;
			int seq;
			int childseq;
			int startatomx;
			int restype;
			int imagingx;
		} _Residue;

		// !entry.ALA.unit.residuesPdbSequenceNumber array int
		//  0
		int ResiduesPdbSequenceNumber_; // число загруженных остатков

		// !entry.ALA.unit.solventcap array dbl
		// -1.000000
		typedef real_t _Solventcap; // don't used

		// !entry.ALA.unit.velocities table  dbl x  dbl y  dbl z
		//  0.0 0.0 0.0
		typedef struct { real_t x, y, z; } _Velocity; // don't used

		/// residue name (key of sort)
		_Name name;

		/// external bonds [0] - prev, [1] - next residues
		_Connect connect_[2];

		/// box
		_BoundBox boundbox_;

		/// atoms
		std::vector<_Atom> atoms;
		std::vector<_Rotamer> rotamers;
		std::vector<_Bond> bonds; // bonds
		std::vector<_Position> positions; // x,y,z

		Params() {}
		Params(index_type name__) : name(name__) {}
			// need for residome sorting & search functions

		bool has_connects() const
		{
			if (bonds.size() || connect_[0] || connect_[1]) return true;
			return false;
		}

		void clear()
		{
			atoms.clear();
			rotamers.clear();
			boundbox_.clear();
			connect_[0] = 0; connect_[1] = 0;
			bonds.clear();
			name = "";
			positions.clear();
		}

		bool load(std::ifstream &file);
		unsigned count(_I2T<ATOM_>) const { return atoms.size(); }

	};

	#define TEMPLATE_HEADER
	#define TEMPLATE_ARG     RESIDUE_, RESIDOME_AMBER_

	TEMPLATE_HEADER
	inline bool Params<TEMPLATE_ARG>
	::load(std::ifstream &file)
	{
		char line[MAX_STRING_LEN]; // field to read the string
		char buf1[16], buf2[16]; // tmp arrays
		clear();

		_Atom atom;
		_Bond bond;
		_Position position;
	#ifdef FILE_FORMAT_DEBUG
		{
			std::string msg = _S("[DEBUG] new residue reading start...");
			PRINT_MESSAGE(msg);
		}
	#endif

		// used that we are positioned on first atom line
		file.getline(line, MAX_STRING_LEN);

		while (true)
		{
			std::istringstream s(make_string(line));
			s >> buf1 >> buf2 >> atom.typex >> atom.resx >> atom.flags >> atom.seq >> atom.elmnt >> atom.chg;
			atom.chg = IntCharge(atom.chg);

			atom.name = trim(_S(buf1), " \"");
			atom.type = trim(_S(buf2), " \"");
			atoms.push_back(atom);
			file.getline(line, MAX_STRING_LEN);
			if (::strncmp(line, "!entry", 6) == 0) break;
		}

		{
			std::string ctrl(line);
			if (ctrl.find("unit.rotamers", 7) != std::string::npos)
			{
				_Rotamer rotamer;
				while (true)
				{
					file.getline(line, MAX_STRING_LEN);
					if (::strncmp(line, "!entry", 6) == 0) break;

					std::istringstream is(line);
					is >> buf1 >> buf2 >> rotamer.atomnum;
					rotamer.axis.first = trim(_S(buf1), " \"");
					rotamer.axis.second = trim(_S(buf2), " \"");
					for (int i=0; i<rotamer.atomnum; i++)
					{
						is >> buf1;
						rotamer.atom.push_back(fstring(trim(_S(buf1), " \"")));
					}
					rotamers.push_back(rotamer);
					rotamer.atom.clear();
				}
			}
		}

	#ifdef FILE_FORMAT_DEBUG
		{
			std::string ctrl(line);
			std::string::size_type pos = ctrl.find("unit.atomspertinfo", 7);
				// skip "!entry.                           ^^^
			if (pos == std::string::npos)
			{
				std::string msg = _S("bad file format\n");
				msg += _S("    must be :") + _S("*.unit.atomspertinfo.*\n");
				msg += _S("    here is :") + _S(line);
				PRINT_ERR(msg);
			}
		}
	#endif
		int skip_count = (int)atoms.size() + 1;
		for (int i=0; i<skip_count; i++) file.getline(line, MAX_STRING_LEN);
			// skip atomspertinfo table + !..boundbox header
	#ifdef FILE_FORMAT_DEBUG
		{
			std::string ctrl(line);
			std::string::size_type pos = ctrl.find("unit.boundbox", 7);
			if (pos == std::string::npos)
			{
				std::string msg = _S("bad file format\n");
				msg += _S("    must be :") + _S("*.unit.boundbox.*\n");
				msg += _S("    here is :") + _S(line);
				PRINT_ERR(msg);
			}
		}
	#endif
		{
			file.getline(line, MAX_STRING_LEN);
			{
				std::istringstream s(make_string(line));
				s >> boundbox_.flag;
			}
			file.getline(line, MAX_STRING_LEN);
			{
				std::istringstream s(make_string(line));
				s >> boundbox_.angle;
			}
			file.getline(line, MAX_STRING_LEN);
			{
				std::istringstream s(make_string(line));
				s >> boundbox_.translation[0];
			}
			file.getline(line, MAX_STRING_LEN);
			{
				std::istringstream s(make_string(line));
				s >> boundbox_.translation[1];
			}
			file.getline(line, MAX_STRING_LEN);
			{
				std::istringstream s(make_string(line));
				s >> boundbox_.translation[2];
			}
			file.getline(line, MAX_STRING_LEN);
		}
	#ifdef FILE_FORMAT_DEBUG
		{
			std::string ctrl(line);
			std::string::size_type pos = ctrl.find("unit.childsequence", 7);
			if (pos == std::string::npos)
			{
				std::string msg = _S("bad file format\n");
				msg += _S("    must be :") + _S("*.unit.childsequence.*\n");
				msg += _S("    here is :") + _S(line);
				PRINT_ERR(msg);
			}
		}
	#endif
		for (int i=0; i<2; i++) file.getline(line, MAX_STRING_LEN);
			// skip childsequence(1) + !..connect header
	#ifdef FILE_FORMAT_DEBUG
		{
			std::string ctrl(line);
			std::string::size_type pos = ctrl.find("unit.connect", 7);
			if (pos == std::string::npos)
			{
				std::string msg = _S("bad file format\n");
				msg += _S("    must be :") + _S("*.unit.connect.*\n");
				msg += _S("    here is :") + _S(line);
				PRINT_ERR(msg);
			}
		}
	#endif
		file.getline(line, MAX_STRING_LEN);
		::sscanf(line, "%d", &connect_[0]);
		file.getline(line, MAX_STRING_LEN);
		::sscanf(line, "%d", &connect_[1]);
		file.getline(line, MAX_STRING_LEN); // !..connectivity header
		std::string ctrl(line);
		std::string::size_type pos = ctrl.find("connectivity", 13);
			// skip "!entry.*.unit.                                      ^^^
		if (pos != std::string::npos) // reads connectivity block
		{
			file.getline(line, MAX_STRING_LEN);
			while (true)
			{
				::sscanf(line, "%d %d %d", &bond.atom1x, &bond.atom2x, &bond.flags);
				bonds.push_back(bond);
				file.getline(line, MAX_STRING_LEN);
				if (::strncmp(line, "!entry", 6) == 0) break;
			}
		}
	#ifdef FILE_FORMAT_DEBUG
		{
			std::string ctrl(line);
			std::string::size_type pos = ctrl.find("unit.hierarchy", 7);
			if (pos == std::string::npos)
			{
				std::string msg = _S("bad file format\n");
				msg += _S("    must be :") + _S("*.unit.hierarchy.*\n");
				msg += _S("    here is :") + _S(line);
				PRINT_ERR(msg);
			}
		}
	#endif
		file.getline(line, MAX_STRING_LEN);
		while (::strncmp(&line[0], "!entry", 6) != 0)
			file.getline(line, MAX_STRING_LEN);
			// skip hierarchy table
	#ifdef FILE_FORMAT_DEBUG
		{
			std::string ctrl(line);
			std::string::size_type pos = ctrl.find("unit.name", 7);
			if (pos == std::string::npos)
			{
				std::string msg = _S("bad file format\n");
				msg += _S("    must be :") + _S("*.unit.name.*\n");
				msg += _S("    here is :") + _S(line);
				PRINT_ERR(msg);
			}
		}
	#endif
		file.getline(line, MAX_STRING_LEN);
	#ifdef FILE_FORMAT_DEBUG
		{
			std::string msg = _S("[DEBUG] residue name is :") + _S(line);
			PRINT_MESSAGE(msg);
		}
	#endif
		name = trim(_S(line), " \"");
		to_upper(name); // все остатки имеют имена в верхнем регистре (AMBER поле)

		file.getline(line, MAX_STRING_LEN);
			// skip !..positions table header
	#ifdef FILE_FORMAT_DEBUG
		{
			std::string ctrl(line);
			std::string::size_type pos = ctrl.find("unit.positions", 7);
			if (pos == std::string::npos)
			{
				std::string msg = _S("bad file format\n");
				msg += _S("    must be :") + _S("*.unit.positions.*\n");
				msg += _S("    here is :") + _S(line);
				PRINT_ERR(msg);
			}
		}
	#endif
		{
			file.getline(line, MAX_STRING_LEN);
			while (true)
			{
				std::istringstream s(make_string(line));
				s >> position.x >> position.y >> position.z;
				positions.push_back(position);
				file.getline(line, MAX_STRING_LEN);
				if (::strncmp(line, "!entry", 6) == 0) break;
			}
		}
	#ifdef FILE_FORMAT_DEBUG
		{
			std::string ctrl(line);
			std::string::size_type pos = ctrl.find("unit.residueconnect", 7);
			if (pos == std::string::npos)
			{
				std::string msg = _S("bad file format\n");
				msg += _S("    must be :") + _S("*.unit.residueconnect.*\n");
				msg += _S("    here is :") + _S(line);
				PRINT_ERR(msg);
			}
		}
	#endif
		file.getline(line, MAX_STRING_LEN);
		while (::strncmp(&line[0], "!entry", 6) != 0)
			file.getline(line, MAX_STRING_LEN);
			// skip residueconnect table
	#ifdef FILE_FORMAT_DEBUG
		{
			std::string ctrl(line);
			std::string::size_type pos = ctrl.find("unit.residues", 7);
			if (pos == std::string::npos)
			{
				std::string msg = _S("bad file format\n");
				msg += _S("    must be :") + _S("*.unit.residues.*\n");
				msg += _S("    here is :") + _S(line);
				PRINT_ERR(msg);
			}
		}
	#endif
		file.getline(line, MAX_STRING_LEN);
		while (::strncmp(&line[0], "!entry", 6) != 0)
			file.getline(line, MAX_STRING_LEN);
			// skip residues table
	#ifdef FILE_FORMAT_DEBUG
		{
			std::string ctrl(line);
			std::string::size_type pos = ctrl.find("unit.residuesPdbSequenceNumber", 7);
			if (pos == std::string::npos)
			{
				std::string msg = _S("bad file format\n");
				msg += _S("    must be :") + _S("*.unit.residuesPdbSequenceNumber.*\n");
				msg += _S("    here is :") + _S(line);
				PRINT_ERR(msg);
			}
		}
	#endif
		ResiduesPdbSequenceNumber_ = 0;

		file.getline(line, MAX_STRING_LEN);
		while (::strncmp(&line[0], "!entry", 6) != 0)
		{
			if (::strlen(&line[0]) != 0) ++ResiduesPdbSequenceNumber_;
			file.getline(line, MAX_STRING_LEN);
		}
			// skip residueconnect table
	#ifdef FILE_FORMAT_DEBUG
		{
			std::string ctrl(line);
			std::string::size_type pos = ctrl.find("unit.solventcap", 7);

			if (pos == std::string::npos)
			{
				std::string msg = _S("bad file format\n");
				msg += _S("    must be :") + _S("*.unit.solventcap.*\n");
				msg += _S("    here is :") + _S(line);
				PRINT_ERR(msg);
			}
		}
	#endif
		file.getline(line, MAX_STRING_LEN);
		while (::strncmp(&line[0], "!entry", 6) != 0)
			file.getline(line, MAX_STRING_LEN);
			// skip residueconnect table
	#ifdef FILE_FORMAT_DEBUG
		{
			std::string ctrl(line);
			std::string::size_type pos = ctrl.find("unit.velocities", 7);
			if (pos == std::string::npos)
			{
				std::string msg = _S("bad file format\n");
				msg += _S("    must be :") + _S("*.unit.velocities.*\n");
				msg += _S("    here is :") + _S(line);
				PRINT_ERR(msg);
			}
		}
	#endif
		skip_count = (int)atoms.size() + 1;
			// skip velocities(sz) + next residue header(1)
		for (int i=0; i<skip_count; i++) file.ignore(MAX_STRING_LEN, '\n');
	#ifdef FILE_FORMAT_DEBUG
		{
			std::string msg = _S("[DEBUG] residue reading finish...\n");
			PRINT_MESSAGE(msg);
		}
	#endif
		return true;
	}

	#undef TEMPLATE_HEADER
	#undef TEMPLATE_ARG

	INLINE bool operator<(const Params<RESIDUE_, RESIDOME_AMBER_> &a1,
		const Params<RESIDUE_, RESIDOME_AMBER_> &a2)
	{ return a1.name < a2.name; }

	INLINE bool operator==(const Params<RESIDUE_, RESIDOME_AMBER_> &a1,
		const Params<RESIDUE_, RESIDOME_AMBER_> &a2)
	{ return a1.name == a2.name; }

	/**
	* AMBER residue
	*/
	template <> class Residue_<RESIDOME_AMBER_>
	: public basic_residue_<RESIDOME_AMBER_>
	{
		typedef basic_residue_<RESIDOME_AMBER_>    _Base;
		typedef Params<RESIDUE_, RESIDOME_AMBER_>  _Param;

	public:

		typedef std::string index_type; // residue name

		Residue_() {}
		void init(const _Param &param);
	};

	#define TEMPLATE_HEADER
	#define TEMPLATE_ARG     RESIDOME_AMBER_

	TEMPLATE_HEADER
	inline void Residue_<TEMPLATE_ARG>
	::init(const _Param &param)
	{
		if (param.bonds.size() || param.connect_[0] || param.connect_[1])
			_Base::has_connects_ = true;

		if (param.connect_[0] || param.connect_[1])
			_Base::has_external_connects_ = true;

		_Base::name = param.name;
		_Base::ext_contacts_[0] = param.connect_[0] - 1;
		_Base::ext_contacts_[1] = param.connect_[1] - 1;
			// convert to C++ convention

		_Base::box_ = box_t(0., 0., 0.,
			param.boundbox_.translation[0],
			param.boundbox_.translation[1],
			param.boundbox_.translation[2]
		);
		_Base::box_.moveto(vector_t(0., 0., 0.));
			// сдвигаем, поскольку координаты ящика и атомов в *.lib файле расходятся

		__ResidueAtomdata tmp;
		for (int i=0, sz=(int)param.atoms.size(); i<sz; i++)
		{
			const _Param::_Atom &atom = param.atoms[i];
			tmp.name   = atom.name;
			tmp.fftype = atom.type;
			tmp.resx   = atom.resx;
			tmp.charge = atom.chg;
			tmp.X0[0] = param.positions[i].x;
			tmp.X0[1] = param.positions[i].y;
			tmp.X0[2] = param.positions[i].z;
			tmp.na = 0, tmp.nh = 0;
			for (int i__=0, sz__=(int)param.bonds.size(); i__<sz__; i__++)
			{
				const _Param::_Bond &bond = param.bonds[i__];
				int atom1x__ = bond.atom1x - 1;
				int atom2x__ = bond.atom2x - 1;
					// convert to C++ convention
				if (i != atom1x__ && i != atom2x__) continue;
				int atomx__ = (i == atom1x__) ? atom2x__ : atom1x__;
				std::string name = make_string(param.atoms[atomx__].name);
				if (name[0] == 'H') { tmp.nhid[tmp.nh] = atomx__; tmp.nh++; }
				else { tmp.nid[tmp.na] = atomx__; tmp.na++; }
			}
			_Base::value_type value(_Base::index_type(atom.name), tmp);
			_Base::residue_atomdata_.push_back(value);
		}

		unsigned sz = param.rotamers.size();
		_Base::rotamers_.resize(sz);
		for (unsigned i=0; i<sz; i++)
		{
			const _Param::_Rotamer &rotamer = param.rotamers[i];
			_Base::rotamers_[i].axis = rotamer.axis;
			_Base::rotamers_[i].atom.assign(rotamer.atom.begin(), rotamer.atom.end());
		}

		_Base::chain_count_ = (unsigned) param.ResiduesPdbSequenceNumber_;
	}

	#undef TEMPLATE_HEADER
	#undef TEMPLATE_ARG

}
#endif
