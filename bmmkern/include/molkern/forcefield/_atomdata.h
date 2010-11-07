#ifndef _ATOMDATA__F9ED1116_3BFB_5a4c_6CDD_F74523C50B00__H
#define _ATOMDATA__F9ED1116_3BFB_5a4c_6CDD_F74523C50B00__H

#include "molkern/__moldefs.h"

namespace molkern
{
	using namespace prgkern;

	/// полная упорядоченная (для fstring) таблица имен типов для AMBER & GAFF,
	/// необходимая для упаковки данных в bmm формате
	const fstring fftypes[] =
	{
		fstring("C"),  fstring("F"),  fstring("H"),  fstring("I"),  fstring("K"),
		fstring("N"),  fstring("O"),  fstring("P"),  fstring("S"),  fstring("c"),
		fstring("f"),  fstring("i"),  fstring("n"),  fstring("o"),  fstring("s"),
		fstring("C*"), fstring("N*"), fstring("C0"), fstring("H1"), fstring("c1"),
		fstring("h1"), fstring("n1"), fstring("H2"), fstring("N2"), fstring("O2"),
		fstring("c2"), fstring("h2"), fstring("n2"), fstring("p2"), fstring("s2"),
		fstring("H3"), fstring("N3"), fstring("c3"), fstring("h3"), fstring("n3"),
		fstring("p3"), fstring("H4"), fstring("h4"), fstring("n4"), fstring("p4"),
		fstring("s4"), fstring("H5"), fstring("h5"), fstring("p5"), fstring("s6"),
		fstring("CA"), fstring("HA"), fstring("NA"), fstring("CB"), fstring("IB"),
		fstring("NB"), fstring("CC"), fstring("HC"), fstring("NC"), fstring("CD"),
		fstring("FE"), fstring("MG"), fstring("OH"), fstring("SH"), fstring("CK"),
		fstring("CM"), fstring("IM"), fstring("CN"), fstring("HO"), fstring("HP"),
		fstring("IP"), fstring("LP"), fstring("CQ"), fstring("CR"), fstring("HS"),
		fstring("OS"), fstring("CT"), fstring("NT"), fstring("CU"), fstring("CV"),
		fstring("CW"), fstring("HW"), fstring("OW"), fstring("CY"), fstring("NY"),
		fstring("CZ"), fstring("HZ"), fstring("Na"), fstring("ca"), fstring("ha"),
		fstring("na"), fstring("Rb"), fstring("nb"), fstring("pb"), fstring("cc"),
		fstring("hc"), fstring("nc"), fstring("pc"), fstring("cd"), fstring("nd"),
		fstring("pd"), fstring("ce"), fstring("ne"), fstring("pe"), fstring("cf"),
		fstring("nf"), fstring("pf"), fstring("cg"), fstring("ch"), fstring("nh"),
		fstring("oh"), fstring("sh"), fstring("Li"), fstring("Cl"), fstring("cl"),
		fstring("Zn"), fstring("hn"), fstring("ho"), fstring("no"), fstring("cp"),
		fstring("hp"), fstring("cq"), fstring("Br"), fstring("br"), fstring("Cs"),
		fstring("hs"), fstring("os"), fstring("ss"), fstring("cu"), fstring("cv"),
		fstring("hw"), fstring("ow"), fstring("cx"), fstring("hx"), fstring("px"),
		fstring("sx"), fstring("cy"), fstring("py"), fstring("sy")
	};
	const unsigned fftypes_len = sizeof(fftypes) / sizeof(fstring);

	/**
	*  Извлекает информацию о наличии контакта из подходящей таблицы контакта.
	* @param na,nb индексы атомов
	* @param ta,tb таблицы контактов соответствующих атомов
	* @return true при наличии контакта (от 1-1 до 1-4)
	*/
	INLINE bool has_contact(unsigned na, unsigned nb, unsigned ta, unsigned tb)
	{
	#ifdef USE_FIXED_ATOMS
		if ((na & nb & 0x80000000)) return false;
		na &= 0x7FFFFFFF; nb &= 0x7FFFFFFF;
	#endif

		unsigned diff = (na < nb) ? nb - na : na - nb; // разность индексов,
			// так как таблица построена для относительных индексов

		if (diff < 32) // таблицы не содержат разности более чем 31
		{
			unsigned t = (na < nb) ? ta : tb;
				// Пользуемся таблицей контактов атома с минимальным индексом,
				// так как таблицы контактов содержат информацию только "вперед".
			return (bool)((t >> diff) & 0x00000001);
		}
		return false;
	}

	/**
	 * Выдает порядок валентности атома по его символу. Для неопределенных символом порядок равен 1.
	 * @param c - символ валентности
	 * @return порядок валентности
	 */
	inline float valency_order(char c)
	{
		switch (c)
		{
		case 's' : return 1.f;
		case 'd' : return 2.f;
		case 't' : return 3.f;
		case 'a' : return 1.5f;
		}
		return 1.f;
	}

	/**
	* @brief uncorrected external data storage
	* @note vs. atom (internal data storage corrected for fast calculations)
	*/
	struct Atomdata_
	{
		//______________________________SHARE FIELDS_______________________________

		/// atom forcefield type name
		fstring fftype;
			// atom forcefield name is defined by forcefield
			// for example, atom carbon has name CT (AMBER) or ca (GAFF)

		/// atom name (C, O, H, .. etc.)
		fstring name;

		/// Self IDent or serial number (for HIN, PDB, ... format)
		int sid;
			// this field may be string but the PDB, HIN .. suppose the integer type

		/// cartesian coordinates of atom
		vector_t X; // [angstrom]
			// these values are loaded from file and aren't changed
			// these can be used for restore start atom positions

		/// скорости атомов (они храняться в BMM формате)
		vector_t V; // [angstrom / ps]

		/// charge of atom
		real_t charge;

		/// van der Waals radius of atom
		union { real_t radius; real_t sigma; };
		real_t eps;

		/// mass of atom
		real_t mass;

		/// заряд ядра
		unsigned nuclear;

		//..............................SHARE FIELDS...............................

		//__________________________HIN FORMAT FIELDS______________________________
		/// max bond (neighbour) of atoms = 6
		static const int max_bond = MAX_ATOM_BOND;

		/// number of bonds (neighbour)
		unsigned nbond; // used for HIN format only

		int nid[max_bond]; /// абсолютный идентификатор соседа
		int rnid[max_bond]; /// относительный идентификатор соседа
			// относительный идентификатор для удобства работы с bmm файлами

		/// neighbour atom bond valency
		char nvalency[max_bond];

		//..........................HIN FORMAT FIELDS..............................

		//__________________________PDB FORMAT FIELDS______________________________

		/// PDB atom name (for example, "HD12")
		fstring pdb_name;
			// "HD12" isn't atom name or type name
			// it's atom identifier in given residue

		/// identifier of alternative atom location (' ', 'A', 'B')
		char altloc;

		/// residue or some atom group name ("HIS", "Zn+", "UNK" etc.)
		std::string residue;

		/// identifier of protein chain ("A", "B")
		char chain;

		/// residue sequence number
		int res_seq;

		//..........................PDB FORMAT FIELDS..............................

		unsigned_t connect_data; // упакованная таблица коннектов от 1-1 до 1-4
		unsigned_t insert_data; // тип вставки

		Atomdata_() : fftype("  "), name("  "), sid(0), charge(0.), radius(0.),
			mass(0.), nuclear(0), nbond(0), pdb_name("   "), altloc(' '), residue("UNK"),
			chain(' '), res_seq(1), connect_data(1), insert_data(0)
				// вставили контакт с самим собой для действительных атомов
		{}

		/**
		*  Проверка на то, что атом является водородом. Необходимость проверки связана
		*  с тем, что во многих алгоритмах можно удалить водороды, как атомов увеличивающих
		*  расчетное время, но не играющих важной роли.
		*/
		bool is_hydrogen() const
		{
			return nuclear == 1 // быстрая проверка для сформированной записи
			|| toupper(make_string(pdb_name)[0]) == 'H' // проверка для несформированной записи pdb
			|| toupper(make_string(name)[0]) == 'H'; // проверка для несформированной записи hin
		}

		bool is_pseudo() const
		{
			return toupper(make_string(pdb_name)[0]) == 'E' // проверка для несформированной записи pdb
			|| toupper(make_string(name)[0]) == 'E'; // проверка для несформированной записи hin
		}
	};

	/**
	* @brief convert atomdata into string
	* @note conversion is used for testing
	*/
	inline std::string make_string(const Atomdata_ &atomdata)
	{
		std::string nids = _S("[");
		if (atomdata.nbond > 0) nids += _S(itoa(atomdata.rnid[0]));

		for (unsigned i=1; i<atomdata.nbond; i++)
		{
			nids += _S(",") + _S(itoa(atomdata.rnid[i]));
			}
		nids += _S("]");

		char line[256];
		std::string pdb_name__(make_string(atomdata.pdb_name));
		if (pdb_name__ == _S("")) pdb_name__ = make_string(atomdata.name);
			// output the most detailed type if possible
		std::string fftype__ = make_string(atomdata.fftype);

		std::string residue__(make_string(atomdata.residue));
		if (residue__ == _S("")) residue__ = _S("UNK");

		const char *format3 = "%4d  %-3s %-4s%-4s%c%4d%8.3f%8.3f%8.3f%6.2f  %1d ";
		const char *format4 = "%4d %-4s %-4s%-4s%c%4d%8.3f%8.3f%8.3f%6.2f  %1d ";

		const char *format = format3;
		if (pdb_name__.length() == 4) format = format4;

		::sprintf(line, format, atomdata.sid, pdb_name__.c_str(), fftype__.c_str(),
			residue__.c_str(), atomdata.chain, atomdata.res_seq,
			(real_t)atomdata.X[0], (real_t)atomdata.X[1], (real_t)atomdata.X[2],
			(real_t)atomdata.charge, atomdata.nbond);

		return _S(line) + nids;
	}

	/**
	* @brief convert atomdata into FORMAT_HIN_ string
	*/
	inline std::string make_string(_I2T<FORMAT_HIN_>, const Atomdata_ &atomdata)
	{
		char line[256];
		// atom 1 - C CT - 0.1801 18.35900 26.86900 52.95500 2 2 s 9 s
		const char *format = "atom %4d - %2s %2s - %8.5f %10.5f %10.5f %10.5f %1d ";
		std::string name__ = make_string(atomdata.name);
		std::string type__ = make_string(atomdata.fftype);

		// восстановим ссылки на sid (hin файла) из ссылок на atomdata индексы
		// операция обратная той, что сделана в archetype::build(_I2T<FORMAT_HIN_>
		// при этом не гарантируем точного восстановления номеров, хотя гарантируем
		// точность восстановления топологии связей

		::sprintf(line, format, atomdata.sid + 1, name__.c_str(), type__.c_str(),
			(real_t)atomdata.charge, (real_t)atomdata.X[0], (real_t)atomdata.X[1],
			(real_t)atomdata.X[2], atomdata.nbond);

		std::string msg = _S(line);
		for (unsigned i=0; i<atomdata.nbond; i++)
			msg += itoa(atomdata.nid[i] + 1) + _S(" ") + ctoa(atomdata.nvalency[i]) + _S(" ");

		return msg;
	}

	/**
	* @brief convert atomdata into FORMAT_PDB_ string
	*/
	inline std::string make_string(_I2T<FORMAT_PDB_>, const Atomdata_ &atomdata)
	{
		//ATOM   5233 HD21 LEU C 336      76.946  41.633  52.823  1.00  0.00
		//0....5...10...15...20...25...30...35...40...45...50...55...60...65...
		std::string pdb_name(make_string(atomdata.pdb_name));
		std::string residue(make_string(atomdata.residue));

		if (pdb_name.length() == 0) pdb_name = make_string(atomdata.name);
		std::string outline = "ATOM  ";

		if (atomdata.sid) outline += make_string("%5d ", atomdata.sid);
		else outline += make_string("      ");

		if (pdb_name.length() != 4)
			outline += make_string(" %-3s ", pdb_name.c_str());
		else outline += make_string("%-4s ", pdb_name.c_str());

		unsigned res_seq = atomdata.res_seq % 10000;
			// так как выделено всего 4 позиции, то режем начальные цифры

		outline += make_string("%-4s", residue.c_str());
		outline += make_string("%c", atomdata.chain);
		outline += make_string("%4d    ", res_seq);
		outline += make_string("%8.3f", (real_t)atomdata.X[0]);
		outline += make_string("%8.3f", (real_t)atomdata.X[1]);
		outline += make_string("%8.3f", (real_t)atomdata.X[2]);

		rtrim(outline);
		return outline;
	}

	/**
	* @brief convert atomdata into FORMAT_BMM_ string
	*/
	inline std::string make_string(_I2T<FORMAT_BMM_>, const Atomdata_ &)
	{
		return _S("");
	}

	/**
	* @brief extracts atomdata from FORMAT_HIN_ file stream
	*/
	inline int extract_object(_I2T<FORMAT_HIN_>, Atomdata_ &atomdata, std::ifstream *src)
	{
		// atom 1 - C CT - 0.1801 18.35900 26.86900 52.95500 2 2 s 9 s

		char skip[256]; // buffer to skip unused data
		char name[MAX_ATOM_NAME];

		*src >> skip; if (::strncmp(skip, "atom", 4) != 0)
			return CODE_ERROR; // test file line

		*src >> atomdata.sid >> skip;
		*src >> name; atomdata.name = name;
		*src >> name; atomdata.fftype = name;
		*src >> skip >> atomdata.charge
			>> atomdata.X[0] >> atomdata.X[1] >> atomdata.X[2];
		*src >> atomdata.nbond;
		for (unsigned i=0; i<atomdata.nbond; i++)
			*src >> atomdata.nid[i] >> atomdata.nvalency[i];

		return CODE_SUCCESS;
	}

	/**
	* @brief extracts atomdata from FORMAT_HIN_ file stream
	*/
	inline int extract_object(_I2T<FORMAT_MOL2_>, Atomdata_ &atomdata, std::ifstream *src)
	{
		char skip[256]; // buffer to skip unused data
		char name[MAX_ATOM_NAME];

		// sample of format MOL2 line
		//      1 C1          2.3610    0.7390    0.7030 c         1 <0>         0.9501

		*src >> skip; if (skip[0] == '@')
			return CODE_ERROR; // test file line

		atomdata.sid = ::atoi(skip);
		*src >> name; atomdata.name = name;
		*src >> atomdata.X[0] >> atomdata.X[1] >> atomdata.X[2];
		*src >> name; atomdata.fftype = name;
		*src >> atomdata.res_seq;
		*src >> name; atomdata.residue = name;
		*src >> atomdata.charge;

		return CODE_SUCCESS;
	}

	/**
	* @brief extracts atomdata from FORMAT_PDB_ file stream
	* @note FORMAT OF PDB SECTION : ATOM (!)FORTRAN COLUMNS NUMBERING
	* -------------------------------------------------------------------------
	* COLUMNS    DATA TYPE       FIELD         DEFINITION
	* -------------------------------------------------------------------------
	*  1 -  6    Record name     "ATOM  "
	*  7 - 11    Integer         serial        Atom serial number.
	* 13 - 16    Atom            name          Atom name.
	[!!! MAYBE ERROR !!!] atom name 14-16 otherwise RASMOL has errors
	[ only long atom names used 13-16 positions]

	* 17         Character       altLoc        Alternate location indicator.
	* 18 - 20    Residue name    resName       Residue name.
	[!!! ERROR !!!] residue name include 4 positions

	* 22         Character       chainID       Chain identifier.
	* 23 - 26    Integer         resSeq        Residue sequence number.
	* 27         AChar           iCode         Code for insertion of residues.
	* 31 - 38    Real(8.3)       x             Orthogonal coord for X in A
	* 39 - 46    Real(8.3)       y             Orthogonal coord for Y in A
	* 47 - 54    Real(8.3)       z             Orthogonal coord for Z in A
	* 55 - 60    Real(6.2)       occupancy     Occupancy.
	* 61 - 66    Real(6.2)       tempFactor    Temperature factor.
	* 73 - 76    LString(4)      segID         Segment ident, left-justified.
	* 77 - 78    LString(2)      element       Element symb, right-justified.
	* 79 - 80    LString(2)      charge        Charge on the atom.
	*--------------------------------------------------------------------------
	*/


	/**
	*  Извлечь данные из строки PDB файла по формату поля ATOM.
	*  Предполагается, что контроль правильности строки сделан перед ее вызовом.
	*/
	inline void extract_object(_I2T<FORMAT_PDB_>, Atomdata_ &atomdata, char *line)
	{
		char skip[256]; // buffer to skip unused data

		::strncpy(skip, &line[6], 5); skip[5] = 0; // self ident
		atomdata.sid = std::atoi(skip);

		::strncpy(skip, &line[12], 4); skip[4] = 0; // name
		atomdata.pdb_name = trim(_S(skip, 4));

		atomdata.altloc = line[16]; // altloc

		::strncpy(skip, &line[17], 4); skip[4] = 0; // residue
		atomdata.residue = trim(_S(skip, 4));

		atomdata.chain = line[21]; // chain

		::strncpy(skip, &line[22], 4); skip[4] = 0; // sequence number
		atomdata.res_seq = ::atoi(skip);

		::strncpy(skip, &line[30], 8); skip[8] = 0; // x
		atomdata.X[0] = (real_t)atof(skip);

		::strncpy(skip, &line[38], 8); skip[8] = 0; // y
		atomdata.X[1] = (real_t)atof(skip);

		::strncpy(skip, &line[46], 8); skip[8] = 0; // z
		atomdata.X[2] = (real_t)atof(skip);
	}

	inline int extract_object(_I2T<FORMAT_PDB_>, Atomdata_ &atomdata, std::ifstream *src)
	{
		char line[256]; // буфер для чтения одной строки файла

		do {
			src->getline(line, 256);
			if (src->eof()) return CODE_EOF;
		} while ( !(::strncmp(line, "ATOM", 4) == 0 && ::isspace(line[26])) );
			// игнорируем HETATOM, iCode (текущие ограничения программы)

		extract_object(_I2T<FORMAT_PDB_>(), atomdata, line);
		return CODE_SUCCESS;
	}

	/**
	*  Считываем данные из файла, распаковывает бинарные данные.
	* @param src бинарный файл
	* @param atomdata атомныйе данные
	*/
	inline int extract_object(_I2T<FORMAT_BMM_>, Atomdata_ &atomdata, std::ifstream *src)
	{
		char buffer[BMM_RECORD_LEN];

		src->read(buffer, BMM_RECORD_LEN);
		if (src->eof()) return CODE_EOF;

		unsigned32_t atom = *(unsigned32_t *)buffer;
		atomdata.fftype = fftypes[atom & 0x000000FF];
		atom >>= 8;

		atomdata.nbond = atom & 0x00000007;
		atom >>= 3;
		for (unsigned i=0; i<atomdata.nbond; i++)
		{
			atomdata.rnid[i] = atom & 0x0000001F;
			atom >>= 5;
		}

		atomdata.charge = (real_t)(*(signed16_t*)&buffer[4]) / BMM_QGRID_STEP;
		atomdata.X[0] = (real_t)(*(signed16_t*)&buffer [6]) / BMM_XGRID_STEP;
		atomdata.V[0] = (real_t)(*(signed16_t*)&buffer [8]) / BMM_VGRID_STEP;
		atomdata.X[1] = (real_t)(*(signed16_t*)&buffer[10]) / BMM_XGRID_STEP;
		atomdata.V[1] = (real_t)(*(signed16_t*)&buffer[12]) / BMM_VGRID_STEP;
		atomdata.X[2] = (real_t)(*(signed16_t*)&buffer[14]) / BMM_XGRID_STEP;
		atomdata.V[2] = (real_t)(*(signed16_t*)&buffer[16]) / BMM_VGRID_STEP;

		return CODE_SUCCESS;
	}

	/**
	* @brief convert FORMAT_HIN_ input line into Atom<ATOM_>
	* @param param - forcefield index
	* @note throw exception in case of bad casting
	*/
	inline void make_object(_I2T<FORMAT_HIN_>, Atomdata_ &atomdata,
		const Forcefield_<FORCEFIELD_AMBER_> *forcefield,
		const Residome_<RESIDOME_AMBER_> *residome)
	{
		assert(_NE((void*)forcefield, (void*)NULL));
		assert(_NE((void*)residome, (void*)NULL));

		{
			typedef Params<ATOM_, FORCEFIELD_AMBER_> _Param;
			typedef _Param::index_type _Index;
			int nff = forcefield->get_data_ndx(_I2T<ATOM_>(), _Index(atomdata.fftype));
			if (nff == nill)
			{
				fstring eqv = forcefield->get(_I2T<EQUI_>(), _Index(atomdata.fftype));
					// get vdw equivalent type
				nff = forcefield->get_data_ndx(_I2T<ATOM_>(), _Index(eqv));
				if (nff == nill) { NO_EQUIVALENTS(_Index(eqv)); }
			}
			const _Param *param__ = forcefield->get(_I2T<ATOM_>(), nff);
			atomdata.sigma = param__->sigma;
			atomdata.eps = param__->eps;
		}
		{
			typedef Params<NUCLEAR_, FORCEFIELD_AMBER_> _Param;
			_Param param__;
			_Param::index_type index__(atomdata.fftype);
			forcefield->get_data(&param__, index__);
			atomdata.mass = param__.mass;

			atomdata.nuclear = find_nuclear(_I2T<MASS_>(), atomdata.mass);
			atomdata.name = nuclears[atomdata.nuclear].name;

		#ifdef USE_HEAVY_HYDROGENS
			if (atomdata.nuclear == 1) atomdata.mass *= HEAVY_MASS_FACTOR;
		#endif
		}
	}

	/**
	* @brief convert FORMAT_PDB_ input line into atomdata
	* @note throw exception in case of bad casting
	*/
	inline void make_object(_I2T<FORMAT_PDB_>, Atomdata_ &atomdata,
		const Forcefield_<FORCEFIELD_AMBER_> *forcefield,
		const Residome_<RESIDOME_AMBER_> *residome)
	{
		assert(_NE((void*)forcefield, (void*)NULL));
		assert(_NE((void*)residome, (void*)NULL));

		if (atomdata.pdb_name == fstring("OT1")) atomdata.pdb_name = fstring("O");
		if (atomdata.pdb_name == fstring("OT2")) atomdata.pdb_name = fstring("OXT");
			// convert pdb_name to AMBER agreements
		{
			typedef Residome_<RESIDOME_AMBER_>::residue_type _Residue;
			const _Residue *residue = residome->get_data(atomdata.residue);
			assert(_NE((void*)residue, (void*)0));

			int ndx = residue->get_data_ndx(atomdata.pdb_name);
			const __ResidueAtomdata *residue_atom = residue->get_data(ndx);

			atomdata.fftype = residue_atom->fftype;
			atomdata.charge = residue_atom->charge;
		}

		{
			typedef Params<ATOM_, FORCEFIELD_AMBER_> _Param;
			typedef _Param::index_type _Index;
			int nff = forcefield->get_data_ndx(_I2T<ATOM_>(), _Index(atomdata.fftype));
			if (nff == nill)
			{
				fstring eqv = forcefield->get(_I2T<EQUI_>(), _Index(atomdata.fftype));
					// get vdw equivalent type
				nff = forcefield->get_data_ndx(_I2T<ATOM_>(), _Index(eqv));
				if (nff == nill) { NO_EQUIVALENTS(_Index(eqv)); }
			}
			const _Param *param__ = forcefield->get(_I2T<ATOM_>(), nff);
			atomdata.sigma = param__->sigma;
			atomdata.eps = param__->eps;
		}
		{
			typedef Params<NUCLEAR_, FORCEFIELD_AMBER_> _Param;
			_Param param__;
			_Param::index_type index__(atomdata.fftype);
			forcefield->get_data(&param__, index__);
			atomdata.mass = param__.mass;

			atomdata.nuclear = find_nuclear(_I2T<MASS_>(), atomdata.mass);
			atomdata.name = nuclears[atomdata.nuclear].name;

		#ifdef USE_HEAVY_HYDROGENS
			if (atomdata.nuclear == 1) atomdata.mass *= HEAVY_MASS_FACTOR;
		#endif
		}
	}

	INLINE void make_object(_I2T<FORMAT_MOL2_>, Atomdata_ &atomdata,
		const Forcefield_<FORCEFIELD_AMBER_> *forcefield,
		const Residome_<RESIDOME_AMBER_> *residome)
	{
		make_object(_I2T<FORMAT_HIN_>(), atomdata, forcefield, residome);
	}

	INLINE void make_object(_I2T<FORMAT_BMM_>, Atomdata_ &atomdata,
		const Forcefield_<FORCEFIELD_AMBER_> *forcefield,
		const Residome_<RESIDOME_AMBER_> *residome)
	{
		make_object(_I2T<FORMAT_HIN_>(), atomdata, forcefield, residome);
	}

}
#endif
