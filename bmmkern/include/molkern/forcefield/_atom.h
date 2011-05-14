#ifndef _ATOM__F9ED1116_DBE3_5136_EADB_F745B2F10100__H
#define _ATOM__F9ED1116_DBE3_5136_EADB_F745B2F10100__H

#include "molkern/__moldefs.h"
#include "molkern/forcefield/_nuclear.h"
#include "molkern/forcefield/_atomdata.h"
#include "molkern/forcefield/_forcefield.h"
#include "molkern/forcefield/_forcefield_amber.h"
#include "molkern/forcefield/_residome.h"
#include "molkern/forcefield/_residome_amber.h"

namespace molkern
{
	using namespace prgkern;

	/**
	* @brief convert atom into string
	*/
	template <typename _Atom>
	inline std::string make_string(const _Atom &atom, const Atomdata_ &atomdata)
	{
		//ATOM   5233 HD21 LEU C 336      76.946  41.633  52.823  1.00  0.00
		//0....5...10...15...20...25...30...35...40...45...50...55...60...65...
		char line[256];
		std::string pdb_name(make_string(atomdata.pdb_name));
		std::string residue(make_string(atomdata.residue));
		const char *format3 = "%5d  %-3s %-4s%c%4d    %8.3f%8.3f%8.3f%8.3f";
		const char *format4 = "%5d %-4s %-4s%c%4d    %8.3f%8.3f%8.3f%8.3f";
		const char *format = format3;
		if (pdb_name.length() == 4) format = format4;
		::sprintf(line, format, atomdata.sid, pdb_name.c_str(),
			residue.c_str(), atomdata.chain, atomdata.res_seq,
			(float)atom.X[0], (float)atom.X[1], (float)atom.X[2],
			(float)ExtCharge(atom.charge));
		std::string name__ = make_string(atomdata.name);
		std::string type__ = make_string(atomdata.fftype);
		std::string msg = _S(line) + _S(" ") + name__ + _S(" ") + type__
			+ _S(" ") + itoa(atomdata.nbond);
		if (atomdata.nbond != 0)
		{
			msg += _S("[") + itoa(atomdata.nid[0]);
			for (unsigned i=1; i<atomdata.nbond; i++)
				msg += _S(",") + itoa(atomdata.nid[i]);
			msg += _S("]");
		}

		return msg;
	}

	/**
	*  Конвертирует атом в строку вывода для HIN файла. Данная сигнатура функции
	*  позволяет делать вывод как текстовой, так и бинарной информации, в отличии
	*  от сигнатуры ранних версий, использующих std::string.
	*/
	template <typename _Atom>
	inline const char *make_string(_I2T<FORMAT_HIN_>, char *s,
		const _Atom &atom, const Atomdata_ &atomdata)
	{
		std::string name__ = trim(make_string(atomdata.name));
		if (name__.length() == 0)
		{
			unsigned ndx = find_nuclear(_I2T<MASS_>(), atomdata.mass);
			name__ = make_string(nuclears[ndx].name);
		}
		std::string type__ = make_string(atomdata.fftype);

		// восстановим ссылки на sid (hin файла) из ссылок на atomdata индексы
		// операция обратная той, что сделана в archetype::build(_I2T<FORMAT_HIN_>
		// при этом не гарантируем точного восстановления номеров, хотя гарантируем
		// точность восстановления топологии связей
		::sprintf(s, "atom %5d - %2s %2s - %8.5f %10.5f %10.5f %10.5f  %1u ",
			atomdata.sid, name__.c_str(), type__.c_str(),
			(float)ExtCharge(atom.charge),
			(float)atom.X[0], (float)atom.X[1], (float)atom.X[2],
			atomdata.nbond);

		char s__[12];
		for (unsigned i=0; i<atomdata.nbond; i++)
		{
			sprintf(s__, " %d %c", atomdata.nid[i], atomdata.nvalency[i]);
			strcat(s, s__);
		}
		return s;
	}

	/**
	* @brief convert atom into FORMAT_MOL2_ string
	*/
	template <typename _Atom>
	inline const char *make_string(_I2T<FORMAT_MOL2_>, char *s,
		const _Atom &atom, const Atomdata_ &atomdata)
	{
		std::string name__ = make_string(atomdata.name);
		std::string type__ = make_string(atomdata.fftype);
		std::string residue__ = make_string(atomdata.residue);
		::sprintf(s, "  %5d %-4s     %9.4f %9.4f %9.4f %2s     %4d %-3s      %9.4f",
			atomdata.sid, name__.c_str(),
			(float)atom.X[0], (float)atom.X[1], (float)atom.X[2], type__.c_str(),
			atomdata.res_seq, residue__.c_str(), (float)ExtCharge(atom.charge));
		return s;
	}

	/**
	* @brief convert atom into FORMAT_PDB_ string
	*/
	template <typename _Atom>
	inline const char *make_string(_I2T<FORMAT_PDB_>, char *s,
		const _Atom &atom, const Atomdata_ &atomdata)
	{
		std::string pdb_name(make_string(atomdata.pdb_name));
		std::string residue(make_string(atomdata.residue));

		unsigned res_seq = atomdata.res_seq % 10000;
			// так как выделено всего 4 позиции, то режем начальные цифры

		char s__[128];
		if (pdb_name.length() == 0)
		{
			sprintf(s, "HETATM%5d %-4s UNK     1    ",
				atomdata.sid, make_string(atomdata.name).c_str());
		}
		else
		{
			sprintf(s, "ATOM  %5d ", atomdata.sid);
			if (pdb_name.length() != 4) sprintf(s__, " %-3s ", pdb_name.c_str());
			else sprintf(s__, "%-4s ", pdb_name.c_str());
			strcat(s, s__);
			sprintf(s__, "%-4s%c%4d    ", residue.c_str(), atomdata.chain, res_seq);
			strcat(s, s__);
		}
		sprintf(s__,"%8.3f%8.3f%8.3f", (float)atom.X[0], (float)atom.X[1], (float)atom.X[2]);
		strcat(s, s__);

		return s;
	}

	/**
	*  Формирует сжатую строку байт (18) требуемого формата.
	* @param s место записи данных (должно быть не менее 18 байт)
	* @param atomdata данные по атому
	* @param atom данные по атому
	*/
	template <typename _Atom>
	inline const char *make_string(_I2T<FORMAT_BMM_>, char *s,
		const _Atom &atom, const Atomdata_ &atomdata)
	{
		DECLARE_AS_RESTORING_ROUND_MODE
		_MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);

		unsigned32_t zip_atom = 0;
		unsigned32_t count = 0;
		for (unsigned i=0; i<atomdata.nbond; i++)
		{
			if (atomdata.rnid[i] > 0) // сохраняем только ссылки вперед
			{
				zip_atom = (zip_atom << 5) + (atomdata.rnid[i] & 0x0000001F);
				count++;
			}
		}
		zip_atom = (zip_atom << 3) + (count & 0x00000007); // сохраняем число связей

		std::pair<const fstring*, const fstring*> p = std::equal_range(&fftypes[0],
			&fftypes[0] + fftypes_len, atomdata.fftype);
		assert(_NE(p.first, p.second));
		unsigned n = p.first - &fftypes[0];

		zip_atom = (zip_atom << 8) + (n & 0x000000FF); // сохраняем идентификатор атома
		*(unsigned32_t *)s = zip_atom;


		*(signed16_t*)(s +  4) = (signed16_t)(round(ExtCharge(atomdata.charge) * BMM_QGRID_STEP));
		*(signed16_t*)(s +  6) = (signed16_t)(round(atom.X[0] * BMM_XGRID_STEP));
		*(signed16_t*)(s +  8) = (signed16_t)(round(atom.V[0] * BMM_VGRID_STEP));
		*(signed16_t*)(s + 10) = (signed16_t)(round(atom.X[1] * BMM_XGRID_STEP));
		*(signed16_t*)(s + 12) = (signed16_t)(round(atom.V[1] * BMM_VGRID_STEP));
		*(signed16_t*)(s + 14) = (signed16_t)(round(atom.X[2] * BMM_XGRID_STEP));
		*(signed16_t*)(s + 16) = (signed16_t)(round(atom.V[2] * BMM_VGRID_STEP));

		return s;
	}

}
#endif
