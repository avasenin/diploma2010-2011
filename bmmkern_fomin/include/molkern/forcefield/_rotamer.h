#ifndef _ROTAMER__F9ED1116_DBE3_5136_EADB_F745B2F10101__H
#define _ROTAMER__F9ED1116_DBE3_5136_EADB_F745B2F10101__H

#include <boost/tuple/tuple.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "molkern/__moldefs.h"
#include "molkern/forcefield/_nuclear.h"
#include "molkern/forcefield/_atomdata.h"
#include "molkern/complex/_geom_tool.h"

namespace molkern
{
	using namespace prgkern;

	class Chain_
	{
		unsigned *ndx_; // начало массива идентификаторов атомов цепи
		unsigned ndx_count_; // количество атомов в цепи

	public:

		void set(_I2T<ATOM_>, unsigned *ndx) { ndx_ = ndx; }
		const unsigned *get(_I2T<ATOM_>, unsigned i=0) const { return ndx_ + i; }

		unsigned count(_I2T<ATOM_>) const { return ndx_count_; }
		unsigned &count(_I2T<ATOM_>) { return ndx_count_; }
	};

	struct RotamerConnect_
	{
		unsigned rotamer_from;
		unsigned stick_from;
		unsigned rotamer_to;
		unsigned stick_to;

		RotamerConnect_() {}
		RotamerConnect_(unsigned r1, unsigned s1, unsigned r2, unsigned s2)
		: rotamer_from(r1), stick_from(s1), rotamer_to(r2), stick_to(s2) {}
	};

	class Rotamer_
	{
		typedef index_<3> _Stick;

		std::vector<_Stick> sticks_; // индексы атомов, определяющих угол вращения ротамера
			// ВНИМАНИЕ! В индексы (как вспомогательный атом) могут попасть атомы не из тела ротамера,
			// а связанные со стиком и разделяемые другим ротамером. Поскольку такие атомы вращаются,
			// то нет необходимости специально ограничивать индекс только телом.
		unsigned *ndx_; // индексы атомов ротамера (инициализируется в Archetype_)
			// совместно вращаемые атомы (атомы попавшие на стики от других ротамеров) НЕ включены
			// в список атомов ротамера
		unsigned ndx_count_;

	public:

		typedef _Stick  stick_type;
		DEFINE_VECTOR_ACCESS_FUNCTION(STICK_, _Stick, sticks_);

		/**
		 *  Делает стик из номеров переданных атомов. Все стики ВЫХОДЯТ из ротамера.
		 * @param i1 - номер атома внешнего ротамера
		 * @param i2 - номер атома, принадлежащего ротамеру и лежащего на стике
		 * @param i3 - номер атома ротамера, связанного со стиком
		 */
		void insert_stick(unsigned i1, unsigned i2, unsigned i3)
		{ sticks_.push_back(stick_type(i1, i2, i3)); }

		unsigned get(_I2T<ROTAMER_EXT_ATOM_>, unsigned n) const { return sticks_[n][0]; }
		unsigned get(_I2T<ROTAMER_INT_ATOM_>, unsigned n) const { return sticks_[n][1]; }
		unsigned get(_I2T<ROTAMER_ADD_ATOM_>, unsigned n) const { return sticks_[n][2]; }

		void set(_I2T<ATOM_>, unsigned *ndx) { ndx_ = ndx; }
		const unsigned *get(_I2T<ATOM_>, unsigned i=0) const { return ndx_ + i; }

		unsigned count(_I2T<ATOM_>) const { return ndx_count_; }
		unsigned &count(_I2T<ATOM_>) { return ndx_count_; }

		/**
		 * отладочная функция
		 */
		void print() const
		{
			unsigned ssz = sticks_.size();
			unsigned asz = count(ATOM);

			_S msg = make_string("The rotamer includes %d sticks and %d atoms.\n", ssz, asz);

			for (unsigned i=0; i<ssz; i++)
				msg += make_string("stick %d : %s\n", i, make_string(sticks_[i]).c_str());

			msg += make_string("atoms : ");
			for (unsigned i=0; i<asz; i++) msg += make_string("%d ", ndx_[i]);
			msg += _S("\n");

			PRINT_MESSAGE(msg);
		}
	};

	/**
	 *
	 * @param A
	 * @param B
	 * @param C
	 * @param D
	 * @return
	 */
	INLINE real_t rotamer_angle(const vector_t &A, const vector_t &B,
		const vector_t &C, const vector_t &D)
	{
		vector_t BA = A - B;
		vector_t CD = D - C;
		vector_t BC = C - B;

		vector_t U = vector_product(BA, BC);
		vector_t T = vector_product(CD, BC);

		// определяем угол между плоскостями U и T (угол между ротамерами)
		real_t phi = get_angle(U, T);

		// определяем знак угла между плоскостями U и T (угол между ротамерами)
		vector_t UT = vector_product(U, T);
		if (scalar_product(UT, BC) < 0) phi = -phi;

		return phi;
	}

	template <typename A>
	INLINE real_t rotamer_angle(const A *atoms, const Rotamer_ &rotamer1,
		const Rotamer_ &rotamer2, unsigned s1, unsigned s2)
	{
		unsigned ni = rotamer1.get(ROTAMER_INT_ATOM, s1);
		unsigned na = rotamer1.get(ROTAMER_ADD_ATOM, s1);
		unsigned ni__ = rotamer2.get(ROTAMER_INT_ATOM, s2);
		unsigned na__ = rotamer2.get(ROTAMER_ADD_ATOM, s2);

		return rotamer_angle(atoms.X[na], atoms.X[ni], atoms.X[ni__], atoms.X[na__]);
	}

	/**
	 * Приклеивает второй ротамер к первому, точнее меняет координаты его атомов.
	 * @param a1 - объекты, имеющие координаты (поля X), конечные
	 * @param a2 - объекты, имеющие координаты (поля X), стартовые
	 * @param rotamers - массив ротамеров
	 * @param s1 - индекс стика константного ротамера
	 * @param s2 - индекс стика приклеиваемого ротамера
	 * @param alpha - угол поворота ротамера, при условии, то угол между добавочными атомами = 0
	 */
	template <typename A1, typename A2>
	inline void rotamer_paste(A1 *a1, const A2 *a2, const Rotamer_ &rotamer1,
		const Rotamer_ &rotamer2, unsigned s1, unsigned s2, float alpha)
	{
		unsigned atom_count2 = rotamer2.count(ATOM); // число атомов во втором ротамере
		const unsigned *ndx2 = rotamer2.get(ATOM); // начало массива идентификаторов атомов второго ротамера
		unsigned stick_count2 = rotamer2.count(STICK);

		//----------------------------------------------------------------------------
		//  сдвигаем второй ротамер в пространстве, совмещая с первым по одной точке
		//----------------------------------------------------------------------------
		unsigned ni = rotamer1.get(ROTAMER_INT_ATOM, s1);
		unsigned ne = rotamer1.get(ROTAMER_EXT_ATOM, s1);
		unsigned na = rotamer1.get(ROTAMER_ADD_ATOM, s1);
		unsigned ni__ = rotamer2.get(ROTAMER_INT_ATOM, s2);
		unsigned ne__ = rotamer2.get(ROTAMER_EXT_ATOM, s2);
		unsigned na__ = rotamer2.get(ROTAMER_ADD_ATOM, s2);

		vector_t XI = a1[ni].X;
		vector_t XE = a1[ne].X; // (!) копируется заранее, так X меняется при сдвиге 2-го ротамера
		vector_t XA = a1[na].X;

		// смещаем тело ротамера
		vector_t R = XI - a2[ne__].X; // вектор смещения атомов
		for (unsigned i=0; i<atom_count2; i++)
		{
			unsigned n = ndx2[i];
			a1[n].X = a2[n].X + R; // перетаскиваем атомы из хранилища 2 в новые координаты в хранилище 1
		}
		for (unsigned i=0; i<stick_count2; i++)
		{
			unsigned n = rotamer2.get(ROTAMER_EXT_ATOM, i);
			a1[n].X = a2[n].X + R;
		}

		//----------------------------------------------------------------------------
		//                 совмещаем оси стиков ротамеров
		//----------------------------------------------------------------------------
		vector_t R1 = XE - XI; // вектор оси первого ротамера
		R1.normalize();

		vector_t R1__ = a1[ni__].X - a1[ne__].X; // обращенный вектор оси второго ротамера
		R1__.normalize();

		vector_t W = vector_product(R1__, R1); // проверка угла вращения на нуль
		if (W.length2() > sqr(std::numeric_limits<real_t>::epsilon()))   ////?
		{
			real_t phi = get_angle(R1__, R1);
			Rotator<AXIS_ROTATOR_, real_t> rotator(W, phi, XI);

			// смещаем тело ротамера
			for (unsigned i=0; i<atom_count2; i++)
			{
				unsigned n = ndx2[i];
				a1[n].X = rotator(a1[n].X);
					// мы уже перетащили все атомы из хранилища 2 в новые координаты в хранилище 1
			}
			for (unsigned i=0; i<stick_count2; i++)
			{
				unsigned n = rotamer2.get(ROTAMER_EXT_ATOM, i);
				a1[n].X = rotator(a1[n].X);
			}
		}

		//----------------------------------------------------------------------------
		//                    вращаем вдоль оси ротамера
		//----------------------------------------------------------------------------
		real_t beta = rotamer_angle(a1[na].X, a1[ni].X, a1[ni__].X, a1[na__].X);
			// текущий угол между ротамерами

		real_t angle = alpha - beta; // добавка к уже существующему повороту
		{
			Rotator<AXIS_ROTATOR_, real_t> rotator(R1, angle, XI);

			// смещаем тело ротамера
			for (unsigned i=0; i<atom_count2; i++)
			{
				unsigned n = ndx2[i];
				a1[n].X = rotator(a1[n].X);
			}
			for (unsigned i=0; i<stick_count2; i++)
			{
				unsigned n = rotamer2.get(ROTAMER_EXT_ATOM, i);
				a1[n].X = rotator(a1[n].X);
			}
		}

	}

	/**
	 * @param W - единичный вектор оси вращения
	 * @param O - точка оси
	 * @param F - вектор силы
	 * @param X - положение точки приложения силы
	 * @return значение момента силы
	 */
	INLINE real_t atom_moment(const vector_t &W, const vector_t &O, const vector_t &F,
		const vector_t &X)
	{
		vector_t A = vector_product(W, X - O);
		return scalar_product(F, A);
	}

	template <typename A>
	INLINE real_t rotamer_moment(const A *atoms, const Rotamer_ *rotamer,
		const vector_t &W, const vector_t &O)
	{
		real_t moment = (real_t)0.;
		unsigned atom_count = rotamer->count(ATOM); // число атомов в ротамере
		const unsigned *ndx = rotamer->get(ATOM); // начало массива идентификаторов атомов ротамера

		for (unsigned i=0; i<atom_count; i++)
		{
			unsigned n = ndx[i];
//			real_t mmm = atom_moment(W, O, atoms[n].F, atoms[n].X);
//			std::cout <<"i=" << i << " n=" << n << " mmm=" << mmm << std::endl;
//			std::cout << "F=" << make_string(atoms[n].F) << "X=" << make_string(atoms[n].X) << std::endl;
//			std::cout << "XX=" << make_string(atoms[n].F * 0.3 + atoms[n].X) << std::endl;
//			std::cout << "W=" << make_string(W) << "O=" << make_string(O) << std::endl;

			moment += atom_moment(W, O, atoms[n].F, atoms[n].X);
		}
		return moment;
	}

}
#endif
