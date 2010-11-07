#ifndef _ENSEMBLE__F9ED1116_EDB9_5e17_25FF_F745B15D0101__H
#define _ENSEMBLE__F9ED1116_EDB9_5e17_25FF_F745B15D0101__H

#include "molkern/__moldefs.h"

namespace molkern
{
	using namespace prgkern;

	/**
	*  Расчет текущей кинетической энергии объекта в предположении максвеловского
	*  распределения.
	* @param molecule указатель на молекулу или комплекс молекул
	* @return энергия [KJ/mol]
	*/
	template <typename _Molecule>
	inline _E(real_t) kinetic_energy(_Molecule *molecule)
	{
		typedef typename _Molecule::atom_type          _Atom;
		typedef array_iterator<_Atom, range_iterator>  _Iterator;

		_E(real_t) kenergy = 0.;
		unsigned atom_count = molecule->count(_I2T<ATOM_>());

		_Iterator it = molecule->make_iterator(range_iterator(0));
		_Iterator ite = molecule->make_iterator(range_iterator(atom_count));
		for (; it!=ite; ++it)
		{
			const _Atom &atom = *it;
			kenergy += atom->mass * scalar_product(atom.V, atom.V);
		}
		return (real_t)0.5 * kenergy / atom_count; // нормируем на одну частицу
	}

	/**===========================================================================
	*             Различные типы термодинамических ансамблей.
	*  Предполагается, что существует единый термостат, чтобы облегить код
	*  на уровне представления пользователя. Однако этот термостат имеет свойства,
	*  связанные с типом сбора статистики (NVE, NVT, NVP).
	============================================================================*/
	enum { NVE, NVT, NVP };
	template <int TYPE> class Ensemble_;

	/**
	*  Статистика NVE ансамбля должна хранит только полную энергию системы.
	*  Но для более комфортного использования мы будем хранить в нем, как
	*  потенциальную, так и кинетическую энергии. Накопления постоянно
	*  сбрасываются в консоль для контроля поведения системы.
	*  Статистика по числу частиц и объему системы (в текущей версии) не ведется,
	*  так как они не меняются.
	*/
	template <> class Ensemble_<NVE>
	{
		Average_<real_t> T_; // температура для разных точек фазового пространства
		real_t target_energy_; // константа данного ансамбля

	public:

		Ensemble_() {}
		Ensemble_(real_t target_energy, unsigned ns) : T_(ns), target_energy_(target_energy) {}

		template <typename _Molecule>
		void sample_statistics(_Molecule *molecule)
		{
			T_.push(get(_I2T<TEMPERATURE_>(), molecule));
		}

		/**
		*  Масштабирование скоростей указанных атомов согласно заданной полной энергии.
		*  Данный алгоритм позволяет поддерживть точное значение энергии для выбранной
		*  фазовой точки, которая идет в статистику. Алгоритм не требует параметра,
		*  связанного со временем релаксации (шагом набора статистики), поскольку
		*  масштабирование делается точно и может быть сделано для каждой точки.
		* @param molecule указатель на молекулу или комплекс молекул
		* @return коэффициент lambda
		*/
		template <typename _Molecule>
		std::string scaling(_Molecule *molecule)
		{
			typedef typename _Molecule::atom_type          _Atom;
			typedef array_iterator<_Atom, range_iterator>  _Iterator;

			real_t U = (real_t) molecule->U();
			real_t K = (real_t) kinetic_energy(molecule);
			real_t lambda = (real_t) sqrt((target_energy_ - U) / K);
				// разность энергий может быть покрыта только изменением кинетической
				// энергии, потому используем (E - <E>) / K

			unsigned atom_count = molecule->count(_I2T<ATOM_>());
			_Iterator it = molecule->make_iterator(range_iterator(0));
			_Iterator ite = molecule->make_iterator(range_iterator(atom_count));
			for (; it!=ite; ++it) (*it).V *= lambda;
			K *= sqr(lambda);

			return make_string("K = %12.5e  U = %12.5e  E = %12.5e  T = %6.2f",
				K, U, (float)(K + U), (float)T_());
		}
	};

	/**
	*  Статистика NVT ансамбля должна хранит только температуру.
	*/
	template <> class Ensemble_<NVT>
	{
		Average_<real_t> T_; // температура для разных точек фазового пространства
		real_t target_temperature_; // константа данного ансамбля

	public:

		Ensemble_() {}

		/**
		*  Расчет NVT статистики текущего состояния.
		* @param target_temperature целевая температура
		* @param ns максимальное хранимое число точек для построения статистики
		*/
		Ensemble_(real_t target_temperature, unsigned ns)
		: T_(ns), target_temperature_(target_temperature)
		{
		}

		/**
		*  Расчет NVT статистики текущего состояния.
		* @param molecule указатель на молекулу или комплекс молекул
		*/
		template <typename _Molecule>
		void sample_statistics(_Molecule *molecule)
		{
			T_.push(get(_I2T<TEMPERATURE_>(), molecule));
		}

		/**
		*  Масштабирование скоростей указанных атомов согласно заданной температуре.
		*  Для масштабирования используется метод Берендсена.
		* @param molecule указатель на молекулу или комплекс молекул
		* @param target_temperature целевая температура для масштабирования скоростей
		* @param time_step шаг интегрирования
		* @param relax_time время релаксации (подгоночный параметр метода)
		* @return коэффициент lambda
		*/
		template <typename _Molecule>
		std::string scaling(_Molecule *molecule)
		{
			typedef typename _Molecule::atom_type          _Atom;
			typedef array_iterator<_Atom, range_iterator>  _Iterator;

			real_t T__ = T_(); // усредним статистику по полному набору
			real_t T = get(_I2T<TEMPERATURE_>(), molecule);

			real_t lambda = (real_t) sqrt(1. + (target_temperature_ - T__) / T);
				// полностью покрываем разность кинетических энергий без релаксации

			unsigned atom_count = molecule->count(_I2T<ATOM_>());
			_Iterator it = molecule->make_iterator(range_iterator(0));
			_Iterator ite = molecule->make_iterator(range_iterator(atom_count));
			for (; it!=ite; ++it) (*it).V *= lambda;

			return make_string("T = %7.2f", (float)T__);
		}
	};

	/**
	*  Статистика NVP ансамбля должна хранит только давление.
	*/
	template<> class Ensemble_<NVP>
	{
		Average_<real_t> pressure_;
		real_t target_pressure_;

	public:

		Ensemble_() {}

		Ensemble_(real_t target_pressure, unsigned ns)
		: pressure_(ns), target_pressure_(target_pressure) {}

		/**
		*  Расчет текущего внутреннего вириала W, кинетической энергии K и давления.
		*  Заметим, что давление может быть определено только для системы с
		*  заданным объемом. Расчеты давления показывают, что оно скачет в диапозоне
		*  [-1e+4, 1e+4], то есть не может быть рассчитано точно.
		* @param box ячейка, в которой рассчитывается давление
		* @param molecule указатель на молекулу или комплекс молекул
		*/
		template <typename _Molecule>
		void sample_statistics(_Molecule *molecule)
		{
			pressure_.push(get(_I2T<PRESSURE_>(), molecule));
		}

		template <typename _Molecule>
		std::string scaling(_Molecule *molecule)
		{
			typedef typename _Molecule::atom_type          _Atom;
			typedef array_iterator<_Atom, range_iterator>  _Iterator;

			real_t P__ = pressure_(); // усредним статистику по полному набору
			real_t P = get(_I2T<PRESSURE_>(), molecule);

			real_t lambda = (real_t) sqrt(1. + (target_pressure_ - P__) / P);
				// полностью покрываем разность кинетических энергий без релаксации

			unsigned atom_count = molecule->count(_I2T<ATOM_>());
			_Iterator it = molecule->make_iterator(range_iterator(0));
			_Iterator ite = molecule->make_iterator(range_iterator(atom_count));
			for (; it!=ite; ++it) (*it).V *= lambda;

			return make_string("P = %6.2f", (float)P__);
		}

	};

}
#endif
