#ifndef _THERMOSTAT__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H
#define _THERMOSTAT__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H

#include "molkern/__moldefs.h"
#include "molkern/complex/_ensemble.h"

namespace molkern
{
	using namespace prgkern;

	/**
	* @brief максвеловское (нормальное) распределение по проекциям скоростей
	*/
	template <typename _Real=real_t>
	class maxwell_distribution_
	: public variate_generator_<normal_distribution_<_Real> >
	{
		typedef normal_distribution_<_Real>        _Distribution;
		typedef variate_generator_<_Distribution>  _Base;

		_Real temperature_;
		_Real mass_;

	public:

		/**
		* @brief конструктор максвеловского распределения по проекциям скоростей
		* @param mass масса атома (a.u.m.)
		* @param temperature температура (K)
		*/
		explicit maxwell_distribution_(_Real mass, _Real temperature)
		: _Base(_Distribution(sqrt(KT(temperature)/(mass)))),
			temperature_(temperature), mass_(mass) {}

		_Real temperature() const { return temperature_; }
		_Real mass() const { return mass_; }
	};

	typedef maxwell_distribution_<real_t>  maxwell_distribution;

	/**
	* @brief ансамбль генераторов максвеловских распределений для разных масс атомов
	*/
	template <typename _Real=real_t>
	class maxwell_distribution_ensemble_
	{
		typedef maxwell_distribution_<_Real>        _Distribution;
		typedef std::pair<unsigned, _Distribution>  _Pair;
		typedef std::map<unsigned, _Distribution>   _Ensemble;

		_Ensemble ensemble_; // набор генераторов (mass, distr)
		_Real temperature_; // температура распределений

	public:

		/**
		* @brief конструирование генераторов для атомов с зарядом ядра в диапозоне [1, MAX_NUCLEAR_INDEX)
		* @param temperature температура (K)
		*/
		maxwell_distribution_ensemble_(_Real temperature=300.) : temperature_(temperature) {}

		/**
		* Переустановка температуры генераторов.
		* @param temperature новая температура
		*/
		void set(_I2T<TEMPERATURE_>, _Real temperature)
		{
			temperature_ = temperature;
			ensemble_.clear();
		}

		/**
		* @brief генерация очередного случайного числа для заданного типа атома
		* @param nuclear заряд ядра атома в Менделеевской таблице
		*/
		_Real generate(_Real mass)
		{
			unsigned key = unsigned(mass);

			typename _Ensemble::iterator it = ensemble_.find(key);
			if (it == ensemble_.end())
			{
				ensemble_.insert(_Pair(key, _Distribution(mass, temperature_)));
			}
			return ensemble_[key]->operator()();
		}

		/**
		* @brief генерация набора случайных чисел для заданного типа атома
		* @param mass масса частицы
		* @param n число генераций
		* @param ranv[out] набор случайных чисел
		*/
		void generate(_Real mass, unsigned n, _Real *ranv)
		{
			unsigned key = unsigned(mass);

			typename _Ensemble::iterator it = ensemble_.find(key);
			if (it == ensemble_.end())
			{
				ensemble_.insert(_Pair(key, _Distribution(mass, temperature_)));
			}
			it = ensemble_.find(key);

			for (unsigned i=0; i<n; i++)
			{
				_Distribution &dist = (*it).second;
				ranv[i] = dist();
			}
		}
	};
	typedef maxwell_distribution_ensemble_<>  maxwell_distribution_ensemble;

	/**
	*  Установка начальных скоростей атомов объекта, согласно максвеловскому
	*  распределению в газе при заданной температуре. Здесь предполагается, что
	*  все атомы являются независимыми друг от друга. Данный способ установления
	*  скоростей является стандартным, но он дает полностью нескоррелированное
	*  распределение скоростей в молекуле (клеши в пространстве скоростей).
	*  Лучше пользоваться функциями которые устанавливают распределение для
	*  каждых групп степеней свободы отдельно.
	* @param molecule указатель на молекулу или комплекс молекул
	* @param temperature температура [K]
	*/
	template <typename _Molecule, typename _Iterator>
	INLINE void set(_I2T<TEMPERATURE_>, real_t temperature,
		_Molecule *molecule, _Iterator start, _Iterator end)
	{
		typedef typename _Molecule::atom_type     _Atom;
		typedef array_iterator<_Atom, _Iterator>  iterator;

		maxwell_distribution_ensemble emsemble(temperature);

		iterator it = molecule->make_iterator(start);
		iterator ite = molecule->make_iterator(end);
		for (; it!=ite; ++it)
		{
			emsemble.generate((*it).atomdata->mass, 3, &(*it).V[0]);
		}
	}

	/**
	*  Расчет текущей температуры объекта в предположении максвеловского
	*  распределения.
	*  Степени свободы всей системы как целое игнорируются.
	* @param molecule указатель на молекулу или комплекс молекул
	* @return температура [K]
	*/
	template <typename _Molecule>
	INLINE real_t get(_I2T<TEMPERATURE_>, _Molecule *molecule)
	{
		typedef typename _Molecule::atom_type          _Atom;
		typedef array_iterator<_Atom, range_iterator>  _Iterator;

		_E(real_t) kenergy = 0.;
		unsigned atom_count = molecule->count(ATOM);

		_Iterator it = molecule->make_iterator(range_iterator(0));
		_Iterator ite = molecule->make_iterator(range_iterator(atom_count));
		for (; it!=ite; ++it)
		{
			const _Atom &atom = *it;
			kenergy += atom->mass * scalar_product(atom.V, atom.V);
		}
		kenergy /= 3 * atom_count;
		return Temperature((real_t)kenergy);
	}

	// удаление движения центра масс
	template <typename _Molecule>
	INLINE void remove(_I2T<IMPULSE_>, _Molecule *molecule)
	{
		typedef typename _Molecule::atom_type          _Atom;
		typedef array_iterator<_Atom, range_iterator>  _Iterator;

		unsigned atom_count = molecule->count(ATOM);
		vector_t V = 0; real_t mass = 0.;
		_Iterator it = molecule->make_iterator(range_iterator(0));
		_Iterator ite = molecule->make_iterator(range_iterator(atom_count));
		for (; it!=ite; ++it)
		{
			const _Atom &atom = *it;
			mass += atom->mass;
			V += atom.V * atom->mass;
		}

		V *= (1. / mass);
//		{
//			_S msg = make_string("Deleting the mass center velocity [%10.7e %10.7e %10.7e]", V[0], V[1], V[2]);
//			PRINT_MESSAGE(msg);
//		}

		it = molecule->make_iterator(range_iterator(0));
		for (; it!=ite; ++it)
		{
			_Atom &atom = *it;
			atom.V -= V;
		}
	}

	template <typename _Molecule>
	INLINE real_t get(_I2T<PRESSURE_>, _Molecule *molecule)
	{
		typedef typename _Molecule::atom_type          _Atom;
		typedef array_iterator<_Atom, range_iterator>  _Iterator;

		_E(real_t) kenergy = 0.;
		_E(real_t) virial = 0.;

		unsigned atom_count = molecule->count(ATOM);
		_Iterator it = molecule->make_iterator(range_iterator(0));
		_Iterator ite = molecule->make_iterator(range_iterator(atom_count));
		for (; it!=ite; ++it)
		{
			const _Atom &atom = *it;
			kenergy += atom->mass * scalar_product(atom.V, atom.V);
			virial += scalar_product(atom.X, atom.F);
		}

		box_t box = molecule->get(BOX);
		return (real_t) (ATMOSPHERE_FACTOR * (kenergy + virial) / (3 * box.volume()));
	}

	/**
	*  Данный класс определяет термодинамические свойства объекта.
	*  Формулы использованы из forcefield_review.pdf
	*
	*  Термостат ответственнен:
	*  (1) за установку начальных скоростей атомов при заданной температуре,
	*  (2) определение текущих термодинамических параметров (температура, ..)
	*  (3) поддержание (масштабирование) термодинамических параметров
	*/
	class Thermostat
	{
		typedef Descriptor_<THERMOSTATE_>  _Descriptor;

		_Descriptor param_; // все параметры термостата
		void *ensemble_; // текущий термодинамический ансамбль

	public:

		typedef _Descriptor  descriptor_type;

		template <typename _Molecule>
		Thermostat(const Descriptor_<THERMOSTATE_> &desc, _Molecule *molecule)
		: param_(desc), ensemble_(0)
		{
			unsigned count = molecule->count(ATOM);
			set(TEMPERATURE, desc.temperature, molecule, range_iterator(0), range_iterator(count));
				// приведем в контакт с термостатом и нагреем молекулу

		#define NEW(type, param) \
			ensemble_ = new Ensemble_<type>(param, param_.sampling_time/param_.integration_time + 1); \

			if (param_.ensamble == _S("NVE")) NEW(NVE, molecule->U() + kinetic_energy(molecule));
				// сосчитаем полную энергию и установим ее как целевую константу
			if (param_.ensamble == _S("NVT")) NEW(NVT, param_.temperature)
			if (param_.ensamble == _S("NVP")) NEW(NVP, param_.pressure)
		#undef NEW
		}

		~Thermostat()
		{
		#define DELETE(type) delete (Ensemble_<type>*) ensemble_;
			if (param_.ensamble == _S("NVE")) DELETE(NVE);
			if (param_.ensamble == _S("NVT")) DELETE(NVT);
			if (param_.ensamble == _S("NVP")) DELETE(NVP);
		#undef DELETE
		}

		template <typename _Molecule>
		void sample_statistics(_Molecule *molecule) const
		{
		#define SAMPLE(type) ((Ensemble_<type>*)ensemble_)->sample_statistics(molecule);
			if (param_.ensamble == _S("NVE")) SAMPLE(NVE);
			if (param_.ensamble == _S("NVT")) SAMPLE(NVT);
			if (param_.ensamble == _S("NVP")) SAMPLE(NVP);
		#undef MACRO
		}

		template <typename _Molecule>
		std::string scaling(_Molecule *molecule) const
		{
		#define SCALING(type) return ((Ensemble_<type>*)ensemble_)->scaling(molecule)
			if (param_.ensamble == _S("NVE")) SCALING(NVE);
			if (param_.ensamble == _S("NVT")) SCALING(NVT);
			if (param_.ensamble == _S("NVP")) SCALING(NVP);
		#undef SCALING
			return _S("");
		}
	};

}
#endif
