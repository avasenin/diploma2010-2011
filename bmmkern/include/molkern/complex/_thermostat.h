#ifndef _THERMOSTAT__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H
#define _THERMOSTAT__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H

#include "molkern/__moldefs.h"
#include "molkern/__config.h"
#include "molkern/complex/_ensemble.h"
#include "molkern/complex/_geom_tool.h"


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
		: _Base(_Distribution(sqrt(KT(temperature)/mass))),
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
			unsigned key = unsigned(ExtMass(mass));
				// ключ не может работать правильно для внутренних масс

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
			unsigned key = unsigned(ExtMass(mass));

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
	template <typename _System>
	INLINE void set(_I2T<TEMPERATURE_>, real_t temperature, _System *system)
	{
		typedef typename _System::atom_type            _Atom;
		typedef array_iterator<_Atom, range_iterator>  _Iterator;

		maxwell_distribution_ensemble emsemble(temperature);
		unsigned N = system->count(ATOM);

		_Iterator it = system->make_iterator(range_iterator(0));
		_Iterator ite = system->make_iterator(range_iterator(N));
		for (; it!=ite; ++it)
		{
			emsemble.generate((*it).atomdata->mass, 3, &(*it).V[0]);
				// корректно, для внутренних масс
		}
	}

	/*
	 *  Удаление движения центра масс у всей системы в целом
	 */
	template <typename _System>
	INLINE void remove(_I2T<IMPULSE_>, _System *system)
	{
		typedef typename _System::atom_type            _Atom;
		typedef array_iterator<_Atom, range_iterator>  _Iterator;

		unsigned N = system->count(ATOM);
		vector_t V = 0; real_t mass = 0.;
		_Iterator it = system->make_iterator(range_iterator(0));
		_Iterator ite = system->make_iterator(range_iterator(N));
		for (; it!=ite; ++it)
		{
			const _Atom &atom = *it;
			mass += atom->mass;
			V += atom.V * atom->mass;
		}
		V *= (1. / mass);

		it = system->make_iterator(range_iterator(0));
		for (; it!=ite; ++it)
		{
			_Atom &atom = *it;
			atom.V -= V;
		}
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
	template <typename _Ensemble, typename _System>
	class Thermostat_ : public notified_object_
	{
		mutable _Ensemble *ensemble_;      // микроансамбль для "мгновенного" усреднения
		mutable _System *system_;          // система, с которой общается термостат
		_S ensemble_type_;       // тип ансамбля
		real_t targetP_;         // целевое давление
		real_t targetT_;         // целевая температура
		real_t targetE_;         // целевая энергия

		/*
		 * Термодинамические свойства системы (например, коэффициенты передачи тепла и давления)
		 * правильней было бы разместить в классе COMPLEX, но я не хочу его перегружать термодинамикой.
		 */
		real_t couplingT_coef_;
		real_t couplingP_coef_;

		real_t integration_time_;
		real_t relaxation_time_;

		virtual void process_notification()
		{
			notified_object_::process_notification();

			/*
			 *  Cкалирование при каждом уведомлении. При малом числе циклов динамики скалирование
			 *  делаем по предельному варианту, то есть точное, чтобы адаптировать систему от клешей.
			 *  Далее скалируем ее термодинамическим образом, через постепенную релаксацию.
			 */
			if (counter_ < 100) scaling_(direct_scalingT_coef_(), 1.);
			else scaling();
		}

	public:

		Thermostat_(const Configure *conf, _Ensemble *ensemble, _System *system)
		: notified_object_  (1)
		, ensemble_         (ensemble)
		, system_           (system)
		, couplingT_coef_   (conf->couplingT_coef)
		, couplingP_coef_   (conf->couplingP_coef)
		, integration_time_ (conf->integration_time)
		, relaxation_time_  (conf->relaxation_time)
		{
			Interface_<ENSEMBLE_> interface = conf->get_interface(ENSEMBLE);
			ensemble_type_ = interface.ensemble_type;
			targetP_       = interface.pressure;
			targetT_       = interface.temperature;
			targetE_       = interface.energy + system->get(POTENT_ENERGY);
			set(TEMPERATURE, targetT_, system_); // приведем в контакт с термостатом и нагреем молекулу
			remove(IMPULSE, system_); // удалим случано набранный импульс
		}

		void scaling() const
		{
			if (ensemble_type_ == _S("NVT"))
			{
				scaling_(berendsen_scalingT_coef_(), 1.);
			}
			if (ensemble_type_ == _S("NVP"))
			{
				scaling_(berendsen_scalingT_coef_(), berendsen_scalingP_coef_());
			}
		}

		real_t direct_scalingT_coef_() const
		{
			real_t currentT = ensemble_->T_.top();
			return sqrt(targetT_ / currentT);
		}

		real_t berendsen_scalingT_coef_() const
		{
			real_t currentT = ensemble_->T_.average();
				// использую среднее по интервалу вместо мгновенного, так как оно менее подвержено скачкам

			real_t lambda = 1. - (integration_time_ * couplingT_coef_ / relaxation_time_ )
				* (currentT - targetT_) / currentT;
			if (lambda > 0.) lambda = sqrt(lambda);
			else lambda = 1.;

			return lambda;
		}

		real_t berendsen_scalingP_coef_() const
		{
			real_t currentP = ensemble_->P_.average();
				// использую среднее по интервалу вместо мгновенного, так как оно менее подвержено скачкам

			real_t lambda = 1. + (integration_time_ * couplingP_coef_ / relaxation_time_ )
				* (currentP - targetP_);
			if (lambda > 0.) lambda = pow(lambda, 1./3);
			else lambda = 1.;

			return lambda;
		}

	protected:

		void scaling_(real_t vlambda, real_t xlambda = 1.) const
		{
			typedef typename _System::region_type          _Region;
			typedef typename _System::atom_type            _Atom;
			typedef array_iterator<_Atom, range_iterator>  _Iterator;

			unsigned N = system_->count(_I2T<ATOM_>());
			_Iterator it = system_->make_iterator(range_iterator(0));
			_Iterator ite = system_->make_iterator(range_iterator(N));
			for (; it!=ite; ++it)
			{
				(*it).V *= vlambda;
				(*it).X *= xlambda;
			}

			_Region *region = system_->get(REGION);
			region->rescale(xlambda);
		}

		void print() const
		{
			_S msg = make_string(global_model_time) + _S(" ") + ensemble_->print_statistics()
				+ make_string(current_time() - global_start_time);
			PRINT_MESSAGE(msg);
		}

	};
	typedef Thermostat_<Ensemble, Complex>  Thermostat;

}
#endif
