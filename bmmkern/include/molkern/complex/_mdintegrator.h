#ifndef _MDINTEGRATOR__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H
#define _MDINTEGRATOR__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H

#include "molkern/__moldefs.h"
#include "molkern/__config.h"
#include "molkern/complex/_thermostat.h"

namespace molkern
{
	using namespace prgkern;

	/**
	*  Обменивает импульсы двух частиц при столкновении для 1D случая. Функция
	*  работает для любого потенциала взаимодействия, однако она, принципиально
	*  в силу своей простоты, не может учесть сдвиг фазы при столкновении и,
	*  следовательно, найти время столкновения, что  могло бы быть полезным
	*  в динамике. При свопинге не учитываются координаты объектов, то есть
	*  свопинг выполняется всегда, даже если реально столкновения быть не может.
	*  Например, это случай, когда частица с меньшей скоростью догоняет
	*  частицу, движующуюся с большей скоростью.
	*
	* @param m1,m2 массы частиц
	* @param v1,v2 скорости частиц ([in]/[out] до/после столкновния)
	*/
	INLINE void impulse_swap(real_t m1, real_t &v1, real_t m2, real_t &v2)
	{
		real_t P = v1 * m1 + v2 * m2;
		real_t E = m1 * sqr(v1) + m2 * sqr(v2);

		real_t _1m1 = (real_t) 1. / m1;
		real_t alpha = m2 * (m2 * _1m1 + (real_t) 1.);
		real_t beta  = mult2(P * m2 * _1m1);
		real_t gamma = sqr(P) * _1m1 - E;
		real_t _1__2alpha = (real_t) 1. / mult2(alpha);

		assert(_GE(sqr(beta) - 4 * alpha * gamma, 0.));
			// проверка дискриминанта (никогда не должна выскочить)

		real_t d = (real_t) sqrt(sqr(beta) - 4 * alpha * gamma);
		real_t s = beta - v2 * mult2(alpha);

		if (fabs(s + d) < fabs(s - d)) v2 = (beta - d) * _1__2alpha;
		else v2 = (beta + d) * _1__2alpha;
			// делаем выбор по максимальной разнице,
			// чтобы из двух решений не попасть в исходное

		v1 = (P - m2 * v2) * _1m1;

		assert(_LT(fabs(P - m1*v1 - m2*v2), std::numeric_limits<real_t>::epsilon()));
		assert(_LT(fabs(E - m1*sqr(v1) - m2*sqr(v2)), std::numeric_limits<real_t>::epsilon()));
			// проверка сохранения импульса и энергии
	}

	/**
	*  Обменивает импульсы двух частиц для 3D случая. В обмене не участвуют
	*  компоненты импульса ортогональные к вектору R, соединяющему частицы.
	*  Обмен производится, если расстояние между атомы станет меньше критического,
	*  либо для ситуации когда атомы "проходят" через друг друга. Второй случай
	*  учитывается тем, что проверяется величины с учетом знака.
	*
	* @param R вектор, направленный от 1 частицы ко 2.
	* @param sigma критическое расстояние
	* @param dt шаг по времени
	* @param m1,m2 массы частиц
	* @param v1,v2 скорости частиц ([in]/[out] до/после столкновния)
	*/
	template <unsigned N>
	INLINE void impulse_swap(vdense_<N, real_t> R, real_t sigma, real_t dt,
		real_t m1, vdense_<N, real_t> &v1, real_t m2, vdense_<N, real_t> &v2)
	{
		assert(_GT(m1, 0.)); assert(_GT(m2, 0.));
		assert(_GE(dt, 0.)); assert(_GT(sigma, 0.));
			// функция определена только для положительных времен и масс

		real_t r = R.normalize();
		assert(_GT(r, sigma));
			// в начальном состоянии атомы еще не столкнулись

		real_t s1 = scalar_product(v1, R); // проекция скорости 1 на R
		real_t s2 = scalar_product(v2, R); // проекция скорости 2 на R
		if (r + (s2 - s1) * dt >= sigma) return; // нет столкновения
			// столкновение - это когда частицы сближаются любым образом (сталкиваются,
			// либо одна догоняет другую)

		vdense_<N, real_t> v1__ = v1 - R * s1; // ортогональная компонента 1
		vdense_<N, real_t> v2__ = v2 - R * s2; // ортогональная компонента 2
		impulse_swap(m1, s1, m2, s2); // свопируем проекции на R

		v1 = v1__ + R * s1;
		v2 = v2__ + R * s2;
	}

	/**
	*  Обменивает импульс частицы для 3D случая при столкновении со "стенкой".
	*  Под "стенкой" понимается жестко зажатая частица в некоторой координате.
	*   В обмене не участвуют компоненты импульса ортогональные к вектору R,
	*  соединяющему частицы. Обмен производится, если расстояние между атомом и
	*  стенкой  станет меньше критического, либо для ситуации когда атом "проходит"
	*  через стенку.
	*
	* @param R вектор, направленный от частицы к "стенке"
	* @param sigma критическое расстояние
	* @param dt шаг по времени
	* @param m масса частицы
	* @param v скорость частицы ([in]/[out] до/после столкновния)
	*/
	template <unsigned N>
	INLINE void impulse_swap(vdense_<N, real_t> R, real_t sigma, real_t dt,
		real_t m, vdense_<N, real_t> &v)
	{
		assert(_GT(m, 0.)); assert(_GE(dt, 0.)); assert(_GT(sigma, 0.));
			// функция определена только для положительных времен и масс

		real_t r = R.normalize();
		assert(_GT(r, sigma));
			// в начальном состоянии атомы еще не столкнулись

		real_t s = scalar_product(v, R); // проекция скорости на R
		if (r - s * dt >= sigma) return; // нет столкновения

		v -= R * (s + s);
	}

	enum { LEAP_FROG_};
	template <int INTEGRATOR_TYPE, typename _Ensemble, typename _System>
	class Integrator_;

	template <typename _Ensemble, typename _System>
	class Integrator_<LEAP_FROG_, _Ensemble, _System>
	: public system_notifier_t
	{
		mutable _Ensemble *ensemble_;      // статистика для счета "мгновенных" термодинамических параметров
		mutable _System *system_;          // система, для которой выполняется интегрирование и статистика

		real_t dt_;              // шаг интегрирования (число ps)
		real_t process_time_;    // полное время интегрирования
		real_t max_velocity_;    // максимально разрешенная скорость атомов, чтобы не развалить систему
		real_t M_;               // суммарная масса всех атомов системы

		std::vector<real_t> invmass_, mass_;    // инверсные и обычные массы атомов (1 / mass)
		std::vector<vector_t> X0_, X1_;    // координаты в предыдущей, текущей точке
		std::vector<vector_t> V0_, V1_;    // скорости в предыдущей, текущей точке
		std::vector<vector_t> A0_, A1_;    // ускорения в предыдущей, текущей точке

	public:

		Integrator_() {}

		/**
		*  Инициализация интегратора и установление начальных скоростей.
		* @param molecule комплекс молекул
		* @param dt шаг интегрирования
		* @param average_time время усреднения
		*/
		Integrator_(const Configure *conf, _Ensemble *ensemble, _System *system)
		: system_notifier_t (global_model_time)
		, ensemble_         (ensemble)
		, system_           (system)
		, dt_               (conf->integration_time)
		, process_time_     (conf->process_time)
		, max_velocity_     (conf->max_displacement / dt_)
		{
			typedef typename _System::atom_type            _Atom;
			typedef array_iterator<_Atom, range_iterator>  _Iterator;

			_Iterator it; unsigned i; // постоянно используемые переменные цикла
			unsigned N = system->count(ATOM);

			invmass_.resize(N);
			mass_.resize(N);
			X0_.resize(N); X1_.resize(N);
			V0_.resize(N); V1_.resize(N);
			A0_.resize(N); A1_.resize(N);

			//------------------------------------------------------------------------
			//               инициализация цикла молекулярной динамики
			//------------------------------------------------------------------------
			system->dU__dX(); // предвычисление сил

			_Iterator itb = system->make_iterator(range_iterator(0));
			_Iterator ite = system->make_iterator(range_iterator(N));
			for (it=itb,i=0; it!=ite; ++it, i++)
			{
				const _Atom &atom = *it;
				mass_[i] = atom->mass;
				invmass_[i] = (real_t)1. / mass_[i];
				M_ += mass_[i];
				X0_[i] = atom.X;
				A0_[i] = atom.F * invmass_[i];
				V0_[i] = atom.V - 0.5 * dt_ * A0_[i]; // скорость в промежуточной точке -dt/2
			}
			ensemble_->sample_statistics(); // чтобы занести 0-точку и избежать счета с нулями в термостате

			//------------------------------------------------------------------------
			//               основной цикл молекулярной динамики
			//------------------------------------------------------------------------
			PRINT_MESSAGE("\n      ***** ____ start of molecular dynamics ____ *****\n");

			global_model_time.init(conf->integration_time, "ps"); // обнулим счетчики времени
			global_start_time = current_time(); // фиксируем физическое время старта симуляции
		}

		~Integrator_()
		{
			//------------------------------------------------------------------------
			//              конец основного цикла молекулярной динамики
			//------------------------------------------------------------------------
			PRINT_MESSAGE("\n      ***** _____ end of molecular dynamics _____ *****");
		}

		/**
		*  Выполненение интегрирования уравнений молекулярной динамики.
		* @param thermostat тип термодинамического ансамбля
		* @param system комплекс молекул
		* @param process_time время процесса [fs]
		* @return число выполненых шагов динамики
		*  Возвращается целое число, чтобы точно отсчитывать пройденное время.
		*/
		template <typename _Thermostat>
		void run(const _Thermostat *thermostat)
		{
			typedef typename _System::atom_type            _Atom;
			typedef array_iterator<_Atom, range_iterator>  _Iterator;

			unsigned N = system_->count(ATOM);
			_Atom *atoms = system_->get(ATOM);
			real_t max_velocity2_ = sqr(max_velocity_);

			_Iterator it; // постоянно используемые переменные цикла
			_Iterator itb = system_->make_iterator(range_iterator(0));
			_Iterator ite = system_->make_iterator(range_iterator(N));

			long unsigned cycles_to_be_processed = round(process_time_ / dt_);
			for (long unsigned ti=0; ti<cycles_to_be_processed; ti++)
			{
				vector_t P = 0.;
				// скорости и координаты в следующей точке
				for (unsigned i=0; i<N; i++)
				{
					V1_[i] = V0_[i] + A0_[i] * dt_;
					if (V1_[i].length2() > max_velocity2_)
					{
						// скинем очень большие скорости, что в свою очередь ограничит смещение
						V1_[i].normalize(max_velocity_);
					}
					P += V1_[i] * mass_[i];
				}
				P *= (1. / M_);

				for (unsigned i=0; i<N; i++)
				{
					V1_[i] -= P; // удалим движение центра масс
					X1_[i] = X0_[i] + V1_[i] * dt_;
				}

				// сохраним новые координаты и вычислим силы в новой конфигурации
				system_->read_md_position(&X1_[0], &V1_[0]);
				system_->dU__dX();

				++(*this); // увеличим счетчик событий и уведомим об этом все связанные объекты

				for (unsigned i=0; i<N; i++)
				{
					V1_[i] = atoms[i].V; // сохраним скалированные термостатом скорости
					A1_[i] = atoms[i].F * invmass_[i];
				}

				X1_.swap(X0_);
				V1_.swap(V0_);
				A1_.swap(A0_);
			}
		}
	};

	typedef Integrator_<LEAP_FROG_, Ensemble, Complex> Integrator;

}
#endif
