#ifndef _MDYNAMICS__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H
#define _MDYNAMICS__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H

#include "molkern/__moldefs.h"
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
	template <int INTEGRATOR_TYPE> class Integrator_;

	template <> class Integrator_<LEAP_FROG_>
	{
		unsigned dt_; // шаг интегрирования (число fs)
		unsigned sampling_time_; // шаг сбора статистики (fs)
		std::ofstream file1_; // базовый файл статистики
		std::ofstream file2_; // вспомогательный файл статистики
		std::string filename_; // имя файла статистики
		bool prn_water_; // печатать ли воду в файл статистики

		std::vector<real_t> _1mass_;      // инверсные массы атомов (1 / mass)
		std::vector<vector_t> X0_, X1_; // координаты в предыдущей, текущей точке
		std::vector<vector_t> V0_, V1_; // скорости в предыдущей, текущей точке
		std::vector<vector_t> A0_, A1_; // ускорения в предыдущей, текущей точке

	public:

		Integrator_() {}

		/**
		*  Инициализация интегратора и установление начальных скоростей.
		* @param molecule комплекс молекул
		* @param dt шаг интегрирования
		* @param average_time время усреднения
		*/
		template <typename _Molecule>
		Integrator_(_Molecule *molecule,
			unsigned dt=DEFAULT_INTEGRATION_TIME_STEP,
			unsigned sampling_time=DEFAULT_SAMPLING_TIME_STEP,
			const std::string &statfile=_S(""),
			bool prn_water=true)
		: dt_(dt), sampling_time_(sampling_time), prn_water_(prn_water)
		{
			_S filename = basename(statfile); // удаляем путь и расширение

			// открываем файл статистики, если он задан
			if (filename != _S(""))
			{
				std::string file1_name = statfile + _S(".bmm");
				std::string file2_name = statfile + _S(".bmf");
				file1_.open(file1_name.c_str(), std::ios_base::binary);
				file2_.open(file2_name.c_str(), std::ios_base::binary);
				if (!file1_ || !file2_)
				{
					std::string msg = _S("[ERROR] can't open statistics file ") + statfile;
					PRINT_BREAK(msg);
				}
				filename_ = statfile;

				// сохраним ящик и число атомов
				molecule->write_header(_I2T<FORMAT_BMM_>(), file1_);
			}

			typedef typename _Molecule::atom_type          _Atom;
			typedef array_iterator<_Atom, range_iterator>  _Iterator;

			_Iterator it; unsigned i; // постоянно используемые переменные цикла
			unsigned atom_count = molecule->count(ATOM);

			_1mass_.resize(atom_count);
			X0_.resize(atom_count); X1_.resize(atom_count);
			V0_.resize(atom_count); V1_.resize(atom_count);
			A0_.resize(atom_count); A1_.resize(atom_count);

			//------------------------------------------------------------------------
			//               инициализация цикла молекулярной динамики
			//------------------------------------------------------------------------
			molecule->dU__dX(); // предвычисление сил

			_Iterator itb = molecule->make_iterator(range_iterator(0));
			_Iterator ite = molecule->make_iterator(range_iterator(atom_count));
			for (it=itb,i=0; it!=ite; ++it, i++)
			{
				const _Atom &atom = *it;
				_1mass_[i] = (real_t)1. / atom->mass;
				X0_[i] = atom.X;
				A0_[i] = atom.F * _1mass_[i];
				V0_[i] = atom.V - A0_[i] * ((real_t)0.0005 * dt_);
					// скорость в промежуточной точке -dt/2 + преобразование в ps
			}
			remove(IMPULSE, molecule);

			//------------------------------------------------------------------------
			//               основной цикл молекулярной динамики
			//------------------------------------------------------------------------
			PRINT_MESSAGE("\n      ***** ____ start of molecular dynamics ____ *****\n");
		}

		~Integrator_()
		{
			//------------------------------------------------------------------------
			//              конец основного цикла молекулярной динамики
			//------------------------------------------------------------------------
			PRINT_MESSAGE("\n      ***** _____ end of molecular dynamics _____ *****");

			file1_.close();
			file2_.close();

			if (filename_!=_S("")) SAVED_OK_MESSAGE(filename_ + _S("{.bmm, .bmf}"));
		}

		/**
		*  Выполненение интегрирования уравнений молекулярной динамики.
		* @param thermostat тип термодинамического ансамбля
		* @param molecule комплекс молекул
		* @param process_time время процесса [fs]
		* @return число выполненых шагов динамики
		*  Возвращается целое число, чтобы точно отсчитывать пройденное время.
		*/
		template <typename _Thermostat, typename _Molecule>
		void run(const _Thermostat *thermostat, _Molecule *molecule,
			unsigned process_time)
		{
			typedef typename _Molecule::atom_type          _Atom;
			typedef array_iterator<_Atom, range_iterator>  _Iterator;
			unsigned atom_count = molecule->count(ATOM);

			_Iterator it; unsigned i; // постоянно используемые переменные цикла
			_Iterator itb = molecule->make_iterator(range_iterator(0));
			_Iterator ite = molecule->make_iterator(range_iterator(atom_count));

			global_model_time = 0; unsigned curr_time = 0; // обнулим счетчики времени
			real_t dt = (real_t) 0.001 * dt_; // физическое время в [ps]
			real_t maxv = (real_t) DEFAULT_MAX_X / dt; // максимально разрешенная скорость

			system_time_t start_time = current_time(); // фиксируем физическое время старта симуляции

			TIME_TESTING_START("Full cycle of dynamisc finished...", 1)
			while ((global_model_time += dt_) < process_time)
			{
				unsigned stopped_atoms = 0; // число атомов для которых скорости сброшены

				// скорости и координаты в следующей точке
	//			real_t v2=0, v2m; unsigned c=0;
				for (i=0; i<atom_count; i++)
				{
					V1_[i] = V0_[i] + A0_[i] * dt;
					real_t sp = scalar_product(V1_[i], V1_[i]);
//					if (_1mass_[i] < 0.1) { v2 += sp; c++; if (v2m<sp) v2m=sp; }
					if (sp > sqr(maxv))
					{
						// скинем очень большие скорости, что в свою очередь ограничит смещение
						V1_[i].normalize(maxv);
						stopped_atoms++;
					}
					X1_[i] = X0_[i] + V1_[i] * dt;
				}
//				v2 /= c;
//				v2 = sqrt(v2);
//				v2m = sqrt(v2m);

				// сохраним новые координаты и вычислим силы в новой конфигурации
				molecule->read_md_position(&X1_[0], &V1_[0]);

				molecule->dU__dX();
				thermostat->sample_statistics(molecule);
				//thermostat->scaling(molecule);
				if ( (curr_time += dt_) >= sampling_time_ )
				{
//					thermostat->scaling(molecule);
					// сделаем скейлинг скоростей согласно термодинамическому ансамблю
					remove(IMPULSE, molecule);
					_S msg = _S(" ") + make_string(current_model_time()) + _S("   ")
						+ thermostat->scaling(molecule)
						+ _S("   time =") + make_string(current_time() - start_time)
						; //+ _S(" ") + make_string(function_timer_t<CALC_TIMER_>::last_time);
//					if (stopped_atoms) msg += _S(" *");
//					msg += make_string(" <v>=%e mv=%e", v2, v2m);
					PRINT_MESSAGE(msg);
	//				std::cout << make_string(current_model_time()) << " " <<  make_string(current_time() - start_time) << std::endl;

//					if (file1_.is_open() && file2_.is_open())
//					{
//						std::ofstream::pos_type pos = file1_.tellp();
//						file2_.write((char*)&pos, sizeof(pos));
//							// позиция записи в поток сохраняется во вспомогательном файле
//						molecule->save(FORMAT_BMM, file1_, prn_water_);
//
//						// реальная запись, чтобы можно было оборвать выполнение и
//						// файлы при этом сохранились
//						file2_.flush();
//						file1_.flush();
//					}
					curr_time = 0;
				}
				for (it=itb,i=0; it!=ite; ++it, i++)
				{
					V1_[i] = (*it).V;
					A1_[i] = (*it).F * _1mass_[i];
				}

				X1_.swap(X0_);
				V1_.swap(V0_);
				A1_.swap(A0_);
			}
			TIME_TESTING_FINISH;

#ifdef USE_VERLET_TABLE
			{
//				_S msg = make_string("verlet table rebuilding count = %d", molecule->count(REBUILD));
//				PRINT_MESSAGE(msg);
			}
#endif
//			unsigned md_steps = global_model_time / dt_;
//			_S msg = make_string("md steps = %d", md_steps) + _S(" time per md step = ")
//					+ make_string(double(current_time() - start_time) / md_steps);
//			PRINT_MESSAGE(msg);
		}

	};

	typedef Integrator_<LEAP_FROG_> Integrator;

}
#endif
