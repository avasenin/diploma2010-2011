#ifndef _INTERACTION__F9ED1116_D4D6_5cf5_C20E_AC432B290700__H
#define _INTERACTION__F9ED1116_D4D6_5cf5_C20E_AC432B290700__H

#include "molkern/__moldefs.h"

namespace molkern
{

/*==============================================================================
 *                     ОПИСАНИЕ ОСОБЕННОСТЕЙ РЕАЛИЗАЦИИ
 *  Взаимодействия считаются без существенного накопления каких-либо сумм,
 *  по этой причине вполне достаточную точность обеспечивают типы real_t.
 *
 *  В ранних реализациях функции U(*coul, *vdw), обладали побочным эффектом,
 *  они не только высчитывали парциальные энергии, но производили наколение их
 *  в каких-то элементах. Поскольку накопление требует другой точности, то
 *  нужно определять типы для переменных coul и vdw как шаблонные типы. Однако,
 *  более правильным является уничтожение побочного эффекта, за накопление
 *  пусть отвечает вызывающая программа. Это может потребовать чуть больших
 *  затрат на передачу результата, но для U(), которая не используется в
 *  массивных вычислениях, это приемлимо.
 *
 *  Используются следующие потенциалы:
 *  (1) COUL - кулоновский потенциал 1/r, который сглаживается на расстоянии
 *      cutoff методом сдвига. Он обеспечивает гладкость только энергии.
 *  (2) VDW - ван-дер-ваальс 6-12 потенциал. Так как он спадает быстро, то
 *      сглаживаение к нему НЕ применяется, что позволяет избавится от расчета
 *      константы сдвига для каждой пары
 *  (3) ERF - кулоновский потенциал (ближний и дальний), который масштабируется
 *      умножением на функцию erf. Так как эта функция заставляет потенциал
 *      резко спадать, то сглаживание НЕ применяется.
 *
 *  Потенциалы COUL + VDW или ERF + VDW объединены вместе в одну структуру
 *  данных, так как совместно используют некоторые величины, типа r2, tau2...
 *============================================================================*/

	/**
	* @param TYPE тип взаимодействия (NON_BONDED_)
	* @param N порядок сглаживания на cutoff
	*/
	template <int TYPE> struct Interaction_;

	// типы совместных потенциалов COUL + VDW или ERF + VDW
	enum { C612, E612 };


	template <typename T>
	INLINE T calculate_sigma(T sigma, T sigma__)
	{
	#ifdef USE_HARM_AVERAGE_OF_SIGMA
		return (sigma * sigma__) * 4.f;
	#else
		return sqr(sigma + sigma__);
	#endif
	}

	/**
	*   Совместный COUL + VDW потенциал.
	*/
	template <> struct Interaction_<C612>
	{
	private:

	#ifdef USE_PSEUDO_POTENTIAL_IN_ZERO
		// Функция кулона (1/x) erf(sqrt(pi) x/s) сшивается с 1/x в точке x = s
		// через произведение с полиномом второго порядка P = 1 + o(x)
		//     P(x) = 1 + c1 x/s + c2 (x/s)^2)
		// Важно, что этот полином равен 1 в нуле. Отказ от этого условия
		// приводит к тому, что хотя и в x=s обеспечивается гладкая сшивка,
		// но в точках дальних (например, в 0) полином сильно уходит в сторону.
		//
		// Введя обозначения E = exp(-pi) и R = erf(sqrt(pi)) константы равны:
		//    с1 = 2 * ( (1 + E/R)/R - 1)
		//    с2 = 1 - (1 + 2 E/R)/R
		static const real_t c1__;
		static const real_t c2__;
		static const real_t sqrt_pi;

		static real_t alpha_; ///< значение vdw псевдопотенциала в 0 (в единицах eps)
		static real_t tau2_; ///< квадрат точки сшивки псевдопотенциала tau2 = sqr(r / sigma)
		static real_t gamma_; ///< коэффициент при полиноме V = alpha + gamma * r^12
			///< (в единицах eps)

	#endif // USE_PSEUDO_POTENTIAL_IN_ZERO

		static real_t shift_factor_; // константа сдвига
		static real_t rcutoff_; // радиус взаимодействия

	public:

		/**
		* @param cutoff граница области взаимодействия
		* @param barrier значение барьера потенциала в 0 (в единицах eps)
		*  То есть, если мы хотим, чтобы в нуле потенциал был в N раз больше
		*  чем глубина ямы (eps), то устанавливает barrier = N.
		* @note Важно, алгоритм определен только для barrier > -1
		*/
		Interaction_(real_t cutoff, real_t barrier)
		{
			rcutoff_ = cutoff;
			shift_factor_ = (real_t) 1. / cutoff;
		#ifdef USE_PSEUDO_POTENTIAL_IN_ZERO
			set_barrier(barrier);
				// инициализация барьера и функций сшивки на нем
		#endif
		}

		static real_t interaction_radius() { return rcutoff_; }

	#ifdef USE_PSEUDO_POTENTIAL_IN_ZERO
		/**
		* @note Важно, алгоритм определен только для barrier > -1
		* @param barrier значение барьера потенциала в 0 (в единицах eps)
		*/
		static void set_barrier(real_t barrier)
		{
			assert(_GE(barrier, -1.));

			alpha_ = barrier;
			real_t _tau6 = (real_t) sqrt((8./9.) * alpha_ + 1.);
			_tau6 = (real_t) 0.75 * _tau6 + (real_t) 0.75;
			real_t _tau12 = sqr(_tau6);
			gamma_ = _tau12 * (_tau6 - _tau12);
			tau2_ = (real_t) pow(1. / _tau6, (real_t) 1./3.);
		}
	#endif

		/**
		* @brief энергия парного кулона
		* @param sigma2 sqr(sigma1 + sigma2), нужен для вычисления потенциала в нуле
		* @param charge произведение зарядов атомов
		* @param r2 квадрат расстояния r^2
		* @return вклад кулона в энергию
		*/
		static real_t U(_I2T<COUL_>, real_t sigma2, real_t charge,
			real_t r2, real_t tau2, real_t r_1)
		{
			real_t coul_energy = 0.;

		#ifndef SKIP_COUL_ENERGY

		#ifdef USE_PSEUDO_POTENTIAL_IN_ZERO
			if (tau2 < tau2_)
			{
				real_t tau = (real_t) sqrt(tau2);
				coul_energy = (real_t) 1. + c1__ * tau + c2__ * tau2;

				if (r2 < sqr(std::numeric_limits<real_t>::epsilon())) coul_energy *= (real_t) (2. / sqrt(sigma2));
				else coul_energy *= r_1 * (real_t) erf(sqrt_pi * tau);
			}
			else
		#endif // USE_PSEUDO_POTENTIAL_IN_ZERO
			coul_energy = r_1; // (осторожно!) оператор привязан к #ifdef

			coul_energy = charge * (coul_energy - shift_factor_);
		#endif // SKIP_COUL_ENERGY

			return coul_energy;
		}

		/**
		* @brief энергия парного ван-дер-ваальса
		* @param sigma sigma1 + sigma2
		* @param eps eps1 * eps2 (в реальности хранится epsi <- sqrt(epsi))
		* @param r2 квадрат расстояния r^2
		* @return ван-дер-ваальсовый вклад в энергию
		*/
		static real_t U(_I2T<VDW_>, real_t sigma2, real_t eps,
			real_t r2, real_t tau2, real_t r_1)
		{
			real_t vdw_energy = 0.;

		#ifndef SKIP_VDW_ENERGY
				real_t tau6 = cube(tau2);

		#ifdef USE_PSEUDO_POTENTIAL_IN_ZERO
			if (tau2 < tau2_)
			{
				vdw_energy = eps * (alpha_ + gamma_ * sqr(tau6));
			}
			else
		#endif // USE_PSEUDO_POTENTIAL_IN_ZERO
			{
				real_t tau_6 = (real_t) 1. / tau6;
				vdw_energy = eps * tau_6 * (tau_6 - (real_t) 2.);
			}
		#endif // SKIP_VDW_ENERGY

			return vdw_energy;
		}


		/**
		* @brief энергия и силы парного взаимодействия
		* @param du__dq[out] первая производная (внешний накопитель, значит нельзя обнулять)
		* @param sigma2 sqr(sigma1 + sigma2), нужен для вычисления потенциала в нуле
		* @param charge произведение зарядов атомов
		* @param r2 квадрат расстояния r^2
		* @return энергия
		*/
		static real_t dU__dX(_I2T<COUL_>, real_t *du__dq, real_t sigma2, real_t charge,
			real_t r2, real_t tau2, real_t r_1)
		{
			real_t coul_energy = 0.;

		#ifndef SKIP_COUL_ENERGY
				real_t r_2 = sqr(r_1);

		#ifdef USE_PSEUDO_POTENTIAL_IN_ZERO
			if (tau2 < tau2_)
			{
				real_t tau = (real_t) sqrt(tau2);
				real_t c1 = c1__ * tau;
				real_t c2 = c2__ * tau2;
				real_t exp_ = (real_t) mult2(tau * exp(-M_PI * tau2));
				real_t erf_ = (real_t) erf(sqrt_pi * tau);

				coul_energy = (real_t) 1. + c1 + c2;
				if (r2 < sqr(std::numeric_limits<real_t>::epsilon())) coul_energy *= (real_t)(2. / sqrt(sigma2));
				else coul_energy *= r_1 * erf_;

				*du__dq += charge * r_2 * (exp_ - erf_ + c1 * exp_ + c2 * (exp_ + erf_));
			}
			else
		#endif // USE_PSEUDO_POTENTIAL_IN_ZERO
			{
				coul_energy = r_1;
				*du__dq += -charge * r_2;
			}

			coul_energy = charge * (coul_energy - shift_factor_);
		#endif // SKIP_COUL_ENERGY

			return coul_energy;
		}

		/**
		* @brief энергия и силы парного взаимодействия
		* @param du__dq[out] первая производная
		* @param sigma sigma1 + sigma2
		* @param eps eps1 * eps2 (в реальности хранится epsi <- sqrt(epsi))
		* @param charge произведение зарядов атомов
		* @param r2 квадрат расстояния r^2
		* @return энергия
		*/
		static real_t dU__dX(_I2T<VDW_>, real_t *du__dq, real_t sigma2, real_t eps,
			real_t r2, real_t tau2, real_t r_1)
		{
			real_t vdw_energy = 0.;

		#ifndef SKIP_VDW_ENERGY
			real_t tau6 = cube(tau2);

		#ifdef USE_PSEUDO_POTENTIAL_IN_ZERO
			if (tau2 < tau2_)
			{
				real_t tau12 = sqr(tau6);
				vdw_energy = eps * (alpha_ + gamma_ * tau12);
				*du__dq += (real_t) 12. * eps * gamma_ * tau12 * r_1;
			}
			else
		#endif // USE_PSEUDO_POTENTIAL_IN_ZERO
			{
				real_t tau_6 = (real_t) 1. / tau6;
				real_t eps_tau_6 = eps * tau_6;

				vdw_energy = eps_tau_6 * (tau_6 - (real_t) 2.);
				*du__dq += (real_t) 12. * eps_tau_6 * ((real_t) 1. - tau_6) * r_1;
			}
		#endif // SKIP_VDW_ENERGY

			return vdw_energy;
		}

		static vreal_t dU__dX(_I2T<COUL_>, vreal_t *du__dq, vreal_t sigma2, vreal_t charge,
			vreal_t tau2, vreal_t r_1)
		{
			vreal_t coul_energy = 0.f;

		#ifndef SKIP_COUL_ENERGY
			vreal_t r_2 = sqr(r_1);

		#ifdef USE_PSEUDO_POTENTIAL_IN_ZERO
			if (tau2 < tau2_)
			{
				real_t tau = (real_t) sqrt(tau2);
				real_t c1 = c1__ * tau;
				real_t c2 = c2__ * tau2;
				real_t exp_ = (real_t) mult2(tau * exp(-M_PI * tau2));
				real_t erf_ = (real_t) erf(sqrt_pi * tau);

				coul_energy = (real_t) 1. + c1 + c2;
				if (r2 < sqr(std::numeric_limits<real_t>::epsilon())) coul_energy *= (real_t)(2. / sqrt(sigma2));
				else coul_energy *= r_1 * erf_;

				*du__dq += charge * r_2 * (exp_ - erf_ + c1 * exp_ + c2 * (exp_ + erf_));
			}
			else
		#endif // USE_PSEUDO_POTENTIAL_IN_ZERO
			{
				coul_energy = r_1;
				*du__dq -= charge * r_2;
			}

			coul_energy = charge * (coul_energy - shift_factor_);

		#endif // SKIP_COUL_ENERGY

			return coul_energy;
		}

		static vreal_t dU__dX(_I2T<VDW_>, vreal_t *du__dq, vreal_t sigma2, vreal_t eps,
			vreal_t tau2, vreal_t r_1)
		{
			vreal_t vdw_energy = 0.f;

		#ifndef SKIP_VDW_ENERGY
			vreal_t tau6 = cube(tau2);

		#ifdef USE_PSEUDO_POTENTIAL_IN_ZERO
			if (tau2 < tau2_)
			{
				real_t tau12 = sqr(tau6);
				vdw_energy = eps * (alpha_ + gamma_ * tau12);
				*du__dq += (real_t) 12. * eps * gamma_ * tau12 * r_1;
			}
			else
		#endif // USE_PSEUDO_POTENTIAL_IN_ZERO
			{
				vreal_t tau_6 = vreal_t(1.) / tau6;
				vreal_t eps_tau_6 = eps * tau_6;

				vdw_energy = eps_tau_6 * (tau_6 - vreal_t(2.));
				*du__dq += 12.f * eps_tau_6 * (1.f - tau_6) * r_1;
			}

		#endif // SKIP_VDW_ENERGY

			return vdw_energy;
		}

	};

	template <> struct Interaction_<E612>
	{
	private:

		static real_t alpha_; // параметр, связанный с радиусом разделения потенциала на ближнюю и дальнюю части

	public:

		/**
		* @param split_radius - расстояние разделения потенциала на ближнюю и дальнюю части.
		*  В области (0, split_radius) из дальнего потенциала вырезана кулоновская особенность
		*  в 0 и потенциал близок к константе. В области (split_radius, inf) потенциал
		*  асимптотически приближается к 1/r.
		*/
		Interaction_(real_t split_radius=DEFAULT_COUL_SPLIT_RADIUS)
		{
			alpha_ = (real_t) M_SQRT_PI / split_radius;
		}

		/**
		* @brief энергия парного кулона
		* @param sigma2 sqr(sigma1 + sigma2), нужен для вычисления потенциала в нуле
		* @param charge произведение зарядов атомов
		* @param r2 квадрат расстояния r^2
		* @return вклад кулона в энергию
		*/
		static real_t U(_I2T<COUL_>, real_t sigma2, real_t charge,
			real_t r2, real_t tau2, real_t r_1)
		{
			real_t coul_energy = 0.;

		#ifndef SKIP_COUL_ENERGY

			assert(_GT(r2, (real_t)sqr(std::numeric_limits<real_t>::epsilon())));
				// потенциал определен для r >= M_SAFETY

			coul_energy = charge * ((real_t)1. - (real_t)erf(alpha_ / r_1)) * r_1;

		#endif // SKIP_COUL_ENERGY

			return coul_energy;
		}

		/**
		* @brief энергия парного ван-дер-ваальса
		* @param sigma sigma1 + sigma2
		* @param eps eps1 * eps2 (в реальности хранится epsi <- sqrt(epsi))
		* @param r2 квадрат расстояния r^2
		* @return ван-дер-ваальсовый вклад в энергию
		*/
		static real_t U(_I2T<VDW_>, real_t sigma2, real_t eps,
			real_t r2, real_t tau2, real_t r_1)
		{
			real_t vdw_energy = 0.;

		#ifndef SKIP_VDW_ENERGY
			real_t tau6 = cube(tau2);
			real_t tau_6 = (real_t) 1. / tau6;
			vdw_energy = eps * tau_6 * (tau_6 - (real_t) 2.);
		#endif // SKIP_VDW_ENERGY

			return vdw_energy;
		}


		/**
		* @brief энергия и силы парного взаимодействия
		* @param du__dq[out] первая производная (внешний накопитель, значит нельзя обнулять)
		* @param sigma2 sqr(sigma1 + sigma2), нужен для вычисления потенциала в нуле
		* @param charge произведение зарядов атомов
		* @param r2 квадрат расстояния r^2
		* @return энергия
		*/
		static real_t dU__dX(_I2T<COUL_>, real_t *du__dq, real_t sigma2, real_t charge,
			real_t r2, real_t tau2, real_t r_1)
		{
			real_t coul_energy = 0.;

		#ifndef SKIP_COUL_ENERGY
			assert(_GT(r2, (real_t)sqr(std::numeric_limits<real_t>::epsilon())));
				// потенциал определен для r >= M_SAFETY

			real_t coef = (real_t)2.0 * alpha_ / (real_t)M_SQRT_PI;
			real_t s = coef * (real_t)exp(-sqr(alpha_) * r2);

			real_t erf_ = (real_t)erf(alpha_ / r_1);
			real_t p = ((real_t)1. - erf_) * r_1;
			*du__dq -= charge * (p + s) * r_1;
			coul_energy = charge * p;

		#endif // SKIP_COUL_ENERGY

			return coul_energy;
		}

		/**
		* @brief энергия и силы парного взаимодействия
		* @param du__dq[out] первая производная
		* @param sigma sigma1 + sigma2
		* @param eps eps1 * eps2 (в реальности хранится epsi <- sqrt(epsi))
		* @param charge произведение зарядов атомов
		* @param r2 квадрат расстояния r^2
		* @return энергия
		*/
		static real_t dU__dX(_I2T<VDW_>, real_t *du__dq, real_t sigma2, real_t eps,
			real_t r2, real_t tau2, real_t r_1)
		{
			real_t vdw_energy = 0.;

		#ifndef SKIP_VDW_ENERGY
			real_t tau6 = cube(tau2);

			real_t tau_6 = (real_t) 1. / tau6;
			real_t eps_tau_6 = eps * tau_6;

			vdw_energy = eps_tau_6 * (tau_6 - (real_t) 2.);
			*du__dq += (real_t) 12. * eps_tau_6 * ((real_t) 1. - tau_6) * r_1;
		#endif // SKIP_VDW_ENERGY

			return vdw_energy;
		}

		static vreal_t dU__dX(_I2T<COUL_>, vreal_t *du__dq, vreal_t sigma2, vreal_t charge,
			vreal_t tau2, vreal_t r_1)
		{
			vreal_t coul_energy = 0.f;

		#ifndef SKIP_COUL_ENERGY
			//assert(_GT(r2, (real_t)sqr(std::numeric_limits<real_t>::epsilon())));
				// потенциал определен для r >= M_SAFETY

			vreal_t coef = (vreal_t)(alpha_ * (2./M_SQRT_PI));
			vreal_t w = (vreal_t)alpha_ / r_1;
			vreal_t s = coef * (vreal_t)exp(0-sqr(w));

			vreal_t erf_ = (vreal_t)erf(w);
			vreal_t p = ((vreal_t)1. - erf_) * r_1;
			*du__dq -= charge * (p + s) * r_1;
			coul_energy = charge * p;

		#endif // SKIP_COUL_ENERGY

			return coul_energy;
		}

		static vreal_t dU__dX(_I2T<VDW_>, vreal_t *du__dq, vreal_t sigma2, vreal_t eps,
			vreal_t tau2, vreal_t r_1)
		{
			vreal_t vdw_energy = 0.f;

		#ifndef SKIP_VDW_ENERGY
			vreal_t tau6 = cube(tau2);

		#ifdef USE_PSEUDO_POTENTIAL_IN_ZERO
			if (tau2 < tau2_)
			{
				real_t tau12 = sqr(tau6);
				vdw_energy = eps * (alpha_ + gamma_ * tau12);
				*du__dq += (real_t) 12. * eps * gamma_ * tau12 * r_1;
			}
			else
		#endif // USE_PSEUDO_POTENTIAL_IN_ZERO
			{
				vreal_t tau_6 = vreal_t(1.) / tau6;
				vreal_t eps_tau_6 = eps * tau_6;

				vdw_energy = eps_tau_6 * (tau_6 - vreal_t(2.));
				*du__dq += 12.f * eps_tau_6 * (1.f - tau_6) * r_1;
			}

		#endif // SKIP_VDW_ENERGY

			return vdw_energy;
		}

	};

#ifdef USE_RARE_GAS
	#define CONTROL_CONTACT(is_interaction, atom1, atom2)                                             \
		/* удаление пустых атомов первого массива (во вротом их нет) */                                 \
		is_interaction = multiple(atom1.connect_data > 0, is_interaction);                              \
		is_interaction = multiple(atom2.connect_data > 0, is_interaction);                              \
		if (is_interaction == false) continue;                                                          \

#else
	#define CONTROL_CONTACT(is_interaction, atom1, atom2)                                             \
		{                                                                                               \
			vint_t n = atom1.insert_data - atom2.insert_data;                                             \
			vbool_t f = n < 0; /* условие для коррекции таблиц связности t далее */                       \
			n = abs(n); /* общая операция вместо специализированной для SSE */                            \
			vbool_t f__ = n < 31; /* данные по связности есть только для диапозона [0..30] */             \
			if (f__ == true)                                                                              \
			{ /* есть контакты, удалим взаимодействия с контактами  */                                    \
				vint_t t = multiple(f, atom1.connect_data) + multiple_not(f, atom2.connect_data);           \
				t = multiple(f__, t); /* удаление отсутствующих компонент вектора контактов */              \
				n = multiple(f__, n); /* удаление отсутствующих компонент вектора смещений */               \
				vint_t mask = vint_t(1) << n; /* маска для выделения нужного бита по позиции */             \
				t = t & mask; /* выделение бита контакта */                                                 \
				is_interaction = multiple_not(t > 0, is_interaction);                                       \
				if (is_interaction == false) continue;                                                      \
			}                                                                                             \
			/* удаление пустых атомов первого массива (во вротом их нет) */                               \
			is_interaction = multiple(atom1.connect_data > 0, is_interaction);                            \
			is_interaction = multiple(atom2.connect_data > 0, is_interaction);                            \
			if (is_interaction == false) continue;                                                        \
		}
#endif

	#define LJ_CALCULATE(is_interaction, energy, fx, fy, fz, rmax2, atom1, atom2)                                     \
		vreal_t fx = atom2.x - atom1.x;                                                                 \
		vreal_t fy = atom2.y - atom1.y;                                                                 \
		vreal_t fz = atom2.z - atom1.z;                                                                 \
		region_->make_nearest_image_vector(fx, fy, fz);                                                 \
		vreal_t energy(0.f);                                                                            \
		vbool_t is_interaction(false);                                                                  \
		{                                                                                               \
			vreal_t r2 = sqr(fx) + sqr(fy) + sqr(fz);                                                     \
			is_interaction = r2 < rmax2;                                                                  \
			if (is_interaction == false) continue;                                                        \
			CONTROL_CONTACT(is_interaction, atom1, atom2)                                                 \
			vreal_t qq = atom1.charge * atom2.charge;                                                     \
			vreal_t ss = calculate_sigma(atom1.sigma, atom2.sigma);                                       \
			vreal_t ee = atom1.eps * atom2.eps;                                                           \
			vreal_t tau2 = r2 / ss;                                                                       \
			vreal_t r_1 = vreal_t(1.f) / sqrt(r2);                                                        \
			vreal_t coef_du__dq = 0.f;                                                                    \
			vreal_t coul_energy = _Interaction::dU__dX(_I2T<COUL_>(), &coef_du__dq, ss, qq, tau2, r_1);   \
			vreal_t vdw_energy  = _Interaction::dU__dX(_I2T<VDW_ >(), &coef_du__dq, ss, ee, tau2, r_1);   \
			energy = multiple(is_interaction, coul_energy + vdw_energy);                                  \
			coef_du__dq = multiple(is_interaction, coef_du__dq);                                          \
			coef_du__dq *= r_1;                                                                           \
			fx *= coef_du__dq;                                                                            \
			fy *= coef_du__dq;                                                                            \
			fz *= coef_du__dq;                                                                            \
		}

	#define LJ_INTERACT(is_interaction, rmax2, atom1, atom2)                                          \
		vreal_t fx = atom2.x - atom1.x;                                                                 \
		vreal_t fy = atom2.y - atom1.y;                                                                 \
		vreal_t fz = atom2.z - atom1.z;                                                                 \
		region_.make_nearest_image_vector(fx, fy, fz);                                                  \
		vreal_t r2 = sqr(fx) + sqr(fy) + sqr(fz);                                                       \
		vbool_t is_interaction = r2 < rmax2;                                                            \
		if (is_interaction == false) continue;                                                          \
		CONTROL_CONTACT(is_interaction, atom1, atom2)                                                   \

#ifdef USE_ERF
	typedef Interaction_<E612> Interaction;
#else
	typedef Interaction_<C612> Interaction;
#endif

}
#endif
